import pytest
import torch

from sief_druco.decoders import DECODERS, Decagon

N_SIDE_EFFECTS, D = 963, 16


def build(name, seed=0):
    torch.manual_seed(seed)
    base_rate = torch.rand(N_SIDE_EFFECTS) * 0.1 + 1e-3
    return DECODERS[name](D, base_rate)


def embeddings(n=50, seed=1):
    g = torch.Generator().manual_seed(seed)
    return torch.randn(n, D, generator=g), torch.randn(n, D, generator=g)


@pytest.mark.parametrize("name", DECODERS)
def test_decoder_is_symmetric_and_gives_full_rows(name):
    dec = build(name)
    z_a, z_b = embeddings()
    out = dec(z_a, z_b)
    assert out.shape == (50, N_SIDE_EFFECTS)
    torch.testing.assert_close(out, dec(z_b, z_a), rtol=1e-5, atol=1e-4)


@pytest.mark.parametrize("name", DECODERS)
def test_triples_are_the_matching_entries_of_full_rows(name):
    dec = build(name)
    z_a, z_b = embeddings()
    r = torch.randint(N_SIDE_EFFECTS, (50,))
    expected = dec(z_a, z_b)[torch.arange(50), r]
    torch.testing.assert_close(dec.triples(z_a, z_b, r), expected, rtol=1e-5, atol=1e-4)


def test_decagon_one_matmul_equals_the_broadcast_form():
    dec = build("decagon")
    z_a, z_b = embeddings()
    # z_aᵀ D_r R D_r z_b + b_r, broadcast to [n, n_side_effects, d].
    left = z_a[:, None, :] * dec.D  # [n, r, d]
    right = z_b[:, None, :] * dec.D
    broadcast = torch.einsum("nri,ij,nrj->nr", left, dec.R, right) + dec.b
    torch.testing.assert_close(dec(z_a, z_b), broadcast, rtol=1e-4, atol=1e-3)


def test_decagon_start_values():
    torch.manual_seed(0)
    dec = Decagon(64, torch.full((N_SIDE_EFFECTS,), 0.02))
    assert dec.D.std().item() == pytest.approx(1, rel=0.05)
    assert dec.R0.std().item() == pytest.approx(1 / 8, rel=0.05)
    torch.testing.assert_close(dec.R, dec.R.T)


@pytest.mark.parametrize("name", DECODERS)
def test_bias_starts_at_the_logit_of_the_base_rate(name):
    base_rate = torch.linspace(0.001, 0.2, N_SIDE_EFFECTS)
    dec = DECODERS[name](D, base_rate)
    torch.testing.assert_close(torch.sigmoid(dec.b), base_rate)
