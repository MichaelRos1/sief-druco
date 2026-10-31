# Polypharmacy Side-Effect Burden

Predicting side effects of drugs taken together, and weighting those predictions by how often people actually take them together.

## Language

**Drug pair**:
Exactly two distinct drugs taken together; the unit every prediction, prevalence and burden figure is about.
_Avoid_: Combination, combo, interaction

**Drug set**:
Three or more drugs taken together, such as one survey participant's full medication list; scored only by combining its drug pairs' scores.
_Avoid_: Combination, regimen, cocktail

**Labelled pair**:
A drug pair with at least one side effect in the Decagon data.
_Avoid_: Positive pair, interacting pair, known interaction

**Unlabelled pair**:
A drug pair with no side effect in the Decagon data; this means "not reported", not "safe".
_Avoid_: Negative pair, safe pair, non-interacting pair

**Burden**:
How many people take a drug pair, as a share of the population, times the pair's predicted side-effect signal; it measures exposure to a signal, not cases of harm.
_Avoid_: Risk, harm, incidence

**Novelty candidate**:
A drug pair, or pair × side effect, that is neither known to Decagon nor clinically documented; a hypothesis, not a finding.
_Avoid_: New interaction, undiscovered interaction, novel DDI

**Not in the model**:
A survey drug that reached an RxNorm ingredient outside Decagon's 645 drugs; the model can't score it, and that is reported, not fixed.
_Avoid_: Unmapped, missing drug

**Mapping failed**:
A survey drug that never reached an RxNorm ingredient, or a Decagon drug that never reached RxNorm; a gap in the name matching, not in the model.
_Avoid_: Unmapped, not in Decagon
