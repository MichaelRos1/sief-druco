# Papers

Key papers for the project. PDFs are git-ignored (`papers/*.pdf`); this index is tracked.
Status as of 2026-10-28. "Missing" = not downloadable by script (publisher or Europe PMC returned 403 / a bot-check page, or not open access).

| # | Paper | Why it matters | Link | PDF |
|---|---|---|---|---|
| 01 | Zitnik, Agrawal & Leskovec 2018, Modeling polypharmacy side effects with graph convolutional networks (Decagon) | dataset, RGCN encoder, decoder, evaluation | https://arxiv.org/abs/1802.00543 | ✅ |
| 02 | Tatonetti et al. 2012, Data-driven prediction of drug effects and interactions (TWOSIDES) | how Decagon's labels were made | https://doi.org/10.1126/scitranslmed.3003377 | missing |
| 03 | Schlichtkrull et al. 2018, Modeling relational data with graph convolutional networks | RGCN, basis decomposition | https://arxiv.org/abs/1703.06103 | ✅ |
| 04 | Malone, García-Durán & Niepert 2018, Knowledge graph completion to predict polypharmacy side effects | tensor-factorisation baseline beating Decagon (#2) | https://arxiv.org/abs/1810.09227 | ✅ (arXiv version) |
| 05 | Xu et al. 2019, TIP: tri-graph information propagation | staged propagation encoder (#2) | https://arxiv.org/abs/2001.10516 | ✅ |
| 06 | Elkan & Noto 2008, Learning classifiers from only positive and unlabeled data | PU learning | https://doi.org/10.1145/1401890.1401920 | ✅ (author copy; title not machine-checked) |
| 07 | Bekker & Davis 2018, Beyond the selected completely at random assumption | PU with biased reporting | https://arxiv.org/abs/1809.03207 | ✅ |
| 08 | Dal Pozzolo et al. 2015, Calibrating probability with undersampling for unbalanced classification | negative-sampling correction of `p` | https://doi.org/10.1109/SSCI.2015.33 | ✅ (author copy) |
| 09 | Gottlieb et al. 2015, Ranking adverse drug reactions with crowdsourcing | severity weights (#14, #20) | https://doi.org/10.2196/jmir.3962 | missing |
| 10 | Cheng, Kovács & Barabási 2019, Network-based prediction of drug combinations | why protein-network proximity matters | https://doi.org/10.1038/s41467-019-09186-x | missing |
| 11 | Guo et al. 2017, On calibration of modern neural networks | temperature scaling, ECE (#9) | https://arxiv.org/abs/1706.04599 | ✅ |
| 12 | Szklarczyk et al. 2016, STITCH 5 | flat vs. stereo compound IDs (#12) | https://doi.org/10.1093/nar/gkv1277 | ✅ |
| 13 | Kontsioti et al. 2022, A reference set of clinically relevant adverse drug-drug interactions (CRESCENDDI) | side-effect-level "documented" reference (#13, #19) | https://doi.org/10.1038/s41597-022-01159-y | missing |
| 14 | Xiong et al. 2022, DDInter | pair-level "documented" reference (#13, #19) | https://doi.org/10.1093/nar/gkab880 | missing |
| 15 | DDInter 2.0 (NAR 2025) | newer DDInter release (#13) | https://doi.org/10.1093/nar/gkae726 | missing |
| 16 | Phansalkar et al. 2012, High-priority drug–drug interactions for use in EHRs (ONC list) | positive controls (#13, `7.1`) | https://doi.org/10.1136/amiajnl-2011-000612 | missing |
| 17 | Kee et al. 2024, Are polypharmacy side effects predicted by public data still valid in real-world data? | closest prior work (#10) | https://doi.org/10.1016/j.heliyon.2024.e24620 | missing |
| 18 | Sutherland et al. 2015, Co-prescription trends in a large cohort of subjects predict substantial drug-drug interactions | known DDI × NHANES co-prescription (#10) | https://doi.org/10.1371/journal.pone.0118991 | ✅ |
| 19 | Dagenais et al. 2024, framework to quantify potential impact of pharmacokinetic DDIs using real-world data | prevalence × risk class (#10) | https://doi.org/10.1111/cts.13741 | missing |
| 20 | Hecker et al. 2023, deep learning DDI prediction for patients with multiple sclerosis | predicted DDIs on real medication plans (#10) | https://doi.org/10.3390/pharmaceutics16010003 | missing |
| 21 | Qato et al. 2008, Use of prescription and OTC medications and dietary supplements among older adults in the US | co-prescription prevalence (#10) | https://doi.org/10.1001/jama.2008.892 | missing |
| 22 | Qato et al. 2016, Changes in prescription and OTC medication use among older adults, 2005 vs 2011 | co-prescription prevalence (#10) | https://doi.org/10.1001/jamainternmed.2015.8581 | missing |
| 23 | Sánchez-Valle et al. 2024, Prevalence and differences in the co-administration of drugs known to interact (BMC Med) | known-DDI co-administration in three populations (#10) | https://doi.org/10.1186/s12916-024-03384-1 | missing |
| 24 | Saha et al. 2024, Adverse event classification from co-prescribed drugs (LNCS) | possible prior work, unread lead (#10) | https://doi.org/10.1007/978-3-031-12700-7_35 | missing (paywalled) |

Data documentation (not PDFs): NHANES 2017–18 Prescription Medications, https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/RXQ_RX_J.htm
