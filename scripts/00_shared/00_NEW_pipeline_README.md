# Sex-stratified RA biomarker pipeline — run order

**This pipeline is sex-STRATIFIED only.** Every model is fitted separately
within each sex, contrasting RA vs Control inside that sex. The sex-SPECIFIC
(diagnosis x sex interaction) branch was removed on 2026-07-24 — see
`results/README_GOALS.md` for what went and why.

The candidate step uses the **data-driven disease modules — green + brown — for
BOTH sexes** (|cor_RA| >= 0.5 & p < 1e-8 in the overall RA-vs-control
module–trait analysis), intersected with the sex-stratified DEGs. Outputs go to
`data/processed/` (+ `new/`), `results/tables/`, `results/figures/` (+ `new/`,
`current/`).

## `scripts/00_shared/` — shared upstream, run first

| # | Script | Purpose | Key output |
|---|--------|---------|-----------|
| 01 | `01_packages.R` | install/load all packages | — |
| 01/02 | `01_data_loading.R`, `02_data_loading.R` | GEO series -> raw RDS | `data/raw/GSE*_raw.rds` |
| 03 | `03_normalize_batch.R` (+`_figure`) | qnorm + ComBat | `combined_train.rds`, `fig_combine_pca_combat` |
| 03 | `03_dataset_gene_overlap_venn.R` | gene overlap across datasets | `fig_dataset_gene_overlap_venn` |
| 04 | `04_apply_holdout.R` | 70/30 split | `internal_val_holdout_processed.rds` |
| 05 | `05_dge.R` (+`_figure`) | **sex-stratified** limma DEG + enrichment | `DEG_{all,female,male}_*.csv`, `dge_results.rds` |
| 05b | `05b_dge_sensitivity.R` | ComBat-then-DE vs batch-in-model; interaction screen | `DEG_sensitivity_combat_vs_batch.csv`, `DEG_interaction_summary.csv` |
| 05c | `05c_deconvolution.R` | **CIBERSORT/LM22 + MCP-counter cell fractions; composition-adjusted DE** | `CELL_fractions_*.csv`, `DEG_celladjusted_*.csv`, `cell_fractions.rds` |
| 05d | `05d_interaction_report.R` | **diagnosis x sex interaction, reported in full** | `DEG_interaction_{full,significant,patterns}.csv` |
| 06 | `06_WGCNA.R` (+`_figure`) | modules, eigengenes, kME, preservation | `wgcna_results.rds`, `WGCNA_*.csv` |
| 08 | `08_module_trait_RA_control.R` (+`_figure`) | module–trait RA vs Control, per sex | `module_trait_RAvsControl_{all,female,male}.csv` |
| 09 | `09_disease_module_deg_intersect.R` | (green ∪ brown) ∩ sex-DEG -> candidates | `candidates_{female,male}_disease.csv` |
| 09 | `09_greenmodule_deg_venn.R` | green module ∩ DEG venns | `greenmod_DEG_intersection_*.csv` |
| 10 | `10_MR.R` | **two-sample MR, end to end**: eQTLGen cis-eQTL instruments -> Okada RA outcome (chunked) -> IVW/Egger/median -> within-stratum FDR | `MR_{sex}_TABLE1-4.csv`, `MR_causal_FDR_{sex}.csv`, `FS_input_{female,male}.csv` |
| 10c | `10c_MR_mhc_sensitivity.R` | **MHC-excluded re-run of the whole MR, as a parallel column** (offline, from the 10_MR.R cache) | `MR_MHC_sensitivity_{female,male}.csv`, `_summary.csv`, `_panel_fate.csv`, `FS_input_*_noMHC.csv` |
| 10d | `10d_coloc_panel_genes.R` | **coloc.abf on regional eQTLGen vs Okada stats** for every causal + panel gene; two priors | `COLOC_results.csv`, `COLOC_panel_genes.csv`, `COLOC_summary.csv` |
| 30 | `30_methodology_flowchart.R` | methods flowchart | `fig_methodology_flowchart` |
| — | `eda.R` (+`_figure`) | exploratory QC | — |

> `10b` is the authoritative MR run. `10` and `11` are its cached / resume
> predecessors and produce the same result.

## `scripts/goal2_sex_stratified/` — the panel branch

| # | Script | Purpose | Key output |
|---|--------|---------|-----------|
| 12 | `12_feature_selection.R` | LASSO ∩ RF ∩ SVM-RFE, **within each sex** | `ml_features.rds` (7F / 4M) |
| 13 | `13_roc_validation_train_internal_blood.R` | ROC: train / internal / external blood | `mr_roc_panel_auc.csv` |
| 13b | `13b_panel_auc_celladjusted.R` | **panel vs the composition-only benchmark**: 4 models + LRT, per sex and dataset | `PANEL_auc_celladjusted.csv`, `PANEL_incremental_value_LRT.csv`, `PANEL_auc_celladjusted_summary.csv` |
| 14 | `14_nested_cv.R` | leakage-free nested CV | `mr_nested_cv_summary.csv` |
| 15 | `15_elasticnet_panel.R` | elastic-net over full MR set | `mr_elasticnet_summary.csv` |
| 16 | `16_final_panel.R` | final recommended models | `mr_final_panel_summary.csv` |
| 17 | `17_final_roc_figure.R` | final ROC per sex | `fig_mr_final_roc_{female,male}` |
| 18 | `18_pergene_roc_overlay.R` | per-gene ROC (blood/train) | `fig_mr_pergene_roc_{female,male}` |
| 19 | `19_fig4A_expression_groups.R` | Fig 4A expression by group | `fig_mr_gene_expression_groups` |
| 19 | `19_nested_cv.R` (+`_figure`) | nested-CV pooled ROC | `fig_nested_cv_roc` |
| 20 | `20_fig4BC_logfc_heatmaps.R` | Fig 4B/C logFC heatmaps | `fig_mr_panelB/C` |
| 22 | `22_fig4C_clustered_heatmap.R` | clustered per-sample heatmap | `fig_mr_panelC_clustered_heatmap` |
| 23 | `23_fig4_composite.R` | assemble Fig 4 **A–C** | `fig_mr_NEW_FIG4_composite` |
| 24 | `24_diagnostic_model_validation.R` | nomogram / calibration / DCA / CIC | `fig_diag_validation_{female,male}`, `diag_dca_*` |
| 25 | `25_pathway_enrichment.R` | GO/KEGG on causal genes | `fig_enrich_female`, `enrich_*` |
| 28 | `28_validate_synovium.R` | cross-tissue synovium (GSE89408) | `val_synovium_pergene_*.csv` |
| 29 | `29_crosstissue_summary.R` | synovium panel AUC + concordance | `fig_crosstissue` |
| 30 | `30_synovium_roc.R` | synovium per-gene + panel ROC | `fig_syn_pergene_roc_{female,male}` |
| 31 | `31_allvalidation_roc.R` | panel ROC across all 4 datasets | `fig_allvalidation_roc` |
| 32 | `32_pergene_auc_alltissues.R` | per-gene AUC bars, 4 datasets | `fig_pergene_auc_alltissues` |
| 33 | `33_pergene_roc_alltissues.R` | per-gene ROC, 4 datasets | `fig_pergene_roc_alltissues` |
| 35 | `35_crossancestry_eqtl_mr.R` | cross-ancestry MR (EUR / EAS) | `MR35_crossancestry_*.csv` |
| 36 | `36_crossancestry_eqtl_figure.R` | cross-ancestry figures | `fig_mr35_crossancestry_*` |
| 37 | `37_mr_diagnostic_figures.R` | MR forest / SNP support / funnel / LOO | `FIG_MR_01`–`FIG_MR_04` |

## Removed 2026-07-24, and what replaced it (2026-07-28)

- **Step 12b** `sex_interaction_panel_test` — the diagnosis x sex interaction
  test. **Superseded, not abandoned:** the interaction is now tested properly in
  `05b_dge_sensitivity.R` and reported in full by `05d_interaction_report.R`
  (53 genes at FDR < 0.05 on the pre-ComBat matrix with batch modelled).
- **Step 21** `fig4DE_cibersort` and `26_cibersort` — CIBERSORT Female-RA vs
  Male-RA. Fig 4 is A–C only. **Deconvolution itself is restored** in
  `05c_deconvolution.R`, for a different and more important purpose: measuring
  whether the RA signature and the diagnostic panel are confounded by leukocyte
  composition rather than comparing composition between sexes.
- **Steps 26/27** `gene_correlation_bysex`, `curated_correlations` — between-sex
  gene–gene correlation comparisons. Still removed.

## Robustness layer added 2026-07-28

Five analyses that the primary pipeline assumed rather than tested. Each was
written so its conclusion could be negative, and two of them are:

| Script | Question | Answer |
|---|---|---|
| `10c` | Do the MR causal genes survive excluding the MHC? | **Partly.** 32 -> 14 female, 25 -> 14 male. 3 of 6 genes in each panel are MHC-confounded. |
| `10d` | Do the eQTL and RA signals share a causal variant? | **No.** 0 of 9 panel genes colocalise; 6 of 9 show PP.H3 >= 0.8, i.e. positive evidence of distinct variants. |
| `05c` | Is the DE signature a differential white-cell count? | **Partly.** DEGs 5,131 -> 2,709 (F) and 5,820 -> 1,450 (M) after adjustment; 10 of 12 panel gene-sex pairs survive. |
| `13b` | Does the panel beat a composition-only model? | **Yes, in females.** Train 0.831 vs 0.662; external 1.000 vs 0.429. Internal test is marginal (0.721 vs 0.714). |
| `05d` | Does RA differ transcriptionally between the sexes? | **Yes, at a floor of 53 genes** — but only 7 survive composition adjustment. |

## Notes

- **Step 10 -> 11:** step 10 extracts instruments (OpenGWAS, ~1 s/gene). Its own
  MR call can time out on the outcome fetch; step 11 resumes from the saved
  instruments with a chunked outcome query and writes `FS_input` + per-sex XLSX.
- Headline result: **female 6-gene consensus panel** (C6orf136, ESYT1, GNL1,
  IKZF3, MED1, SMARCC2), nested-CV AUC 0.809 (0.738–0.879). The **male 6-gene
  panel** (ESYT1, HLA-DMA, INPP5B, MED1, SMARCC2, VPS52) is **EXPLORATORY**:
  train n = 38, internal n = 13, external n = 9. Its internal AUC of 1.000 is a
  perfect-separation artefact and is flagged `SEPARATION` in
  `mr_final_panel_summary.csv`. It is not a validated panel.
- The panels overlap by half (ESYT1, MED1, SMARCC2 in both) and carry
  **identical MR estimates by construction**, because both GWAS inputs are
  sex-combined. Sex-stratification enters only upstream, at the DE step.
- Disease modules are **yellow (up in RA) + brown (down in RA)**. Earlier drafts
  said green+brown; that was the variance-filtered configuration and the colour
  names are size-rank labels, not stable identities.
- Figure -> script mapping: `results/FIGURE_PROVENANCE.md` — **not present in
  this tree**; neither is `scripts/build_figure_provenance.py`. Regenerate or
  delete the reference before deposition.
