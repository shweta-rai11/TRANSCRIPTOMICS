# Sex-stratified RA biomarker pipeline — run order

**This pipeline is sex-STRATIFIED.** Every model is fitted separately within each
sex, contrasting RA vs Control inside that sex. The panels are not sex-SPECIFIC:
a gene selected in one sex is not shown to behave differently in the other. The
diagnosis × sex interaction *is* tested — separately, in `05d` — and it is
positive (53 genes), but no panel gene is among them.

Scripts are numbered in run order and named by **workflow stage**, so the
filename says what the step does:

```
FEATURE SELECTION → VENN → MODEL TRAINING → INTERNAL TESTING (BLOOD)
  → EXTERNAL TESTING (BLOOD) → EXTERNAL TESTING (SYNOVIUM)
  → CROSS-TISSUE BIOMARKER DISCOVERY → CROSS-ANCESTRY BIOMARKER
```

Nothing here `source()`s anything else — every script is standalone and reads
from `data/processed/`, so any stage can be re-run on its own.

---

## What this pipeline claims, and what it does not

This is the single most important section for anyone reading the results.

| Claim | Status |
|---|---|
| A sex-stratified **diagnostic** panel discriminates RA from control in women | **Supported.** Nested-CV AUC 0.78–0.81, and it beats a cell-composition-only benchmark on all three blood datasets. |
| The panel is more than a differential white-cell count | **Supported** (`17b`). |
| The panel genes are **causal** for RA | **NOT supported.** No panel gene colocalises with the RA association (`10d`). Genes are described as **MR-prioritised**, never as causal. |
| The **male** panel is a validated diagnostic | **NOT supported.** n = 38 / 13 / 9. Reported as EXPLORATORY throughout. |
| RA differs transcriptionally between the sexes | **Supported at a floor of 53 genes** (`05d`), but only 7 survive cell-composition adjustment. |

The MR step is a **genetically-informed filter** that narrowed ~2,045 candidates
to a few dozen. That is a legitimate and unusual way to build a candidate set. It
is not evidence of causality, and the pipeline no longer says it is.

---

## `scripts/00_shared/` — shared upstream, run first

| # | Script | Purpose | Key output |
|---|--------|---------|-----------|
| 01 | `01_packages.R` | install/load all packages | — |
| 02 | `02_data_loading.R` | GEO series → raw RDS | `data/raw/GSE*_raw.rds` |
| 03 | `03_normalize_batch.R` (+`_figure`) | 70:30 split **first**, then qnorm + ComBat on train only | `combined_train.rds`, `internal_val_holdout.rds` |
| 03 | `03_dataset_gene_overlap_venn.R` | gene overlap across datasets | `fig_dataset_gene_overlap_venn` |
| 04 | `04_apply_holdout.R` | project the sealed hold-out using **frozen** train parameters | `internal_val_holdout_processed.rds` |
| 05 | `05_dge.R` (+`_figure`) | sex-stratified limma DEG + enrichment | `DEG_{all,female,male}_*.csv` |
| 05b | `05b_dge_sensitivity.R` | ComBat-then-DE vs batch-in-model; TREAT thresholds; interaction screen | `DEG_sensitivity_combat_vs_batch.csv` |
| **05c** | `05c_deconvolution.R` | **CIBERSORT/LM22 + MCP-counter; composition-adjusted DE** | `CELL_fractions_*.csv`, `DEG_celladjusted_*.csv` |
| **05d** | `05d_interaction_report.R` | **diagnosis × sex interaction, reported in full** | `DEG_interaction_{full,significant,patterns}.csv` |
| 06 | `06_WGCNA.R` (+`_figure`) | modules, eigengenes, kME, 200-perm preservation | `wgcna_results.rds`, `WGCNA_*.csv` |
| 08 | `08_module_trait_RA_control.R` | module–trait RA vs Control, per sex | `module_trait_RAvsControl_*.csv` |
| 09 | `09_disease_module_deg_intersect.R` (+`09b` venn) | (yellow ∪ brown) ∩ sex-DEG → candidates | `candidates_{female,male}_disease.csv` |
| 10 | `10_MR.R` | two-sample MR: eQTLGen cis-eQTL → Okada RA, within-stratum FDR | `FS_input_{female,male}.csv` |
| **10c** | `10c_MR_mhc_sensitivity.R` | **MHC-excluded re-run of the whole MR** (offline, from the `10_MR.R` cache) | `MR_MHC_sensitivity_*.csv`, `FS_input_*_noMHC.csv` |
| **10d** | `10d_coloc_panel_genes.R` | **colocalisation (coloc.abf, two priors) for every prioritised gene** | `COLOC_results.csv`, `COLOC_panel_genes.csv` |
| **10e** | `10e_coloc_susie_mhc.R` | **multi-causal-variant colocalisation (coloc.susie) for the MHC**, where coloc.abf's assumption fails | `COLOC_SUSIE_mhc.csv`, `COLOC_combined_abf_susie.csv` |
| 30 | `30_methodology_flowchart.R` | methods flowchart | `fig_methodology_flowchart` |

---

## `scripts/goal2_sex_stratified/` — the panel branch, by workflow stage

### Stage 1 — FEATURE SELECTION

| # | Script | Purpose | Key output |
|---|--------|---------|-----------|
| 12 | `12_feature_selection.R` | LASSO ∩ RF ∩ SVM-RFE within each sex, on the MR-prioritised set | `ml_features.rds` — **6F / 6M** |
| 12b | `12b_feature_selection_noMHC.R` | identical procedure on the **MHC-free** set | `ml_features_noMHC.rds` — **4F / 5M** |

### Stage 2 — VENN DIAGRAM

| # | Script | Purpose | Key output |
|---|--------|---------|-----------|
| 13 | `13_feature_selection_venn.R` | three-selector overlap per sex, both candidate sets, plus female-vs-male membership | `FIG_G2_01`–`FIG_G2_04`, `FS_venn_membership.csv` |

### Stage 3 — MODEL TRAINING (female and male separately)

| # | Script | Purpose | Key output |
|---|--------|---------|-----------|
| 14 | `14_model_training_nested_cv.R` | leakage-free nested CV; selection redone in every fold | `mr_nested_cv_summary.csv` |
| 15 | `15_model_training_elasticnet.R` | elastic net over the full prioritised set | `mr_elasticnet_summary.csv` |
| 16 | `16_model_training_final_panel.R` | locked final models, **with evidence tiers** | `mr_final_panel_summary.csv` |
| 16b | `16b_model_training_final_panel_noMHC.R` | **primary vs MHC-free, head-to-head under nested CV** | `PANEL_primary_vs_noMHC_{nestedcv,delong,performance}.csv` |
| 16c | `16c_model_training_nested_cv_transcriptomewide.R` | transcriptome-wide LASSO baseline (no MR prior) | `fig_nested_cv_roc` |
| **16d** | `16d_nested_cv_reconciliation.R` | **the ONE authoritative nested-CV table** — every candidate-set x selector variant, one process, one seed policy | `NESTED_CV_AUTHORITATIVE.csv`, `NESTED_CV_legacy_reconciliation.csv` |

### Stage 4 + 5 — INTERNAL TESTING (BLOOD) and EXTERNAL TESTING (BLOOD)

| # | Script | Purpose | Key output |
|---|--------|---------|-----------|
| 17 | `17_testing_blood_internal_external.R` | ROC across train / internal hold-out / external blood | `mr_roc_panel_auc.csv` |
| 17b | `17b_testing_blood_celladjusted.R` | **panel vs a composition-only benchmark**, both panels, 4 models + LRT | `PANEL_auc_celladjusted*.csv` |
| 18 | `18_testing_blood_pergene_roc.R` | per-gene ROC overlays | `fig_mr_pergene_roc_*` |
| 19 | `19_testing_blood_clinical_utility.R` | nomogram, calibration, decision curve, clinical impact | `fig_diag_validation_*` |

### Stage 6 — EXTERNAL TESTING (SYNOVIUM)

| # | Script | Purpose | Key output |
|---|--------|---------|-----------|
| 20 | `20_testing_synovium_external.R` | cross-tissue validation in GSE89408 | `val_synovium_pergene_*.csv` |
| 21 | `21_testing_synovium_roc.R` | synovium per-gene and panel ROC | `fig_syn_pergene_roc_*` |

### Stage 7 — CROSS-TISSUE BIOMARKER DISCOVERY

| # | Script | Purpose | Key output |
|---|--------|---------|-----------|
| 22 | `22_crosstissue_biomarker_discovery.R` | blood↔synovium concordance and panel AUC | `fig_crosstissue` |
| 23 | `23_crosstissue_roc_all_datasets.R` | panel ROC across all four datasets | `fig_allvalidation_roc` |
| 24 | `24_crosstissue_pergene_auc.R` | per-gene AUC bars, four datasets | `fig_pergene_auc_alltissues` |
| 25 | `25_crosstissue_pergene_roc.R` | per-gene ROC, four datasets | `fig_pergene_roc_alltissues` |

### Stage 8 — CROSS-ANCESTRY BIOMARKER

| # | Script | Purpose | Key output |
|---|--------|---------|-----------|
| 26 | `26_crossancestry_biomarker_mr.R` | European eQTL exposures vs East-Asian RA outcome | `MR35_crossancestry_*.csv` |
| 27 | `27_crossancestry_biomarker_figure.R` | cross-ancestry figures | `fig_mr35_crossancestry_*` |

> Cross-ancestry MR uses European exposures against an East-Asian outcome. That
> is exploratory support, **not** like-for-like replication, and must be
> described as such.

### Supporting — figures and enrichment

| # | Script | Purpose |
|---|--------|---------|
| 28–33 | `28_figure_final_roc.R`, `29_figure_nested_cv_roc.R`, `30_figure_fig4A_expression_groups.R`, `31_figure_fig4BC_logfc_heatmaps.R`, `32_figure_fig4C_clustered_heatmap.R`, `33_figure_fig4_composite.R` | Figure 4 components and composite |
| 34 | `34_pathway_enrichment.R` | GO/KEGG on the prioritised genes |
| 35 | `35_figure_mr_diagnostics.R` | MR forest, SNP support, funnel, leave-one-out |

---

## The robustness layer — why this pipeline is defensible

Five analyses test assumptions that biomarker pipelines of this kind normally
make silently. **Each was written so that its answer could be negative, and two
of them are.** That is the point: a result that survives a test designed to break
it is worth more than one that was never tested.

| Script | Question it could have failed | Answer |
|---|---|---|
| `10c` | Do the prioritised genes survive excluding the MHC? | **Partly.** 32→14 female, 25→14 male. 3 of 6 genes in each panel were MHC-affected. |
| `10d` | Do the eQTL and RA signals share a causal variant? | **No.** 0 of 9 panel genes colocalise. Outside the MHC, 9 genes show PP.H3 ≥ 0.8. **The causal claim was withdrawn on this evidence.** |
| `05c` | Is the DE signature a differential white-cell count? | **Partly.** DEGs 5,131→2,709 (F), 5,820→1,450 (M). 10 of 12 panel gene-sex pairs survive. |
| `17b` | Does the panel beat a composition-only model? | **Yes.** The MHC-free female panel wins on **all three** blood datasets. |
| `05d` | Does RA differ transcriptionally between the sexes? | **Yes, at a floor of 53 genes** — but only 7 survive composition adjustment, and no panel gene is among them. |
| `16b` | Does removing the MHC cost diagnostic performance? | **No.** Nested CV 0.805→0.781 (F) and 0.796→0.734 (M); DeLong p = 0.35 / 0.36. |

### The headline consequence

**The MHC-free panel is the one to carry forward.** It is statistically
indistinguishable from the original under nested CV, it is free of HLA linkage
disequilibrium, and it is *better behaved* against the composition benchmark —
beating cell fractions on all three blood datasets where the original panel was
marginal on the internal test (Δ AUC 0.074 vs 0.007), and showing no
perfect-separation artefact externally.

| Panel | Female genes | Male genes |
|---|---|---|
| Primary (MHC retained) | C6orf136, ESYT1, GNL1, IKZF3, MED1, SMARCC2 | ESYT1, HLA-DMA, INPP5B, MED1, SMARCC2, VPS52 |
| **MHC-free (recommended)** | **CDC37, IKZF3, MED1, SMARCC2** | **MED1, NCOA5, PHF19, SMARCC2, TAB1** |

Report both. Do not quietly substitute one for the other.

---

## Reporting rules (non-negotiable)

1. **"MR-prioritised", never "MR-causal".** No gene survives colocalisation.
0. **Cite `NESTED_CV_AUTHORITATIVE.csv` for every nested-CV figure.** Nothing else. The recommended row is `noMHC` x `consensus`.
2. Every panel AUC is reported with its **composition-only benchmark** and the
   LRT of panel+composition against composition alone.
3. The **male arm is EXPLORATORY** in every table, figure and sentence. Never
   quote an interval flagged `SEPARATION` — at n = 13 a perfect AUC is an
   artefact of 42 case-control pairs, not a performance estimate.
4. Report **53** interaction genes from the pre-ComBat model, not the inflated
   270 from the ComBat matrix, and state that only 7 survive composition
   adjustment.
5. A gene may be called **cell-intrinsic** only if it survives `05c` adjustment.
6. Disease modules are **yellow (up) + brown (down)**. Colour names are
   size-rank labels, not stable identities.
7. Inside the MHC, coloc.abf's single-causal-variant assumption is violated, so
   PP.H3 there is **not** clean evidence of distinct variants — the correct
   statement is that MR and coloc are *both* unreliable in that region.

## Known gaps

- `results/README_GOALS.md`, `results/FIGURE_PROVENANCE.md`,
  `scripts/build_figure_provenance.py` and `scripts/AUDIT_independent_checks.R`
  are referenced in places but **do not exist** in this tree.
- `data/raw` is a **symlink** outside the project; the tree is not self-contained.
- No medication or disease-activity adjustment.
- MHC colocalisation is only partly settled: `coloc.susie` (`10e`) resolved 3 of
  14 genes; the other 11 need in-sample LD, which is not released by eQTLGen or
  Okada.
