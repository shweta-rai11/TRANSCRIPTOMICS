# Results — verification report

_Generated 2026-09-11 by `scripts/verify_results_numbers.py`. Do not hand-edit — re-run the script to regenerate. Every number below is read directly out of a CSV in `results/tables/` (or a run log in `results/logs/`); none is retyped by hand._

**333 claims across 15 sections**, organised by the canonical §2.1–§2.15 scheme in `results/METHODS_00_INDEX.md`. Each claim shows: the value, its claim ID (also a row in `results/RESULTS_PROVENANCE.tsv`), the exact source CSV/log it was read from, the R script that produced that source file, and the derivation (which column/filter was applied).

## Contents

- [§2.1 Datasets](#21-datasets) (13 claims)
- [§2.2 Data pre-processing, leakage-safe partition, normalisation and batch correction](#22-data-pre-processing-leakage-safe-partition-normalisation-and-batch-correction) (8 claims)
- [§2.3 Differential gene expression analysis](#23-differential-gene-expression-analysis) (48 claims)
- [§2.4 Co-expression network analysis (WGCNA)](#24-co-expression-network-analysis-wgcna) (35 claims)
- [§2.5 Candidate gene identification (disease module ∩ sex-stratified DEG)](#25-candidate-gene-identification-disease-module-∩-sex-stratified-deg) (11 claims)
- [§2.6 Mendelian randomisation, incl. MHC sensitivity analysis](#26-mendelian-randomisation-incl-mhc-sensitivity-analysis) (23 claims)
- [§2.7 Bayesian colocalisation (coloc.abf, coloc.susie)](#27-bayesian-colocalisation-colocabf-colocsusie) (18 claims)
- [§2.8 Sex-stratified feature selection (LASSO / RF / SVM-RFE consensus)](#28-sex-stratified-feature-selection-lasso--rf--svm-rfe-consensus) (11 claims)
- [§2.9 Diagnostic model development and evaluation](#29-diagnostic-model-development-and-evaluation) (72 claims)
- [§2.10 Diagnosis-by-sex interaction testing](#210-diagnosis-by-sex-interaction-testing) (13 claims)
- [§2.11 Cross-tissue evaluation (synovium)](#211-cross-tissue-evaluation-synovium) (16 claims)
- [§2.12 Cross-ancestry evaluation](#212-cross-ancestry-evaluation) (2 claims)
- [§2.13 Functional enrichment](#213-functional-enrichment) (24 claims)
- [§2.14 Immune deconvolution and composition-adjusted expression](#214-immune-deconvolution-and-composition-adjusted-expression) (23 claims)
- [§2.15 Nomogram construction and clinical evaluation](#215-nomogram-construction-and-clinical-evaluation) (16 claims)

---

## §2.1 Datasets

- **Training set size:** 183
  <br>*[R-001]* — source: `results/tables/combined_cohort_summary.csv` (sum of Freq) · script: `scripts/00_shared/03_normalize_batch.R`
- **Training RA:** 103
  <br>*[R-002]* — source: `results/tables/combined_cohort_summary.csv` (sum Freq where group==RA) · script: `scripts/00_shared/03_normalize_batch.R`
- **Training HC:** 80
  <br>*[R-003]* — source: `results/tables/combined_cohort_summary.csv` (sum Freq where group==HC) · script: `scripts/00_shared/03_normalize_batch.R`
- **Training female:** 145
  <br>*[R-004]* — source: `results/tables/combined_cohort_summary.csv` (sum Freq where sex==F) · script: `scripts/00_shared/03_normalize_batch.R`
- **Training male:** 38
  <br>*[R-005]* — source: `results/tables/combined_cohort_summary.csv` (sum Freq where sex==M) · script: `scripts/00_shared/03_normalize_batch.R`
- **Training samples from GSE110169:** 111
  <br>*[R-006.GSE110169]* — source: `results/tables/combined_cohort_summary.csv` (sum Freq where dataset==GSE110169) · script: `scripts/00_shared/03_normalize_batch.R`
- **Training samples from GSE93272:** 72
  <br>*[R-006.GSE93272]* — source: `results/tables/combined_cohort_summary.csv` (sum Freq where dataset==GSE93272) · script: `scripts/00_shared/03_normalize_batch.R`
- **Hold-out size:** 74
  <br>*[R-010]* — source: `results/tables/internal_val_holdout_meta.csv` (row count) · script: `scripts/00_shared/03_normalize_batch.R`
- **Hold-out female RA:** 36
  <br>*[R-011]* — source: `results/tables/internal_val_holdout_meta.csv` (count group==RA & sex==F) · script: `scripts/00_shared/03_normalize_batch.R`
- **Hold-out female HC:** 25
  <br>*[R-012]* — source: `results/tables/internal_val_holdout_meta.csv` (count group==HC & sex==F) · script: `scripts/00_shared/03_normalize_batch.R`
- **Hold-out male RA:** 6
  <br>*[R-013]* — source: `results/tables/internal_val_holdout_meta.csv` (count group==RA & sex==M) · script: `scripts/00_shared/03_normalize_batch.R`
- **Hold-out male HC:** 7
  <br>*[R-014]* — source: `results/tables/internal_val_holdout_meta.csv` (count group==HC & sex==M) · script: `scripts/00_shared/03_normalize_batch.R`
- **Combined discovery cohort (train + hold-out):** 257
  <br>*[R-015]* — source: `results/tables/combined_cohort_summary.csv + results/tables/internal_val_holdout_meta.csv` (sum of both) · script: `scripts/00_shared/03_normalize_batch.R`

## §2.2 Data pre-processing, leakage-safe partition, normalisation and batch correction

- **Genes after probe collapse, GSE93272:** 20848
  <br>*[R-020.GSE93272]* — source: `results/tables/dataset_gene_overlap.csv` (column n_genes) · script: `scripts/00_shared/03_dataset_gene_overlap_venn.R`
- **Genes after probe collapse, GSE110169:** 19041
  <br>*[R-020.GSE110169]* — source: `results/tables/dataset_gene_overlap.csv` (column n_genes) · script: `scripts/00_shared/03_dataset_gene_overlap_venn.R`
- **Genes after probe collapse, common:** 15763
  <br>*[R-020.common]* — source: `results/tables/dataset_gene_overlap.csv` (column n_genes) · script: `scripts/00_shared/03_dataset_gene_overlap_venn.R`
- **SD of per-sample medians, before_qnorm:** 0.574
  <br>*[R-030.before_qnorm]* — source: `results/tables/normalization_diagnostics.csv` (row stage==before_qnorm, column median_sd) · script: `scripts/00_shared/03_normalize_batch.R`
- **SD of per-sample IQRs, before_qnorm:** 0.799
  <br>*[R-031.before_qnorm]* — source: `results/tables/normalization_diagnostics.csv` (row stage==before_qnorm, column iqr_sd) · script: `scripts/00_shared/03_normalize_batch.R`
- **SD of per-sample medians, after_qnorm:** 0.000334
  <br>*[R-030.after_qnorm]* — source: `results/tables/normalization_diagnostics.csv` (row stage==after_qnorm, column median_sd) · script: `scripts/00_shared/03_normalize_batch.R`
- **SD of per-sample IQRs, after_qnorm:** 0.000457
  <br>*[R-031.after_qnorm]* — source: `results/tables/normalization_diagnostics.csv` (row stage==after_qnorm, column iqr_sd) · script: `scripts/00_shared/03_normalize_batch.R`
- **Number of ComBat batches resolved:** 6
  <br>*[R-035]* — source: `results/logs/05d_interaction.log` (literal 'batches 6' in the cohort summary line) · script: `scripts/00_shared/05d_interaction_report.R`

## §2.3 Differential gene expression analysis

- **Significant DEGs, All:** 6422
  <br>*[R-040.All]* — source: `results/tables/DEG_summary.csv` (row All, col significant) · script: `scripts/00_shared/05_dge.R`
- **Up in RA, All:** 2760
  <br>*[R-041.All]* — source: `results/tables/DEG_summary.csv` (row All, col up_in_RA) · script: `scripts/00_shared/05_dge.R`
- **Down in RA, All:** 3662
  <br>*[R-042.All]* — source: `results/tables/DEG_summary.csv` (row All, col down_in_RA) · script: `scripts/00_shared/05_dge.R`
- **n / RA / HC, All:** 183 / 103 / 80
  <br>*[R-043.All]* — source: `results/tables/DEG_summary.csv` (row All, cols RA+HC) · script: `scripts/00_shared/05_dge.R`
- **% of transcriptome significant, All:** 40.7%
  <br>*[R-044.All]* — source: `results/tables/DEG_summary.csv` (significant / genes_tested) · script: `scripts/00_shared/05_dge.R`
- **Significant DEGs, Female:** 5131
  <br>*[R-040.Female]* — source: `results/tables/DEG_summary.csv` (row Female, col significant) · script: `scripts/00_shared/05_dge.R`
- **Up in RA, Female:** 2238
  <br>*[R-041.Female]* — source: `results/tables/DEG_summary.csv` (row Female, col up_in_RA) · script: `scripts/00_shared/05_dge.R`
- **Down in RA, Female:** 2893
  <br>*[R-042.Female]* — source: `results/tables/DEG_summary.csv` (row Female, col down_in_RA) · script: `scripts/00_shared/05_dge.R`
- **n / RA / HC, Female:** 145 / 86 / 59
  <br>*[R-043.Female]* — source: `results/tables/DEG_summary.csv` (row Female, cols RA+HC) · script: `scripts/00_shared/05_dge.R`
- **% of transcriptome significant, Female:** 32.6%
  <br>*[R-044.Female]* — source: `results/tables/DEG_summary.csv` (significant / genes_tested) · script: `scripts/00_shared/05_dge.R`
- **Significant DEGs, Male:** 5820
  <br>*[R-040.Male]* — source: `results/tables/DEG_summary.csv` (row Male, col significant) · script: `scripts/00_shared/05_dge.R`
- **Up in RA, Male:** 2510
  <br>*[R-041.Male]* — source: `results/tables/DEG_summary.csv` (row Male, col up_in_RA) · script: `scripts/00_shared/05_dge.R`
- **Down in RA, Male:** 3310
  <br>*[R-042.Male]* — source: `results/tables/DEG_summary.csv` (row Male, col down_in_RA) · script: `scripts/00_shared/05_dge.R`
- **n / RA / HC, Male:** 38 / 17 / 21
  <br>*[R-043.Male]* — source: `results/tables/DEG_summary.csv` (row Male, cols RA+HC) · script: `scripts/00_shared/05_dge.R`
- **% of transcriptome significant, Male:** 36.9%
  <br>*[R-044.Male]* — source: `results/tables/DEG_summary.csv` (significant / genes_tested) · script: `scripts/00_shared/05_dge.R`
- **Array weight median [min-max], All:** 1.042 [0.380-2.055]
  <br>*[R-050.All]* — source: `results/tables/DEG_array_weights.csv` (median/min/max of weight where comparison==All) · script: `scripts/00_shared/05_dge.R`
- **Arrays with weight < 0.5, All:** 8
  <br>*[R-051.All]* — source: `results/tables/DEG_array_weights.csv` (count weight<0.5 where comparison==All) · script: `scripts/00_shared/05_dge.R`
- **Array weight median [min-max], Female:** 1.060 [0.380-2.051]
  <br>*[R-050.Female]* — source: `results/tables/DEG_array_weights.csv` (median/min/max of weight where comparison==Female) · script: `scripts/00_shared/05_dge.R`
- **Arrays with weight < 0.5, Female:** 5
  <br>*[R-051.Female]* — source: `results/tables/DEG_array_weights.csv` (count weight<0.5 where comparison==Female) · script: `scripts/00_shared/05_dge.R`
- **Array weight median [min-max], Male:** 1.101 [0.338-1.702]
  <br>*[R-050.Male]* — source: `results/tables/DEG_array_weights.csv` (median/min/max of weight where comparison==Male) · script: `scripts/00_shared/05_dge.R`
- **Arrays with weight < 0.5, Male:** 3
  <br>*[R-051.Male]* — source: `results/tables/DEG_array_weights.csv` (count weight<0.5 where comparison==Male) · script: `scripts/00_shared/05_dge.R`
- **Median array weight, GSE110169 (All contrast):** 1.006
  <br>*[R-052.GSE110169]* — source: `results/tables/DEG_array_weights.csv` (median weight where comparison==All & dataset==GSE110169) · script: `scripts/00_shared/05_dge.R`
- **Median array weight, GSE93272 (All contrast):** 1.087
  <br>*[R-052.GSE93272]* — source: `results/tables/DEG_array_weights.csv` (median weight where comparison==All & dataset==GSE93272) · script: `scripts/00_shared/05_dge.R`
- **Top DEGs by FDR, female:** HMGB2 (+0.67, FDR 2.7e-12); MAGED1 (-0.39, FDR 2.7e-12); KDM1A (-0.29, FDR 3.3e-12); C1GALT1C1 (+0.65, FDR 4.4e-12); CSGALNACT2 (+0.46, FDR 2.6e-11); S100A8 (+0.50, FDR 5.1e-11)
  <br>*[R-060.female]* — source: `results/tables/DEG_female_significant.csv` (6 smallest adj.P.Val) · script: `scripts/00_shared/05_dge.R`
- **Largest |logFC| among significant, female:** ARG1 (+1.36); CLEC4D (+1.17); BCL2A1 (+1.16); COX7B (+1.02); DEFA4 (+1.00); SCOC (+0.97)
  <br>*[R-061.female]* — source: `results/tables/DEG_female_significant.csv` (6 largest |logFC|) · script: `scripts/00_shared/05_dge.R`
- **Median |logFC| among significant, female:** 0.192
  <br>*[R-062.female]* — source: `results/tables/DEG_female_significant.csv` (median of |logFC|) · script: `scripts/00_shared/05_dge.R`
- **Minimum attainable FDR, female:** 2.7e-12
  <br>*[R-063.female]* — source: `results/tables/DEG_female_significant.csv` (min adj.P.Val) · script: `scripts/00_shared/05_dge.R`
- **Top DEGs by FDR, male:** STARD3NL (+0.70, FDR 3.3e-08); TBC1D15 (+1.01, FDR 3.3e-08); MIER1 (+0.65, FDR 3.3e-08); TTC33 (+1.30, FDR 3.3e-08); TRIM23 (+1.21, FDR 3.4e-08); SMARCA4 (-0.56, FDR 3.4e-08)
  <br>*[R-060.male]* — source: `results/tables/DEG_male_significant.csv` (6 smallest adj.P.Val) · script: `scripts/00_shared/05_dge.R`
- **Largest |logFC| among significant, male:** BCL2A1 (+1.84); COX7B (+1.70); RPL22L1 (+1.62); COMMD6 (+1.60); SCOC (+1.58); ZNF117 (+1.58)
  <br>*[R-061.male]* — source: `results/tables/DEG_male_significant.csv` (6 largest |logFC|) · script: `scripts/00_shared/05_dge.R`
- **Median |logFC| among significant, male:** 0.337
  <br>*[R-062.male]* — source: `results/tables/DEG_male_significant.csv` (median of |logFC|) · script: `scripts/00_shared/05_dge.R`
- **Minimum attainable FDR, male:** 3.3e-08
  <br>*[R-063.male]* — source: `results/tables/DEG_male_significant.csv` (min adj.P.Val) · script: `scripts/00_shared/05_dge.R`
- **DEGs significant in both sexes:** 3857
  <br>*[R-070]* — source: `results/tables/DEG_female_significant.csv + results/tables/DEG_male_significant.csv` (set intersection on gene) · script: `scripts/00_shared/05_dge.R`
- **Female-list-only DEGs:** 1274
  <br>*[R-071]* — source: `results/tables/DEG_female_significant.csv + results/tables/DEG_male_significant.csv` (set difference) · script: `scripts/00_shared/05_dge.R`
- **Male-list-only DEGs:** 1963
  <br>*[R-072]* — source: `results/tables/DEG_female_significant.csv + results/tables/DEG_male_significant.csv` (set difference) · script: `scripts/00_shared/05_dge.R`
- **Union of the two DEG lists:** 7094
  <br>*[R-073]* — source: `results/tables/DEG_female_significant.csv + results/tables/DEG_male_significant.csv` (set union) · script: `scripts/00_shared/05_dge.R`
- **Shared DEGs concordant in direction:** 3854/3857 (99.9%)
  <br>*[R-074]* — source: `results/tables/DEG_female_significant.csv + results/tables/DEG_male_significant.csv` (sign(logFC) agreement) · script: `scripts/00_shared/05_dge.R`
- **TREAT vs post-hoc filter, All at 1.07x:** 6422 -> 3331
  <br>*[R-080.All.1.07x]* — source: `results/tables/DEG_treat_sensitivity.csv` (cols posthoc_filter_n and treat_n) · script: `scripts/00_shared/05_dge.R`
- **TREAT vs post-hoc filter, All at 1.20x:** 1942 -> 464
  <br>*[R-080.All.1.20x]* — source: `results/tables/DEG_treat_sensitivity.csv` (cols posthoc_filter_n and treat_n) · script: `scripts/00_shared/05_dge.R`
- **TREAT vs post-hoc filter, All at 1.41x:** 321 -> 22
  <br>*[R-080.All.1.41x]* — source: `results/tables/DEG_treat_sensitivity.csv` (cols posthoc_filter_n and treat_n) · script: `scripts/00_shared/05_dge.R`
- **TREAT vs post-hoc filter, Female at 1.07x:** 5131 -> 1771
  <br>*[R-080.Female.1.07x]* — source: `results/tables/DEG_treat_sensitivity.csv` (cols posthoc_filter_n and treat_n) · script: `scripts/00_shared/05_dge.R`
- **TREAT vs post-hoc filter, Female at 1.20x:** 1300 -> 86
  <br>*[R-080.Female.1.20x]* — source: `results/tables/DEG_treat_sensitivity.csv` (cols posthoc_filter_n and treat_n) · script: `scripts/00_shared/05_dge.R`
- **TREAT vs post-hoc filter, Female at 1.41x:** 163 -> 3
  <br>*[R-080.Female.1.41x]* — source: `results/tables/DEG_treat_sensitivity.csv` (cols posthoc_filter_n and treat_n) · script: `scripts/00_shared/05_dge.R`
- **TREAT vs post-hoc filter, Male at 1.07x:** 5820 -> 3493
  <br>*[R-080.Male.1.07x]* — source: `results/tables/DEG_treat_sensitivity.csv` (cols posthoc_filter_n and treat_n) · script: `scripts/00_shared/05_dge.R`
- **TREAT vs post-hoc filter, Male at 1.20x:** 4126 -> 745
  <br>*[R-080.Male.1.20x]* — source: `results/tables/DEG_treat_sensitivity.csv` (cols posthoc_filter_n and treat_n) · script: `scripts/00_shared/05_dge.R`
- **TREAT vs post-hoc filter, Male at 1.41x:** 1225 -> 48
  <br>*[R-080.Male.1.41x]* — source: `results/tables/DEG_treat_sensitivity.csv` (cols posthoc_filter_n and treat_n) · script: `scripts/00_shared/05_dge.R`
- **ComBat-then-DE vs batch-in-model, All:** 6422 vs 6580 (ratio 1.025)
  <br>*[R-085.All]* — source: `results/tables/DEG_sensitivity_combat_vs_batch.csv` (cols combat_then_DE, batch_in_model, ratio) · script: `scripts/00_shared/05b_dge_sensitivity.R`
- **ComBat-then-DE vs batch-in-model, Female:** 5131 vs 4780 (ratio 0.932)
  <br>*[R-085.Female]* — source: `results/tables/DEG_sensitivity_combat_vs_batch.csv` (cols combat_then_DE, batch_in_model, ratio) · script: `scripts/00_shared/05b_dge_sensitivity.R`
- **ComBat-then-DE vs batch-in-model, Male:** 5820 vs 5612 (ratio 0.964)
  <br>*[R-085.Male]* — source: `results/tables/DEG_sensitivity_combat_vs_batch.csv` (cols combat_then_DE, batch_in_model, ratio) · script: `scripts/00_shared/05b_dge_sensitivity.R`

## §2.4 Co-expression network analysis (WGCNA)

_Methodology: see [`results/METHODS_2.4_WGCNA_expanded.md`](../results/METHODS_2.4_WGCNA_expanded.md) for the full specification of this step._

- **Signed R2 / mean k at power 3:** R2 0.873, mean k 2240.3
  <br>*[R-090.b3]* — source: `results/tables/WGCNA_01_soft_threshold.csv` (row power==3) · script: `scripts/00_shared/06_WGCNA.R`
- **Signed R2 / mean k at power 6:** R2 0.861, mean k 446.6
  <br>*[R-090.b6]* — source: `results/tables/WGCNA_01_soft_threshold.csv` (row power==6) · script: `scripts/00_shared/06_WGCNA.R`
- **Signed R2 / mean k at power 9:** R2 0.897, mean k 125.0
  <br>*[R-090.b9]* — source: `results/tables/WGCNA_01_soft_threshold.csv` (row power==9) · script: `scripts/00_shared/06_WGCNA.R`
- **Signed R2 / mean k at power 10:** R2 0.903, mean k 87.8
  <br>*[R-090.b10]* — source: `results/tables/WGCNA_01_soft_threshold.csv` (row power==10) · script: `scripts/00_shared/06_WGCNA.R`
- **Signed R2 / mean k at power 12:** R2 0.901, mean k 47.4
  <br>*[R-090.b12]* — source: `results/tables/WGCNA_01_soft_threshold.csv` (row power==12) · script: `scripts/00_shared/06_WGCNA.R`
- **Signed R2 / mean k at power 14:** R2 0.908, mean k 28.2
  <br>*[R-090.b14]* — source: `results/tables/WGCNA_01_soft_threshold.csv` (row power==14) · script: `scripts/00_shared/06_WGCNA.R`
- **Signed R2 / mean k at power 20:** R2 0.914, mean k 8.5
  <br>*[R-090.b20]* — source: `results/tables/WGCNA_01_soft_threshold.csv` (row power==20) · script: `scripts/00_shared/06_WGCNA.R`
- **Modules detected (excluding grey):** 12
  <br>*[R-100]* — source: `results/tables/WGCNA_02_module_trait.csv` (row count, module!=grey) · script: `scripts/00_shared/06_WGCNA.R`
- **Genes unassigned (grey):** 6019
  <br>*[R-101]* — source: `results/tables/WGCNA_02_module_trait.csv` (row module==grey, col size) · script: `scripts/00_shared/06_WGCNA.R`
- **Module sizes:** turquoise 2584; blue 2554; brown 1878; yellow 708; green 639; red 525; black 295; pink 172; magenta 127; purple 120; greenyellow 108; tan 34
  <br>*[R-102]* — source: `results/tables/WGCNA_02_module_trait.csv` (cols module,size sorted desc) · script: `scripts/00_shared/06_WGCNA.R`
- **cor(ME,RA) and p, yellow:** +0.556 (p 1.91e-15)
  <br>*[R-103.yellow]* — source: `results/tables/WGCNA_02_module_trait.csv` (row yellow, cols cor_RA,p_RA) · script: `scripts/00_shared/06_WGCNA.R`
- **cor(ME,RA) and p, turquoise:** +0.369 (p 6.09e-07)
  <br>*[R-103.turquoise]* — source: `results/tables/WGCNA_02_module_trait.csv` (row turquoise, cols cor_RA,p_RA) · script: `scripts/00_shared/06_WGCNA.R`
- **cor(ME,RA) and p, green:** -0.362 (p 9.59e-07)
  <br>*[R-103.green]* — source: `results/tables/WGCNA_02_module_trait.csv` (row green, cols cor_RA,p_RA) · script: `scripts/00_shared/06_WGCNA.R`
- **cor(ME,RA) and p, brown:** -0.590 (p 1.29e-17)
  <br>*[R-103.brown]* — source: `results/tables/WGCNA_02_module_trait.csv` (row brown, cols cor_RA,p_RA) · script: `scripts/00_shared/06_WGCNA.R`
- **Disease modules selected:** yellow (n=708, r=+0.556, p=1.91e-15, UP in RA); brown (n=1878, r=-0.590, p=1.29e-17, DOWN in RA)
  <br>*[R-110]* — source: `results/tables/WGCNA_03_disease_modules.csv` (all rows) · script: `scripts/00_shared/06_WGCNA.R`
- **Disease-module gene background:** 2586
  <br>*[R-111]* — source: `results/tables/WGCNA_03_disease_modules.csv` (sum of size) · script: `scripts/00_shared/06_WGCNA.R`
- **cor(ME,Male) for disease modules:** yellow -0.176 (p 0.0206); brown +0.190 (p 0.0122)
  <br>*[R-112]* — source: `results/tables/WGCNA_03_disease_modules.csv` (cols cor_Male,p_Male) · script: `scripts/00_shared/06_WGCNA.R`
- **cor(ME,RA) yellow in All samples (n=173):** +0.556 (p 1.91e-15)
  <br>*[R-115.yellow.All_samples]* — source: `results/tables/module_trait_RAvsControl_ALLSTRATA.csv` (row module==yellow & stratum==All samples) · script: `scripts/00_shared/08_module_trait_RA_control.R`
- **cor(ME,RA) brown in All samples (n=173):** -0.590 (p 1.29e-17)
  <br>*[R-115.brown.All_samples]* — source: `results/tables/module_trait_RAvsControl_ALLSTRATA.csv` (row module==brown & stratum==All samples) · script: `scripts/00_shared/08_module_trait_RA_control.R`
- **cor(ME,RA) yellow in Female (n=138):** +0.485 (p 1.69e-09)
  <br>*[R-115.yellow.Female]* — source: `results/tables/module_trait_RAvsControl_ALLSTRATA.csv` (row module==yellow & stratum==Female) · script: `scripts/00_shared/08_module_trait_RA_control.R`
- **cor(ME,RA) brown in Female (n=138):** -0.515 (p 1.00e-10)
  <br>*[R-115.brown.Female]* — source: `results/tables/module_trait_RAvsControl_ALLSTRATA.csv` (row module==brown & stratum==Female) · script: `scripts/00_shared/08_module_trait_RA_control.R`
- **cor(ME,RA) yellow in Male (n=35):** +0.772 (p 5.51e-08)
  <br>*[R-115.yellow.Male]* — source: `results/tables/module_trait_RAvsControl_ALLSTRATA.csv` (row module==yellow & stratum==Male) · script: `scripts/00_shared/08_module_trait_RA_control.R`
- **cor(ME,RA) brown in Male (n=35):** -0.813 (p 2.88e-09)
  <br>*[R-115.brown.Male]* — source: `results/tables/module_trait_RAvsControl_ALLSTRATA.csv` (row module==brown & stratum==Male) · script: `scripts/00_shared/08_module_trait_RA_control.R`
- **Hub genes in yellow module:** 217
  <br>*[R-120.yellow]* — source: `results/tables/WGCNA_06_disease_module_hubs.csv` (count is_hub==TRUE & module==yellow) · script: `scripts/00_shared/06_WGCNA.R`
- **Top-connectivity hubs, yellow:** ZNF267 (kME 0.946, GS +0.534, k 117.9); ZFYVE16 (kME 0.953, GS +0.490, k 114.8); SP3 (kME 0.946, GS +0.520, k 113.6); FAR1 (kME 0.951, GS +0.519, k 111.9); SNX13 (kME 0.937, GS +0.526, k 111.8); SPOPL (kME 0.955, GS +0.505, k 110.7)
  <br>*[R-121.yellow]* — source: `results/tables/WGCNA_06_disease_module_hubs.csv` (6 largest connectivity among hubs) · script: `scripts/00_shared/06_WGCNA.R`
- **Hub genes in brown module:** 144
  <br>*[R-120.brown]* — source: `results/tables/WGCNA_06_disease_module_hubs.csv` (count is_hub==TRUE & module==brown) · script: `scripts/00_shared/06_WGCNA.R`
- **Top-connectivity hubs, brown:** KHSRP (kME 0.911, GS -0.564, k 130.1); GPI (kME 0.907, GS -0.510, k 124.5); SF3B3 (kME 0.894, GS -0.567, k 122.0); CLSTN1 (kME 0.896, GS -0.556, k 119.1); XRCC6 (kME 0.879, GS -0.489, k 118.4); EXOSC10 (kME 0.902, GS -0.549, k 117.9)
  <br>*[R-121.brown]* — source: `results/tables/WGCNA_06_disease_module_hubs.csv` (6 largest connectivity among hubs) · script: `scripts/00_shared/06_WGCNA.R`
- **Preservation Zsummary, gold (random benchmark):** 34.19
  <br>*[R-130]* — source: `results/tables/WGCNA_10_module_preservation.csv` (row module==gold) · script: `scripts/00_shared/06_WGCNA.R`
- **Preservation Zsummary, yellow:** 48.14
  <br>*[R-131.yellow]* — source: `results/tables/WGCNA_10_module_preservation.csv` (row module==yellow) · script: `scripts/00_shared/06_WGCNA.R`
- **yellow Zsummary minus gold benchmark:** +13.95
  <br>*[R-132.yellow]* — source: `results/tables/WGCNA_10_module_preservation.csv` (Zsummary(module) - Zsummary(gold)) · script: `scripts/00_shared/06_WGCNA.R`
- **Preservation Zsummary, brown:** 36.95
  <br>*[R-131.brown]* — source: `results/tables/WGCNA_10_module_preservation.csv` (row module==brown) · script: `scripts/00_shared/06_WGCNA.R`
- **brown Zsummary minus gold benchmark:** +2.76
  <br>*[R-132.brown]* — source: `results/tables/WGCNA_10_module_preservation.csv` (Zsummary(module) - Zsummary(gold)) · script: `scripts/00_shared/06_WGCNA.R`
- **Sample-outlier cut height:** 69.32
  <br>*[R-140]* — source: `results/logs/06_WGCNA_run3.log` (literal in QC block) · script: `scripts/00_shared/06_WGCNA.R`
- **Sample outliers removed:** 10
  <br>*[R-141]* — source: `results/logs/06_WGCNA_run3.log` (literal in QC block) · script: `scripts/00_shared/06_WGCNA.R`
- **Network matrix after outlier removal:** 173 samples x 15763 genes
  <br>*[R-142]* — source: `results/logs/06_WGCNA_run3.log` (literal 'final matrix' line) · script: `scripts/00_shared/06_WGCNA.R`

## §2.5 Candidate gene identification (disease module ∩ sex-stratified DEG)

- **Candidates (Female) = disease modules INTERSECT Female DEGs:** 2045
  <br>*[R-150.Female]* — source: `results/tables/candidate_summary.csv` (col n_candidates) · script: `scripts/00_shared/09_disease_module_deg_intersect.R`
- **Candidate split by module, Female:** yellow=595 brown=1450
  <br>*[R-151.Female]* — source: `results/tables/candidate_summary.csv` (col per_module) · script: `scripts/00_shared/09_disease_module_deg_intersect.R`
- **DEGs lost to variance filter, Female:** 0 (0%)
  <br>*[R-152.Female]* — source: `results/tables/candidate_summary.csv` (cols DEGs_lost_to_variance_filter, pct_DEG_lost) · script: `scripts/00_shared/09_disease_module_deg_intersect.R`
- **Candidates (Male) = disease modules INTERSECT Male DEGs:** 2079
  <br>*[R-150.Male]* — source: `results/tables/candidate_summary.csv` (col n_candidates) · script: `scripts/00_shared/09_disease_module_deg_intersect.R`
- **Candidate split by module, Male:** yellow=585 brown=1494
  <br>*[R-151.Male]* — source: `results/tables/candidate_summary.csv` (col per_module) · script: `scripts/00_shared/09_disease_module_deg_intersect.R`
- **DEGs lost to variance filter, Male:** 0 (0%)
  <br>*[R-152.Male]* — source: `results/tables/candidate_summary.csv` (cols DEGs_lost_to_variance_filter, pct_DEG_lost) · script: `scripts/00_shared/09_disease_module_deg_intersect.R`
- **Directional consistency (consistent / inconsistent overlaps):** Female-yellow: 595/0; Female-brown: 1450/0; Male-yellow: 585/0; Male-brown: 1494/0
  <br>*[R-155]* — source: `results/tables/diseasemod_DEG_direction_summary.csv` (cols n_consistent, n_inconsistent) · script: `scripts/00_shared/09b_disease_module_deg_venn.R`
- **Candidates shared by both sexes:** 1773
  <br>*[R-156]* — source: `results/tables/WGCNA_11_candidates_{female,male}.csv` (set intersection) · script: `scripts/00_shared/06_WGCNA.R`
- **Female-list-only candidates:** 272
  <br>*[R-157]* — source: `results/tables/WGCNA_11_candidates_{female,male}.csv` (set difference) · script: `scripts/00_shared/06_WGCNA.R`
- **Male-list-only candidates:** 306
  <br>*[R-158]* — source: `results/tables/WGCNA_11_candidates_{female,male}.csv` (set difference) · script: `scripts/00_shared/06_WGCNA.R`
- **Union of candidates carried into MR:** 2351
  <br>*[R-159]* — source: `results/tables/WGCNA_11_candidates_{female,male}.csv` (set union) · script: `scripts/00_shared/06_WGCNA.R`

## §2.6 Mendelian randomisation, incl. MHC sensitivity analysis

_Methodology: see [`results/METHODS_2.6_mendelian_randomisation.md`](../results/METHODS_2.6_mendelian_randomisation.md) for the full specification of this step._

- **Genes with an MR estimate, Female:** 1477
  <br>*[R-160.Female]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (col genes_tested_primary) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **MR-prioritised at within-stratum FDR<0.05, Female:** 32
  <br>*[R-161.Female]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (col causal_FDR05_primary) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **Surviving MHC exclusion, Female:** 14
  <br>*[R-162.Female]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (col causal_FDR05_noMHC) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **MHC verdicts (robust/untestable/MHC-dep/FDR-rank), Female:** 14/14/0/4
  <br>*[R-163.Female]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (cols causal_robust, causal_untestable_noMHC, causal_MHC_dependent, causal_FDR_rank_only) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **Genes tested after MHC exclusion, Female:** 1448
  <br>*[R-164.Female]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (col genes_tested_noMHC) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **Panel-gene MHC verdicts (robust/untestable/dep/rank), Female:** 3/2/0/1 of 6
  <br>*[R-165.Female]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (panel_* columns) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **Genes with an MR estimate, Male:** 1478
  <br>*[R-160.Male]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (col genes_tested_primary) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **MR-prioritised at within-stratum FDR<0.05, Male:** 25
  <br>*[R-161.Male]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (col causal_FDR05_primary) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **Surviving MHC exclusion, Male:** 14
  <br>*[R-162.Male]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (col causal_FDR05_noMHC) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **MHC verdicts (robust/untestable/MHC-dep/FDR-rank), Male:** 14/10/0/1
  <br>*[R-163.Male]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (cols causal_robust, causal_untestable_noMHC, causal_MHC_dependent, causal_FDR_rank_only) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **Genes tested after MHC exclusion, Male:** 1456
  <br>*[R-164.Male]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (col genes_tested_noMHC) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **Panel-gene MHC verdicts (robust/untestable/dep/rank), Male:** 3/2/0/1 of 6
  <br>*[R-165.Male]* — source: `results/tables/MR_MHC_sensitivity_summary.csv` (panel_* columns) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **FS input gene count, female:** 32
  <br>*[R-170.female]* — source: `results/tables/FS_input_female.csv` (row count) · script: `scripts/00_shared/10_MR.R`
- **Strongest MR associations, female:** HLA-DRB1 OR 0.337, FDR 1.5e-247, 1 SNP, Wald ratio; WDR46 OR 4.531, FDR 3.8e-39, 1 SNP, Wald ratio; AIF1 OR 0.556, FDR 1.2e-21, 2 SNP, Inverse variance weighted; VPS52 OR 0.414, FDR 1.5e-21, 1 SNP, Wald ratio; GNL1 OR 4.629, FDR 2.7e-15, 1 SNP, Wald ratio; VARS2 OR 1.203, FDR 2.7e-14, 3 SNP, Inverse variance weighted
  <br>*[R-171.female]* — source: `results/tables/FS_input_female.csv` (6 smallest MR_pval) · script: `scripts/00_shared/10_MR.R`
- **Single-instrument genes among FS input, female:** 22/32
  <br>*[R-172.female]* — source: `results/tables/FS_input_female.csv` (count nSNP==1) · script: `scripts/00_shared/10_MR.R`
- **MHC-free FS input gene count, female:** 14
  <br>*[R-173.female]* — source: `results/tables/FS_input_female_noMHC.csv` (row count) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **FS input gene count, male:** 25
  <br>*[R-170.male]* — source: `results/tables/FS_input_male.csv` (row count) · script: `scripts/00_shared/10_MR.R`
- **Strongest MR associations, male:** HLA-DRB1 OR 0.337, FDR 1.5e-247, 1 SNP, Wald ratio; WDR46 OR 4.531, FDR 3.8e-39, 1 SNP, Wald ratio; AIF1 OR 0.556, FDR 1.2e-21, 2 SNP, Inverse variance weighted; VPS52 OR 0.414, FDR 1.5e-21, 1 SNP, Wald ratio; VARS2 OR 1.203, FDR 3.2e-14, 3 SNP, Inverse variance weighted; AP4B1 OR 1.424, FDR 3.2e-10, 1 SNP, Wald ratio
  <br>*[R-171.male]* — source: `results/tables/FS_input_male.csv` (6 smallest MR_pval) · script: `scripts/00_shared/10_MR.R`
- **Single-instrument genes among FS input, male:** 19/25
  <br>*[R-172.male]* — source: `results/tables/FS_input_male.csv` (count nSNP==1) · script: `scripts/00_shared/10_MR.R`
- **MHC-free FS input gene count, male:** 14
  <br>*[R-173.male]* — source: `results/tables/FS_input_male_noMHC.csv` (row count) · script: `scripts/00_shared/10c_MR_mhc_sensitivity.R`
- **MR-prioritised genes shared by both strata:** 24
  <br>*[R-175]* — source: `results/tables/FS_input_{female,male}.csv` (set intersection) · script: `scripts/00_shared/10_MR.R`
- **Female-list-only MR-prioritised genes:** 8 (FCGR2B, GNL1, HLA-DPA1, IKZF4, LAX1, LSM2, RETSAT, ZBTB9)
  <br>*[R-176]* — source: `results/tables/FS_input_{female,male}.csv` (set difference) · script: `scripts/00_shared/10_MR.R`
- **Male-list-only MR-prioritised genes:** 1 (TAB1)
  <br>*[R-177]* — source: `results/tables/FS_input_{female,male}.csv` (set difference) · script: `scripts/00_shared/10_MR.R`

## §2.7 Bayesian colocalisation (coloc.abf, coloc.susie)

_Methodology: see [`results/METHODS_2.7_colocalisation.md`](../results/METHODS_2.7_colocalisation.md) for the full specification of this step._

- **Colocalisation tally: All causal genes:** n=33, colocalised=0, prior-fragile=2, suggestive=4, distinct=9, MHC-unreliable=13, inconclusive=5
  <br>*[R-180.All_causal_genes]* — source: `results/tables/COLOC_summary.csv` (all columns of that row) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Colocalisation tally: MHC genes:** n=14, colocalised=0, prior-fragile=0, suggestive=0, distinct=0, MHC-unreliable=13, inconclusive=1
  <br>*[R-180.MHC_genes]* — source: `results/tables/COLOC_summary.csv` (all columns of that row) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Colocalisation tally: non-MHC genes:** n=19, colocalised=0, prior-fragile=2, suggestive=4, distinct=9, MHC-unreliable=0, inconclusive=4
  <br>*[R-180.non-MHC_genes]* — source: `results/tables/COLOC_summary.csv` (all columns of that row) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Colocalisation tally: Final panels (F union M):** n=9, colocalised=0, prior-fragile=0, suggestive=1, distinct=2, MHC-unreliable=4, inconclusive=2
  <br>*[R-180.Final_panels_(F_union_M)]* — source: `results/tables/COLOC_summary.csv` (all columns of that row) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Coloc posteriors, IKZF3:** PP.H3 0.226, PP.H4 0.774, PP.H4@p12=1e-6 0.255, 811 SNPs, MHC=FALSE
  <br>*[R-181.IKZF3]* — source: `results/tables/COLOC_panel_genes.csv` (cols PP3, PP4, PP4_cons, nsnp_coloc, MHC_gene) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Coloc posteriors, MED1:** PP.H3 0.552, PP.H4 0.431, PP.H4@p12=1e-6 0.070, 722 SNPs, MHC=FALSE
  <br>*[R-181.MED1]* — source: `results/tables/COLOC_panel_genes.csv` (cols PP3, PP4, PP4_cons, nsnp_coloc, MHC_gene) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Coloc posteriors, SMARCC2:** PP.H3 0.709, PP.H4 0.275, PP.H4@p12=1e-6 0.037, 365 SNPs, MHC=FALSE
  <br>*[R-181.SMARCC2]* — source: `results/tables/COLOC_panel_genes.csv` (cols PP3, PP4, PP4_cons, nsnp_coloc, MHC_gene) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Coloc posteriors, INPP5B:** PP.H3 0.929, PP.H4 0.069, PP.H4@p12=1e-6 0.007, 888 SNPs, MHC=FALSE
  <br>*[R-181.INPP5B]* — source: `results/tables/COLOC_panel_genes.csv` (cols PP3, PP4, PP4_cons, nsnp_coloc, MHC_gene) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Coloc posteriors, ESYT1:** PP.H3 0.912, PP.H4 0.068, PP.H4@p12=1e-6 0.007, 317 SNPs, MHC=FALSE
  <br>*[R-181.ESYT1]* — source: `results/tables/COLOC_panel_genes.csv` (cols PP3, PP4, PP4_cons, nsnp_coloc, MHC_gene) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Coloc posteriors, GNL1:** PP.H3 1.000, PP.H4 0.000, PP.H4@p12=1e-6 0.000, 1885 SNPs, MHC=TRUE
  <br>*[R-181.GNL1]* — source: `results/tables/COLOC_panel_genes.csv` (cols PP3, PP4, PP4_cons, nsnp_coloc, MHC_gene) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Coloc posteriors, C6orf136:** PP.H3 1.000, PP.H4 0.000, PP.H4@p12=1e-6 0.000, 1911 SNPs, MHC=TRUE
  <br>*[R-181.C6orf136]* — source: `results/tables/COLOC_panel_genes.csv` (cols PP3, PP4, PP4_cons, nsnp_coloc, MHC_gene) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Coloc posteriors, VPS52:** PP.H3 1.000, PP.H4 0.000, PP.H4@p12=1e-6 0.000, 2182 SNPs, MHC=TRUE
  <br>*[R-181.VPS52]* — source: `results/tables/COLOC_panel_genes.csv` (cols PP3, PP4, PP4_cons, nsnp_coloc, MHC_gene) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **Coloc posteriors, HLA-DMA:** PP.H3 1.000, PP.H4 0.000, PP.H4@p12=1e-6 0.000, 4193 SNPs, MHC=TRUE
  <br>*[R-181.HLA-DMA]* — source: `results/tables/COLOC_panel_genes.csv` (cols PP3, PP4, PP4_cons, nsnp_coloc, MHC_gene) · script: `scripts/00_shared/10d_coloc_panel_genes.R`
- **coloc.susie: MHC genes attempted / resolved:** 14 attempted, 3 resolved
  <br>*[R-190]* — source: `results/tables/COLOC_SUSIE_mhc.csv` (count rows; status=='ok') · script: `scripts/00_shared/10e_coloc_susie_mhc.R`
- **SuSiE result, WDR46:** 1 eQTL sets, 9 RA sets, best PP.H4 0.918 — COLOCALISED with an RA signal (multi-variant model)
  <br>*[R-191.WDR46]* — source: `results/tables/COLOC_SUSIE_mhc.csv` (cols n_cs_eqtl, n_cs_gwas, best_PP_H4, susie_verdict) · script: `scripts/00_shared/10e_coloc_susie_mhc.R`
- **SuSiE result, HLA-DMA:** 3 eQTL sets, 9 RA sets, best PP.H4 0.323 — INCONCLUSIVE under the multi-variant model
  <br>*[R-191.HLA-DMA]* — source: `results/tables/COLOC_SUSIE_mhc.csv` (cols n_cs_eqtl, n_cs_gwas, best_PP_H4, susie_verdict) · script: `scripts/00_shared/10e_coloc_susie_mhc.R`
- **SuSiE result, GNL1:** 1 eQTL sets, 10 RA sets, best PP.H4 0.000 — DISTINCT VARIANTS (valid: multiple signals modelled)
  <br>*[R-191.GNL1]* — source: `results/tables/COLOC_SUSIE_mhc.csv` (cols n_cs_eqtl, n_cs_gwas, best_PP_H4, susie_verdict) · script: `scripts/00_shared/10e_coloc_susie_mhc.R`
- **Independent RA credible sets found in MHC regions:** 9-10
  <br>*[R-192]* — source: `results/tables/COLOC_SUSIE_mhc.csv` (range of n_cs_gwas (non-zero)) · script: `scripts/00_shared/10e_coloc_susie_mhc.R`

## §2.8 Sex-stratified feature selection (LASSO / RF / SVM-RFE consensus)

_Methodology: see [`results/METHODS_2.8_feature_selection.md`](../results/METHODS_2.8_feature_selection.md) for the full specification of this step._

- **Selector yields, Female (primary):** LASSO 7, RF 10, SVM-RFE 8, consensus 6
  <br>*[R-200.Female]* — source: `results/tables/mr_fs_summary.csv` (cols n_lasso, n_rf, n_svmrfe, n_consensus) · script: `scripts/goal2_sex_stratified/12_feature_selection.R`
- **Primary consensus panel, Female:** C6orf136; ESYT1; GNL1; IKZF3; MED1; SMARCC2
  <br>*[R-201.Female]* — source: `results/tables/mr_fs_summary.csv` (col consensus_genes) · script: `scripts/goal2_sex_stratified/12_feature_selection.R`
- **Tuned hyperparameters, Female:** lambda.min 0.0434, mtry 2, SVM cost 0.01
  <br>*[R-202.Female]* — source: `results/tables/mr_fs_summary.csv` (tuned_* columns) · script: `scripts/goal2_sex_stratified/12_feature_selection.R`
- **Selector yields, Male (primary):** LASSO 10, RF 10, SVM-RFE 11, consensus 6
  <br>*[R-200.Male]* — source: `results/tables/mr_fs_summary.csv` (cols n_lasso, n_rf, n_svmrfe, n_consensus) · script: `scripts/goal2_sex_stratified/12_feature_selection.R`
- **Primary consensus panel, Male:** ESYT1; HLA-DMA; INPP5B; MED1; SMARCC2; VPS52
  <br>*[R-201.Male]* — source: `results/tables/mr_fs_summary.csv` (col consensus_genes) · script: `scripts/goal2_sex_stratified/12_feature_selection.R`
- **Tuned hyperparameters, Male:** lambda.min 0.0204, mtry 1, SVM cost 0.01
  <br>*[R-202.Male]* — source: `results/tables/mr_fs_summary.csv` (tuned_* columns) · script: `scripts/goal2_sex_stratified/12_feature_selection.R`
- **MHC-free consensus panel, Female:** 4 genes: CDC37; IKZF3; MED1; SMARCC2
  <br>*[R-205.Female]* — source: `results/tables/mr_fs_summary_noMHC.csv` (cols n_consensus, consensus_genes) · script: `scripts/goal2_sex_stratified/12b_feature_selection_noMHC.R`
- **MHC-free consensus panel, Male:** 5 genes: MED1; NCOA5; PHF19; SMARCC2; TAB1
  <br>*[R-205.Male]* — source: `results/tables/mr_fs_summary_noMHC.csv` (cols n_consensus, consensus_genes) · script: `scripts/goal2_sex_stratified/12b_feature_selection_noMHC.R`
- **Panel membership change: DROPPED:** C6orf136(F); ESYT1(F); GNL1(F); ESYT1(M); HLA-DMA(M); INPP5B(M); VPS52(M)
  <br>*[R-206.DROPPED]* — source: `results/tables/PANEL_primary_vs_noMHC_membership.csv` (rows with status starting DROPPED) · script: `scripts/goal2_sex_stratified/12b_feature_selection_noMHC.R`
- **Panel membership change: NEW:** CDC37(F); NCOA5(M); PHF19(M); TAB1(M)
  <br>*[R-206.NEW]* — source: `results/tables/PANEL_primary_vs_noMHC_membership.csv` (rows with status starting NEW) · script: `scripts/goal2_sex_stratified/12b_feature_selection_noMHC.R`
- **Panel membership change: retained:** IKZF3(F); MED1(F); SMARCC2(F); MED1(M); SMARCC2(M)
  <br>*[R-206.retained]* — source: `results/tables/PANEL_primary_vs_noMHC_membership.csv` (rows with status starting retained) · script: `scripts/goal2_sex_stratified/12b_feature_selection_noMHC.R`

## §2.9 Diagnostic model development and evaluation

_Methodology: see [`results/METHODS_2.9_diagnostic_model.md`](../results/METHODS_2.9_diagnostic_model.md) for the full specification of this step._

- **Nested CV AUC — Female, noMHC, consensus:** 0.798 (0.725-0.871), per-repeat SD 0.006, median genes 3, n=145, recommended=TRUE
  <br>*[R-210.Female.noMHC.consensus]* — source: `results/tables/NESTED_CV_AUTHORITATIVE.csv` (one row of the 8-row grid) · script: `scripts/goal2_sex_stratified/16d_nested_cv_reconciliation.R`
- **Nested CV AUC — Female, noMHC, elasticnet:** 0.791 (0.716-0.866), per-repeat SD 0.003, median genes 7, n=145, recommended=FALSE
  <br>*[R-210.Female.noMHC.elasticnet]* — source: `results/tables/NESTED_CV_AUTHORITATIVE.csv` (one row of the 8-row grid) · script: `scripts/goal2_sex_stratified/16d_nested_cv_reconciliation.R`
- **Nested CV AUC — Female, primary, consensus:** 0.816 (0.747-0.884), per-repeat SD 0.01, median genes 5, n=145, recommended=FALSE
  <br>*[R-210.Female.primary.consensus]* — source: `results/tables/NESTED_CV_AUTHORITATIVE.csv` (one row of the 8-row grid) · script: `scripts/goal2_sex_stratified/16d_nested_cv_reconciliation.R`
- **Nested CV AUC — Female, primary, elasticnet:** 0.801 (0.73-0.873), per-repeat SD 0.009, median genes 12, n=145, recommended=FALSE
  <br>*[R-210.Female.primary.elasticnet]* — source: `results/tables/NESTED_CV_AUTHORITATIVE.csv` (one row of the 8-row grid) · script: `scripts/goal2_sex_stratified/16d_nested_cv_reconciliation.R`
- **Nested CV AUC — Male, noMHC, consensus:** 0.924 (0.841-1), per-repeat SD 0.056, median genes 3, n=38, recommended=TRUE
  <br>*[R-210.Male.noMHC.consensus]* — source: `results/tables/NESTED_CV_AUTHORITATIVE.csv` (one row of the 8-row grid) · script: `scripts/goal2_sex_stratified/16d_nested_cv_reconciliation.R`
- **Nested CV AUC — Male, noMHC, elasticnet:** 0.961 (0.909-1), per-repeat SD 0.028, median genes 13, n=38, recommended=FALSE
  <br>*[R-210.Male.noMHC.elasticnet]* — source: `results/tables/NESTED_CV_AUTHORITATIVE.csv` (one row of the 8-row grid) · script: `scripts/goal2_sex_stratified/16d_nested_cv_reconciliation.R`
- **Nested CV AUC — Male, primary, consensus:** 0.896 (0.8-0.993), per-repeat SD 0.06, median genes 3, n=38, recommended=FALSE
  <br>*[R-210.Male.primary.consensus]* — source: `results/tables/NESTED_CV_AUTHORITATIVE.csv` (one row of the 8-row grid) · script: `scripts/goal2_sex_stratified/16d_nested_cv_reconciliation.R`
- **Nested CV AUC — Male, primary, elasticnet:** 0.958 (0.902-1), per-repeat SD 0.015, median genes 20, n=38, recommended=FALSE
  <br>*[R-210.Male.primary.elasticnet]* — source: `results/tables/NESTED_CV_AUTHORITATIVE.csv` (one row of the 8-row grid) · script: `scripts/goal2_sex_stratified/16d_nested_cv_reconciliation.R`
- **LOCKED-TRANSFER AUC — Female, primary panel, Train (apparent):** 0.869 (0.813-0.925) [n=145]
  <br>*[R-220.Female.primary.Train_(apparent)]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **LOCKED-TRANSFER AUC — Female, primary panel, Internal test:** 0.823 (0.719-0.928) [n=61]
  <br>*[R-220.Female.primary.Internal_test]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **LOCKED-TRANSFER AUC — Female, primary panel, External blood:** 0.957 (0.886-1.000) [n=24]
  <br>*[R-220.Female.primary.External_blood]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **LOCKED-TRANSFER AUC — Male, primary panel, Train (apparent):** 1.000 (1.000-1.000) [n=38] SEPARATION
  <br>*[R-220.Male.primary.Train_(apparent)]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **LOCKED-TRANSFER AUC — Male, primary panel, Internal test:** 1.000 (1.000-1.000) [n=13] SEPARATION
  <br>*[R-220.Male.primary.Internal_test]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **LOCKED-TRANSFER AUC — Male, primary panel, External blood:** 0.750 (0.399-1.000) [n=9]
  <br>*[R-220.Male.primary.External_blood]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **LOCKED-TRANSFER AUC — Female, noMHC panel, Train (apparent):** 0.845 (0.782-0.908) [n=145]
  <br>*[R-220.Female.noMHC.Train_(apparent)]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **LOCKED-TRANSFER AUC — Female, noMHC panel, Internal test:** 0.806 (0.692-0.919) [n=61]
  <br>*[R-220.Female.noMHC.Internal_test]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **LOCKED-TRANSFER AUC — Female, noMHC panel, External blood:** 0.743 (0.535-0.951) [n=24]
  <br>*[R-220.Female.noMHC.External_blood]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **LOCKED-TRANSFER AUC — Male, noMHC panel, Train (apparent):** 1.000 (1.000-1.000) [n=38] SEPARATION
  <br>*[R-220.Male.noMHC.Train_(apparent)]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **LOCKED-TRANSFER AUC — Male, noMHC panel, Internal test:** 0.917 (0.750-1.000) [n=13]
  <br>*[R-220.Male.noMHC.Internal_test]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **LOCKED-TRANSFER AUC — Male, noMHC panel, External blood:** 0.700 (0.400-1.000) [n=9]
  <br>*[R-220.Male.noMHC.External_blood]* — source: `results/tables/PANEL_primary_vs_noMHC_performance.csv` (col reported) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **DeLong, primary vs MHC-free nested CV, Female:** 0.816 vs 0.798, p=0.322
  <br>*[R-225.Female]* — source: `results/tables/PANEL_primary_vs_noMHC_delong.csv` (cols AUC_primary, AUC_noMHC, delong_p) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **DeLong, primary vs MHC-free nested CV, Male:** 0.896 vs 0.924, p=0.582
  <br>*[R-225.Male]* — source: `results/tables/PANEL_primary_vs_noMHC_delong.csv` (cols AUC_primary, AUC_noMHC, delong_p) · script: `scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R`
- **Apparent / flat / nested and optimism, Female:** apparent 0.869, flat 0.801, nested 0.816, optimism 0.054
  <br>*[R-230.Female]* — source: `results/tables/mr_nested_cv_summary.csv` (cols apparent_AUC, flat_CV_AUC, nested_CV_AUC, optimism_*) · script: `scripts/goal2_sex_stratified/14_model_training_nested_cv.R`
- **Apparent / flat / nested and optimism, Male:** apparent 1, flat 0.821, nested 0.896, optimism 0.104
  <br>*[R-230.Male]* — source: `results/tables/mr_nested_cv_summary.csv` (cols apparent_AUC, flat_CV_AUC, nested_CV_AUC, optimism_*) · script: `scripts/goal2_sex_stratified/14_model_training_nested_cv.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Female, primary, Train:** panel 0.831 vs composition 0.662 (delta 0.169), LRT p=1.21e-10, n=145, separation=FALSE
  <br>*[R-240.Female.primary.Train]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Female, primary, Internal test:** panel 0.721 vs composition 0.714 (delta 0.007), LRT p=0.00142, n=61, separation=FALSE
  <br>*[R-240.Female.primary.Internal_test]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Female, primary, External blood:** panel 1 vs composition 0.429 (delta 0.571), LRT p=2.42e-05, n=24, separation=TRUE
  <br>*[R-240.Female.primary.External_blood]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Male, primary, Train:** panel 0.933 vs composition 0.857 (delta 0.076), LRT p=5.21e-05, n=38, separation=TRUE
  <br>*[R-240.Male.primary.Train]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Male, primary, Internal test:** panel 1 vs composition 0.798 (delta 0.202), LRT p=1, n=13, separation=TRUE
  <br>*[R-240.Male.primary.Internal_test]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Male, primary, External blood:** panel 0.725 vs composition 0.575 (delta 0.15), LRT p=0.118, n=9, separation=TRUE
  <br>*[R-240.Male.primary.External_blood]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Female, noMHC, Train:** panel 0.812 vs composition 0.662 (delta 0.15), LRT p=4.08e-09, n=145, separation=FALSE
  <br>*[R-240.Female.noMHC.Train]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Female, noMHC, Internal test:** panel 0.789 vs composition 0.714 (delta 0.074), LRT p=0.0322, n=61, separation=FALSE
  <br>*[R-240.Female.noMHC.Internal_test]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Female, noMHC, External blood:** panel 0.843 vs composition 0.429 (delta 0.414), LRT p=0.000803, n=24, separation=FALSE
  <br>*[R-240.Female.noMHC.External_blood]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Male, noMHC, Train:** panel 0.888 vs composition 0.857 (delta 0.031), LRT p=1.97e-05, n=38, separation=TRUE
  <br>*[R-240.Male.noMHC.Train]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Male, noMHC, Internal test:** panel 1 vs composition 0.798 (delta 0.202), LRT p=1, n=13, separation=TRUE
  <br>*[R-240.Male.noMHC.Internal_test]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **WITHIN-DATASET RESAMPLED panel vs composition — Male, noMHC, External blood:** panel 0.7 vs composition 0.575 (delta 0.125), LRT p=0.118, n=9, separation=TRUE
  <br>*[R-240.Male.noMHC.External_blood]* — source: `results/tables/PANEL_incremental_value_LRT.csv` (cols AUC_panel, AUC_composition, delta_*, LRT_p_*) · script: `scripts/goal2_sex_stratified/17b_testing_blood_celladjusted.R`
- **Per-gene AUC (train-fixed orientation) — C6orf136, Female, Train:** 0.784 (0.707-0.862), concordant=TRUE
  <br>*[R-250.Female.Train.C6orf136]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — ESYT1, Female, Train:** 0.81 (0.741-0.879), concordant=TRUE
  <br>*[R-250.Female.Train.ESYT1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — GNL1, Female, Train:** 0.765 (0.689-0.842), concordant=TRUE
  <br>*[R-250.Female.Train.GNL1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — IKZF3, Female, Train:** 0.738 (0.658-0.819), concordant=TRUE
  <br>*[R-250.Female.Train.IKZF3]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — MED1, Female, Train:** 0.793 (0.718-0.868), concordant=TRUE
  <br>*[R-250.Female.Train.MED1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — SMARCC2, Female, Train:** 0.795 (0.721-0.868), concordant=TRUE
  <br>*[R-250.Female.Train.SMARCC2]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — C6orf136, Female, Internal test:** 0.804 (0.691-0.918), concordant=TRUE
  <br>*[R-250.Female.Internal_test.C6orf136]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — ESYT1, Female, Internal test:** 0.8 (0.691-0.909), concordant=TRUE
  <br>*[R-250.Female.Internal_test.ESYT1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — GNL1, Female, Internal test:** 0.593 (0.443-0.743), concordant=TRUE
  <br>*[R-250.Female.Internal_test.GNL1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — IKZF3, Female, Internal test:** 0.66 (0.518-0.802), concordant=TRUE
  <br>*[R-250.Female.Internal_test.IKZF3]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — MED1, Female, Internal test:** 0.713 (0.582-0.845), concordant=TRUE
  <br>*[R-250.Female.Internal_test.MED1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — SMARCC2, Female, Internal test:** 0.837 (0.731-0.942), concordant=TRUE
  <br>*[R-250.Female.Internal_test.SMARCC2]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — C6orf136, Female, External blood:** 0.964 (0.89-1), concordant=TRUE
  <br>*[R-250.Female.External_blood.C6orf136]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — ESYT1, Female, External blood:** 0.836 (0.663-1), concordant=TRUE
  <br>*[R-250.Female.External_blood.ESYT1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — GNL1, Female, External blood:** 0.7 (0.479-0.921), concordant=TRUE
  <br>*[R-250.Female.External_blood.GNL1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — IKZF3, Female, External blood:** 0.521 (0.277-0.765), concordant=TRUE
  <br>*[R-250.Female.External_blood.IKZF3]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — MED1, Female, External blood:** 0.586 (0.337-0.835), concordant=TRUE
  <br>*[R-250.Female.External_blood.MED1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — SMARCC2, Female, External blood:** 0.657 (0.432-0.882), concordant=TRUE
  <br>*[R-250.Female.External_blood.SMARCC2]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — ESYT1, Male, Train:** 0.938 (0.863-1), concordant=TRUE
  <br>*[R-250.Male.Train.ESYT1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — HLA-DMA, Male, Train:** 0.927 (0.849-1), concordant=TRUE
  <br>*[R-250.Male.Train.HLA-DMA]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — INPP5B, Male, Train:** 0.812 (0.678-0.946), concordant=TRUE
  <br>*[R-250.Male.Train.INPP5B]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — MED1, Male, Train:** 0.863 (0.727-0.999), concordant=TRUE
  <br>*[R-250.Male.Train.MED1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — SMARCC2, Male, Train:** 0.938 (0.869-1), concordant=TRUE
  <br>*[R-250.Male.Train.SMARCC2]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — VPS52, Male, Train:** 0.877 (0.771-0.982), concordant=TRUE
  <br>*[R-250.Male.Train.VPS52]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — ESYT1, Male, Internal test:** 1 (1-1), concordant=TRUE
  <br>*[R-250.Male.Internal_test.ESYT1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — HLA-DMA, Male, Internal test:** 0.952 (0.81-1), concordant=TRUE
  <br>*[R-250.Male.Internal_test.HLA-DMA]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — INPP5B, Male, Internal test:** 0.952 (0.81-1), concordant=TRUE
  <br>*[R-250.Male.Internal_test.INPP5B]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — MED1, Male, Internal test:** 0.976 (0.857-1), concordant=TRUE
  <br>*[R-250.Male.Internal_test.MED1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — SMARCC2, Male, Internal test:** 1 (1-1), concordant=TRUE
  <br>*[R-250.Male.Internal_test.SMARCC2]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — VPS52, Male, Internal test:** 1 (1-1), concordant=TRUE
  <br>*[R-250.Male.Internal_test.VPS52]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — ESYT1, Male, External blood:** 0.7 (0.25-1), concordant=TRUE
  <br>*[R-250.Male.External_blood.ESYT1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — HLA-DMA, Male, External blood:** 0.975 (0.85-1), concordant=TRUE
  <br>*[R-250.Male.External_blood.HLA-DMA]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — INPP5B, Male, External blood:** 0.55 (0-1), concordant=TRUE
  <br>*[R-250.Male.External_blood.INPP5B]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — MED1, Male, External blood:** 0.35 (0-0.8), concordant=FALSE
  <br>*[R-250.Male.External_blood.MED1]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — SMARCC2, Male, External blood:** 0.6 (0.1-1), concordant=TRUE
  <br>*[R-250.Male.External_blood.SMARCC2]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`
- **Per-gene AUC (train-fixed orientation) — VPS52, Male, External blood:** 0.35 (0-0.8), concordant=FALSE
  <br>*[R-250.Male.External_blood.VPS52]* — source: `results/tables/mr_roc_pergene_auc.csv` (cols AUC, AUC_lo, AUC_hi, concordant) · script: `scripts/goal2_sex_stratified/17_testing_blood_internal_external.R`

## §2.10 Diagnosis-by-sex interaction testing

- **Interaction genes at FDR<0.05, model = qnorm + batch_full:** 53
  <br>*[R-260.qnorm_+_batch_full]* — source: `results/tables/DEG_interaction_summary.csv` (col n_interaction_FDR05) · script: `scripts/00_shared/05b_dge_sensitivity.R`
- **Interaction genes at FDR<0.05, model = ComBat:** 270
  <br>*[R-260.ComBat]* — source: `results/tables/DEG_interaction_summary.csv` (col n_interaction_FDR05) · script: `scripts/00_shared/05b_dge_sensitivity.R`
- **Interaction pattern count: MALE-restricted:** 35
  <br>*[R-261.MALE-restricted]* — source: `results/tables/DEG_interaction_patterns.csv` (col N) · script: `scripts/00_shared/05d_interaction_report.R`
- **Interaction pattern count: MAGNITUDE difference:** 12
  <br>*[R-261.MAGNITUDE_difference]* — source: `results/tables/DEG_interaction_patterns.csv` (col N) · script: `scripts/00_shared/05d_interaction_report.R`
- **Interaction pattern count: OPPOSITE direction:** 4
  <br>*[R-261.OPPOSITE_direction]* — source: `results/tables/DEG_interaction_patterns.csv` (col N) · script: `scripts/00_shared/05d_interaction_report.R`
- **Interaction pattern count: neither sex significant:** 1
  <br>*[R-261.neither_sex_significant]* — source: `results/tables/DEG_interaction_patterns.csv` (col N) · script: `scripts/00_shared/05d_interaction_report.R`
- **Interaction pattern count: FEMALE-restricted:** 1
  <br>*[R-261.FEMALE-restricted]* — source: `results/tables/DEG_interaction_patterns.csv` (col N) · script: `scripts/00_shared/05d_interaction_report.R`
- **Sex-differential genes (primary model):** 53
  <br>*[R-262]* — source: `results/tables/DEG_interaction_significant.csv` (row count) · script: `scripts/00_shared/05d_interaction_report.R`
- **Sex-differential genes surviving composition adjustment:** 7/53 (ZNF800, BCLAF1, MAP4K5, UBXN4, RIF1, RAB2A, KDM5D)
  <br>*[R-263]* — source: `results/tables/DEG_interaction_significant.csv` (count interaction_FDR_celladj < 0.05) · script: `scripts/00_shared/05d_interaction_report.R`
- **Leading sex-differential genes:** ZNF800 int +0.750 FDR 0.0024 [F -0.104 (0.276), M +0.727 (2.1e-05)]; BCLAF1 int +0.796 FDR 0.0024 [F -0.128 (0.186), M +0.868 (1.18e-05)]; MAP4K5 int +1.057 FDR 0.0024 [F -0.144 (0.314), M +1.053 (6.3e-06)]; UBXN4 int +0.517 FDR 0.0052 [F +0.007 (0.944), M +0.620 (2.04e-05)]; RIF1 int +0.829 FDR 0.0056 [F -0.047 (0.722), M +1.012 (3.11e-06)]; RAB2A int +0.419 FDR 0.0067 [F +0.091 (0.097), M +0.534 (1.44e-05)]
  <br>*[R-264]* — source: `results/tables/DEG_interaction_significant.csv` (6 smallest interaction_FDR) · script: `scripts/00_shared/05d_interaction_report.R`
- **PTPN22 interaction:** int +0.585 FDR 0.0382; F +0.171 (FDR 0.0457); M +0.893 (FDR 1.18e-06)
  <br>*[R-265]* — source: `results/tables/DEG_interaction_significant.csv` (row gene==PTPN22) · script: `scripts/00_shared/05d_interaction_report.R`
- **Sex-chromosome-linked genes among the sex-differential set:** KDM5D
  <br>*[R-266]* — source: `results/tables/DEG_interaction_significant.csv` (membership test against a Y-linked/XIST gene list) · script: `scripts/00_shared/05d_interaction_report.R`
- **Panel genes among the sex-differential set:** NONE
  <br>*[R-267]* — source: `results/tables/DEG_interaction_significant.csv + results/tables/mr_fs_summary*.csv` (set intersection) · script: `scripts/00_shared/05d_interaction_report.R`, `scripts/goal2_sex_stratified/12_feature_selection.R`, `scripts/goal2_sex_stratified/12b_feature_selection_noMHC.R`

## §2.11 Cross-tissue evaluation (synovium)

_Methodology: see [`results/METHODS_2.11_crosstissue.md`](../results/METHODS_2.11_crosstissue.md) for the full specification of this step._

- **Synovium per-gene, C6orf136 (female):** log2FC -1.121 (FDR 2.14e-21), concordant=TRUE, AUC best-direction 0.906, AUC train-oriented 0.906
  <br>*[R-270.female.C6orf136]* — source: `results/tables/val_synovium_pergene_female.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Synovium per-gene, ESYT1 (female):** log2FC -0.215 (FDR 3.68e-02), concordant=TRUE, AUC best-direction 0.511, AUC train-oriented 0.511
  <br>*[R-270.female.ESYT1]* — source: `results/tables/val_synovium_pergene_female.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Synovium per-gene, GNL1 (female):** log2FC -1.235 (FDR 4.98e-28), concordant=TRUE, AUC best-direction 0.927, AUC train-oriented 0.927
  <br>*[R-270.female.GNL1]* — source: `results/tables/val_synovium_pergene_female.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Synovium per-gene, IKZF3 (female):** log2FC +2.532 (FDR 5.04e-04), concordant=FALSE, AUC best-direction 0.901, AUC train-oriented 0.099
  <br>*[R-270.female.IKZF3]* — source: `results/tables/val_synovium_pergene_female.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Synovium per-gene, MED1 (female):** log2FC +0.493 (FDR 3.30e-10), concordant=FALSE, AUC best-direction 0.768, AUC train-oriented 0.232
  <br>*[R-270.female.MED1]* — source: `results/tables/val_synovium_pergene_female.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Synovium per-gene, SMARCC2 (female):** log2FC -1.105 (FDR 1.13e-32), concordant=TRUE, AUC best-direction 0.996, AUC train-oriented 0.996
  <br>*[R-270.female.SMARCC2]* — source: `results/tables/val_synovium_pergene_female.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Panel genes concordant in direction blood->synovium, female:** 4/6
  <br>*[R-271.female]* — source: `results/tables/val_synovium_pergene_female.csv` (count concordant==TRUE) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Synovium per-gene, ESYT1 (male):** log2FC -0.215 (FDR 3.68e-02), concordant=TRUE, AUC best-direction 0.696, AUC train-oriented 0.696
  <br>*[R-270.male.ESYT1]* — source: `results/tables/val_synovium_pergene_male.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Synovium per-gene, HLA-DMA (male):** log2FC +1.240 (FDR 7.16e-09), concordant=FALSE, AUC best-direction 0.759, AUC train-oriented 0.241
  <br>*[R-270.male.HLA-DMA]* — source: `results/tables/val_synovium_pergene_male.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Synovium per-gene, INPP5B (male):** log2FC -0.197 (FDR 4.35e-03), concordant=TRUE, AUC best-direction 0.575, AUC train-oriented 0.575
  <br>*[R-270.male.INPP5B]* — source: `results/tables/val_synovium_pergene_male.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Synovium per-gene, MED1 (male):** log2FC +0.493 (FDR 3.30e-10), concordant=FALSE, AUC best-direction 0.672, AUC train-oriented 0.328
  <br>*[R-270.male.MED1]* — source: `results/tables/val_synovium_pergene_male.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Synovium per-gene, SMARCC2 (male):** log2FC -1.105 (FDR 1.13e-32), concordant=TRUE, AUC best-direction 0.983, AUC train-oriented 0.983
  <br>*[R-270.male.SMARCC2]* — source: `results/tables/val_synovium_pergene_male.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Synovium per-gene, VPS52 (male):** log2FC -0.388 (FDR 1.51e-08), concordant=TRUE, AUC best-direction 0.839, AUC train-oriented 0.839
  <br>*[R-270.male.VPS52]* — source: `results/tables/val_synovium_pergene_male.csv` (cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **Panel genes concordant in direction blood->synovium, male:** 4/6
  <br>*[R-271.male]* — source: `results/tables/val_synovium_pergene_male.csv` (count concordant==TRUE) · script: `scripts/goal2_sex_stratified/20_testing_synovium_external.R`
- **WITHIN-SYNOVIUM refit panel AUC, Female:** apparent 1, 10-fold CV 0.986 (0.969-1), n=120 (RA 106, Normal 14)
  <br>*[R-275.Female]* — source: `results/tables/crosstissue_panel_auc.csv` (cols apparent_AUC, CV_AUC, CV_lo, CV_hi, n, n_RA, n_Normal) · script: `scripts/goal2_sex_stratified/22_crosstissue_biomarker_discovery.R`
- **WITHIN-SYNOVIUM refit panel AUC, Male:** apparent 0.994, 10-fold CV 0.776 (0.573-0.978), n=60 (RA 46, Normal 14)
  <br>*[R-275.Male]* — source: `results/tables/crosstissue_panel_auc.csv` (cols apparent_AUC, CV_AUC, CV_lo, CV_hi, n, n_RA, n_Normal) · script: `scripts/goal2_sex_stratified/22_crosstissue_biomarker_discovery.R`

## §2.12 Cross-ancestry evaluation

_Methodology: see [`results/METHODS_2.12_crossancestry.md`](../results/METHODS_2.12_crossancestry.md) for the full specification of this step._

- **Cross-ancestry classification, female:** tested 32: EUR-replicated 23, EAS-transferable 7, EAS-untestable 2, EUR-only 5
  <br>*[R-280.female]* — source: `results/tables/MR35_crossancestry_summary.csv` (all columns of that row) · script: `scripts/goal2_sex_stratified/26_crossancestry_biomarker_mr.R`
- **Cross-ancestry classification, male:** tested 25: EUR-replicated 18, EAS-transferable 6, EAS-untestable 2, EUR-only 3
  <br>*[R-280.male]* — source: `results/tables/MR35_crossancestry_summary.csv` (all columns of that row) · script: `scripts/goal2_sex_stratified/26_crossancestry_biomarker_mr.R`

## §2.13 Functional enrichment

_Methodology: see [`results/METHODS_2.13_functional_enrichment.md`](../results/METHODS_2.13_functional_enrichment.md) for the full specification of this step._

- **Enriched terms, female DEGs:** 1408
  <br>*[R-290.female]* — source: `results/tables/dge_enriched_terms_female.csv` (row count) · script: `scripts/00_shared/05_dge.R`
- **Enriched terms, male DEGs:** 1513
  <br>*[R-290.male]* — source: `results/tables/dge_enriched_terms_male.csv` (row count) · script: `scripts/00_shared/05_dge.R`
- **Enriched terms, all DEGs:** 1659
  <br>*[R-290.all]* — source: `results/tables/dge_enriched_terms_all.csv` (row count) · script: `scripts/00_shared/05_dge.R`
- **GSEA KEGG, All:** 19 significant, 18 suppressed, 1 activated
  <br>*[R-295.All]* — source: `results/tables/gsea_kegg.csv` (rows group==All, sign of NES) · script: `scripts/00_shared/05_dge.R`
- **Most suppressed KEGG sets, All:** ATP-dependent chromatin remodeling (NES -1.97, FDR 1.0e-03); Carbon metabolism (NES -1.92, FDR 1.0e-03); Primary immunodeficiency (NES -1.87, FDR 2.6e-02); Biosynthesis of amino acids (NES -1.86, FDR 1.0e-02)
  <br>*[R-296.All]* — source: `results/tables/gsea_kegg.csv` (4 most negative NES) · script: `scripts/00_shared/05_dge.R`
- **Activated KEGG sets, All:** Ribosome (NES +1.89)
  <br>*[R-297.All]* — source: `results/tables/gsea_kegg.csv` (rows with NES>0) · script: `scripts/00_shared/05_dge.R`
- **GSEA KEGG, Female:** 17 significant, 16 suppressed, 1 activated
  <br>*[R-295.Female]* — source: `results/tables/gsea_kegg.csv` (rows group==Female, sign of NES) · script: `scripts/00_shared/05_dge.R`
- **Most suppressed KEGG sets, Female:** Primary immunodeficiency (NES -2.23, FDR 3.6e-04); ATP-dependent chromatin remodeling (NES -2.02, FDR 3.6e-04); Spliceosome (NES -1.90, FDR 8.6e-04); Biosynthesis of amino acids (NES -1.80, FDR 2.2e-02)
  <br>*[R-296.Female]* — source: `results/tables/gsea_kegg.csv` (4 most negative NES) · script: `scripts/00_shared/05_dge.R`
- **Activated KEGG sets, Female:** Ribosome (NES +1.89)
  <br>*[R-297.Female]* — source: `results/tables/gsea_kegg.csv` (rows with NES>0) · script: `scripts/00_shared/05_dge.R`
- **GSEA KEGG, Male:** 22 significant, 21 suppressed, 1 activated
  <br>*[R-295.Male]* — source: `results/tables/gsea_kegg.csv` (rows group==Male, sign of NES) · script: `scripts/00_shared/05_dge.R`
- **Most suppressed KEGG sets, Male:** Carbon metabolism (NES -2.15, FDR 5.0e-05); Citrate cycle (TCA cycle) (NES -1.99, FDR 1.7e-02); Fc gamma R-mediated phagosome formation (NES -1.95, FDR 1.7e-03); Other glycan degradation (NES -1.89, FDR 3.3e-02)
  <br>*[R-296.Male]* — source: `results/tables/gsea_kegg.csv` (4 most negative NES) · script: `scripts/00_shared/05_dge.R`
- **Activated KEGG sets, Male:** Ribosome (NES +1.82)
  <br>*[R-297.Male]* — source: `results/tables/gsea_kegg.csv` (rows with NES>0) · script: `scripts/00_shared/05_dge.R`
- **GO:BP terms enriched among MR-prioritised genes, female:** 46
  <br>*[R-305.female]* — source: `results/tables/enrich_GO_BP_female.csv` (row count) · script: `scripts/goal2_sex_stratified/34_pathway_enrichment.R`
- **Top GO:BP terms, MR-prioritised female genes:** antigen processing and presentation of exogenous peptide antigen via MHC class II (4/31, FDR 1.2e-04); antigen processing and presentation of peptide antigen via MHC class II (4/31, FDR 1.2e-04); antigen processing and presentation of peptide or polysaccharide antigen via MHC class II (4/31, FDR 1.2e-04); regulation of lymphocyte activation (9/31, FDR 1.2e-04); antigen processing and presentation of exogenous peptide antigen (4/31, FDR 1.8e-04)
  <br>*[R-306.female]* — source: `results/tables/enrich_GO_BP_female.csv` (first 5 rows as written (sorted by p.adjust)) · script: `scripts/goal2_sex_stratified/34_pathway_enrichment.R`
- **GO:BP terms enriched among MR-prioritised genes, male:** 45
  <br>*[R-305.male]* — source: `results/tables/enrich_GO_BP_male.csv` (row count) · script: `scripts/goal2_sex_stratified/34_pathway_enrichment.R`
- **Top GO:BP terms, MR-prioritised male genes:** positive regulation of T cell activation (5/24, FDR 1.8e-02); positive regulation of leukocyte cell-cell adhesion (5/24, FDR 1.8e-02); regulation of lymphocyte activation (6/24, FDR 1.8e-02); positive regulation of lymphocyte activation (5/24, FDR 1.8e-02); positive regulation of cell-cell adhesion (5/24, FDR 1.8e-02)
  <br>*[R-306.male]* — source: `results/tables/enrich_GO_BP_male.csv` (first 5 rows as written (sorted by p.adjust)) · script: `scripts/goal2_sex_stratified/34_pathway_enrichment.R`
- **KEGG pathways enriched among MR-prioritised genes, female:** 25
  <br>*[R-307.female]* — source: `results/tables/enrich_KEGG_female.csv` (row count) · script: `scripts/goal2_sex_stratified/34_pathway_enrichment.R`
- **Top KEGG pathways, MR-prioritised female genes:** Staphylococcus aureus infection (4/19, FDR 1.4e-03); Asthma (3/19, FDR 1.4e-03); Allograft rejection (3/19, FDR 1.6e-03); Graft-versus-host disease (3/19, FDR 1.6e-03); Type I diabetes mellitus (3/19, FDR 1.7e-03)
  <br>*[R-308.female]* — source: `results/tables/enrich_KEGG_female.csv` (first 5 rows as written (sorted by p.adjust)) · script: `scripts/goal2_sex_stratified/34_pathway_enrichment.R`
- **KEGG pathways enriched among MR-prioritised genes, male:** 14
  <br>*[R-307.male]* — source: `results/tables/enrich_KEGG_male.csv` (row count) · script: `scripts/goal2_sex_stratified/34_pathway_enrichment.R`
- **Top KEGG pathways, MR-prioritised male genes:** Herpes simplex virus 1 infection (4/16, FDR 1.6e-02); Leishmaniasis (3/16, FDR 1.6e-02); Toxoplasmosis (3/16, FDR 2.7e-02); Asthma (2/16, FDR 2.7e-02); Allograft rejection (2/16, FDR 3.0e-02)
  <br>*[R-308.male]* — source: `results/tables/enrich_KEGG_male.csv` (first 5 rows as written (sorted by p.adjust)) · script: `scripts/goal2_sex_stratified/34_pathway_enrichment.R`
- **GO terms enriched in disease-module genes:** 218
  <br>*[R-300]* — source: `results/tables/WGCNA_08_disease_GO.csv` (row count) · script: `scripts/00_shared/06_WGCNA.R`
- **KEGG pathways enriched in disease-module genes:** 16
  <br>*[R-301]* — source: `results/tables/WGCNA_09_disease_KEGG.csv` (row count) · script: `scripts/00_shared/06_WGCNA.R`
- **Top disease-module KEGG pathways:** Spliceosome (43/1224, FDR 3.7e-05); Citrate cycle (TCA cycle) (17/1224, FDR 1.9e-04); Polycomb repressive complex (28/1224, FDR 5.9e-03); ATP-dependent chromatin remodeling (29/1224, FDR 6.1e-03); 2-Oxocarboxylic acid metabolism (15/1224, FDR 6.1e-03); Nucleocytoplasmic transport (31/1224, FDR 6.1e-03)
  <br>*[R-302]* — source: `results/tables/WGCNA_09_disease_KEGG.csv` (first 6 rows as written) · script: `scripts/00_shared/06_WGCNA.R`
- **Top disease-module GO terms:** mRNA processing (118/2342, FDR 1.5e-07); ribonucleoprotein complex biogenesis (111/2342, FDR 1.5e-07); RNA splicing, via transesterification reactions (83/2342, FDR 1.8e-07); RNA splicing, via transesterification reactions with bulged adenosine as nucleophile (82/2342, FDR 1.8e-07); mRNA splicing, via spliceosome (82/2342, FDR 1.8e-07)
  <br>*[R-303]* — source: `results/tables/WGCNA_08_disease_GO.csv` (first 5 rows as written) · script: `scripts/00_shared/06_WGCNA.R`

## §2.14 Immune deconvolution and composition-adjusted expression

_Methodology: see [`results/METHODS_2.14_deconvolution.md`](../results/METHODS_2.14_deconvolution.md) for the full specification of this step._

- **Cell fraction RA vs HC, Female T_cells_CD8:** HC 0.1622 -> RA 0.0931 (diff -0.0691, FDR 1.61e-08)
  <br>*[R-310.Female.T_cells_CD8]* — source: `results/tables/CELL_fraction_group_tests.csv` (rows with FDR<0.05) · script: `scripts/00_shared/05c_deconvolution.R`
- **Cell fraction RA vs HC, Female T_cells_gamma_delta:** HC 0.0079 -> RA 0.0207 (diff 0.0128, FDR 1.73e-02)
  <br>*[R-310.Female.T_cells_gamma_delta]* — source: `results/tables/CELL_fraction_group_tests.csv` (rows with FDR<0.05) · script: `scripts/00_shared/05c_deconvolution.R`
- **Cell fraction RA vs HC, Female B_cells_naive:** HC 0.0457 -> RA 0.0335 (diff -0.0121, FDR 4.65e-02)
  <br>*[R-310.Female.B_cells_naive]* — source: `results/tables/CELL_fraction_group_tests.csv` (rows with FDR<0.05) · script: `scripts/00_shared/05c_deconvolution.R`
- **Cell fraction RA vs HC, Male T_cells_gamma_delta:** HC 0.0039 -> RA 0.0288 (diff 0.0249, FDR 1.45e-02)
  <br>*[R-310.Male.T_cells_gamma_delta]* — source: `results/tables/CELL_fraction_group_tests.csv` (rows with FDR<0.05) · script: `scripts/00_shared/05c_deconvolution.R`
- **Cell fraction RA vs HC, Male Eosinophils:** HC 0 -> RA 0.0045 (diff 0.0045, FDR 1.45e-02)
  <br>*[R-310.Male.Eosinophils]* — source: `results/tables/CELL_fraction_group_tests.csv` (rows with FDR<0.05) · script: `scripts/00_shared/05c_deconvolution.R`
- **Cell fraction RA vs HC, Male T_cells_CD8:** HC 0.1542 -> RA 0.0905 (diff -0.0637, FDR 1.83e-02)
  <br>*[R-310.Male.T_cells_CD8]* — source: `results/tables/CELL_fraction_group_tests.csv` (rows with FDR<0.05) · script: `scripts/00_shared/05c_deconvolution.R`
- **Cell fraction RA vs HC, Male T_cells_CD4_memory_resting:** HC 0 -> RA 0.0088 (diff 0.0088, FDR 1.83e-02)
  <br>*[R-310.Male.T_cells_CD4_memory_resting]* — source: `results/tables/CELL_fraction_group_tests.csv` (rows with FDR<0.05) · script: `scripts/00_shared/05c_deconvolution.R`
- **Cell subsets differing by disease status at FDR<0.05:** 3 female, 4 male
  <br>*[R-311]* — source: `results/tables/CELL_fraction_group_tests.csv` (count FDR<0.05 by sex) · script: `scripts/00_shared/05c_deconvolution.R`
- **Composition-adjusted DEG count, Female:** 5131 -> 2709 (49.4% retained), Spearman rho 0.889
  <br>*[R-315.Female]* — source: `results/tables/DEG_celladjusted_summary.csv` (cols n_DEG_unadjusted, n_DEG_adjusted, pct_retained, logFC_spearman) · script: `scripts/00_shared/05c_deconvolution.R`
- **Composition-adjusted DEG count, Male:** 5820 -> 1450 (23.8% retained), Spearman rho 0.842
  <br>*[R-315.Male]* — source: `results/tables/DEG_celladjusted_summary.csv` (cols n_DEG_unadjusted, n_DEG_adjusted, pct_retained, logFC_spearman) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Female ESYT1:** logFC -0.381 -> -0.28 (73.5% retained), retained
  <br>*[R-320.Female.ESYT1]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Female SMARCC2:** logFC -0.336 -> -0.208 (61.8% retained), retained
  <br>*[R-320.Female.SMARCC2]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Female MED1:** logFC -0.188 -> -0.141 (75.2% retained), retained
  <br>*[R-320.Female.MED1]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Female IKZF3:** logFC -0.428 -> -0.34 (79.5% retained), retained
  <br>*[R-320.Female.IKZF3]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Female GNL1:** logFC -0.216 -> -0.159 (73.4% retained), retained
  <br>*[R-320.Female.GNL1]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Female C6orf136:** logFC -0.19 -> -0.125 (65.6% retained), retained
  <br>*[R-320.Female.C6orf136]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Male ESYT1:** logFC -0.484 -> -0.401 (82.8% retained), retained
  <br>*[R-320.Male.ESYT1]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Male MED1:** logFC -0.237 -> -0.292 (123% retained), retained
  <br>*[R-320.Male.MED1]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Male SMARCC2:** logFC -0.532 -> -0.386 (72.7% retained), retained
  <br>*[R-320.Male.SMARCC2]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Male HLA-DMA:** logFC -0.394 -> -0.314 (79.7% retained), retained
  <br>*[R-320.Male.HLA-DMA]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Male VPS52:** logFC -0.286 -> -0.199 (69.6% retained), lost to composition
  <br>*[R-320.Male.VPS52]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene composition robustness, Male INPP5B:** logFC -0.387 -> -0.219 (56.7% retained), lost to composition
  <br>*[R-320.Male.INPP5B]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (cols logFC_unadj, logFC_adj, pct_logFC_retained, status) · script: `scripts/00_shared/05c_deconvolution.R`
- **Panel gene-sex pairs retained after composition adjustment:** 10/12
  <br>*[R-321]* — source: `results/tables/CELL_panel_gene_adjustment.csv` (count status=='retained') · script: `scripts/00_shared/05c_deconvolution.R`

## §2.15 Nomogram construction and clinical evaluation

_Methodology: see [`results/METHODS_2.15_clinical_utility_nomogram.md`](../results/METHODS_2.15_clinical_utility_nomogram.md) for the full specification of this step._

> **Provenance note.** This section had no executed outputs anywhere in the repository before this report was built — `results/tables/diag_dca_{female,male}.csv` and `results/figures/new/fig_diag_validation_{female,male}.png/pdf` did not exist. `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R` was run once (inputs: `data/processed/combined_train.rds`, `data/processed/new/ml_features.rds`, both already present) to generate them. The numbers below are read from that run's output, not asserted.

- **Nomogram panel fit on the training cohort, female:** 6 genes, n=145 (RA=86, HC=59)
  <br>*[R-330.female]* — source: `results/logs/19_testing_blood_clinical_utility.log` (literal 'wrote fig_diag_validation_female.png' line) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Calibration error (200-rep bootstrap bias-correction), female:** MAE 0.036, MSE 0.00155, 90th-percentile absolute error 0.053
  <br>*[R-331.female]* — source: `results/logs/19_testing_blood_clinical_utility.log` (literal rms::calibrate() console block preceding the wrote-line) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Net benefit at threshold 0.10, female:** panel 0.5563 vs treat-all 0.5479 vs treat-none 0.0000
  <br>*[R-332.female.0.10]* — source: `results/tables/diag_dca_female.csv` (row nearest threshold==0.10) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Net benefit at threshold 0.20, female:** panel 0.5069 vs treat-all 0.4914 vs treat-none 0.0000
  <br>*[R-332.female.0.20]* — source: `results/tables/diag_dca_female.csv` (row nearest threshold==0.20) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Net benefit at threshold 0.30, female:** panel 0.4759 vs treat-all 0.4187 vs treat-none 0.0000
  <br>*[R-332.female.0.30]* — source: `results/tables/diag_dca_female.csv` (row nearest threshold==0.30) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Net benefit at threshold 0.50, female:** panel 0.3655 vs treat-all 0.1862 vs treat-none 0.0000
  <br>*[R-332.female.0.50]* — source: `results/tables/diag_dca_female.csv` (row nearest threshold==0.50) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Threshold range where the panel beats treat-all, female:** 0.01-0.99
  <br>*[R-333.female]* — source: `results/tables/diag_dca_female.csv` (thresholds where NB_panel > NB_all) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Threshold range with positive net benefit, female:** 0.01-0.99
  <br>*[R-334.female]* — source: `results/tables/diag_dca_female.csv` (thresholds where NB_panel > 0) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Nomogram panel fit on the training cohort, male:** 6 genes, n=38 (RA=17, HC=21)
  <br>*[R-330.male]* — source: `results/logs/19_testing_blood_clinical_utility.log` (literal 'wrote fig_diag_validation_male.png' line) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Calibration error (200-rep bootstrap bias-correction), male:** MAE 0.112, MSE 0.01286, 90th-percentile absolute error 0.127
  <br>*[R-331.male]* — source: `results/logs/19_testing_blood_clinical_utility.log` (literal rms::calibrate() console block preceding the wrote-line) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Net benefit at threshold 0.10, male:** panel 0.4094 vs treat-all 0.3860 vs treat-none 0.0000
  <br>*[R-332.male.0.10]* — source: `results/tables/diag_dca_male.csv` (row nearest threshold==0.10) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Net benefit at threshold 0.20, male:** panel 0.4013 vs treat-all 0.3092 vs treat-none 0.0000
  <br>*[R-332.male.0.20]* — source: `results/tables/diag_dca_male.csv` (row nearest threshold==0.20) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Net benefit at threshold 0.30, male:** panel 0.4023 vs treat-all 0.2105 vs treat-none 0.0000
  <br>*[R-332.male.0.30]* — source: `results/tables/diag_dca_male.csv` (row nearest threshold==0.30) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Net benefit at threshold 0.50, male:** panel 0.4211 vs treat-all -0.1053 vs treat-none 0.0000
  <br>*[R-332.male.0.50]* — source: `results/tables/diag_dca_male.csv` (row nearest threshold==0.50) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Threshold range where the panel beats treat-all, male:** 0.01-0.99
  <br>*[R-333.male]* — source: `results/tables/diag_dca_male.csv` (thresholds where NB_panel > NB_all) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`
- **Threshold range with positive net benefit, male:** 0.01-0.99
  <br>*[R-334.male]* — source: `results/tables/diag_dca_male.csv` (thresholds where NB_panel > 0) · script: `scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R`

---

## How to verify any number above

1. Find the claim's `[R-xxx]` tag.
2. Open `results/RESULTS_PROVENANCE.tsv` (same tag, one row) or the `results/tables/*.csv` file named in the source line directly above.
3. Apply the stated derivation (a column read, a filter, a set operation) by eye — every derivation here is a single filter/aggregation, not a multi-step calculation.
4. To re-derive from scratch, re-run the script named after `script:` and then `python3 scripts/verify_results_numbers.py`.

