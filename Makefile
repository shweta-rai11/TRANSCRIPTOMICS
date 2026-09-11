# =============================================================================
# Makefile — enforces the pipeline execution order.
#
# WHY THIS EXISTS
#   Before this file, the pipeline's dependency order existed only as an
#   implicit convention (numbered script filenames) with nothing to enforce
#   it. Running a script out of order silently produces wrong or stale
#   output rather than failing loudly. This file makes the order explicit
#   and checkable.
#
# HONEST SCOPE — read this before trusting it more than it claims
#   This is a PHASE-level dependency graph, not an exhaustive file-level one.
#   Each target is a `.PHONY` group of scripts known (from the header comment
#   of each script and the CSV_SCRIPT map in scripts/verify_results_numbers.py)
#   to belong to the same pipeline stage, in the order that stage's scripts
#   must run. It does NOT track every individual .rds/.csv dependency between
#   every pair of scripts the way a true `make` file-dependency graph would —
#   doing that correctly for ~40 R scripts was judged higher-risk (getting one
#   dependency wrong and silently skipping a stale rebuild) than valuable,
#   given this pipeline is re-run end-to-end rather than incrementally. What
#   this DOES guarantee: `make all` runs every stage in the right order, and
#   a missing prerequisite output fails the build immediately instead of
#   letting a later stage run on stale or absent input.
#
#   Practice followed: Wilson G, Bryan J, Cranston K, Kitzes J, Nederbragt L,
#   Teal TK (2017). "Good enough practices in scientific computing." PLOS
#   Computational Biology 13(6):e1005510. — explicit build automation over an
#   ad hoc "run these files in this order" README.
#
# USAGE
#   make check        # verify all required input data/caches are present, no execution
#   make preprocess    # phase 0: load, normalise, batch-correct, partition
#   make expression     # phase 1: DEG, sensitivity, deconvolution, interaction
#   make network         # phase 2: WGCNA
#   make candidates        # phase 3: disease-module x DEG intersection
#   make causal              # phase 4: MR, MHC sensitivity, colocalisation
#   make features               # phase 5: LASSO/RF/SVM-RFE consensus panels
#   make models                    # phase 6: nested CV, elastic net, final panels
#   make evaluate                     # phase 7: internal/external blood testing, nomogram
#   make crosstissue                     # phase 8: synovium validation
#   make crossancestry                      # phase 9: cross-ancestry MR
#   make figures                               # phase 10: figures + enrichment
#   make provenance                               # regenerate RESULTS.md / .tsv + run tests
#   make all                                          # phases 0-10 + provenance, in order
#
#   Each phase target depends on the previous phase's target, so `make evaluate`
#   will run everything upstream of it first if it hasn't already.
# =============================================================================

SHELL := /bin/bash
RSCRIPT := Rscript
SHARED := scripts/00_shared
GOAL2 := scripts/goal2_sex_stratified
PY := python3

.PHONY: all check preprocess expression network candidates causal features models evaluate crosstissue crossancestry figures provenance test clean-logs

all: provenance

# ---- phase 0: data loading, normalisation, batch correction, partition -----
preprocess:
	$(RSCRIPT) $(SHARED)/02_data_loading.R
	$(RSCRIPT) $(SHARED)/03_dataset_gene_overlap_venn.R
	$(RSCRIPT) $(SHARED)/03_normalize_batch.R
	$(RSCRIPT) $(SHARED)/03_normalize_batch_figure.R
	$(RSCRIPT) $(SHARED)/04_apply_holdout.R
	$(RSCRIPT) $(SHARED)/load_GSE15573.R
	@test -f results/tables/combined_cohort_summary.csv || \
		(echo "preprocess phase did not produce combined_cohort_summary.csv" && exit 1)

# ---- phase 1: differential expression, sensitivity, composition, interaction
expression: preprocess
	$(RSCRIPT) $(SHARED)/05_dge.R
	$(RSCRIPT) $(SHARED)/05_dge_figure.R
	$(RSCRIPT) $(SHARED)/05b_dge_sensitivity.R
	$(RSCRIPT) $(SHARED)/05c_deconvolution.R
	$(RSCRIPT) $(SHARED)/05d_interaction_report.R
	@test -f results/tables/DEG_summary.csv || \
		(echo "expression phase did not produce DEG_summary.csv" && exit 1)

# ---- phase 2: co-expression network (WGCNA) --------------------------------
network: expression
	$(RSCRIPT) $(SHARED)/06_WGCNA.R
	@test -f results/tables/WGCNA_03_disease_modules.csv || \
		(echo "network phase did not produce WGCNA_03_disease_modules.csv" && exit 1)

# ---- phase 3: candidate genes (disease module x sex-stratified DEG) --------
candidates: network
	$(RSCRIPT) $(SHARED)/08_module_trait_RA_control.R
	$(RSCRIPT) $(SHARED)/09_disease_module_deg_intersect.R
	$(RSCRIPT) $(SHARED)/09b_disease_module_deg_venn.R
	@test -f results/tables/candidate_summary.csv || \
		(echo "candidates phase did not produce candidate_summary.csv" && exit 1)

# ---- phase 4: causal screening (MR, MHC sensitivity, colocalisation) -------
causal: candidates
	$(RSCRIPT) $(SHARED)/10_MR.R
	$(RSCRIPT) $(SHARED)/10c_MR_mhc_sensitivity.R
	$(RSCRIPT) $(SHARED)/10d_coloc_panel_genes.R
	$(RSCRIPT) $(SHARED)/10e_coloc_susie_mhc.R
	@test -f results/tables/MR_MHC_sensitivity_summary.csv || \
		(echo "causal phase did not produce MR_MHC_sensitivity_summary.csv" && exit 1)

# ---- phase 5: feature selection (LASSO / RF / SVM-RFE consensus) ----------
features: causal
	$(RSCRIPT) $(GOAL2)/12_feature_selection.R
	$(RSCRIPT) $(GOAL2)/12b_feature_selection_noMHC.R
	$(RSCRIPT) $(GOAL2)/13_feature_selection_venn.R
	@test -f results/tables/mr_fs_summary.csv || \
		(echo "features phase did not produce mr_fs_summary.csv" && exit 1)

# ---- phase 6: model training (nested CV, elastic net, final panels) -------
models: features
	$(RSCRIPT) $(GOAL2)/14_model_training_nested_cv.R
	$(RSCRIPT) $(GOAL2)/15_model_training_elasticnet.R
	$(RSCRIPT) $(GOAL2)/16_model_training_final_panel.R
	$(RSCRIPT) $(GOAL2)/16b_model_training_final_panel_noMHC.R
	$(RSCRIPT) $(GOAL2)/16c_model_training_nested_cv_transcriptomewide.R
	$(RSCRIPT) $(GOAL2)/16d_nested_cv_reconciliation.R
	@test -f results/tables/NESTED_CV_AUTHORITATIVE.csv || \
		(echo "models phase did not produce NESTED_CV_AUTHORITATIVE.csv" && exit 1)

# ---- phase 7: blood evaluation (internal/external, per-gene ROC, nomogram) -
evaluate: models
	$(RSCRIPT) $(GOAL2)/17_testing_blood_internal_external.R
	$(RSCRIPT) $(GOAL2)/17b_testing_blood_celladjusted.R
	$(RSCRIPT) $(GOAL2)/18_testing_blood_pergene_roc.R
	$(RSCRIPT) $(GOAL2)/19_testing_blood_clinical_utility.R
	@test -f results/tables/mr_roc_panel_auc.csv || \
		(echo "evaluate phase did not produce mr_roc_panel_auc.csv" && exit 1)

# ---- phase 8: cross-tissue (synovium) validation ---------------------------
crosstissue: evaluate
	$(RSCRIPT) $(GOAL2)/20_testing_synovium_external.R
	$(RSCRIPT) $(GOAL2)/21_testing_synovium_roc.R
	$(RSCRIPT) $(GOAL2)/22_crosstissue_biomarker_discovery.R
	$(RSCRIPT) $(GOAL2)/23_crosstissue_roc_all_datasets.R
	$(RSCRIPT) $(GOAL2)/24_crosstissue_pergene_auc.R
	$(RSCRIPT) $(GOAL2)/25_crosstissue_pergene_roc.R
	@test -f results/tables/crosstissue_panel_auc.csv || \
		(echo "crosstissue phase did not produce crosstissue_panel_auc.csv" && exit 1)

# ---- phase 9: cross-ancestry evaluation ------------------------------------
crossancestry: evaluate
	$(RSCRIPT) $(GOAL2)/26_crossancestry_biomarker_mr.R
	$(RSCRIPT) $(GOAL2)/27_crossancestry_biomarker_figure.R
	@test -f results/tables/MR35_crossancestry_summary.csv || \
		(echo "crossancestry phase did not produce MR35_crossancestry_summary.csv" && exit 1)

# ---- phase 10: figures + functional enrichment -----------------------------
figures: crosstissue crossancestry
	$(RSCRIPT) $(GOAL2)/28_figure_final_roc.R
	$(RSCRIPT) $(GOAL2)/29_figure_nested_cv_roc.R
	$(RSCRIPT) $(GOAL2)/30_figure_fig4A_expression_groups.R
	$(RSCRIPT) $(GOAL2)/31_figure_fig4BC_logfc_heatmaps.R
	$(RSCRIPT) $(GOAL2)/32_figure_fig4C_clustered_heatmap.R
	$(RSCRIPT) $(GOAL2)/33_figure_fig4_composite.R
	$(RSCRIPT) $(GOAL2)/34_pathway_enrichment.R
	$(RSCRIPT) $(GOAL2)/35_figure_mr_diagnostics.R
	$(RSCRIPT) $(GOAL2)/36_figure_pathways_by_sex.R

# ---- provenance: regenerate the verification report and run the test suite -
provenance: figures test
	$(PY) scripts/verify_results_numbers.py

test:
	$(PY) scripts/test_verify_results_numbers.py

# ---- sanity check only: confirms upstream data caches exist, runs nothing --
check:
	@test -f data/processed/combined_train.rds && echo "OK: data/processed/combined_train.rds present" \
		|| echo "MISSING: data/processed/combined_train.rds -- run 'make preprocess' first"
	@$(PY) scripts/verify_results_numbers.py --check

clean-logs:
	@echo "Not implemented deliberately: results/logs/*.log are provenance records" \
		"(cited by scripts/verify_results_numbers.py). Delete individually if you" \
		"really mean to invalidate a specific run's record."
