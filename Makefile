# Makefile - enforces the pipeline execution order as phase-level .PHONY targets, not a file-level dependency graph.
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
PREP := scripts/00_preprocess
EXPR := scripts/01_expression
NET := scripts/02_network
CAND := scripts/03_candidates
CAUSAL := scripts/04_causal
FEAT := scripts/05_features
MODELS := scripts/06_models
EVAL := scripts/07_evaluate
CROSST := scripts/08_crosstissue
CROSSA := scripts/09_crossancestry
GFIG := scripts/10_figures
PY := python3

.PHONY: all check preprocess expression network candidates causal features models evaluate crosstissue crossancestry figures provenance test clean-logs

all: provenance

# ---- phase 0: data loading, normalisation, batch correction, partition -----
preprocess:
	$(RSCRIPT) $(PREP)/02_data_loading.R
	$(RSCRIPT) $(PREP)/03_dataset_gene_overlap_venn.R
	$(RSCRIPT) $(PREP)/03_normalize_batch.R
	$(RSCRIPT) $(PREP)/03_normalize_batch_figure.R
	$(RSCRIPT) $(PREP)/04_apply_holdout.R
	$(RSCRIPT) $(PREP)/load_GSE15573.R
	@test -f results/tables/combined_cohort_summary.csv || \
		(echo "preprocess phase did not produce combined_cohort_summary.csv" && exit 1)

# ---- phase 1: differential expression, sensitivity, composition, interaction
expression: preprocess
	$(RSCRIPT) $(EXPR)/05_dge.R
	$(RSCRIPT) $(EXPR)/05_dge_figure.R
	$(RSCRIPT) $(EXPR)/05b_dge_sensitivity.R
	$(RSCRIPT) $(EXPR)/05c_deconvolution.R
	$(RSCRIPT) $(EXPR)/05d_interaction_report.R
	@test -f results/tables/DEG_summary.csv || \
		(echo "expression phase did not produce DEG_summary.csv" && exit 1)

# ---- phase 2: co-expression network (WGCNA) --------------------------------
network: expression
	$(RSCRIPT) $(NET)/06_WGCNA.R
	@test -f results/tables/WGCNA_03_disease_modules.csv || \
		(echo "network phase did not produce WGCNA_03_disease_modules.csv" && exit 1)

# ---- phase 3: candidate genes (disease module x sex-stratified DEG) --------
candidates: network
	$(RSCRIPT) $(CAND)/08_module_trait_RA_control.R
	$(RSCRIPT) $(CAND)/09_disease_module_deg_intersect.R
	$(RSCRIPT) $(CAND)/09b_disease_module_deg_venn.R
	@test -f results/tables/candidate_summary.csv || \
		(echo "candidates phase did not produce candidate_summary.csv" && exit 1)

# ---- phase 4: causal screening (MR, MHC sensitivity, colocalisation) -------
causal: candidates
	$(RSCRIPT) $(CAUSAL)/10_MR.R
	$(RSCRIPT) $(CAUSAL)/10c_MR_mhc_sensitivity.R
	$(RSCRIPT) $(CAUSAL)/10d_coloc_panel_genes.R
	$(RSCRIPT) $(CAUSAL)/10e_coloc_susie_mhc.R
	@test -f results/tables/MR_MHC_sensitivity_summary.csv || \
		(echo "causal phase did not produce MR_MHC_sensitivity_summary.csv" && exit 1)

# ---- phase 5: feature selection (LASSO / RF / SVM-RFE consensus) ----------
features: causal
	$(RSCRIPT) $(FEAT)/12_feature_selection.R
	$(RSCRIPT) $(FEAT)/12b_feature_selection_noMHC.R
	$(RSCRIPT) $(FEAT)/13_feature_selection_venn.R
	@test -f results/tables/mr_fs_summary.csv || \
		(echo "features phase did not produce mr_fs_summary.csv" && exit 1)

# ---- phase 6: model training (nested CV, elastic net, final panels) -------
models: features
	$(RSCRIPT) $(MODELS)/14_model_training_nested_cv.R
	$(RSCRIPT) $(MODELS)/15_model_training_elasticnet.R
	$(RSCRIPT) $(MODELS)/16_model_training_final_panel.R
	$(RSCRIPT) $(MODELS)/16b_model_training_final_panel_noMHC.R
	$(RSCRIPT) $(MODELS)/16c_model_training_nested_cv_transcriptomewide.R
	$(RSCRIPT) $(MODELS)/16d_nested_cv_reconciliation.R
	$(RSCRIPT) $(MODELS)/37_model_training_ml_algorithms.R
	@test -f results/tables/NESTED_CV_AUTHORITATIVE.csv || \
		(echo "models phase did not produce NESTED_CV_AUTHORITATIVE.csv" && exit 1)
	@test -f results/tables/ML_hyperparameter_tuning.csv || \
		(echo "models phase did not produce ML_hyperparameter_tuning.csv" && exit 1)

# ---- phase 7: blood evaluation (internal/external, per-gene ROC, nomogram) -
evaluate: models
	$(RSCRIPT) $(EVAL)/17_testing_blood_internal_external.R
	$(RSCRIPT) $(EVAL)/17b_testing_blood_celladjusted.R
	$(RSCRIPT) $(EVAL)/18_testing_blood_pergene_roc.R
	$(RSCRIPT) $(EVAL)/19_testing_blood_clinical_utility.R
	@test -f results/tables/mr_roc_panel_auc.csv || \
		(echo "evaluate phase did not produce mr_roc_panel_auc.csv" && exit 1)

# ---- phase 8: cross-tissue (synovium) validation ---------------------------
crosstissue: evaluate
	$(RSCRIPT) $(CROSST)/20_testing_synovium_external.R
	$(RSCRIPT) $(CROSST)/21_testing_synovium_roc.R
	$(RSCRIPT) $(CROSST)/22_crosstissue_biomarker_discovery.R
	$(RSCRIPT) $(CROSST)/23_crosstissue_roc_all_datasets.R
	$(RSCRIPT) $(CROSST)/24_crosstissue_pergene_auc.R
	$(RSCRIPT) $(CROSST)/25_crosstissue_pergene_roc.R
	$(RSCRIPT) $(CROSST)/38_testing_ml_algorithms_sametissue_crosstissue.R
	@test -f results/tables/crosstissue_panel_auc.csv || \
		(echo "crosstissue phase did not produce crosstissue_panel_auc.csv" && exit 1)
	@test -f results/tables/ML_performance_crosstissue.csv || \
		(echo "crosstissue phase did not produce ML_performance_crosstissue.csv" && exit 1)

# ---- phase 9: cross-ancestry evaluation ------------------------------------
crossancestry: evaluate
	$(RSCRIPT) $(CROSSA)/26_crossancestry_biomarker_mr.R
	$(RSCRIPT) $(CROSSA)/27_crossancestry_biomarker_figure.R
	@test -f results/tables/MR35_crossancestry_summary.csv || \
		(echo "crossancestry phase did not produce MR35_crossancestry_summary.csv" && exit 1)

# ---- phase 10: figures + functional enrichment -----------------------------
figures: crosstissue crossancestry
	$(RSCRIPT) $(GFIG)/28_figure_final_roc.R
	$(RSCRIPT) $(GFIG)/29_figure_nested_cv_roc.R
	$(RSCRIPT) $(GFIG)/30_figure_fig4A_expression_groups.R
	$(RSCRIPT) $(GFIG)/31_figure_fig4BC_logfc_heatmaps.R
	$(RSCRIPT) $(GFIG)/32_figure_fig4C_clustered_heatmap.R
	$(RSCRIPT) $(GFIG)/33_figure_fig4_composite.R
	$(RSCRIPT) $(GFIG)/34_pathway_enrichment.R
	$(RSCRIPT) $(GFIG)/35_figure_mr_diagnostics.R
	$(RSCRIPT) $(GFIG)/36_figure_pathways_by_sex.R
	$(RSCRIPT) $(GFIG)/39_figure_ml_algorithm_roc.R

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
