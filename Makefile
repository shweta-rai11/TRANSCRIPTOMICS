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
PREP := 01_DATA_PREPROCESSING/SCRIPT
EXPR := 02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT
NET := 03_WGCNA/SCRIPT
CAND := 04_CANDIDATE_GENE/SCRIPT
CAUSAL := 05_MENDELIAN_RANDOMISATION/SCRIPT
FEAT := 06_FEATURE_SELECTION/SCRIPT
MODELS := 07_MACHINE_LEARNING/SCRIPT
EVAL := 08_MODEL_EVALUATION/SCRIPT
CROSST := 09_CROSS_TISSUE_SYNOVIUM/SCRIPT
CROSSA := 10_CROSS_ANCESTRAL/SCRIPT
GFIG := 15_SUMMARY_FIGURES/SCRIPT
UTIL := UTILITIES/SCRIPT
PY := python3

.PHONY: all check preprocess expression network candidates causal features models evaluate crosstissue crossancestry figures provenance test clean-logs

all: provenance

# ---- phase 0: data loading, normalisation, batch correction, partition -----
preprocess:
	$(RSCRIPT) $(PREP)/02_data_loading.R
	$(RSCRIPT) $(PREP)/03_dataset_gene_overlap_venn.R
	$(RSCRIPT) $(PREP)/04_normalize_batch.R
	$(RSCRIPT) $(PREP)/05_normalize_batch_figure.R
	$(RSCRIPT) $(PREP)/11_apply_holdout.R
	$(RSCRIPT) $(PREP)/12_load_GSE15573.R
	@test -f results/tables/combined_cohort_summary.csv || \
		(echo "preprocess phase did not produce combined_cohort_summary.csv" && exit 1)

# ---- phase 1: differential expression, sensitivity, composition, interaction
expression: preprocess
	$(RSCRIPT) $(EXPR)/01_dge.R
	$(RSCRIPT) $(EXPR)/02_dge_figure.R
	$(RSCRIPT) $(EXPR)/03_dge_sensitivity.R
	$(RSCRIPT) $(EXPR)/04_deconvolution.R
	$(RSCRIPT) $(EXPR)/06_interaction_report.R
	@test -f results/tables/DEG_summary.csv || \
		(echo "expression phase did not produce DEG_summary.csv" && exit 1)

# ---- phase 2: co-expression network (WGCNA) --------------------------------
network: expression
	$(RSCRIPT) $(NET)/01_WGCNA.R
	@test -f results/tables/WGCNA_03_disease_modules.csv || \
		(echo "network phase did not produce WGCNA_03_disease_modules.csv" && exit 1)

# ---- phase 3: candidate genes (disease module x sex-stratified DEG) --------
candidates: network
	$(RSCRIPT) $(CAND)/01_module_trait_RA_control.R
	$(RSCRIPT) $(CAND)/02_disease_module_deg_intersect.R
	$(RSCRIPT) $(CAND)/03_disease_module_deg_venn.R
	@test -f results/tables/candidate_summary.csv || \
		(echo "candidates phase did not produce candidate_summary.csv" && exit 1)

# ---- phase 4: causal screening (MR, MHC sensitivity, colocalisation) -------
causal: candidates
	$(RSCRIPT) $(CAUSAL)/01_MR.R
	$(RSCRIPT) $(CAUSAL)/02_MR_mhc_sensitivity.R
	$(RSCRIPT) $(CAUSAL)/03_coloc_panel_genes.R
	$(RSCRIPT) $(CAUSAL)/04_coloc_susie_mhc.R
	@test -f results/tables/MR_MHC_sensitivity_summary.csv || \
		(echo "causal phase did not produce MR_MHC_sensitivity_summary.csv" && exit 1)

# ---- phase 5: feature selection (LASSO / RF / SVM-RFE consensus) ----------
features: causal
	$(RSCRIPT) $(FEAT)/01_feature_selection.R
	$(RSCRIPT) $(FEAT)/02_feature_selection_noMHC.R
	$(RSCRIPT) $(FEAT)/05_feature_selection_venn.R
	@test -f results/tables/mr_fs_summary.csv || \
		(echo "features phase did not produce mr_fs_summary.csv" && exit 1)

# ---- phase 6: model training (nested CV, elastic net, final panels) -------
# 02_model_training_elasticnet.R is deferred to the evaluate phase below: it
# reads results/tables/mr_roc_panel_auc.csv, which is only written by
# 08_MODEL_EVALUATION/01_testing_blood_internal_external.R (a later phase).
models: features
	$(RSCRIPT) $(MODELS)/01_model_training_nested_cv.R
	$(RSCRIPT) $(MODELS)/03_model_training_final_panel.R
	$(RSCRIPT) $(MODELS)/04_model_training_final_panel_noMHC.R
	$(RSCRIPT) $(MODELS)/05_model_training_nested_cv_transcriptomewide.R
	$(RSCRIPT) $(MODELS)/06_nested_cv_reconciliation.R
	$(RSCRIPT) $(MODELS)/07_model_training_ml_algorithms.R
	@test -f results/tables/NESTED_CV_AUTHORITATIVE.csv || \
		(echo "models phase did not produce NESTED_CV_AUTHORITATIVE.csv" && exit 1)
	@test -f results/tables/ML_hyperparameter_tuning.csv || \
		(echo "models phase did not produce ML_hyperparameter_tuning.csv" && exit 1)

# ---- phase 7: blood evaluation (internal/external, per-gene ROC, nomogram) -
evaluate: models
	$(RSCRIPT) $(EVAL)/01_testing_blood_internal_external.R
	$(RSCRIPT) $(MODELS)/02_model_training_elasticnet.R
	$(RSCRIPT) $(EVAL)/02_testing_blood_celladjusted.R
	$(RSCRIPT) $(EVAL)/03_testing_blood_pergene_roc.R
	$(RSCRIPT) $(EVAL)/04_testing_blood_clinical_utility.R
	@test -f results/tables/mr_roc_panel_auc.csv || \
		(echo "evaluate phase did not produce mr_roc_panel_auc.csv" && exit 1)

# ---- phase 8: cross-tissue (synovium) validation ---------------------------
crosstissue: evaluate
	$(RSCRIPT) $(CROSST)/01_testing_synovium_external.R
	$(RSCRIPT) $(CROSST)/02_testing_synovium_roc.R
	$(RSCRIPT) $(CROSST)/03_crosstissue_biomarker_discovery.R
	$(RSCRIPT) $(CROSST)/04_crosstissue_roc_all_datasets.R
	$(RSCRIPT) $(CROSST)/05_crosstissue_pergene_auc.R
	$(RSCRIPT) $(CROSST)/06_crosstissue_pergene_roc.R
	$(RSCRIPT) $(CROSST)/08_testing_ml_algorithms_sametissue_crosstissue.R
	@test -f results/tables/crosstissue_panel_auc.csv || \
		(echo "crosstissue phase did not produce crosstissue_panel_auc.csv" && exit 1)
	@test -f results/tables/ML_performance_crosstissue.csv || \
		(echo "crosstissue phase did not produce ML_performance_crosstissue.csv" && exit 1)

# ---- phase 9: cross-ancestry evaluation ------------------------------------
crossancestry: evaluate
	$(RSCRIPT) $(CROSSA)/01_crossancestry_biomarker_mr.R
	$(RSCRIPT) $(CROSSA)/02_crossancestry_biomarker_figure.R
	@test -f results/tables/MR35_crossancestry_summary.csv || \
		(echo "crossancestry phase did not produce MR35_crossancestry_summary.csv" && exit 1)

# ---- phase 10: figures + functional enrichment -----------------------------
figures: crosstissue crossancestry
	$(RSCRIPT) $(GFIG)/04_figure_final_roc.R
	$(RSCRIPT) $(GFIG)/05_figure_nested_cv_roc.R
	$(RSCRIPT) $(GFIG)/06_figure_fig4A_expression_groups.R
	$(RSCRIPT) $(GFIG)/07_figure_fig4BC_logfc_heatmaps.R
	$(RSCRIPT) $(GFIG)/08_figure_fig4C_clustered_heatmap.R
	$(RSCRIPT) $(GFIG)/09_figure_fig4_composite.R
	$(RSCRIPT) $(GFIG)/10_pathway_enrichment.R
	$(RSCRIPT) $(GFIG)/11_figure_mr_diagnostics.R
	$(RSCRIPT) $(GFIG)/12_figure_pathways_by_sex.R
	$(RSCRIPT) $(GFIG)/13_figure_ml_algorithm_roc.R

# ---- provenance: regenerate the verification report and run the test suite -
provenance: figures test
	$(PY) $(UTIL)/02_verify_results_numbers.py

test:
	$(PY) $(UTIL)/03_test_verify_results_numbers.py

# ---- sanity check only: confirms upstream data caches exist, runs nothing --
check:
	@test -f data/processed/combined_train.rds && echo "OK: data/processed/combined_train.rds present" \
		|| echo "MISSING: data/processed/combined_train.rds -- run 'make preprocess' first"
	@$(PY) $(UTIL)/02_verify_results_numbers.py --check

clean-logs:
	@echo "Not implemented deliberately: results/logs/*.log are provenance records" \
		"(cited by $(UTIL)/02_verify_results_numbers.py). Delete individually if you" \
		"really mean to invalidate a specific run's record."
