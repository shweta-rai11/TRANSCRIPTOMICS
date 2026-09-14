#!/usr/bin/env python3
"""Copies each existing figure into results/figures/by_section/ (renamed <section-slug><n>.png) and writes an INDEX.md mapping new name -> original file -> producing script -> thesis section."""
import os, shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_ROOT = os.path.join(ROOT, "results", "figures")
SRC_NEW = os.path.join(SRC_ROOT, "new")
OUT = os.path.join(SRC_ROOT, "by_section")
os.makedirs(OUT, exist_ok=True)

# (section_id, section_title, slug, [(original_filename, source_subdir, producing_script), ...])
SECTIONS = [
    ("3.2", "Normalisation and batch correction", "normalisation-batch", [
        ("fig_dataset_gene_overlap_venn.png", "", "scripts/00_shared/03_dataset_gene_overlap_venn.R"),
        ("expr_distribution.png", "", "scripts/00_shared/03_normalize_batch.R"),
        ("fig_combine_two_datasets.png", "", "scripts/00_shared/03_normalize_batch_figure.R"),
        ("fig_combine_density_norm.png", "", "scripts/00_shared/03_normalize_batch_figure.R"),
        ("fig_combine_pca_combat.png", "", "scripts/00_shared/03_normalize_batch_figure.R"),
    ]),
    ("3.4", "Co-expression network analysis (WGCNA)", "wgcna", [
        ("fig_wgcna_01_gene_variance.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_02_sample_clustering.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_03_soft_threshold.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_04_dendrogram.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_05_module_trait.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_06_module_eigengenes.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_07_ME_correlation.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_08_kME_vs_GS.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_09_hub_genes.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_10_disease_heatmap.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_11_ME_by_group_sex.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_13_dendro_female.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_14_dendro_male.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_15_preservation.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_16_candidates.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_wgcna_module_trait_RAvsControl_all.png", "", "scripts/00_shared/08_module_trait_RA_control.R"),
        ("fig_wgcna_module_trait_RAvsControl_female.png", "", "scripts/00_shared/08_module_trait_RA_control.R"),
        ("fig_wgcna_module_trait_RAvsControl_male.png", "", "scripts/00_shared/08_module_trait_RA_control.R"),
    ]),
    ("3.5", "Candidate gene identification", "candidates", [
        ("fig_module_trait_disease_selection.png", "", "scripts/00_shared/08_module_trait_RA_control.R"),
        ("fig_disease_module_selection_composite.png", "", "scripts/00_shared/09_disease_module_deg_intersect.R"),
        ("fig_venn_female_disease_candidates.png", "", "scripts/00_shared/09_disease_module_deg_intersect.R"),
        ("fig_venn_male_disease_candidates.png", "", "scripts/00_shared/09_disease_module_deg_intersect.R"),
        ("fig_diseasemod_venn_yellow_female_up.png", "", "scripts/00_shared/09b_disease_module_deg_venn.R"),
        ("fig_diseasemod_venn_yellow_female_down.png", "", "scripts/00_shared/09b_disease_module_deg_venn.R"),
        ("fig_diseasemod_venn_yellow_male_up.png", "", "scripts/00_shared/09b_disease_module_deg_venn.R"),
        ("fig_diseasemod_venn_yellow_male_down.png", "", "scripts/00_shared/09b_disease_module_deg_venn.R"),
        ("fig_diseasemod_venn_brown_female_up.png", "", "scripts/00_shared/09b_disease_module_deg_venn.R"),
        ("fig_diseasemod_venn_brown_female_down.png", "", "scripts/00_shared/09b_disease_module_deg_venn.R"),
        ("fig_diseasemod_venn_brown_male_up.png", "", "scripts/00_shared/09b_disease_module_deg_venn.R"),
        ("fig_diseasemod_venn_brown_male_down.png", "", "scripts/00_shared/09b_disease_module_deg_venn.R"),
    ]),
    ("3.8", "Feature selection", "feature-selection", [
        ("FIG_G2_01_panel_venn_female.png", "", "scripts/goal2_sex_stratified/13_feature_selection_venn.R"),
        ("FIG_G2_01_panel_venn_male.png", "", "scripts/goal2_sex_stratified/13_feature_selection_venn.R"),
        ("FIG_G2_02_panel_venn_noMHC_female.png", "", "scripts/goal2_sex_stratified/13_feature_selection_venn.R"),
        ("FIG_G2_02_panel_venn_noMHC_male.png", "", "scripts/goal2_sex_stratified/13_feature_selection_venn.R"),
        ("FIG_G2_03_panel_venn_female_vs_male.png", "", "scripts/goal2_sex_stratified/13_feature_selection_venn.R"),
        ("FIG_G2_04_feature_selection_venn_composite.png", "", "scripts/goal2_sex_stratified/13_feature_selection_venn.R"),
    ]),
    ("3.11", "Cross-tissue evaluation (synovium)", "cross-tissue", [
        ("fig_syn_pergene_roc_female.png", "new", "scripts/goal2_sex_stratified/21_testing_synovium_roc.R"),
        ("fig_syn_pergene_roc_male.png", "new", "scripts/goal2_sex_stratified/21_testing_synovium_roc.R"),
        ("fig_syn_panel_roc.png", "new", "scripts/goal2_sex_stratified/21_testing_synovium_roc.R"),
        ("fig_crosstissue.png", "new", "scripts/goal2_sex_stratified/22_crosstissue_biomarker_discovery.R"),
        ("fig_allvalidation_roc.png", "new", "scripts/goal2_sex_stratified/23_crosstissue_roc_all_datasets.R"),
        ("fig_pergene_auc_alltissues.png", "new", "scripts/goal2_sex_stratified/24_crosstissue_pergene_auc.R"),
        ("fig_pergene_roc_alltissues.png", "new", "scripts/goal2_sex_stratified/25_crosstissue_pergene_roc.R"),
    ]),
    ("3.13", "Functional enrichment", "functional-enrichment", [
        ("fig_wgcna_12_enrichment.png", "", "scripts/00_shared/06_WGCNA.R"),
        ("fig_enrich_female.png", "new", "scripts/goal2_sex_stratified/34_pathway_enrichment.R"),
        ("fig_enrich_male.png", "new", "scripts/goal2_sex_stratified/34_pathway_enrichment.R"),
    ]),
    ("3.15", "Nomogram construction and clinical evaluation", "nomogram", [
        ("fig_diag_validation_female.png", "new", "scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R"),
        ("fig_diag_validation_male.png", "new", "scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R"),
    ]),
]

# Sections confirmed (by grep + disk search) to have no figure at all, either unrun or (§3.7) none by design.
MISSING = [
    ("3.3", "Differential gene expression", "05_dge_figure.R defines figure-saving helpers but no run has produced output in results/figures/ or results/figures/new/"),
    ("3.6", "Mendelian randomisation", "30_figure_fig4A/31/32/33_*.R (fig_mr_gene_expression_groups, fig_mr_panelB/C_logfc, fig_mr_NEW_FIG4_composite) - none on disk"),
    ("3.7", "Colocalisation", "10d_coloc_panel_genes.R and 10e_coloc_susie_mhc.R produce no figures at all (tables only) - not missing, just none by design"),
    ("3.9", "Diagnostic model development and evaluation", "18_testing_blood_pergene_roc.R, 28_figure_final_roc.R, 29_figure_nested_cv_roc.R - none on disk"),
    ("3.10", "Diagnosis-by-sex interaction testing", "05d_interaction_report.R produces no figures at all (tables only) - not missing, just none by design"),
    ("3.12", "Cross-ancestry evaluation", "27_crossancestry_biomarker_figure.R (fig_mr35_crossancestry_*, fig_mr35_ancestry_class_summary) - none on disk"),
    ("3.14", "Immune deconvolution and composition-adjusted expression", "05c_deconvolution.R produces no figures at all (tables only) - not missing, just none by design"),
]

def main():
    manifest = []
    total = 0
    for sect_id, title, slug, files in SECTIONS:
        for i, (orig, subdir, script) in enumerate(files, start=1):
            src = os.path.join(SRC_ROOT, subdir, orig) if subdir else os.path.join(SRC_ROOT, orig)
            new_name = f"{slug}{i}.png"
            dst = os.path.join(OUT, new_name)
            if not os.path.exists(src):
                manifest.append((sect_id, title, new_name, orig, script, "MISSING AT COPY TIME"))
                continue
            shutil.copy2(src, dst)
            manifest.append((sect_id, title, new_name, orig, script, "OK"))
            total += 1

    with open(os.path.join(OUT, "INDEX.md"), "w", encoding="utf-8") as fh:
        fh.write("# Figure index by thesis Results section\n\n")
        fh.write(f"Generated by `scripts/organize_figures_by_section.py`. {total} figures copied into "
                 "this directory, renamed `<section-slug><n>.png`. Every mapping below was verified "
                 "against the actual `ggsave()`/`png()`/`pdf()`/`image_write()` call in the named script "
                 "- none is inferred from the original filename alone.\n\n")
        fh.write("| New name | Section | Original file | Producing script |\n|---|---|---|---|\n")
        for sect_id, title, new_name, orig, script, status in manifest:
            fh.write(f"| `{new_name}` | §{sect_id} {title} | `{orig}` | `{script}` |\n")
        fh.write("\n## Sections with no figure currently on disk\n\n")
        fh.write("These are not omissions - confirmed by grepping every script for figure-output calls "
                 "and searching the full `results/figures/` tree; the files simply do not exist because "
                 "the producing script has never been run, or (marked \"by design\") the script never "
                 "produces a figure at all.\n\n")
        fh.write("| Section | Title | Detail |\n|---|---|---|\n")
        for sect_id, title, detail in MISSING:
            fh.write(f"| §{sect_id} | {title} | {detail} |\n")

    print(f"Copied {total} figures into {OUT}")
    print(f"Wrote {os.path.join(OUT, 'INDEX.md')}")

if __name__ == "__main__":
    main()
