#!/usr/bin/env python3
"""Provenance harness for the Results chapter: re-derives each [R-xxx] claim from the per-layer TABLE/*.csv files (falling back to results/tables/ and results/logs/ when a live pipeline run is present) and writes results/RESULTS_PROVENANCE.tsv. Run with `python3 UTILITIES/SCRIPT/02_verify_results_numbers.py` (optionally a claim id, or --check)."""
import csv, glob, os, re, sys, datetime, statistics as st

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TAB  = os.path.join(ROOT, "results", "tables")
LOG  = os.path.join(ROOT, "results", "logs")
OUT  = os.path.join(ROOT, "results", "RESULTS_PROVENANCE.tsv")
STAMP = datetime.date.today().isoformat()

CLAIMS = []           # (id, description, value, source, derivation, section)
MISSING = []
CUR_SECTION = "unassigned"

def section(label):
    global CUR_SECTION
    CUR_SECTION = label

def locate(fn):
    """Return the path of a source CSV: results/tables/ first, else the lowest-numbered <layer>/TABLE/ copy."""
    p = os.path.join(TAB, fn)
    if os.path.exists(p):
        return p
    hits = sorted(glob.glob(os.path.join(ROOT, "[0-9][0-9]_*", "TABLE", fn)))
    return hits[0] if hits else None

def rows(fn):
    p = locate(fn)
    if p is None:
        MISSING.append(fn); return []
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))

def logtext(fn):
    p = os.path.join(LOG, fn)
    if not os.path.exists(p):
        MISSING.append(fn); return ""
    with open(p, encoding="utf-8", errors="replace") as fh:
        return fh.read()

def add(cid, desc, val, src, how):
    CLAIMS.append((cid, desc, str(val), src, how, CUR_SECTION))

def f(x):
    return float(x)

# =================================================================== §2.1 datasets
section("2.1 Datasets")
r = rows("combined_cohort_summary.csv")
if r:
    tot = sum(int(x["Freq"]) for x in r)
    add("R-001", "Training set size", tot, "combined_cohort_summary.csv", "sum of Freq")
    add("R-002", "Training RA", sum(int(x["Freq"]) for x in r if x["group"] == "RA"),
        "combined_cohort_summary.csv", "sum Freq where group==RA")
    add("R-003", "Training HC", sum(int(x["Freq"]) for x in r if x["group"] == "HC"),
        "combined_cohort_summary.csv", "sum Freq where group==HC")
    add("R-004", "Training female", sum(int(x["Freq"]) for x in r if x["sex"] == "F"),
        "combined_cohort_summary.csv", "sum Freq where sex==F")
    add("R-005", "Training male", sum(int(x["Freq"]) for x in r if x["sex"] == "M"),
        "combined_cohort_summary.csv", "sum Freq where sex==M")
    for ds in sorted(set(x["dataset"] for x in r)):
        add(f"R-006.{ds}", f"Training samples from {ds}",
            sum(int(x["Freq"]) for x in r if x["dataset"] == ds),
            "combined_cohort_summary.csv", f"sum Freq where dataset=={ds}")

h = rows("internal_val_holdout_meta.csv")
if h:
    add("R-010", "Hold-out size", len(h), "internal_val_holdout_meta.csv", "row count")
    add("R-011", "Hold-out female RA", sum(1 for x in h if x["group"]=="RA" and x["sex"]=="F"),
        "internal_val_holdout_meta.csv", "count group==RA & sex==F")
    add("R-012", "Hold-out female HC", sum(1 for x in h if x["group"]=="HC" and x["sex"]=="F"),
        "internal_val_holdout_meta.csv", "count group==HC & sex==F")
    add("R-013", "Hold-out male RA", sum(1 for x in h if x["group"]=="RA" and x["sex"]=="M"),
        "internal_val_holdout_meta.csv", "count group==RA & sex==M")
    add("R-014", "Hold-out male HC", sum(1 for x in h if x["group"]=="HC" and x["sex"]=="M"),
        "internal_val_holdout_meta.csv", "count group==HC & sex==M")
if r and h:
    add("R-015", "Combined discovery cohort (train + hold-out)",
        sum(int(x["Freq"]) for x in r) + len(h),
        "combined_cohort_summary.csv + internal_val_holdout_meta.csv", "sum of both")

# =================================================================== §2.2 pre-processing / normalisation
section("2.2 Data pre-processing, leakage-safe partition, normalisation and batch correction")
g = rows("dataset_gene_overlap.csv")
for x in g:
    add(f"R-020.{x['dataset']}", f"Genes after probe collapse, {x['dataset']}", x["n_genes"],
        "dataset_gene_overlap.csv", "column n_genes")

n = rows("normalization_diagnostics.csv")
for x in n:
    s = x["stage"]
    add(f"R-030.{s}", f"SD of per-sample medians, {s}", f"{f(x['median_sd']):.3g}",
        "normalization_diagnostics.csv", f"row stage=={s}, column median_sd")
    add(f"R-031.{s}", f"SD of per-sample IQRs, {s}", f"{f(x['iqr_sd']):.3g}",
        "normalization_diagnostics.csv", f"row stage=={s}, column iqr_sd")
t = logtext("05d_interaction.log")
if "batches 6" in t:
    add("R-035", "Number of ComBat batches resolved", 6, "logs/05d_interaction.log",
        "literal 'batches 6' in the cohort summary line")

# =================================================================== §2.3 differential expression
section("2.3 Differential gene expression analysis")
d = rows("DEG_summary.csv")
for x in d:
    c = x["comparison"]
    add(f"R-040.{c}", f"Significant DEGs, {c}", x["significant"], "DEG_summary.csv", f"row {c}, col significant")
    add(f"R-041.{c}", f"Up in RA, {c}", x["up_in_RA"], "DEG_summary.csv", f"row {c}, col up_in_RA")
    add(f"R-042.{c}", f"Down in RA, {c}", x["down_in_RA"], "DEG_summary.csv", f"row {c}, col down_in_RA")
    add(f"R-043.{c}", f"n / RA / HC, {c}", f"{int(x['RA'])+int(x['HC'])} / {x['RA']} / {x['HC']}",
        "DEG_summary.csv", f"row {c}, cols RA+HC")
    add(f"R-044.{c}", f"% of transcriptome significant, {c}",
        f"{100*int(x['significant'])/int(x['genes_tested']):.1f}%", "DEG_summary.csv",
        "significant / genes_tested")

aw = rows("DEG_array_weights.csv")
if aw:
    for comp in ["All", "Female", "Male"]:
        w = [f(x["weight"]) for x in aw if x["comparison"] == comp]
        if w:
            add(f"R-050.{comp}", f"Array weight median [min-max], {comp}",
                f"{st.median(w):.3f} [{min(w):.3f}-{max(w):.3f}]", "DEG_array_weights.csv",
                f"median/min/max of weight where comparison=={comp}")
            add(f"R-051.{comp}", f"Arrays with weight < 0.5, {comp}", sum(1 for i in w if i < 0.5),
                "DEG_array_weights.csv", f"count weight<0.5 where comparison=={comp}")
    for ds in sorted(set(x["dataset"] for x in aw)):
        w = [f(x["weight"]) for x in aw if x["comparison"] == "All" and x["dataset"] == ds]
        if w:
            add(f"R-052.{ds}", f"Median array weight, {ds} (All contrast)", f"{st.median(w):.3f}",
                "DEG_array_weights.csv", f"median weight where comparison==All & dataset=={ds}")

for sx in ["female", "male"]:
    s = rows(f"DEG_{sx}_significant.csv")
    if s:
        byfdr = sorted(s, key=lambda x: f(x["adj.P.Val"]))[:6]
        add(f"R-060.{sx}", f"Top DEGs by FDR, {sx}",
            "; ".join(f"{x['gene']} ({f(x['logFC']):+.2f}, FDR {f(x['adj.P.Val']):.1e})" for x in byfdr),
            f"DEG_{sx}_significant.csv", "6 smallest adj.P.Val")
        byfc = sorted(s, key=lambda x: -abs(f(x["logFC"])))[:6]
        add(f"R-061.{sx}", f"Largest |logFC| among significant, {sx}",
            "; ".join(f"{x['gene']} ({f(x['logFC']):+.2f})" for x in byfc),
            f"DEG_{sx}_significant.csv", "6 largest |logFC|")
        add(f"R-062.{sx}", f"Median |logFC| among significant, {sx}",
            f"{st.median([abs(f(x['logFC'])) for x in s]):.3f}",
            f"DEG_{sx}_significant.csv", "median of |logFC|")
        add(f"R-063.{sx}", f"Minimum attainable FDR, {sx}",
            f"{min(f(x['adj.P.Val']) for x in s):.1e}",
            f"DEG_{sx}_significant.csv", "min adj.P.Val")

fs_, ms_ = rows("DEG_female_significant.csv"), rows("DEG_male_significant.csv")
if fs_ and ms_:
    fg = {x["gene"]: f(x["logFC"]) for x in fs_}
    mg = {x["gene"]: f(x["logFC"]) for x in ms_}
    shared = set(fg) & set(mg)
    conc = sum(1 for g_ in shared if (fg[g_] > 0) == (mg[g_] > 0))
    add("R-070", "DEGs significant in both sexes", len(shared),
        "DEG_female_significant.csv + DEG_male_significant.csv", "set intersection on gene")
    add("R-071", "Female-list-only DEGs", len(set(fg) - set(mg)), "same", "set difference")
    add("R-072", "Male-list-only DEGs", len(set(mg) - set(fg)), "same", "set difference")
    add("R-073", "Union of the two DEG lists", len(set(fg) | set(mg)), "same", "set union")
    add("R-074", "Shared DEGs concordant in direction",
        f"{conc}/{len(shared)} ({100*conc/len(shared):.1f}%)", "same", "sign(logFC) agreement")

tr = rows("DEG_treat_sensitivity.csv")
for x in tr:
    add(f"R-080.{x['comparison']}.{x['fold_change']}",
        f"TREAT vs post-hoc filter, {x['comparison']} at {x['fold_change']}",
        f"{x['posthoc_filter_n']} -> {x['treat_n']}", "DEG_treat_sensitivity.csv",
        "cols posthoc_filter_n and treat_n")

cb = rows("DEG_sensitivity_combat_vs_batch.csv")
for x in cb:
    add(f"R-085.{x['comparison']}", f"ComBat-then-DE vs batch-in-model, {x['comparison']}",
        f"{x['combat_then_DE']} vs {x['batch_in_model']} (ratio {x['ratio']})",
        "DEG_sensitivity_combat_vs_batch.csv", "cols combat_then_DE, batch_in_model, ratio")

# =================================================================== §2.4 WGCNA
section("2.4 Co-expression network analysis (WGCNA)")
sft = rows("WGCNA_01_soft_threshold.csv")
for x in sft:
    if x["power"] in ("3", "6", "9", "10", "12", "14", "20"):
        add(f"R-090.b{x['power']}", f"Signed R2 / mean k at power {x['power']}",
            f"R2 {f(x['signed_R2']):.3f}, mean k {f(x['mean_k']):.1f}",
            "WGCNA_01_soft_threshold.csv", f"row power=={x['power']}")

mt = rows("WGCNA_02_module_trait.csv")
if mt:
    nm = [x for x in mt if x["module"] != "grey"]
    add("R-100", "Modules detected (excluding grey)", len(nm), "WGCNA_02_module_trait.csv", "row count, module!=grey")
    add("R-101", "Genes unassigned (grey)", [x["size"] for x in mt if x["module"] == "grey"][0],
        "WGCNA_02_module_trait.csv", "row module==grey, col size")
    add("R-102", "Module sizes",
        "; ".join(f"{x['module']} {x['size']}" for x in sorted(nm, key=lambda y: -int(y["size"]))),
        "WGCNA_02_module_trait.csv", "cols module,size sorted desc")
    for x in mt:
        if x["module"] in ("yellow", "brown", "turquoise", "green"):
            add(f"R-103.{x['module']}", f"cor(ME,RA) and p, {x['module']}",
                f"{f(x['cor_RA']):+.3f} (p {f(x['p_RA']):.2e})", "WGCNA_02_module_trait.csv",
                f"row {x['module']}, cols cor_RA,p_RA")

dm = rows("WGCNA_03_disease_modules.csv")
if dm:
    add("R-110", "Disease modules selected",
        "; ".join(f"{x['module']} (n={x['size']}, r={f(x['cor_RA']):+.3f}, p={f(x['p_RA']):.2e}, {x['direction']})" for x in dm),
        "WGCNA_03_disease_modules.csv", "all rows")
    add("R-111", "Disease-module gene background", sum(int(x["size"]) for x in dm),
        "WGCNA_03_disease_modules.csv", "sum of size")
    add("R-112", "cor(ME,Male) for disease modules",
        "; ".join(f"{x['module']} {f(x['cor_Male']):+.3f} (p {f(x['p_Male']):.3g})" for x in dm),
        "WGCNA_03_disease_modules.csv", "cols cor_Male,p_Male")

strata = rows("module_trait_RAvsControl_ALLSTRATA.csv")
for x in strata:
    if x["module"] in ("yellow", "brown"):
        add(f"R-115.{x['module']}.{x['stratum'].replace(' ','_')}",
            f"cor(ME,RA) {x['module']} in {x['stratum']} (n={x['n_samples']})",
            f"{f(x['cor_RA']):+.3f} (p {f(x['p_RA']):.2e})",
            "module_trait_RAvsControl_ALLSTRATA.csv", f"row module=={x['module']} & stratum=={x['stratum']}")

hb = rows("WGCNA_06_disease_module_hubs.csv")
if hb:
    hubs = [x for x in hb if x["is_hub"].strip().upper() in ("TRUE", "T", "1")]
    for m in ("yellow", "brown"):
        hm = [x for x in hubs if x["module"] == m]
        add(f"R-120.{m}", f"Hub genes in {m} module", len(hm),
            "WGCNA_06_disease_module_hubs.csv", f"count is_hub==TRUE & module=={m}")
        top = sorted(hm, key=lambda y: -f(y["connectivity"]))[:6]
        add(f"R-121.{m}", f"Top-connectivity hubs, {m}",
            "; ".join(f"{x['gene']} (kME {f(x['kME']):.3f}, GS {f(x['GS_RA']):+.3f}, k {f(x['connectivity']):.1f})" for x in top),
            "WGCNA_06_disease_module_hubs.csv", "6 largest connectivity among hubs")

pres = rows("WGCNA_10_module_preservation.csv")
if pres:
    gold = [f(x["Zsummary"]) for x in pres if x["module"] == "gold"]
    add("R-130", "Preservation Zsummary, gold (random benchmark)", f"{gold[0]:.2f}" if gold else "NA",
        "WGCNA_10_module_preservation.csv", "row module==gold")
    for m in ("yellow", "brown"):
        z = [f(x["Zsummary"]) for x in pres if x["module"] == m]
        if z:
            add(f"R-131.{m}", f"Preservation Zsummary, {m}", f"{z[0]:.2f}",
                "WGCNA_10_module_preservation.csv", f"row module=={m}")
            if gold:
                add(f"R-132.{m}", f"{m} Zsummary minus gold benchmark", f"{z[0]-gold[0]:+.2f}",
                    "WGCNA_10_module_preservation.csv", "Zsummary(module) - Zsummary(gold)")

w = logtext("06_WGCNA_run3.log")
if "outlier cut height : 69.32" in w:
    add("R-140", "Sample-outlier cut height", "69.32", "logs/06_WGCNA_run3.log", "literal in QC block")
if "sample outliers    : 10" in w:
    add("R-141", "Sample outliers removed", 10, "logs/06_WGCNA_run3.log", "literal in QC block")
if "final matrix   : 173 samples x 15763 genes" in w:
    add("R-142", "Network matrix after outlier removal", "173 samples x 15763 genes",
        "logs/06_WGCNA_run3.log", "literal 'final matrix' line")

# =================================================================== §2.5 candidates
section("2.5 Candidate gene identification (disease module ∩ sex-stratified DEG)")
cs = rows("candidate_summary.csv")
for x in cs:
    add(f"R-150.{x['sex']}", f"Candidates ({x['sex']}) = disease modules INTERSECT {x['sex']} DEGs",
        x["n_candidates"], "candidate_summary.csv", "col n_candidates")
    add(f"R-151.{x['sex']}", f"Candidate split by module, {x['sex']}", x["per_module"],
        "candidate_summary.csv", "col per_module")
    add(f"R-152.{x['sex']}", f"DEGs lost to variance filter, {x['sex']}",
        f"{x['DEGs_lost_to_variance_filter']} ({x['pct_DEG_lost']}%)",
        "candidate_summary.csv", "cols DEGs_lost_to_variance_filter, pct_DEG_lost")
ds_ = rows("diseasemod_DEG_direction_summary.csv")
if ds_:
    add("R-155", "Directional consistency (consistent / inconsistent overlaps)",
        "; ".join(f"{x['sex']}-{x['module']}: {x['n_consistent']}/{x['n_inconsistent']}" for x in ds_),
        "diseasemod_DEG_direction_summary.csv", "cols n_consistent, n_inconsistent")

cf, cm = rows("WGCNA_11_candidates_female.csv"), rows("WGCNA_11_candidates_male.csv")
if cf and cm:
    sf = set(x["gene"] for x in cf); sm = set(x["gene"] for x in cm)
    add("R-156", "Candidates shared by both sexes", len(sf & sm),
        "WGCNA_11_candidates_{female,male}.csv", "set intersection")
    add("R-157", "Female-list-only candidates", len(sf - sm), "same", "set difference")
    add("R-158", "Male-list-only candidates", len(sm - sf), "same", "set difference")
    add("R-159", "Union of candidates carried into MR", len(sf | sm), "same", "set union")

# =================================================================== §2.6 MR
section("2.6 Mendelian randomisation, incl. MHC sensitivity analysis")
mh = rows("MR_MHC_sensitivity_summary.csv")
for x in mh:
    s = x["sex"]
    add(f"R-160.{s}", f"Genes with an MR estimate, {s}", x["genes_tested_primary"],
        "MR_MHC_sensitivity_summary.csv", "col genes_tested_primary")
    add(f"R-161.{s}", f"MR-prioritised at within-stratum FDR<0.05, {s}", x["causal_FDR05_primary"],
        "MR_MHC_sensitivity_summary.csv", "col causal_FDR05_primary")
    add(f"R-162.{s}", f"Surviving MHC exclusion, {s}", x["causal_FDR05_noMHC"],
        "MR_MHC_sensitivity_summary.csv", "col causal_FDR05_noMHC")
    add(f"R-163.{s}", f"MHC verdicts (robust/untestable/MHC-dep/FDR-rank), {s}",
        f"{x['causal_robust']}/{x['causal_untestable_noMHC']}/{x['causal_MHC_dependent']}/{x['causal_FDR_rank_only']}",
        "MR_MHC_sensitivity_summary.csv", "cols causal_robust, causal_untestable_noMHC, causal_MHC_dependent, causal_FDR_rank_only")
    add(f"R-164.{s}", f"Genes tested after MHC exclusion, {s}", x["genes_tested_noMHC"],
        "MR_MHC_sensitivity_summary.csv", "col genes_tested_noMHC")
    add(f"R-165.{s}", f"Panel-gene MHC verdicts (robust/untestable/dep/rank), {s}",
        f"{x['panel_robust']}/{x['panel_untestable']}/{x['panel_MHC_dependent']}/{x['panel_FDR_rank_only']} of {x['panel_n']}",
        "MR_MHC_sensitivity_summary.csv", "panel_* columns")

for sx in ("female", "male"):
    fsx = rows(f"FS_input_{sx}.csv")
    if fsx:
        add(f"R-170.{sx}", f"FS input gene count, {sx}", len(fsx), f"FS_input_{sx}.csv", "row count")
        top = sorted(fsx, key=lambda x: f(x["MR_pval"]))[:6]
        add(f"R-171.{sx}", f"Strongest MR associations, {sx}",
            "; ".join(f"{x['gene']} OR {f(x['MR_OR']):.3f}, FDR {f(x['MR_FDR_stratum']):.1e}, {x['nSNP']} SNP, {x['method']}" for x in top),
            f"FS_input_{sx}.csv", "6 smallest MR_pval")
        add(f"R-172.{sx}", f"Single-instrument genes among FS input, {sx}",
            f"{sum(1 for x in fsx if x['nSNP']=='1')}/{len(fsx)}", f"FS_input_{sx}.csv", "count nSNP==1")
    nom = rows(f"FS_input_{sx}_noMHC.csv")
    if nom:
        add(f"R-173.{sx}", f"MHC-free FS input gene count, {sx}", len(nom),
            f"FS_input_{sx}_noMHC.csv", "row count")

ff, fm = rows("FS_input_female.csv"), rows("FS_input_male.csv")
if ff and fm:
    a, b = set(x["gene"] for x in ff), set(x["gene"] for x in fm)
    add("R-175", "MR-prioritised genes shared by both strata", len(a & b),
        "FS_input_{female,male}.csv", "set intersection")
    add("R-176", "Female-list-only MR-prioritised genes",
        f"{len(a-b)} ({', '.join(sorted(a-b))})", "same", "set difference")
    add("R-177", "Male-list-only MR-prioritised genes",
        f"{len(b-a)} ({', '.join(sorted(b-a))})", "same", "set difference")

# =================================================================== §2.7 colocalisation
section("2.7 Bayesian colocalisation (coloc.abf, coloc.susie)")
cl = rows("COLOC_summary.csv")
for x in cl:
    add(f"R-180.{x['set'].replace(' ','_')}", f"Colocalisation tally: {x['set']}",
        f"n={x['n_genes']}, colocalised={x['colocalised']}, prior-fragile={x['coloc_prior_fragile']}, "
        f"suggestive={x['suggestive']}, distinct={x['distinct_variants']}, MHC-unreliable={x['mhc_unreliable']}, "
        f"inconclusive={x['inconclusive']}", "COLOC_summary.csv", "all columns of that row")

cp = rows("COLOC_panel_genes.csv")
for x in cp:
    add(f"R-181.{x['gene']}", f"Coloc posteriors, {x['gene']}",
        f"PP.H3 {f(x['PP3']):.3f}, PP.H4 {f(x['PP4']):.3f}, PP.H4@p12=1e-6 {f(x['PP4_cons']):.3f}, "
        f"{x['nsnp_coloc']} SNPs, MHC={x['MHC_gene']}",
        "COLOC_panel_genes.csv", "cols PP3, PP4, PP4_cons, nsnp_coloc, MHC_gene")

su = rows("COLOC_SUSIE_mhc.csv")
if su:
    ok = [x for x in su if x["status"] == "ok"]
    add("R-190", "coloc.susie: MHC genes attempted / resolved", f"{len(su)} attempted, {len(ok)} resolved",
        "COLOC_SUSIE_mhc.csv", "count rows; status=='ok'")
    for x in ok:
        add(f"R-191.{x['gene']}", f"SuSiE result, {x['gene']}",
            f"{x['n_cs_eqtl']} eQTL sets, {x['n_cs_gwas']} RA sets, best PP.H4 {f(x['best_PP_H4']):.3f} - {x['susie_verdict']}",
            "COLOC_SUSIE_mhc.csv", "cols n_cs_eqtl, n_cs_gwas, best_PP_H4, susie_verdict")
    rng = sorted(set(int(x["n_cs_gwas"]) for x in su if x["n_cs_gwas"] not in ("", "0")))
    if rng:
        add("R-192", "Independent RA credible sets found in MHC regions",
            f"{min(rng)}-{max(rng)}", "COLOC_SUSIE_mhc.csv", "range of n_cs_gwas (non-zero)")

# =================================================================== §2.8 feature selection
section("2.8 Sex-stratified feature selection (LASSO / RF / SVM-RFE consensus)")
fsum = rows("mr_fs_summary.csv")
for x in fsum:
    add(f"R-200.{x['sex']}", f"Selector yields, {x['sex']} (primary)",
        f"LASSO {x['n_lasso']}, RF {x['n_rf']}, SVM-RFE {x['n_svmrfe']}, consensus {x['n_consensus']}",
        "mr_fs_summary.csv", "cols n_lasso, n_rf, n_svmrfe, n_consensus")
    add(f"R-201.{x['sex']}", f"Primary consensus panel, {x['sex']}", x["consensus_genes"],
        "mr_fs_summary.csv", "col consensus_genes")
    add(f"R-202.{x['sex']}", f"Tuned hyperparameters, {x['sex']}",
        f"lambda.min {x['tuned_lasso_lambda_min']}, mtry {x['tuned_rf_mtry']}, SVM cost {x['tuned_svm_cost']}",
        "mr_fs_summary.csv", "tuned_* columns")
fnom = rows("mr_fs_summary_noMHC.csv")
for x in fnom:
    add(f"R-205.{x['sex']}", f"MHC-free consensus panel, {x['sex']}",
        f"{x['n_consensus']} genes: {x['consensus_genes']}", "mr_fs_summary_noMHC.csv",
        "cols n_consensus, consensus_genes")

mem = rows("PANEL_primary_vs_noMHC_membership.csv")
if mem:
    for st_ in ("DROPPED", "NEW", "retained"):
        sel = [x for x in mem if x["status"].startswith(st_)]
        if sel:
            add(f"R-206.{st_}", f"Panel membership change: {st_}",
                "; ".join(f"{x['gene']}({x['sex'][0]})" for x in sel),
                "PANEL_primary_vs_noMHC_membership.csv", f"rows with status starting {st_}")

# =================================================================== §2.9 diagnostic model
section("2.9 Diagnostic model development and evaluation")
nc = rows("NESTED_CV_AUTHORITATIVE.csv")
for x in nc:
    add(f"R-210.{x['sex']}.{x['candidate_set']}.{x['selector']}",
        f"Nested CV AUC - {x['sex']}, {x['candidate_set']}, {x['selector']}",
        f"{x['nested_AUC']} ({x['CI_lo']}-{x['CI_hi']}), per-repeat SD {x['per_repeat_sd']}, "
        f"median genes {x['median_genes_used']}, n={x['n']}, recommended={x['recommended']}",
        "NESTED_CV_AUTHORITATIVE.csv", "one row of the 8-row grid")

perf = rows("PANEL_primary_vs_noMHC_performance.csv")
for x in perf:
    add(f"R-220.{x['sex']}.{x['panel']}.{x['dataset'].replace(' ','_')}",
        f"LOCKED-TRANSFER AUC - {x['sex']}, {x['panel']} panel, {x['dataset']}",
        f"{x['reported']}", "PANEL_primary_vs_noMHC_performance.csv", "col reported")

dl = rows("PANEL_primary_vs_noMHC_delong.csv")
for x in dl:
    add(f"R-225.{x['sex']}", f"DeLong, primary vs MHC-free nested CV, {x['sex']}",
        f"{x['AUC_primary']} vs {x['AUC_noMHC']}, p={x['delong_p']}",
        "PANEL_primary_vs_noMHC_delong.csv", "cols AUC_primary, AUC_noMHC, delong_p")

opt = rows("mr_nested_cv_summary.csv")
for x in opt:
    add(f"R-230.{x['sex']}", f"Apparent / flat / nested and optimism, {x['sex']}",
        f"apparent {x['apparent_AUC']}, flat {x['flat_CV_AUC']}, nested {x['nested_CV_AUC']}, "
        f"optimism {x['optimism_apparent_minus_nested']}",
        "mr_nested_cv_summary.csv", "cols apparent_AUC, flat_CV_AUC, nested_CV_AUC, optimism_*")

lrt = rows("PANEL_incremental_value_LRT.csv")
for x in lrt:
    add(f"R-240.{x['sex']}.{x['panel']}.{x['dataset'].replace(' ','_')}",
        f"WITHIN-DATASET RESAMPLED panel vs composition - {x['sex']}, {x['panel']}, {x['dataset']}",
        f"panel {x['AUC_panel']} vs composition {x['AUC_composition']} (delta {x['delta_panel_minus_composition']}), "
        f"LRT p={x['LRT_p_panel_beyond_composition']}, n={x['n']}, separation={x['separation_warning']}",
        "PANEL_incremental_value_LRT.csv", "cols AUC_panel, AUC_composition, delta_*, LRT_p_*")

pg = rows("mr_roc_pergene_auc.csv")
for x in pg:
    add(f"R-250.{x['sex']}.{x['dataset'].replace(' ','_')}.{x['gene']}",
        f"Per-gene AUC (train-fixed orientation) - {x['gene']}, {x['sex']}, {x['dataset']}",
        f"{x['AUC']} ({x['AUC_lo']}-{x['AUC_hi']}), concordant={x['concordant']}",
        "mr_roc_pergene_auc.csv", "cols AUC, AUC_lo, AUC_hi, concordant")

# =================================================================== §2.10 interaction
section("2.10 Diagnosis-by-sex interaction testing")
isum = rows("DEG_interaction_summary.csv")
for x in isum:
    add(f"R-260.{x['model'].replace(' ','_')}", f"Interaction genes at FDR<0.05, model = {x['model']}",
        x["n_interaction_FDR05"], "DEG_interaction_summary.csv", "col n_interaction_FDR05")
ip = rows("DEG_interaction_patterns.csv")
for x in ip:
    add(f"R-261.{x['pattern'].replace(' ','_')}", f"Interaction pattern count: {x['pattern']}",
        x["N"], "DEG_interaction_patterns.csv", "col N")
isig = rows("DEG_interaction_significant.csv")
if isig:
    add("R-262", "Sex-differential genes (primary model)", len(isig),
        "DEG_interaction_significant.csv", "row count")
    sur = [x for x in isig if x["interaction_FDR_celladj"] not in ("", "NA")
           and f(x["interaction_FDR_celladj"]) < 0.05]
    add("R-263", "Sex-differential genes surviving composition adjustment",
        f"{len(sur)}/{len(isig)} ({', '.join(x['gene'] for x in sur)})",
        "DEG_interaction_significant.csv", "count interaction_FDR_celladj < 0.05")
    top = sorted(isig, key=lambda x: f(x["interaction_FDR"]))[:6]
    add("R-264", "Leading sex-differential genes",
        "; ".join(f"{x['gene']} int {f(x['interaction_logFC']):+.3f} FDR {f(x['interaction_FDR']):.4f} "
                  f"[F {f(x['female_logFC']):+.3f} ({f(x['female_FDR']):.3g}), "
                  f"M {f(x['male_logFC']):+.3f} ({f(x['male_FDR']):.3g})]" for x in top),
        "DEG_interaction_significant.csv", "6 smallest interaction_FDR")
    p22 = [x for x in isig if x["gene"] == "PTPN22"]
    if p22:
        x = p22[0]
        add("R-265", "PTPN22 interaction",
            f"int {f(x['interaction_logFC']):+.3f} FDR {f(x['interaction_FDR']):.4f}; "
            f"F {f(x['female_logFC']):+.3f} (FDR {f(x['female_FDR']):.3g}); "
            f"M {f(x['male_logFC']):+.3f} (FDR {f(x['male_FDR']):.2e})",
            "DEG_interaction_significant.csv", "row gene==PTPN22")
    ylink = {"KDM5D","RPS4Y1","DDX3Y","USP9Y","UTY","EIF1AY","NLGN4Y","TXLNGY","ZFY","XIST"}
    hits = [x["gene"] for x in isig if x["gene"] in ylink]
    add("R-266", "Sex-chromosome-linked genes among the sex-differential set",
        ", ".join(hits) if hits else "none", "DEG_interaction_significant.csv",
        "membership test against a Y-linked/XIST gene list")
    panel = set()
    for t_ in (rows("mr_fs_summary.csv") + rows("mr_fs_summary_noMHC.csv")):
        panel |= {g.strip() for g in t_["consensus_genes"].split(";")}
    inpanel = [x["gene"] for x in isig if x["gene"] in panel]
    add("R-267", "Panel genes among the sex-differential set",
        ", ".join(inpanel) if inpanel else "NONE",
        "DEG_interaction_significant.csv + mr_fs_summary*.csv", "set intersection")

# =================================================================== §2.11 cross-tissue
section("2.11 Cross-tissue evaluation (synovium)")
for sx in ("female", "male"):
    v = rows(f"val_synovium_pergene_{sx}.csv")
    for x in v:
        if x.get("present") == "TRUE":
            add(f"R-270.{sx}.{x['gene']}", f"Synovium per-gene, {x['gene']} ({sx})",
                f"log2FC {f(x['syn_log2FC']):+.3f} (FDR {f(x['syn_adjP']):.2e}), concordant={x['concordant']}, "
                f"AUC best-direction {f(x['auc_sex_bestdir']):.3f}, AUC train-oriented {f(x['auc_sex_trainorient']):.3f}",
                f"val_synovium_pergene_{sx}.csv", "cols syn_log2FC, syn_adjP, concordant, auc_sex_bestdir, auc_sex_trainorient")
    if v:
        conc = sum(1 for x in v if x.get("concordant") == "TRUE")
        add(f"R-271.{sx}", f"Panel genes concordant in direction blood->synovium, {sx}",
            f"{conc}/{len([x for x in v if x.get('present')=='TRUE'])}",
            f"val_synovium_pergene_{sx}.csv", "count concordant==TRUE")

ct = rows("crosstissue_panel_auc.csv")
for x in ct:
    add(f"R-275.{x['sex']}", f"WITHIN-SYNOVIUM refit panel AUC, {x['sex']}",
        f"apparent {x['apparent_AUC']}, 10-fold CV {x['CV_AUC']} ({x['CV_lo']}-{x['CV_hi']}), "
        f"n={x['n']} (RA {x['n_RA']}, Normal {x['n_Normal']})",
        "crosstissue_panel_auc.csv", "cols apparent_AUC, CV_AUC, CV_lo, CV_hi, n, n_RA, n_Normal")

# =================================================================== §2.12 cross-ancestry
section("2.12 Cross-ancestry evaluation")
ca = rows("MR35_crossancestry_summary.csv")
for x in ca:
    add(f"R-280.{x['sex']}", f"Cross-ancestry classification, {x['sex']}",
        f"tested {x['n_causal']}: EUR-replicated {x['n_replicated_EUR']}, EAS-transferable {x['n_transferable_EAS']}, "
        f"EAS-untestable {x['n_untestable_EAS']}, EUR-only {x['n_EUR_only']}",
        "MR35_crossancestry_summary.csv", "all columns of that row")

# =================================================================== §2.13 enrichment
section("2.13 Functional enrichment")
for lab, fn in (("female", "dge_enriched_terms_female.csv"),
                ("male", "dge_enriched_terms_male.csv"),
                ("all", "dge_enriched_terms_all.csv")):
    e = rows(fn)
    if e:
        add(f"R-290.{lab}", f"Enriched terms, {lab} DEGs", len(e), fn, "row count")

gs = rows("gsea_kegg.csv")
if gs:
    for grp in sorted(set(x["group"] for x in gs)):
        sub = [x for x in gs if x["group"] == grp]
        neg = sum(1 for x in sub if f(x["NES"]) < 0)
        add(f"R-295.{grp}", f"GSEA KEGG, {grp}",
            f"{len(sub)} significant, {neg} suppressed, {len(sub)-neg} activated",
            "gsea_kegg.csv", f"rows group=={grp}, sign of NES")
        low = sorted(sub, key=lambda x: f(x["NES"]))[:4]
        add(f"R-296.{grp}", f"Most suppressed KEGG sets, {grp}",
            "; ".join(f"{x['Description']} (NES {f(x['NES']):+.2f}, FDR {f(x['p.adjust']):.1e})" for x in low),
            "gsea_kegg.csv", "4 most negative NES")
        pos = [x for x in sub if f(x["NES"]) > 0]
        add(f"R-297.{grp}", f"Activated KEGG sets, {grp}",
            "; ".join(f"{x['Description']} (NES {f(x['NES']):+.2f})" for x in pos) or "none",
            "gsea_kegg.csv", "rows with NES>0")

for lab, fn in (("female", "enrich_GO_BP_female.csv"), ("male", "enrich_GO_BP_male.csv")):
    e = rows(fn)
    if e:
        add(f"R-305.{lab}", f"GO:BP terms enriched among MR-prioritised genes, {lab}", len(e), fn, "row count")
        top = e[:5]
        add(f"R-306.{lab}", f"Top GO:BP terms, MR-prioritised {lab} genes",
            "; ".join(f"{x['Description']} ({x['GeneRatio']}, FDR {f(x['p.adjust']):.1e})" for x in top),
            fn, "first 5 rows as written (sorted by p.adjust)")
for lab, fn in (("female", "enrich_KEGG_female.csv"), ("male", "enrich_KEGG_male.csv")):
    e = rows(fn)
    if e:
        add(f"R-307.{lab}", f"KEGG pathways enriched among MR-prioritised genes, {lab}", len(e), fn, "row count")
        top = e[:5]
        add(f"R-308.{lab}", f"Top KEGG pathways, MR-prioritised {lab} genes",
            "; ".join(f"{x['Description']} ({x['GeneRatio']}, FDR {f(x['p.adjust']):.1e})" for x in top),
            fn, "first 5 rows as written (sorted by p.adjust)")

go = rows("WGCNA_08_disease_GO.csv")
kg = rows("WGCNA_09_disease_KEGG.csv")
if go: add("R-300", "GO terms enriched in disease-module genes", len(go), "WGCNA_08_disease_GO.csv", "row count")
if kg:
    add("R-301", "KEGG pathways enriched in disease-module genes", len(kg), "WGCNA_09_disease_KEGG.csv", "row count")
    add("R-302", "Top disease-module KEGG pathways",
        "; ".join(f"{x['Description']} ({x['GeneRatio']}, FDR {f(x['p.adjust']):.1e})" for x in kg[:6]),
        "WGCNA_09_disease_KEGG.csv", "first 6 rows as written")
if go:
    add("R-303", "Top disease-module GO terms",
        "; ".join(f"{x['Description']} ({x['GeneRatio']}, FDR {f(x['p.adjust']):.1e})" for x in go[:5]),
        "WGCNA_08_disease_GO.csv", "first 5 rows as written")

# =================================================================== §2.14 deconvolution
section("2.14 Immune deconvolution and composition-adjusted expression")
cg = rows("CELL_fraction_group_tests.csv")
sig = [x for x in cg if x["FDR"] not in ("", "NA") and f(x["FDR"]) < 0.05]
for x in sig:
    add(f"R-310.{x['sex']}.{x['cell']}", f"Cell fraction RA vs HC, {x['sex']} {x['cell']}",
        f"HC {x['mean_HC']} -> RA {x['mean_RA']} (diff {x['diff']}, FDR {f(x['FDR']):.2e})",
        "CELL_fraction_group_tests.csv", "rows with FDR<0.05")
add("R-311", "Cell subsets differing by disease status at FDR<0.05",
    f"{sum(1 for x in sig if x['sex']=='Female')} female, {sum(1 for x in sig if x['sex']=='Male')} male",
    "CELL_fraction_group_tests.csv", "count FDR<0.05 by sex")

ca_ = rows("DEG_celladjusted_summary.csv")
for x in ca_:
    add(f"R-315.{x['sex']}", f"Composition-adjusted DEG count, {x['sex']}",
        f"{x['n_DEG_unadjusted']} -> {x['n_DEG_adjusted']} ({x['pct_retained']}% retained), "
        f"Spearman rho {x['logFC_spearman']}",
        "DEG_celladjusted_summary.csv", "cols n_DEG_unadjusted, n_DEG_adjusted, pct_retained, logFC_spearman")

pa = rows("CELL_panel_gene_adjustment.csv")
for x in pa:
    add(f"R-320.{x['sex']}.{x['gene']}", f"Panel gene composition robustness, {x['sex']} {x['gene']}",
        f"logFC {x['logFC_unadj']} -> {x['logFC_adj']} ({x['pct_logFC_retained']}% retained), {x['status']}",
        "CELL_panel_gene_adjustment.csv", "cols logFC_unadj, logFC_adj, pct_logFC_retained, status")
if pa:
    add("R-321", "Panel gene-sex pairs retained after composition adjustment",
        f"{sum(1 for x in pa if x['status']=='retained')}/{len(pa)}",
        "CELL_panel_gene_adjustment.csv", "count status=='retained'")

# =================================================================== §2.15 nomogram / clinical utility
section("2.15 Nomogram construction and clinical evaluation")
# NOTE: this section's source files did not exist previously; 19_testing_blood_clinical_utility.R was run once to produce them.
diagtxt = logtext("19_testing_blood_clinical_utility.log")

def _clinical_block(sx):
    return re.search(
        r"n=(\d+)\s+Mean absolute error=([\d.]+)\s+Mean squared error=([\d.]+)"
        r"(?:(?!wrote fig_diag_validation).)*?"
        r"0\.9 Quantile of absolute error=([\d.]+)"
        r"(?:(?!wrote fig_diag_validation).)*?"
        r"wrote fig_diag_validation_" + sx + r"\.png\s+\(panel (\d+) genes, n=(\d+), RA=(\d+)\)",
        diagtxt, re.S)

for sx in ("female", "male"):
    m = _clinical_block(sx)
    if m:
        n, mae, mse, q90, ngenes, n2, ra = m.groups()
        add(f"R-330.{sx}", f"Nomogram panel fit on the training cohort, {sx}",
            f"{ngenes} genes, n={n2} (RA={ra}, HC={int(n2)-int(ra)})",
            "logs/19_testing_blood_clinical_utility.log",
            f"literal 'wrote fig_diag_validation_{sx}.png' line")
        add(f"R-331.{sx}", f"Calibration error (200-rep bootstrap bias-correction), {sx}",
            f"MAE {mae}, MSE {mse}, 90th-percentile absolute error {q90}",
            "logs/19_testing_blood_clinical_utility.log",
            "literal rms::calibrate() console block preceding the wrote-line")

    dca = rows(f"diag_dca_{sx}.csv")
    if dca:
        for thr_want in ("0.10", "0.20", "0.30", "0.50"):
            row = min(dca, key=lambda x: abs(f(x["threshold"]) - f(thr_want)))
            add(f"R-332.{sx}.{thr_want}", f"Net benefit at threshold {thr_want}, {sx}",
                f"panel {f(row['NB_panel']):.4f} vs treat-all {f(row['NB_all']):.4f} vs treat-none 0.0000",
                f"diag_dca_{sx}.csv", f"row nearest threshold=={thr_want}")
        beats_all = [f(x["threshold"]) for x in dca if f(x["NB_panel"]) > f(x["NB_all"])]
        positive  = [f(x["threshold"]) for x in dca if f(x["NB_panel"]) > 0]
        add(f"R-333.{sx}", f"Threshold range where the panel beats treat-all, {sx}",
            f"{min(beats_all):.2f}-{max(beats_all):.2f}" if beats_all else "never",
            f"diag_dca_{sx}.csv", "thresholds where NB_panel > NB_all")
        add(f"R-334.{sx}", f"Threshold range with positive net benefit, {sx}",
            f"{min(positive):.2f}-{max(positive):.2f}" if positive else "never",
            f"diag_dca_{sx}.csv", "thresholds where NB_panel > 0")

# =================================================================== provenance: csv/log -> producing R script
CSV_SCRIPT = {
    "combined_cohort_summary.csv": "01_DATA_PREPROCESSING/SCRIPT/04_normalize_batch.R",
    "internal_val_holdout_meta.csv": "01_DATA_PREPROCESSING/SCRIPT/04_normalize_batch.R",
    "normalization_diagnostics.csv": "01_DATA_PREPROCESSING/SCRIPT/04_normalize_batch.R",
    "dataset_gene_overlap.csv": "01_DATA_PREPROCESSING/SCRIPT/03_dataset_gene_overlap_venn.R",
    "05d_interaction.log": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/06_interaction_report.R",
    "DEG_summary.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/01_dge.R",
    "DEG_array_weights.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/01_dge.R",
    "DEG_female_significant.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/01_dge.R",
    "DEG_male_significant.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/01_dge.R",
    "DEG_treat_sensitivity.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/01_dge.R",
    "DEG_sensitivity_combat_vs_batch.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/03_dge_sensitivity.R",
    "dge_enriched_terms_female.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/01_dge.R",
    "dge_enriched_terms_male.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/01_dge.R",
    "dge_enriched_terms_all.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/01_dge.R",
    "gsea_kegg.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/01_dge.R",
    "WGCNA_01_soft_threshold.csv": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "WGCNA_02_module_trait.csv": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "WGCNA_03_disease_modules.csv": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "WGCNA_05_gene_module_assignment.csv": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "WGCNA_06_disease_module_hubs.csv": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "WGCNA_08_disease_GO.csv": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "WGCNA_09_disease_KEGG.csv": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "WGCNA_10_module_preservation.csv": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "WGCNA_11_candidates_female.csv": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "WGCNA_11_candidates_male.csv": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "WGCNA_12_candidate_summary.csv": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "06_WGCNA_run3.log": "03_WGCNA/SCRIPT/01_WGCNA.R",
    "module_trait_RAvsControl_ALLSTRATA.csv": "04_CANDIDATE_GENE/SCRIPT/01_module_trait_RA_control.R",
    "module_trait_RAvsControl_all.csv": "04_CANDIDATE_GENE/SCRIPT/01_module_trait_RA_control.R",
    "module_trait_RAvsControl_female.csv": "04_CANDIDATE_GENE/SCRIPT/01_module_trait_RA_control.R",
    "module_trait_RAvsControl_male.csv": "04_CANDIDATE_GENE/SCRIPT/01_module_trait_RA_control.R",
    "candidate_summary.csv": "04_CANDIDATE_GENE/SCRIPT/02_disease_module_deg_intersect.R",
    "candidates_female_disease.csv": "04_CANDIDATE_GENE/SCRIPT/02_disease_module_deg_intersect.R",
    "candidates_male_disease.csv": "04_CANDIDATE_GENE/SCRIPT/02_disease_module_deg_intersect.R",
    "disease_module_selection.csv": "04_CANDIDATE_GENE/SCRIPT/02_disease_module_deg_intersect.R",
    "diseasemod_DEG_direction_summary.csv": "04_CANDIDATE_GENE/SCRIPT/03_disease_module_deg_venn.R",
    "diseasemod_DEG_intersection_female.csv": "04_CANDIDATE_GENE/SCRIPT/03_disease_module_deg_venn.R",
    "diseasemod_DEG_intersection_male.csv": "04_CANDIDATE_GENE/SCRIPT/03_disease_module_deg_venn.R",
    "MR_MHC_sensitivity_summary.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/02_MR_mhc_sensitivity.R",
    "MR_MHC_sensitivity_female.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/02_MR_mhc_sensitivity.R",
    "MR_MHC_sensitivity_male.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/02_MR_mhc_sensitivity.R",
    "MR_MHC_sensitivity_panel_fate.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/02_MR_mhc_sensitivity.R",
    "FS_input_female.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "FS_input_male.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "FS_input_female_noMHC.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/02_MR_mhc_sensitivity.R",
    "FS_input_male_noMHC.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/02_MR_mhc_sensitivity.R",
    "MR_female_TABLE1_instruments.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "MR_male_TABLE1_instruments.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "MR_female_TABLE2_results_allmethods.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "MR_male_TABLE2_results_allmethods.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "MR_female_TABLE3_pleiotropy.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "MR_male_TABLE3_pleiotropy.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "MR_female_TABLE4_heterogeneity.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "MR_male_TABLE4_heterogeneity.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "MR_female_primary_okada.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "MR_male_primary_okada.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/01_MR.R",
    "COLOC_results.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/03_coloc_panel_genes.R",
    "COLOC_panel_genes.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/03_coloc_panel_genes.R",
    "COLOC_summary.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/03_coloc_panel_genes.R",
    "COLOC_SUSIE_mhc.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/04_coloc_susie_mhc.R",
    "COLOC_combined_abf_susie.csv": "05_MENDELIAN_RANDOMISATION/SCRIPT/04_coloc_susie_mhc.R",
    "mr_fs_summary.csv": "06_FEATURE_SELECTION/SCRIPT/01_feature_selection.R",
    "mr_fs_consensus_female.csv": "06_FEATURE_SELECTION/SCRIPT/01_feature_selection.R",
    "mr_fs_consensus_male.csv": "06_FEATURE_SELECTION/SCRIPT/01_feature_selection.R",
    "mr_fs_selected_bymethod_female.csv": "06_FEATURE_SELECTION/SCRIPT/01_feature_selection.R",
    "mr_fs_selected_bymethod_male.csv": "06_FEATURE_SELECTION/SCRIPT/01_feature_selection.R",
    "mr_fs_summary_noMHC.csv": "06_FEATURE_SELECTION/SCRIPT/02_feature_selection_noMHC.R",
    "PANEL_primary_vs_noMHC_membership.csv": "06_FEATURE_SELECTION/SCRIPT/02_feature_selection_noMHC.R",
    "FS_venn_membership.csv": "06_FEATURE_SELECTION/SCRIPT/05_feature_selection_venn.R",
    "mr_nested_cv_stability_female.csv": "07_MACHINE_LEARNING/SCRIPT/01_model_training_nested_cv.R",
    "mr_nested_cv_stability_male.csv": "07_MACHINE_LEARNING/SCRIPT/01_model_training_nested_cv.R",
    "mr_nested_cv_summary.csv": "07_MACHINE_LEARNING/SCRIPT/01_model_training_nested_cv.R",
    "mr_elasticnet_summary.csv": "07_MACHINE_LEARNING/SCRIPT/02_model_training_elasticnet.R",
    "mr_final_coefs_bysex.csv": "07_MACHINE_LEARNING/SCRIPT/03_model_training_final_panel.R",
    "mr_final_panel_summary.csv": "07_MACHINE_LEARNING/SCRIPT/03_model_training_final_panel.R",
    "PANEL_primary_vs_noMHC_nestedcv.csv": "07_MACHINE_LEARNING/SCRIPT/04_model_training_final_panel_noMHC.R",
    "PANEL_primary_vs_noMHC_delong.csv": "07_MACHINE_LEARNING/SCRIPT/04_model_training_final_panel_noMHC.R",
    "PANEL_primary_vs_noMHC_performance.csv": "07_MACHINE_LEARNING/SCRIPT/04_model_training_final_panel_noMHC.R",
    "nested_cv_stability_female.csv": "07_MACHINE_LEARNING/SCRIPT/05_model_training_nested_cv_transcriptomewide.R",
    "nested_cv_stability_male.csv": "07_MACHINE_LEARNING/SCRIPT/05_model_training_nested_cv_transcriptomewide.R",
    "nested_cv_summary.csv": "07_MACHINE_LEARNING/SCRIPT/05_model_training_nested_cv_transcriptomewide.R",
    "NESTED_CV_AUTHORITATIVE.csv": "07_MACHINE_LEARNING/SCRIPT/06_nested_cv_reconciliation.R",
    "NESTED_CV_legacy_reconciliation.csv": "07_MACHINE_LEARNING/SCRIPT/06_nested_cv_reconciliation.R",
    "mr_roc_panel_auc.csv": "08_MODEL_EVALUATION/SCRIPT/01_testing_blood_internal_external.R",
    "mr_roc_pergene_auc.csv": "08_MODEL_EVALUATION/SCRIPT/01_testing_blood_internal_external.R",
    "PANEL_auc_celladjusted.csv": "08_MODEL_EVALUATION/SCRIPT/02_testing_blood_celladjusted.R",
    "PANEL_incremental_value_LRT.csv": "08_MODEL_EVALUATION/SCRIPT/02_testing_blood_celladjusted.R",
    "PANEL_auc_celladjusted_summary.csv": "08_MODEL_EVALUATION/SCRIPT/02_testing_blood_celladjusted.R",
    "mr_pergene_train_auc.csv": "08_MODEL_EVALUATION/SCRIPT/03_testing_blood_pergene_roc.R",
    "diag_dca_female.csv": "08_MODEL_EVALUATION/SCRIPT/04_testing_blood_clinical_utility.R",
    "diag_dca_male.csv": "08_MODEL_EVALUATION/SCRIPT/04_testing_blood_clinical_utility.R",
    "19_testing_blood_clinical_utility.log": "08_MODEL_EVALUATION/SCRIPT/04_testing_blood_clinical_utility.R",
    "val_synovium_pergene_female.csv": "09_CROSS_TISSUE_SYNOVIUM/SCRIPT/01_testing_synovium_external.R",
    "val_synovium_pergene_male.csv": "09_CROSS_TISSUE_SYNOVIUM/SCRIPT/01_testing_synovium_external.R",
    "crosstissue_panel_auc.csv": "09_CROSS_TISSUE_SYNOVIUM/SCRIPT/03_crosstissue_biomarker_discovery.R",
    "pergene_auc_alltissues.csv": "09_CROSS_TISSUE_SYNOVIUM/SCRIPT/05_crosstissue_pergene_auc.R",
    "MR35_crossancestry_female.csv": "10_CROSS_ANCESTRAL/SCRIPT/01_crossancestry_biomarker_mr.R",
    "MR35_crossancestry_male.csv": "10_CROSS_ANCESTRAL/SCRIPT/01_crossancestry_biomarker_mr.R",
    "MR35_crossancestry_summary.csv": "10_CROSS_ANCESTRAL/SCRIPT/01_crossancestry_biomarker_mr.R",
    "MR35_instrument_transferability_female.csv": "10_CROSS_ANCESTRAL/SCRIPT/01_crossancestry_biomarker_mr.R",
    "MR35_instrument_transferability_male.csv": "10_CROSS_ANCESTRAL/SCRIPT/01_crossancestry_biomarker_mr.R",
    "CELL_fractions_train.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/04_deconvolution.R",
    "CELL_fractions_holdout.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/04_deconvolution.R",
    "CELL_fractions_external.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/04_deconvolution.R",
    "CELL_mcpcounter_train.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/04_deconvolution.R",
    "CELL_mcpcounter_external.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/04_deconvolution.R",
    "CELL_fraction_group_tests.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/04_deconvolution.R",
    "CELL_composition_pca.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/04_deconvolution.R",
    "CELL_panel_gene_adjustment.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/04_deconvolution.R",
    "DEG_celladjusted_female.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/04_deconvolution.R",
    "DEG_celladjusted_male.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/04_deconvolution.R",
    "DEG_celladjusted_summary.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/04_deconvolution.R",
    "DEG_interaction_summary.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/03_dge_sensitivity.R",
    "DEG_interaction_top.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/03_dge_sensitivity.R",
    "DEG_interaction_model_comparison.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/06_interaction_report.R",
    "DEG_interaction_full.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/06_interaction_report.R",
    "DEG_interaction_significant.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/06_interaction_report.R",
    "DEG_interaction_patterns.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/06_interaction_report.R",
    "DEG_interaction_enrichment.csv": "02_DIFFERENTIAL_GENE_EXPRESSION/SCRIPT/06_interaction_report.R",
    "enrich_GO_BP_female.csv": "12_FUNCTIONAL_ENRICHMENT/SCRIPT/10_pathway_enrichment.R",
    "enrich_GO_BP_male.csv": "12_FUNCTIONAL_ENRICHMENT/SCRIPT/10_pathway_enrichment.R",
    "enrich_KEGG_female.csv": "12_FUNCTIONAL_ENRICHMENT/SCRIPT/10_pathway_enrichment.R",
    "enrich_KEGG_male.csv": "12_FUNCTIONAL_ENRICHMENT/SCRIPT/10_pathway_enrichment.R",
}

def scripts_for(src):
    """Resolve a claim's `source` string (possibly ' + '-joined, brace-expanded,
    or wildcarded) to the list of R scripts that write the underlying file(s)."""
    names = []
    for chunk in src.split(" + "):
        chunk = chunk.strip()
        m = re.search(r"\{([^}]+)\}", chunk)
        if m:
            for opt in m.group(1).split(","):
                names.append(chunk.replace("{" + m.group(1) + "}", opt))
        elif chunk.endswith("*.csv"):
            prefix = chunk[:-5]
            names.extend(k for k in CSV_SCRIPT if k.startswith(prefix))
        else:
            names.append(chunk)
    out = []
    for nm in names:
        nm = nm[len("logs/"):] if nm.startswith("logs/") else nm
        s = CSV_SCRIPT.get(nm)
        if s and s not in out:
            out.append(s)
    return out or [f"(script not mapped for {src})"]

def display_source(src):
    parts = [p.strip() for p in src.split(" + ")]
    disp = []
    for p in parts:
        if p.startswith("logs/"):
            disp.append(f"results/{p}"); continue
        loc = locate(p)
        disp.append(os.path.relpath(loc, ROOT) if loc else f"results/tables/{p}")
    return " + ".join(disp)

def resolve_same(claims):
    """Replace src=='same' with the nearest preceding real source in append order."""
    out, last = [], None
    for cid, desc, val, src, how, sect in claims:
        if src == "same":
            src = last
        else:
            last = src
        out.append((cid, desc, val, src, how, sect))
    return out

# ---------------------------------------------------------------- output
def main():
    args = [a for a in sys.argv[1:]]
    if "--check" in args:
        if MISSING:
            print("MISSING SOURCE FILES:", ", ".join(sorted(set(MISSING)))); sys.exit(1)
        print(f"OK - {len(CLAIMS)} claims derived, no missing sources."); sys.exit(0)

    resolved = resolve_same(CLAIMS)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["claim_id", "description", "value", "source_file", "derivation", "section", "script"])
        for cid, desc, val, src, how, sect in resolved:
            w.writerow([cid, desc, val, src, how, sect, "; ".join(scripts_for(src))])


    sel = [a for a in args if a.startswith("R-")]
    show = [c for c in CLAIMS if any(c[0].startswith(s) for s in sel)] if sel else CLAIMS
    for c in show:
        print(f"{c[0]:<34} {c[2]}")
        print(f"{'':<34} src: {c[3]}  |  {c[4]}  |  §{c[5]}")
    print(f"\n{len(CLAIMS)} claims written to {os.path.relpath(OUT, ROOT)}")
    if MISSING:
        print("WARNING - missing source files:", ", ".join(sorted(set(MISSING))))

if __name__ == "__main__":
    main()
