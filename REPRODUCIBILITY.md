# Reproducibility statement

**Short answer: the pipeline is reproducible from the cached artefacts, but NOT
byte-identical if re-run from scratch against live databases.** The parts that
differ, and why, are stated explicitly below. Nothing here is hidden.

> **Scope change (2026-07-24).** The sex-SPECIFIC (diagnosis x sex interaction)
> *panel* branch was removed. The panels remain **sex-stratified only**.
> `results/README_GOALS.md` and `results/FIGURE_PROVENANCE.md` are referenced
> throughout this repository but **are not present in this tree** — regenerate
> or remove the references before deposition.
>
> **Robustness layer (2026-07-28).** Five analyses were added that test
> assumptions the pipeline previously made silently: MHC exclusion (`10c`),
> colocalisation (`10d`), cell-composition deconvolution and adjusted DE
> (`05c`), panel-versus-composition benchmarking (`13b`) and the full sex
> interaction report (`05d`). Two returned negative results that change what
> this project may claim — see Section 5 and `results/RESULTS_ROBUSTNESS.md`.
> The headline change: **no panel gene colocalises with the RA signal, so no
> causal claim survives.** The diagnostic claim does survive, in females.

Environment used to produce all reported results: see `ENVIRONMENT.txt`
(R 4.4.2, aarch64-apple-darwin20; limma 3.62.2, WGCNA 1.74, glmnet 5.0,
TwoSampleMR 0.7.8, clusterProfiler 4.14.6, edgeR 4.4.2, pROC 1.19.0.1).

---

## 1. Fully deterministic — WILL reproduce exactly

These steps depend only on local files and seeded randomness. Re-running them on
the same machine reproduces the reported numbers exactly.

| Step | Why deterministic |
|---|---|
| Preprocessing, ComBat, 70:30 split | `set.seed(70)`, local raw files |
| Differential expression (limma) | closed-form, no randomness |
| WGCNA modules | `randomSeed = 1234`, fixed soft power |
| Feature selection (LASSO / RF / SVM-RFE) | `set.seed(1234)` before every stochastic call |
| Nested CV | outer folds seeded `1000 + repeat` |
| Elastic net | seeded `2000 + repeat` |
| Synovium validation, ROC/AUC | deterministic given inputs |
| **MHC-excluded MR sensitivity** | `10c` reads only the `10_MR.R` cache; no network at all |
| **Cell deconvolution + adjusted DE** | `05c` — CIBERSORT is seeded; LM22 ships with IOBR |
| **Panel vs composition benchmark** | `13b` — folds seeded 1234; bootstrap CIs seeded |
| **Sex-interaction report** | `05d` — limma, closed-form |
| **Nested-CV reconciliation** | `16d` — one engine, one seed policy; reproduces `14` exactly |

## 2. NOT guaranteed to reproduce — live database dependencies

| Step | Script | Live resource | Risk |
|---|---|---|---|
| **MR instrument extraction** | `10_MR.R` | eQTLGen + Okada via **OpenGWAS API**; LD clumping uses a **server-side reference panel** | **HIGH** |
| Cross-ancestry MR | `35` | OpenGWAS (Stahl, Biobank Japan) | **HIGH** |
| **Colocalisation regional stats** | `10d` | OpenGWAS regional queries against eQTLGen + `ieu-a-832` | **MEDIUM** — cached to `coloc_regions.rds` on first run; every later run is offline and exact |
| **coloc.susie LD matrices** | `10e` | OpenGWAS `ld_matrix` (1000G EUR reference) | **MEDIUM** — out-of-sample LD; results are indicative, not definitive |
| KEGG enrichment / GSEA | `05`, `06`, `15`, `25` | KEGG REST API (updated monthly) | MEDIUM |
| STRING PPI network | `06` | string-db.org | LOW |
| GEO download | `01`, `02` | GEO (static; already cached locally) | LOW |

### Why the MR step is the real risk
`10_MR.R` (with `CFG$fresh_extract = TRUE`) re-extracts every instrument live. If
eQTLGen is updated, if the LD reference panel changes, or if a transient network
error causes a gene's instruments to be skipped, then:

    different instruments -> different prioritised genes (FS_input) -> different panels

Because feature selection consumes `FS_input`, a change there propagates all the
way to the reported panels (6/6 primary, 4/5 MHC-free). **This is the single point where a
re-run could legitimately produce different biomarkers.**

## 3. What pins the reported results

Every network-derived result is **cached locally**, so the published numbers can
be regenerated without touching the internet:

| Cached artefact | Pins |
|---|---|
| `data/processed/new/MR_instruments.rds` (285K) | the exact instrument set (4,932 SNPs / 1,980 genes) |
| `data/processed/new/MR_primary_objects.rds` (552K) | harmonised SNP data + all MR estimates |
| `results/tables/FS_input_{female,male}.csv` | the 32 / 25 MR-prioritised genes |
| `results/tables/FS_input_{female,male}_noMHC.csv` | the 14 / 14 MHC-free prioritised genes |
| `data/processed/new/coloc_regions.rds` | the regional eQTL/GWAS extracts behind every coloc posterior |
| `data/processed/new/cell_fractions.rds` | CIBERSORT/MCP-counter fractions for all three blood datasets |
| `data/processed/new/MR_mhc_sensitivity_objects.rds` | the MHC-excluded re-run |
| `data/raw/GSE93272_raw.rds`, `GSE110169_raw.rds` | discovery cohorts |
| `data/raw/GSE15573_raw.rds` | external blood |
| `data/raw/GSE89408_counts.txt.gz` | synovium |

> **`data/raw` is a SYMLINK** to a sibling project directory
> (`../Research_Q2_TRANSCRIPTOMICS_sexstratified/data/raw`), so this tree is
> **not self-contained** and cannot be archived as-is. Resolve the symlink
> before deposition.

**Recommended reproduction route:** start from the cached MR objects and
`FS_input_*.csv` and run the downstream chain (feature selection onward). That
path is fully deterministic and reproduces every reported panel and AUC.

**Full-rebuild route** (re-extracting MR from OpenGWAS) is the honest scientific
replication test, but should be expected to give *similar*, not identical,
prioritised-gene lists.

## 4. Independent verification

An independent audit previously recomputed the critical claims from first
principles (24 checks, 24 PASS), covering leakage, sample alignment, group
coding, limma vs plain t-test (Spearman rho = 1.000), panel AUC reproduction
(0.855 female, 1.000 male apparent), MR instrument F >= 10, OR = exp(beta), no
duplicate SNPs, and PTPN22 up-regulation in RA.

> **The audit script itself is missing from this repository** — the doc referred
> to `scripts/AUDIT_independent_checks.R`, which does not exist here (it was not
> removed by the 2026-07-24 cleanup; it was already absent). The reported PASS
> counts therefore cannot currently be re-verified. Rewrite it before deposition.
> One of its 24 checks ("interaction coefficients reproduce by plain `lm`") is in
> any case now out of scope, since the interaction branch has been removed.

## 5. Known limitations of the analysis itself (not reproducibility)

Five of the limitations previously listed here were assumptions rather than
findings. As of 2026-07-28 each has been TESTED, and where the test came back
negative that is recorded below rather than softened. See
`results/RESULTS_ROBUSTNESS.md` for the full write-up.

**Tested, and the result constrains what may be claimed:**

- **The MHC.** 14 of 32 female and 10 of 25 male prioritised genes are instrumented
  from inside the extended MHC. Re-running the entire MR with MHC instruments
  excluded (`10c_MR_mhc_sensitivity.R`) leaves **14 robust prioritised genes in each
  sex**. Three of the six genes in each final panel do not survive: female
  GNL1 and C6orf136 and male VPS52 and HLA-DMA become untestable (no non-MHC
  instrument exists), and ESYT1 loses FDR significance in both sexes through
  Benjamini-Hochberg re-ranking with its point estimate unchanged.
- **Colocalisation.** `10d_coloc_panel_genes.R` runs coloc.abf on regional
  eQTLGen and Okada summary statistics for all 33 testable prioritised genes.
  **No panel gene colocalises.** Six of nine show PP.H3 >= 0.8 — positive
  evidence that the eQTL and the RA association are driven by *different* causal
  variants, which invalidates the cis-MR estimate for those genes. All four MHC
  panel genes have PP.H3 = 1.000. IKZF3 is the best case at PP.H4 = 0.774, and
  even that collapses to 0.255 under the conservative p12 = 1e-6 prior.
  **Consequence: no gene in either panel may be described as causal.** The MR
  step is retained as a genetically-informed *filter* on the candidate space,
  which is a defensible use of it; it is no longer evidence of causality.
- **Cell composition.** CIBERSORT/LM22 (`05c_deconvolution.R`) shows composition
  differs sharply by disease status: female CD8 T-cell fraction 0.162 -> 0.093
  (p = 8e-10), with gamma-delta T-cell, eosinophil and naive B-cell shifts.
  Adjusting the DE model for composition PCs cuts DEGs from 5,131 to 2,709 in
  women and 5,820 to 1,450 in men. **10 of 12 panel gene-sex pairs survive**,
  retaining 57-83 % of their logFC; male VPS52 and INPP5B do not.
- **Panel vs a white-cell count.** `17b_testing_blood_celladjusted.R` benchmarks the
  panel against a composition-only model. In females the panel wins on train
  (0.831 vs 0.662) and external blood (1.000 vs 0.429), and the composition-
  residualised panel still reaches AUC 0.779 train / 0.954 external. **On the
  internal test the margin is 0.007 (0.721 vs 0.714)** — statistically additive
  by LRT (p = 0.0014) but of no practical size. The panel is therefore not
  merely a differential white-cell count, but that claim rests on two of three
  datasets, not three.
- **The sex interaction.** *The previous statement that "No interaction test is
  performed" was wrong* — `05b_dge_sensitivity.R` has performed one since
  2026-07-27. `05d_interaction_report.R` now reports it in full: **53 genes at
  FDR < 0.05** (pre-ComBat matrix, batch modelled explicitly). 35 are
  male-restricted, 12 magnitude differences, 4 opposite-direction, 1
  female-restricted. PTPN22 shows a magnitude difference (female +0.171, male
  +0.893). **Only 7 of the 53 survive composition adjustment**, and **no panel
  gene is among them** — so the panels remain sex-STRATIFIED, not sex-specific.

**Still true, still untested:**

- Male stratum n = 38 (17 RA / 21 HC); internal hold-out n = 13; external blood
  n = 9. The male arm is now labelled `EXPLORATORY (underpowered)` in
  `mr_final_panel_summary.csv` and its internal AUC of 1.000 (1.000-1.000) is
  flagged `SEPARATION`: at n = 13 there are 42 case-control pairs, so perfect
  separation is unremarkable and the interval is degenerate, not precise.
  **No diagnostic claim rests on the male panel.**
- 58 % of female prioritised genes are instrumented by a **single** cis-eQTL SNP, so
  funnel / leave-one-out / MR-Egger are only possible for 11 genes per sex.
- Neither medication (GSE93272 patients are largely treated) nor disease
  activity is modelled. A panel trained on treated prevalent cases carries an
  unmeasured treatment confound.
- Cross-ancestry MR uses European eQTL exposures against an East-Asian outcome:
  exploratory support, not like-for-like replication.
- Effect sizes are small. Formal TREAT testing at a 1.2x threshold reduces the
  female DEG count from 5,131 to 86; at 1.41x, to 3. The |log2FC| > 0.1 gate
  used throughout is a 7 % change.

## 6. To make a future re-run bit-exact

Not yet done; recommended before deposition:

1. `renv::init()` to lock package versions (`renv.lock`).
2. Cache the KEGG annotation locally instead of calling the REST API.
3. Archive the OpenGWAS query results as a versioned release alongside the code.
4. Rewrite `scripts/AUDIT_independent_checks.R` (see section 4).
5. Add a generator for `FIG_G2_01_panel_venn_{female,male}.png`, currently the
   only figures with no script behind them (`results/FIGURE_PROVENANCE.md`).
