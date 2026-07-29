# Sub-chapter 2 — Sex-stratified transcriptomic biomarkers for RA diagnosis
## Methodology (2.1–2.5) and Results — revised to match the executed pipeline (steps 01–09)

> **Provenance of every number in this document.** All figures are taken from the
> outputs actually written by `scripts/00_shared/01`–`09b` on 27 July 2026, not from
> earlier drafts. Sources: `results/tables/combined_cohort_summary.csv`,
> `normalization_diagnostics.csv`, `DEG_summary.csv`, `DEG_array_weights.csv`,
> `dge_enriched_terms_*.csv`, `gsea_kegg.csv`, `WGCNA_01`–`WGCNA_09*.csv`,
> `data/processed/dge_results.rds`, `data/processed/wgcna_results.rds`.
>
> **Two results are marked `[PENDING]`**: module preservation (WGCNA step 13) and the
> candidate-gene counts (step 14 / script 09). `06_WGCNA.R` was re-run on 27 July with
> the variance filter disabled and was still executing the 200-permutation
> `modulePreservation` call at the time of writing. Everything up to and including
> disease-module enrichment is final.

---

# PART 1 — METHODOLOGY

## 2. Methodology

### 2.1 Datasets

Two whole-blood microarray series formed the discovery cohort: **GSE93272**
(Affymetrix GPL570) and **GSE110169** (Affymetrix GPL13667). Two further series were
reserved for external validation: **GSE15573** (peripheral blood mononuclear cells)
and **GSE89408** (RNA sequencing of synovial tissue).

Both discovery series were deposited already normalised, but by different algorithms:
GSE93272 was processed with frozen robust multi-array analysis (fRMA) and GSE110169
with robust multi-array average (RMA). Because fRMA normalises each array against a
fixed external reference while RMA normalises within the processed batch, the two
series differ systematically in per-sample dynamic range. This difference was measured
rather than assumed, and the diagnostics are reported in the Results (Section 3.2).

### 2.2 Data pre-processing and batch correction

All analyses were performed in R 4.4.2. Series were imported as `ExpressionSet`
objects and the accompanying sample metadata harmonised into a common schema: sex was
standardised to a single-character code (F/M) and disease status recoded to rheumatoid
arthritis (RA), healthy control (HC) or systemic lupus erythematosus (SLE).

**Probe-to-gene summarisation and quality control.** Because the two cohorts were
measured on different arrays, expression was summarised from probe level to gene level
to create a shared feature space. Probes with no detected signal in any sample, probes
mapping ambiguously to multiple symbols (the `///` delimiter) and probes without a
valid gene symbol were removed. Remaining probes were collapsed to unique gene symbols
with the `collapseRows` function of the WGCNA package under the MaxMean rule, which
retains for each gene the probe with the greatest mean expression (Miller et al.,
2011). Following the standard GEO/limma protocol, each matrix was checked for log₂
scale and log₂-transformed only where the data were found to be on a linear scale
(99th percentile > 100); neither series required transformation.

Samples were then restricted to participants with a known group label and a known sex;
SLE and sex-unknown samples were removed. In GSE93272, where several participants
contributed repeated measurements, a single baseline sample per participant was
retained.

The two cleaned matrices were reduced to their common genes — the intersection of the
gene symbols measured on both platforms — and combined by column-binding into a single
merged matrix.

**Leakage-safe partition.** The combined data were divided into a 70 % training set and
a 30 % internal-validation hold-out **before any cross-sample normalisation or batch
correction** (`caret::createDataPartition`, `set.seed(70)`). The split was stratified
by the combination of dataset, group and sex so that RA/HC balance, sex balance and
study representation were preserved in both partitions. The hold-out was sealed in its
pre-normalisation, pre-batch-correction state and was not revisited until the biomarker
panel and classifier had been fully locked (Ambroise & McLachlan, 2002; Simon et al.,
2003; Kaufman et al., 2012).

**Normalisation.** The training set was quantile-normalised
(`normalizeBetweenArrays`, `method = "quantile"`, limma), aligning per-sample intensity
distributions to a common reference distribution learned from the training samples
only. The decision to normalise was made on a measured criterion — the between-sample
spread of per-sample medians — rather than applied unconditionally; the before/after
diagnostics are reported in Section 3.2.

**Batch correction.** Residual technical variance attributable to study and platform
origin was removed from the training set with the ComBat empirical-Bayes method (sva;
Johnson et al., 2007), using parametric priors (`par.prior = TRUE`) on the
quantile-normalised training matrix. The batch variable was defined at the
**study-and-internal-batch** level, with an automated fallback to study-level batch had
the finer definition failed. The biological signal of interest was protected by
supplying a covariate model matrix `mod = ~ group + sex`, so that RA-versus-HC and
sex-associated expression differences were retained rather than removed as batch
effects (Nygaard et al., 2016). All normalisation and ComBat parameters were estimated
on the 70 % training set alone.

**Projection of the internal hold-out.** At validation the frozen hold-out was
projected onto the training scale without re-estimating any parameter. Quantile
normalisation was applied against the fixed training reference distribution
(`normalize.quantiles.use.target`, preprocessCore), after which the training ComBat
standardisation and empirical-Bayes batch parameters were re-applied out of sample.
The covariate design matrix used for the projection contains **sex** — which is known
for any prospective sample — but sets the diagnosis contrast to the training reference
level for every hold-out sample, so no sample's own outcome label contributes to its
own normalised values. This matters both statistically and practically: the ComBat
standardisation term does not cancel across the batch adjustment, so including true
labels would shift each sample in proportion to its own diagnosis; and a pipeline that
required a diagnosis in order to normalise a sample could not be applied to a
prospective patient, whose diagnosis is by definition unknown at the point of assay.
The batch location-and-scale correction is applied in full while the RA-versus-control
difference is left in the data for the locked model to detect. The fidelity of the
frozen procedure was verified by reconstructing the training matrix from the stored
parameters and comparing it with the saved ComBat output.

### 2.3 Differential gene expression analysis

Differential expression was performed on the 70 % training set only
(183 samples × 15,763 common genes), using the quantile-normalised,
ComBat-corrected matrix. Differentially expressed genes (DEGs) between RA and HC were
identified with linear models and empirical-Bayes moderation of the gene-wise variance
in limma (Smyth, 2004; Ritchie et al., 2015). A fixed seed (`set.seed(1234)`) was used
throughout.

**Empirical array quality weights.** Per-array weights were estimated by REML and
applied within the linear-model fit (`arrayWeights`; Ritchie et al., 2006), so that
arrays whose residuals are consistently large across the whole transcriptome are
down-weighted rather than either excluded by an arbitrary quality rule or allowed to
contribute equally. This step was motivated by three properties of the design: the
two-platform merge (GPL570 + GPL13667), variable RNA integrity in GSE93272 (RIN
7.0–9.0), and the small male stratum, in which a single poor-quality array has high
leverage on the estimated fold change. Array weights act **per array** and are
orthogonal to the empirical-Bayes moderation, which borrows strength **per gene**
across the transcriptome.

**Model and contrast.** For each comparison a design matrix `~ group` was constructed
with HC as the reference level, so that the `groupRA` coefficient estimates the log₂
fold change of RA relative to HC. Gene-wise models were fitted with `lmFit`, moderated
with `eBayes`, and moderated *t* statistics, *p* values and log₂ fold changes extracted
with `topTable`. A gene was declared significantly differentially expressed when it met
**both** an effect-size and a statistical criterion: |log₂ fold change| > 0.1 (a
1.07-fold change) **and** a Benjamini–Hochberg false-discovery rate below 0.05
(Benjamini & Hochberg, 1995). Direction was labelled "up in RA" (log₂FC > 0) or
"down in RA" (log₂FC < 0).

**Three configurations at identical thresholds** were fitted: a female-only contrast, a
male-only contrast and a pooled contrast using all training samples. Each is an
RA-versus-HC comparison computed **within** the stratum concerned. No gene is described
as sex-specific on the basis of appearing in one stratum's list and not the other's; the
smaller male stratum can fail to reach significance from reduced power alone (Gelman &
Stern, 2006).

**Diagnosis-by-sex interaction.** Because a set difference between two separately
thresholded lists is a statement about power rather than about biology, the interaction
was tested explicitly as a fourth configuration. The model `~ group * sex + batch_full`
was fitted to the **pre-ComBat quantile-normalised** matrix, and the `groupRA:sexM`
coefficient tested at FDR < 0.05. The pre-ComBat matrix is used deliberately: ComBat was
run with `mod = ~ group + sex`, which protects the two main effects but not their
interaction, so an interaction fitted on the ComBat output is both partly absorbed and
never debited the degrees of freedom ComBat consumed (Nygaard et al., 2016). Both
matrices were fitted and the counts compared; the pre-ComBat result is reported. Each
interaction gene is additionally reported with the two within-sex effects that generate
it, fitted on the same matrix and design, and classified as opposite-direction,
magnitude-difference, or sex-restricted.

**Cell-composition estimation and adjustment.** Whole blood is a mixture, so a
differential expression can reflect either cell-intrinsic regulation or a change in the
proportions of the cells being averaged. Leukocyte composition was estimated with
CIBERSORT and the LM22 22-subset signature matrix (Newman et al., 2015), run in array
mode on linear-scale intensities, and cross-checked with the marker-based MCP-counter
(Becht et al., 2016). Because fractions are compositional and sum to one, they were
centred-log-ratio transformed (Aitchison, 1982) before principal component analysis, and
the first three components — 62 % of composition variance — were added to the within-sex
design as covariates. Differences in composition between RA and control were tested
within each sex by Wilcoxon rank-sum test with Benjamini–Hochberg correction.

The adjusted model is **not** presented as a bias correction. Leukocyte composition in
RA is plausibly a consequence of the disease and therefore a mediator, and adjusting for
a mediator removes part of the true effect rather than a confound. The unadjusted model
is accordingly reported as the estimate of the **total** RA-associated signal — the
relevant quantity for a diagnostic biomarker — and the adjusted model as the estimate of
the signal **not explained by composition**, which is the only basis for a
cell-intrinsic or mechanistic claim about any individual gene.

**Fold-change sensitivity analysis.** The primary rule filters on |log₂FC| *after*
testing against H₀: log₂FC = 0. As a sensitivity analysis the same contrasts were
refitted with `treat`, which places the fold-change threshold inside the null
hypothesis (H₀: |true log₂FC| ≤ τ) and is therefore the statistically correct — and
deliberately conservative — way to demand a minimum effect size (McCarthy & Smyth,
2009). Thresholds τ = 0.1, log₂(1.2) and 0.5 were evaluated and the counts compared
with the post-hoc-filtered counts.

### 2.4 Co-expression network analysis (WGCNA)

The differential-expression analysis of Section 2.3 tests each transcript in isolation,
treating the 15,763 genes as 15,763 independent hypotheses. Genes act through
co-regulated programmes, however, and a coordinated shift distributed across hundreds of
functionally related transcripts may leave every individual gene short of significance
while the programme itself is strongly disease-associated. A gene-by-gene analysis
cannot distinguish a transcript dysregulated as part of a coherent inflammatory
programme from one that moves in isolation, although the two carry very different
evidential weight as candidate biomarkers.

Weighted gene co-expression network analysis (WGCNA) reverses the unit of analysis
(Zhang & Horvath, 2005; Langfelder & Horvath, 2008). Every gene is a node, the weight of
the edge between two genes is a function of their expression correlation across samples,
and the network is partitioned into **modules** of densely interconnected genes. Each
module is summarised by its **module eigengene**, the first principal component of its
member genes' expression (Langfelder & Horvath, 2007), and disease association is tested
at the level of the eigengene rather than the gene. Three properties make this
appropriate here. It is **unsupervised** — modules are defined from the correlation
structure alone, with diagnosis entering only afterwards, so a module that emerges
without knowledge of the phenotype and then proves strongly associated with it is
stronger evidence than a gene set assembled by selecting on that phenotype. It supplies
a **second, independent line of evidence** for candidate selection, module membership
being a network-level property and differential expression a gene-level one
(Section 2.5). And it permits a **sample-size-controlled test** of whether network
structure differs between the sexes (Langfelder et al., 2011), which a comparison of DEG
lists between strata of unequal size cannot provide.

Co-expression is nonetheless a correlational, undirected construct: modules identify
genes that vary together and do not establish regulatory direction, causal ordering or
physical interaction. The network supplies candidates; causal inference is carried by
the Mendelian randomisation that follows.

The network was built on the 70 % training set only (quantile-normalised,
ComBat-corrected), so no information from the sealed hold-out can enter module
definition, disease-module selection or candidate identification. Analyses used R 4.4.2
with WGCNA 1.74 and clusterProfiler 4.14.6, with a fixed seed throughout
(`set.seed(1234)`, and `randomSeed = 1234` passed to `blockwiseModules`). All network
parameters are declared in a single configuration block, and the analysis cache is keyed
to a hash of every network-relevant parameter together with a fingerprint of the input
matrix, so a parameter or data change forces recomputation and a stale network cannot be
silently reused.

**Data-integrity screening.** Genes and samples were first screened with
`goodSamplesGenes` (Langfelder & Horvath, 2008), which removes genes of zero variance or
excessive missingness and samples with too many missing entries. This is a precondition
rather than a filter: a zero-variance gene has an undefined correlation and would
propagate missing values through the entire adjacency matrix. All 15,763 genes and all
183 samples passed, so this step removed nothing and the screen is reported as evidence
that the preprocessing of Section 2.2 left no degenerate features.

**Sample-outlier detection and removal.** Following the procedure of the WGCNA authors'
tutorial (Langfelder & Horvath, 2008), samples were clustered by average-linkage
hierarchical clustering on Euclidean distance, and a data-driven outlier cut height was
derived as the mean plus three standard deviations of the merge heights. Arrays falling
outside the principal cluster at that height were **removed**, not merely reported.

Removal rather than inspection is necessary because a weighted co-expression network is
constructed entirely from between-sample correlations (Zhang & Horvath, 2005), so a
single aberrant array perturbs every gene pair simultaneously and generates
co-expression structure that reflects the artefact rather than the biology.
Co-expression analysis is considerably more sensitive to this failure than differential
expression, in which a poor array inflates the residual variance of one gene at a time
and is further absorbed by the empirical array quality weights of Section 2.3 (Ritchie
et al., 2006); detecting outliers and then retaining them is therefore not a
conservative option but an uncorrected source of spurious modules.

Average linkage was used because it avoids the chaining behaviour of single linkage,
which absorbs isolated arrays into the principal cluster, while not imposing the
equal-sized, spherical clusters assumed by Ward's method (Ward, 1963) — neither property
being appropriate when the target is a small number of isolated arrays against one large
group. The cut height was derived from the merge-height distribution rather than read
off the dendrogram as in the WGCNA tutorial, so that the criterion is a computed
property of the data and transfers unchanged to a new dataset; three standard deviations
is the conventional choice and is the conservative direction, removing fewer arrays than
a two-standard-deviation rule would. Alternative outlier criteria were considered and
not adopted: standardised connectivity (Z.k < −2.5; Oldham et al., 2012) measures the
same construct through a different metric and was not run; principal-component
inspection supplies no threshold and uses only the leading components; and raw-intensity
diagnostics such as `arrayQualityMetrics` (Kauffmann et al., 2009) or the RLE and NUSE
statistics (Brettschneider et al., 2008) no longer apply to a matrix that has already
been quantile-normalised (Bolstad et al., 2003) and ComBat-corrected (Johnson et al.,
2007).

The removed set was instead cross-checked against an independent criterion computed from
an entirely different quantity — the limma REML array quality weights of Section 2.3
(Ritchie et al., 2006), which derive from gene-wise residual variance rather than
between-sample distance. The convergence of two unrelated quality criteria on the same
arrays establishes that the removal identifies genuinely poor-quality data rather than
an arbitrary trim (Section 3.5).

**No variance pre-filter.** All 15,763 genes surviving the integrity screen entered the
network; a variance filter was *not* applied. This departs from common practice: the
WGCNA package FAQ permits filtering by mean or variance (Langfelder & Horvath, 2008),
and published protocols treat it as a standard preprocessing step, in some cases
retaining only the top 5 % of genes by variance (Nguyen & Zeng, 2025). The departure is
deliberate, and it follows from the role the network plays in this particular study.

Where the module is itself the endpoint of the analysis, a variance filter discards
genes that would have contributed little to module structure and costs the result
little. Here the network is not the endpoint. Its sole function is to supply the
disease-module background that is intersected with each sex's DEG list to define the
candidate genes carried into Mendelian randomisation (Section 2.5). A gene removed
before module detection can never be assigned to a disease module, can never enter that
intersection, and can therefore never be tested for a causal relationship with RA. The
filter would operate as a silent and permanent veto over the study's final output,
exercised on marginal variability — a criterion bearing no relationship to the question
being asked.

The cost was measured rather than assumed. An earlier configuration of this pipeline
applied a 40 % variance filter, which removed 6,305 of the 15,763 genes and with them
**14.9 % of the 5,131 female DEGs and 20.9 % of the 5,820 male DEGs**. The loss is
asymmetric, and it falls in the direction least affordable to this study: the male
stratum comprises 38 of the 183 training samples and 17 RA cases, and is already the
weaker of the two in statistical power (Section 3.1), so a filter removing a fifth of
its DEGs erodes precisely the stratum this chapter exists to characterise. A
95th-percentile filter of the kind used in the protocol literature would have been more
severe still. No computational consideration required the filter in compensation: all
15,763 genes across 173 samples were accommodated in a single block
(`maxBlockSize = 20000`), so modules were detected globally rather than block-wise.

The alternative filtering strategies were rejected on the same or stronger grounds. A
mean-expression filter addresses the count-driven mean–variance dependence
characteristic of RNA-seq data (Law et al., 2014), which is absent from a
quantile-normalised, batch-corrected microarray matrix. Retaining a fixed number of the
most variable genes replaces a percentile with an equally arbitrary round number while
inheriting the same structural objection. Filtering to differentially expressed genes is
explicitly forbidden by the WGCNA FAQ (Langfelder & Horvath, 2008) and would render the
analysis circular, since genes selected for their association with RA would then be used
to form modules whose correlation with RA is reported as a finding — the form of
selection bias described by Ambroise & McLachlan (2002) and Simon et al. (2003).

Low-information genes are not thereby ignored; they are excluded softly, by the network
itself. Soft-thresholding at the chosen power drives weak correlations toward zero
(Zhang & Horvath, 2005), so near-invariant genes acquire negligible connectivity and are
assigned to the unassigned
(grey) module rather than distorting genuine ones. The distinction is that this
exclusion is applied *after* the correlation structure is known and on the basis of
connectivity, whereas a variance filter is applied *before* it and on the basis of
marginal variability alone: a gene of modest variance that genuinely co-expresses with a
disease module survives the former but cannot survive the latter.

**Correlation coefficient.** Pearson correlation was used (`corType = "pearson"`).

The WGCNA FAQ recommends the biweight midcorrelation (`bicor`, `maxPOutliers = 0.05`)
in preference to Pearson, because it is resistant to outlying observations *within* a
gene (Langfelder & Horvath, 2008). A small number of extreme expression values can
inflate a Pearson correlation and thereby create edges — and ultimately module
membership — that reflect a handful of arrays rather than a consistent co-expression
relationship.

`bicor` was implemented and evaluated, and was reverted on computational grounds.
Whereas Pearson correlation has a closed-form solution, `bicor` estimates its weights
iteratively, requiring additional working copies of the data at roughly three times the
computational cost. This interacts poorly with the decision to retain all 15,763 genes:
a single 15,763 × 15,763 double-precision matrix occupies approximately 2 GB, and the
topological-overlap calculation holds the adjacency matrix and the TOM in memory
simultaneously together with intermediates. On the analysis machine (16 GB) the `bicor`
run entered swap and did not complete, whereas the Pearson network completed in a single
block.

This is a computational constraint and is not advanced as a claim that Pearson is
preferable to `bicor`. It is also, and should be acknowledged as, a consequence of the
gene-filtering decision above: on the available hardware, retaining every gene and using
the FAQ-recommended correlation were mutually exclusive. Of the two FAQ-permitted
departures, the one preserving every gene's eligibility for the candidate intersection
was chosen, because a gene removed by a variance filter is lost permanently from the
study's output (Section 2.5), whereas the robustness forgone by using Pearson is partly
recoverable by other means.

Robustness to outliers was accordingly addressed at the **sample** level instead. The
ten outlying arrays identified above were removed, and that removal was independently
corroborated against the limma array quality weights. Because an aberrant array
contributes an extreme value to every gene simultaneously, deleting it removes the most
common source of within-gene outliers at source: an array-driven artefact cannot survive
the deletion of the array that produced it.

*Limitation.* Sample-level removal does not address within-gene outliers arising from
causes other than global array quality — for example a probe behaving erratically in a
subset of individuals. Such values remain capable of inflating the Pearson correlations
involving that gene and could in principle strengthen or create its module membership.
Three considerations bound this risk without eliminating it: quantile normalisation
constrains the range of extreme values; the disease modules are large, so no single
gene's correlations determine the module eigengene; and hub genes are defined by |kME|
and |GS| jointly, so a gene whose apparent centrality rested on a few extreme
observations would still require a genuine correlation with disease status to be
reported. The residual risk is nonetheless real and is recorded rather than discounted.

**Network type and topological overlap.** A **signed** network was constructed
(`networkType = "signed"`, `TOMType = "signed"`). An unsigned network uses |cor|, so a
strongly *anti*-correlated gene pair is treated as strongly connected and the two genes
may be placed in the same module. For a disease-programme analysis this is undesirable:
a transcript up-regulated in RA and one down-regulated in RA are biologically opposed,
and a module containing both has no coherent direction, no interpretable eigengene, and
cannot be described as "up in RA" or "down in RA" — precisely the description required
by the directional consistency check of Section 2.5. A signed network preserves the sign
of the correlation, so modules are directionally homogeneous. Signed-hybrid, which sets
negative correlations to zero rather than mapping them to low adjacency, behaves
similarly for module detection; fully signed was preferred for its monotone, continuous
adjacency transformation.

Modules were defined from the topological overlap matrix (TOM) rather than from the
adjacency matrix directly (Zhang & Horvath, 2005). The TOM measures the extent to which
two genes share their network neighbours in addition to being connected to each other,
so it reflects shared regulatory context rather than pairwise association alone. This
suppresses spurious edges arising from noise in individual correlations, since a
spurious pair is unlikely to share a neighbourhood, and yields more reproducible modules
than clustering on correlation directly.

**Soft-thresholding power.** Scale-free topology fit and mean connectivity were
evaluated across candidate powers (1–10, then 12–20 in steps of two) with
`pickSoftThreshold` for a signed network. Soft-thresholding raises the correlation to a
power β so that strong correlations are preserved while weak ones are driven toward zero
— a continuous alternative to a hard correlation cut-off, which would discard the
magnitude information distinguishing a correlation of 0.9 from one of 0.5 (Zhang &
Horvath, 2005). β is chosen so that the resulting network approximates scale-free
topology.

`pickSoftThreshold` returned an estimate of **β = 3**, the first power reaching
R² ≥ 0.85 (R² = 0.873). **This power was rejected on connectivity grounds.** The
estimate is defined as the *first* power clearing the R² threshold rather than the best
one, and at β = 3 the mean connectivity is 2,240 — a network in which the average gene
is connected to one in seven of all others, far too dense for meaningful module
resolution. The signed R² curve plateaus at approximately 0.90 from β = 9 onwards while
mean connectivity falls steeply, so powers in the plateau achieve equivalent scale-free
fit at a far better conditioned connectivity. **β = 12 was used**, where R² = 0.901,
slope = −1.78 and mean connectivity = 47.4. The power was therefore selected on the
joint criterion of scale-free fit **and** network connectivity, with the full
fit-versus-connectivity table reported at every candidate power (Section 3.5) so that
the choice can be audited.

The WGCNA FAQ recommends β = 12 for a signed network with more than 40 samples, which
coincides with the value used. This coincidence is reported as *corroboration only* and
is deliberately not offered as the justification, because the FAQ table is formally
conditional on the scale-free fit **failing**. Here the fit did not fail — an estimate
was returned — so the table does not formally apply, and citing it as the primary reason
would misapply the authors' own guidance.

**Module detection.** Modules were detected in a single block with `blockwiseModules` at
β = 12, using the dynamic tree cut algorithm (Langfelder, Zhang & Horvath, 2008). The
parameters and the reason for each are given in Table 2.1.

**Table 2.1 — Module-detection parameters, with justification.**

| Parameter | Value | Why this value, and why not another |
|---|---|---|
| `maxBlockSize` | 20000 | Exceeds the 15,763 genes, so the network is built in **one block**. The default (5,000) would split the genes into four independently clustered blocks and modules could not span them — an arbitrary partition of the gene space imposed by memory management rather than biology. Raised from 16,000 when the variance filter was removed. |
| `minModuleSize` | 30 | The package default. Smaller values (10–20) fragment the dendrogram into many small modules whose eigengenes are unstable and whose enrichment is uninterpretable; larger values (50+) merge distinct programmes. Left at the default rather than tuned, so that module count was not adjusted toward a preferred outcome. |
| `mergeCutHeight` | 0.25 | Merges modules whose eigengenes correlate above 0.75. Two modules at that correlation are largely redundant summaries of the same variation, and retaining both would inflate the number of module–trait tests without adding information. A stricter value (0.15, r > 0.85) leaves near-duplicate modules; a looser one (0.35, r > 0.65) merges genuinely distinct programmes. |
| `deepSplit` | 2 (default) | Controls dynamic tree cut sensitivity on a 0–4 scale. Left at the default; not tuned. |
| `reassignThreshold` | 0 | Disables post-hoc reassignment of genes between modules on the basis of kME significance, so membership is determined by one criterion (topological overlap) rather than two and cannot drift under a second, correlation-based rule. |
| `pamRespectsDendro` | FALSE | Allows PAM-stage assignment of otherwise unassigned genes on topological overlap without being constrained by dendrogram branch boundaries; the setting used in the authors' tutorial. |
| `numericLabels` | TRUE | Modules are produced as integers and converted to colours only for display. |
| `randomSeed` | 1234 | Fixes the stochastic component so module assignment is exactly reproducible. |

Each module was summarised by its module eigengene, and eigengenes were ordered by
hierarchical clustering (`orderMEs`; Langfelder & Horvath, 2007). The grey module is not
a module — it collects genes not co-expressed with any coherent group — and is excluded
from disease-module selection and from all downstream analysis.

**Module–trait association.** A binary trait matrix encoded disease status (RA, HC),
sex, the four group-by-sex combinations and age. Each trait was correlated with every
module eigengene by Pearson correlation with Student asymptotic *p* values
(`corPvalueStudent`; Langfelder & Horvath, 2008). The reported *p* values are
unadjusted; no correction was applied because the disease-module threshold used for
selection (*p* < 1 × 10⁻⁸, below) is several orders of magnitude more stringent than a
Bonferroni correction over the tests performed, so the selection decision is unaffected
by the choice of adjustment. The same association was additionally recomputed **within each
sex stratum separately** (all samples, female only, male only), so that the disease
association of every module can be read inside each sex. The female and male columns
are separate within-sex contrasts and are not compared with one another: a module
correlating more strongly with RA in men than in women is not evidence of a sex
difference, because the male stratum is far smaller and its correlations correspondingly
noisier and more extreme (Gelman & Stern, 2006).

**Disease-module selection is data-driven, not colour-based.** Disease-associated
modules were selected from the overall, all-sample module–trait analysis by a
pre-specified rule: **|cor(ME, RA)| ≥ 0.5 and association p < 1 × 10⁻⁸**, with the grey
(unassigned) module excluded.

A rule was used rather than a ranking. Describing modules as ranked by their correlation
with RA and then selecting the top ones leaves the number selected undetermined and the
threshold implicit, so the selection can be adjusted after the result is seen; an
explicit, pre-declared threshold fixes the decision before the outcome is known and can
be checked by a reader against the reported correlation table. The two components of the
rule act on different failure modes and both must be met: the correlation cut-off
requires a module to explain a substantial share of disease-related expression variance
rather than merely reaching significance, while the *p* threshold guards against a
spurious module–trait relationship.

Modules meeting the rule in either direction were retained, so that both up-regulated
and down-regulated disease programmes are captured; restricting to positive correlations
would discard half of the disease signal and bias the candidate set toward genes
up-regulated in RA. This is the same criterion applied at the candidate-selection stage
(Section 2.5), and the identical module set is carried forward for both sexes so that no
sex is assigned a different co-expression background — a deliberate correction of an
earlier design in which one module was assigned to each sex, a procedure that guarantees
a difference between the sexes by construction.

It is emphasised that **WGCNA colour names are not stable identifiers.** `labels2colors`
assigns colours by module size rank, so the same biology is renamed whenever the gene
filter, the soft power or the sample set changes. Colour names are therefore recorded
as an *output* of the analysis and are never used as an *input*: no downstream script
refers to a colour literal, and every module is selected by its correlation with RA.

**Gene-level importance.** Within each selected disease module, gene importance was
quantified by three complementary measures (Zhang & Horvath, 2005; Langfelder &
Horvath, 2008): module membership (kME, the correlation between a gene and its module
eigengene, i.e. how representative the gene is of the module's programme); gene
significance (GS, the correlation between a gene and RA status, i.e. its individual
disease association); and intramodular connectivity (row sums of the signed adjacency at
β = 12, computed *within* the module, i.e. how central the gene is to the module's
topology). Connectivity was computed within the module rather than across the whole
network, so that hub status reflects centrality in the disease programme rather than
global connectivity.

Hub genes were defined as genes with **|kME| > 0.8 and |GS| > 0.2**. Requiring both
means a hub must be central to its module *and* disease-associated: kME alone would
select genes central to a module irrespective of their relationship to RA, and GS alone
would select disease-associated genes irrespective of network position, which is merely
a weaker restatement of the differential-expression analysis. The kME threshold of 0.8
is the conventional value and identifies genes tracking the eigengene closely; the GS
threshold of 0.2 is deliberately permissive, because GS is a whole-cohort correlation
with a binary trait and its magnitude is bounded well below that of kME. These analyses
were performed separately for each selected disease module, so that the up-regulated and
down-regulated programmes receive equivalent treatment.

**Functional enrichment.** Gene Ontology enrichment across all three ontologies
(`enrichGO`, `ont = "ALL"`) and KEGG pathway enrichment (`enrichKEGG`,
`organism = "hsa"`, following identifier conversion with `bitr`) were applied to the
pooled disease-module gene set with clusterProfiler (Yu et al., 2012; Wu et al., 2021),
using Benjamini–Hochberg correction (Benjamini & Hochberg, 1995) and p- and q-value
cut-offs of 0.05.

The background universe was restricted to **all genes present in the network** rather
than to the whole genome or the whole annotation. This is the correct comparator: the
disease-module genes were drawn from the genes that entered the network, so those are
the genes that *could* have been selected. A genome-wide universe would test the disease
modules against genes that were never eligible, and would report as enrichment any bias
in the platforms' probe coverage or in the two-platform gene intersection.

**Sex-stratified networks and module preservation.** Two further networks were
constructed independently on the female and the male samples at the same soft-threshold
power. Whether co-expression structure differs between sexes was then tested formally
with `modulePreservation` (Langfelder et al., 2011), using the female data as the
reference network and the male data as the test network, with **200 permutations**.

A permutation statistic was used rather than a direct comparison because preservation
Z-summary statistics are constructed by permuting sample labels within the test set, so
the null distribution is generated at the test set's own size. This controls for the
fact that the male stratum is much smaller than the female stratum. A direct comparison
of male and female module structure, or of male and female DEG lists, does not control
for this and would report reduced power in the male stratum as a biological sex
difference.

The reference module assignment was taken from the **combined-network** modules rather
than from the female-network modules. This is essential: female-network modules are a
different gene partition from the disease modules even where colour names coincide, so
using them would answer the question "are female-network modules preserved in men?"
rather than the question actually required, "are the *disease* modules preserved between
sexes?". A preservation Z-summary below 2 indicates non-preserved structure, between 2
and 10 moderate preservation, and above 10 strong preservation.

Preservation was assessed in one direction only, with the female network as reference
and the male network as test, because the male stratum is too small to serve as a stable
reference network and the question of interest is whether structure established in the
larger stratum recurs in the smaller one. That the comparison is therefore not symmetric
is recorded as a limitation.

This is a permutation-based, sample-size-controlled **network stability** statistic; it
is not a per-gene test of whether the RA effect differs between sexes, and no such
interaction test is performed in this chapter. No gene is described as sex-specific on
the basis of appearing in one stratum's list and not the other's.

**Relationship to the differential-expression analysis.** The differential-expression
contrasts of Section 2.3 were fitted on all 183 training samples, whereas the network
was constructed on the 173 samples surviving array-outlier removal. Because the
candidate set of Section 2.5 intersects a gene list derived from the first analysis with
one derived from the second, the discrepancy requires justification rather than mere
disclosure.

**The two analyses did not ignore array quality differently; they controlled it by the
mechanism each method provides.** limma incorporates array quality *inside* the
estimator, through the REML empirical array weights of Section 2.3 (Ritchie et al.,
2006), which down-weight poor arrays continuously in proportion to their residual
variance rather than deleting them. WGCNA offers no equivalent: `blockwiseModules`
admits no per-sample weight, and every sample contributes equally to every pairwise
correlation, so deletion is the only available lever. Applying weighting where weighting
is available and deletion where it is not is a method-appropriate choice, not an
inconsistency.

**The difference between the two treatments is of degree, not of kind.** All ten arrays
removed from the network carry a limma array weight below 1.0, and six of them fall
below 0.5 — six of only eight such arrays in the entire training set. The median weight
of the removed set is 0.482 against 1.079 for the retained samples, a 2.2-fold
difference. These arrays were therefore already heavily discounted in the
differential-expression analysis; they were not contributing to the DEG results as
full-weight observations. The network deletes what limma down-weights.

**Deletion was not adopted for the differential-expression analysis because its cost
runs in the opposite direction.** Statistical power in a linear model depends directly
on the number of samples per group, and the male stratum comprises only 38 samples and
17 RA cases (Section 3.1); removing ten arrays would erode the very contrast this
chapter is least able to afford to weaken. Weighting retains the information in a
marginal array while discounting its unreliable component. A co-expression network, by
contrast, loses little between 173 and 183 samples — correlation estimates are stable at
either size, and the WGCNA FAQ's minimum is 15 to 20 samples — while being acutely
sensitive to a single aberrant array, which perturbs every gene pair simultaneously. The
trade-off between sample count and array quality therefore resolves differently for the
two methods, and it was resolved separately for each.

**The removal does not alter cohort composition materially.** The proportion of RA cases
changes from 56.3 % (103/183) to 54.9 % (95/173) and the proportion of male participants
from 20.8 % to 20.2 %. The removed set contains eight RA and two HC samples, an
imbalance that does not depart significantly from the cohort ratio (hypergeometric
*p* = 0.108). The network is therefore built on the same population as the
differential-expression analysis, not on a differently constituted subset.

**What the intersection requires.** The operation of Section 2.5 is a set intersection
over gene identities within a shared feature space of 15,763 common genes. Its validity
requires that each list be individually valid on the samples from which it was derived,
not that the two sample sets coincide. What the discrepancy affects is the evidential
basis of each list, and the number of samples entering each analysis is therefore
reported explicitly wherever results are given.

### 2.5 Candidate gene identification (disease module ∩ sex-stratified DEG)

Network-level and gene-level evidence were then combined to define the candidate genes
carried into Mendelian randomisation. The purpose of the integration is to retain only
genes supported concurrently by two independent lines of evidence: **membership of a
disease-associated co-expression module**, which indicates participation in a
coordinated, systemically dysregulated programme, and **significant differential
expression within a given sex**, which indicates an individual disease association.

The disease-module gene background was formed by pooling the genes of the modules
selected by the rule of Section 2.4 (|cor(ME, RA)| ≥ 0.5 and p < 1 × 10⁻⁸ in the
overall, all-sample module–trait analysis). The conservative *p*-value threshold guards
against spurious module–trait relationships while the correlation cut-off ensures that
retained modules capture a substantial share of disease-related expression variance.
Modules changing in either direction were eligible, so the background contains both the
up-regulated and the down-regulated disease programme.

This shared background was then intersected with each sex's significant DEG list
separately:

- **Female candidates** = disease-module genes ∩ female DEGs
- **Male candidates** = disease-module genes ∩ male DEGs

Because the module background is identical for both sexes, **all sex resolution enters
only through the sex-stratified expression contrast**. Any downstream difference between
the two candidate lists is therefore attributable to expression biology rather than to
an analytical asymmetry in which sex was assigned which module. This is a deliberate
correction of an earlier design in which one module was assigned to each sex, a
procedure that guarantees a difference between the sexes by construction.

Genes appearing in one sex's list but not the other are **not** thereby shown to differ
between sexes: a gene may fail significance in one stratum through reduced statistical
power alone, particularly in the smaller male stratum. No claim of sex-specificity is
made on the basis of such a set difference, and no diagnosis-by-sex interaction test is
performed.

Each retained candidate was annotated with its co-expression module, the direction of
dysregulation in RA, the log₂ fold change and FDR-adjusted significance in the relevant
sex-stratified contrast, the module membership kME (how central the gene is within its
module) and the gene significance for RA. Overlap was visualised as a two-set Venn
diagram per sex.

**Directional consistency check.** As a check on the coherence of the intersection, each
disease module's overlap with the up-in-RA and the down-in-RA DEG lists was tabulated
separately within each sex. For a module with cor(ME, RA) > 0 the biologically
consistent overlap is with the up-in-RA DEGs, and for a module with cor(ME, RA) < 0 the
expectation reverses; a well-behaved disease module shows a consistent overlap
substantially larger than its inconsistent one. The expected direction is derived from
each module's own eigengene correlation rather than assumed.

---

# PART 2 — RESULTS

## 3. Sex-stratified transcriptomic signatures for RA

### 3.1 Whole-blood cohort description

After probe-to-gene summarisation, GSE93272 yielded 20,848 genes from 42,894 probes and
GSE110169 yielded 19,041 genes from 47,920 probes. Both series were confirmed to be on
the log₂ scale (99th percentiles 11.86 and 11.64 respectively) and neither required
transformation. Intersecting the two feature spaces gave **15,763 common genes**, the
analysis feature space throughout this chapter.

After removal of SLE samples, sex-unknown samples and repeated measurements, GSE93272
contributed 101 samples (66 RA, 35 HC; 87 female, 14 male) and GSE110169 contributed
156 samples (79 RA, 77 HC; 119 female, 37 male), giving a combined discovery cohort of
**257 samples**.

The stratified 70:30 partition assigned **183 samples (71.2 %) to training** and
**74 samples (28.8 %) to the sealed internal hold-out**. The training set comprised 103
RA and 80 HC participants, of whom 145 were female and 38 male; the hold-out comprised
42 RA and 32 HC, of whom 61 were female and 13 male. Composition by dataset, group and
sex is given in Table 3.1.

**Table 3.1 — Training-set composition (n = 183).**

| Sex | Dataset | HC | RA | Total |
|---|---|---|---|---|
| Female | GSE110169 | 38 | 46 | 84 |
| Female | GSE93272 | 21 | 40 | 61 |
| Male | GSE110169 | 17 | 10 | 27 |
| Male | GSE93272 | 4 | 7 | 11 |
| **Total** | | **80** | **103** | **183** |

The male stratum is 20.8 % of the training set (38 of 183) and contains 17 RA cases.
This imbalance is a property of the source series — both were recruited from
predominantly female RA populations, consistent with the disease's sex ratio — and it
constrains every male-stratum result reported below. It is stated here once and
referred to throughout rather than repeated at each result.

### 3.2 Normalisation and batch correction

The two series differed in per-sample dynamic range exactly as their different origin
normalisations predict. Mean per-sample medians were 5.43 (GSE110169) and 6.48
(GSE93272), and the between-sample spread of medians differed six-fold between the
series (SD 0.327 vs 0.053). Across the merged training matrix the standard deviation of
per-sample medians was 0.574 and of per-sample IQRs 0.799 — sufficient distributional
heterogeneity to trigger the pre-specified normalisation criterion.

Quantile normalisation reduced the standard deviation of per-sample medians from
**0.574 to 3.3 × 10⁻⁴** and of per-sample IQRs from **0.799 to 4.6 × 10⁻⁴**, i.e. by
more than three orders of magnitude, confirming that the merged distributions were
aligned (Figure: `fig_combine_density_norm`).

ComBat was then applied at the study-and-internal-batch level, resolving **six batches**
(GSE110169 batches 1–4 with 24, 35, 35 and 17 samples; GSE93272 batches 1–2 with 33 and
39 samples), with `mod = ~ group + sex` protecting the diagnosis and sex signals.
Principal-component analysis before and after correction shows the separation by study
collapsing while the RA/HC structure is retained (Figure: `fig_combine_pca_combat`).

The frozen-parameter projection of the hold-out was verified by reconstructing the
training matrix from the stored ComBat parameters; maximum absolute reconstruction
error was **1.78 × 10⁻¹⁵**, i.e. at floating-point precision, confirming that the
projection applies exactly the transformation learned on the training data.

### 3.3 Whole-blood differential gene expression

All three contrasts were fitted on the 183 training samples across 15,763 genes at
|log₂FC| > 0.1 and FDR < 0.05 (Table 3.2).

**Table 3.2 — Differential expression summary.**

| Contrast | n | RA | HC | Significant | Up in RA | Down in RA | % of transcriptome | Median \|log₂FC\| |
|---|---|---|---|---|---|---|---|---|
| All | 183 | 103 | 80 | 6,422 | 2,760 | 3,662 | 40.7 % | 0.203 |
| Female | 145 | 86 | 59 | 5,131 | 2,238 | 2,893 | 32.6 % | 0.192 |
| Male | 38 | 17 | 21 | 5,820 | 2,510 | 3,310 | 36.9 % | 0.337 |

Empirical array quality weights had a median of 1.042 and ranged from 0.338 to 2.055;
**eight arrays received a weight below 0.5** and were correspondingly down-weighted
rather than excluded (four in each source series). Median weights were comparable
between series (GSE93272 1.087, GSE110169 1.006), indicating that array quality was not
confounded with study of origin.

#### Female RA versus female control

The female contrast identified **5,131 DEGs (2,238 up, 2,893 down in RA)**. The most
significant genes were *HMGB2* (log₂FC +0.67, FDR 2.7 × 10⁻¹²), *MAGED1* (−0.39),
*KDM1A* (−0.29), *C1GALT1C1* (+0.65) and *CSGALNACT2* (+0.46). Among genes passing
significance, the largest effects were innate-immune and granulocyte transcripts:
*ARG1* (+1.36), *CLEC4D* (+1.17), *BCL2A1* (+1.16), *COX7B* (+1.02), *DEFA4* (+1.00)
and *S100A8* (+0.50, FDR 5.1 × 10⁻¹¹).

#### Male RA versus male control

The male contrast identified **5,820 DEGs (2,510 up, 3,310 down in RA)**. Leading genes
by significance were *STARD3NL* (+0.70), *TBC1D15* (+1.01), *MIER1* (+0.65), *TTC33*
(+1.30), *TRIM23* (+1.21) and *SMARCA4* (−0.56), all at FDR ≈ 3 × 10⁻⁸. The largest
effects were *BCL2A1* (+1.84), *COX7B* (+1.70), *RPL22L1* (+1.62), *COMMD6* (+1.60) and
*SCOC* (+1.59).

**The male stratum yields more DEGs than the female stratum despite having a quarter of
the samples, and this must not be read as a stronger male signal.** The median absolute
fold change among significant genes is 0.337 in men against 0.192 in women — a 1.8-fold
inflation. With n = 38 the per-gene fold-change estimate is far noisier, so more genes
clear the |log₂FC| > 0.1 effect-size gate even though the evidence supporting each is
weaker; the minimum attainable FDR in the male contrast (3.3 × 10⁻⁸) is four orders of
magnitude larger than in the female contrast (2.7 × 10⁻¹²). The male DEG count reflects
estimator variance, not effect magnitude.

#### Shared female and male blood DEGs

Of the two lists, **3,857 genes are significant in both sexes**, 1,274 in women only and
1,963 in men only, giving a union of 7,094 genes. Concordance among the shared genes is
near-complete: **3,854 of 3,857 (99.9 %) change in the same direction in both sexes**,
and the genome-wide correlation between the female and male log₂ fold changes across all
15,763 genes is **r = 0.842**.

This is the central descriptive finding of the differential-expression stage. The two
sexes share a large, strongly concordant RA transcriptional response; the apparent
divergence between the lists is a thresholding artefact of unequal power, not evidence
of two different diseases. Consistent with the analysis plan, no gene is designated
sex-specific on the basis of appearing in one list and not the other.

#### The diagnosis × sex interaction

A set difference between two separately-thresholded lists is a statement about power,
not about biology: the female stratum has 145 samples and the male stratum 38, so the
lists would differ even if the underlying effects were identical. The only test that
speaks directly to the question is the **diagnosis × sex interaction**, whose null
hypothesis is that the RA-versus-control effect on a gene is the same in both sexes.
That test was fitted on the pre-ComBat quantile-normalised matrix with batch modelled
explicitly (`~ group * sex + batch_full`), and identified **53 genes at FDR < 0.05**
(`05d_interaction_report.R`; Table 3.3a).

The pre-ComBat model is reported rather than the ComBat-corrected one, which returns
270 genes. ComBat was run with `mod = ~ group + sex`, which protects the two main
effects but **not their interaction**; any group-by-sex structure correlated with batch
is therefore partly absorbed and partly redistributed, and the downstream model is never
debited the degrees of freedom ComBat consumed (Nygaard et al., 2016). The 270-gene
count is inflated and is retained only as a sensitivity figure.

**Table 3.3a — Interaction patterns among the 53 genes.**

| Pattern | n | Meaning |
|---|---|---|
| Male-restricted | 35 | Significant in men, null in women |
| Magnitude difference | 12 | Same direction, significantly different size |
| Opposite direction | 4 | Significant in both, opposite signs |
| Female-restricted | 1 | Significant in women, null in men |
| Neither sex alone | 1 | Interaction significant, neither within-sex effect is |

The most interpretable result is ***PTPN22***, an established RA susceptibility gene,
up-regulated in both sexes but roughly five times more strongly in men (female log₂FC
+0.171, FDR 0.046; male +0.893, FDR 1.2 × 10⁻⁶; interaction FDR 0.038). The leading
male-restricted genes are *MAP4K5* (+1.053 in men, −0.144 in women), *BCLAF1*, *RIF1*
and *ZNF800*.

Four qualifications must accompany this result, and none of them is optional:

1. **Power.** Interaction tests carry roughly a quarter of the power of main-effect
   tests at the same sample size, and the male stratum contributes 17 RA cases.
   Fifty-three is a **floor**, not an estimate of how many sex-differential genes exist.
2. **The direction of the imbalance.** Thirty-five of the 53 are male-restricted despite
   men contributing a quarter of the samples. This is the same phenomenon described
   above for the DEG counts — male effect estimates are noisier and the median |log₂FC|
   among male DEGs is 1.8 times the female value. Some of these are large true effects
   and some are estimator variance; the class cannot be presented as a biological
   finding without that caveat.
3. **Cell composition.** Only **7 of the 53 survive** adjustment for leukocyte
   composition (Section 3.3a). Most of the apparent sex-differential signal is
   compositional rather than cell-intrinsic.
4. **This does not make the diagnostic panels sex-specific.** No panel gene appears
   among the 53. The panels are sex-STRATIFIED — fitted separately within each sex —
   and this analysis does not retrospectively confer sex-specificity on genes selected
   by a different procedure.

Subject to those four points, the interaction test provides the project's only direct
evidence that the RA blood transcriptome differs between the sexes, and it is positive.

### 3.3a Leukocyte composition and the composition-adjusted signature

Every sample analysed here is whole blood or PBMC, so the measured transcriptome is a
weighted average over circulating leukocyte populations. A gene can therefore appear
differentially expressed either because it is regulated differently *within* a cell type
in RA, or simply because the *proportion* of the cell type expressing it has changed.
The leading DEGs reported above — *ARG1*, *DEFA4*, *S100A8*, *CLEC4D*, *BCL2A1* — are
canonical neutrophil-granule and myeloid-activation transcripts, exactly what a shift in
composition would produce. Composition was therefore estimated with CIBERSORT and the
LM22 signature matrix, cross-checked against MCP-counter (`05c_deconvolution.R`).

Composition differs markedly by disease status within each sex, following the classic
lymphopenia pattern. In women the CD8 T-cell fraction falls from 0.162 in controls to
0.093 in RA (FDR 1.6 × 10⁻⁸), with smaller shifts in γδ T cells (+0.013, FDR 0.017) and
naive B cells (−0.012, FDR 0.047); in men the same CD8 depletion is present
(0.154 → 0.091, FDR 0.018) alongside γδ T-cell and eosinophil increases.

Refitting the within-sex contrasts with three composition principal components (centred
log-ratio space, 62 % of composition variance) added to the design reduces the DEG
counts substantially:

**Table 3.3b — Differential expression before and after composition adjustment.**

| Contrast | DEGs unadjusted | DEGs adjusted | Retained | Spearman ρ (log₂FC) |
|---|---|---|---|---|
| Female | 5,131 | **2,709** | 49.4 % | 0.889 |
| Male | 5,820 | **1,450** | 23.8 % | 0.842 |

Roughly half the female signature and three-quarters of the male signature is
attributable to differences in blood composition rather than to cell-intrinsic
regulation. The high rank correlation between adjusted and unadjusted fold changes
(ρ = 0.89 and 0.84) shows that adjustment attenuates the signature rather than
reorganising it.

**This adjustment must not be read as a bias correction.** Leukocyte composition in RA
is plausibly a *consequence* of the disease — a mediator on the path from RA to
neutrophilia and lymphopenia to the blood transcriptome. Adjusting for a mediator does
not produce a less-confounded estimate of the total disease effect; it removes part of
the real effect. The unadjusted model therefore estimates the **total** RA-associated
signal, which is the appropriate quantity for a **diagnostic** biomarker, and the
adjusted model estimates the signal **not explained by composition**, which is the only
basis on which a **cell-intrinsic or mechanistic** claim may be made about a gene. Both
are reported throughout, and no gene that fails adjustment is described mechanistically.

### 3.4 Functional annotation of the differentially expressed genes

Over-representation analysis returned **1,408 enriched terms for the female DEGs**
(1,295 GO, 113 KEGG), **1,513 for the male DEGs** (1,393 GO, 120 KEGG) and **1,659 for
the pooled contrast** (1,516 GO, 143 KEGG). The female and male term sets overlap
substantially (1,045 shared terms; 1,001 terms shared by all three contrasts of a
2,062-term union), reinforcing the conclusion of Section 3.3 that the two strata index
largely the same biology.

Leading GO biological processes were RNA splicing, Golgi vesicle transport,
ribonucleoprotein complex biogenesis, macroautophagy and regulation of autophagy in
women, and macroautophagy, Golgi vesicle transport, regulation of autophagy, regulation
of cellular catabolic process and proteasome-mediated ubiquitin-dependent protein
catabolism in men. Leading KEGG pathways were human T-cell leukaemia virus 1 infection,
Polycomb repressive complex, *Salmonella* infection, Epstein–Barr virus infection and
protein processing in the endoplasmic reticulum (female); and *Salmonella* infection,
protein processing in the endoplasmic reticulum, Polycomb repressive complex,
endocytosis and insulin signalling (male).

Ranked KEGG gene-set enrichment analysis, which uses the full ranked transcriptome
rather than a thresholded list, gave a strikingly directional result: **16 of 17
significant pathways in women and 21 of 22 in men were suppressed in RA**. The strongest
suppressed sets were primary immunodeficiency (NES −2.23), ATP-dependent chromatin
remodelling (−2.02) and spliceosome (−1.90) in women, and carbon metabolism (−2.15),
citrate cycle (−1.99), Fcγ-receptor-mediated phagosome formation (−1.95) and B-cell
receptor signalling (−1.84) in men. The single activated set in the female and pooled
analyses was **ribosome** (NES +1.89).

The suppression of spliceosome and chromatin-remodelling programmes anticipates the
co-expression result of Section 3.5, where the module most strongly *negatively*
correlated with RA is dominated by exactly these processes — two analytically
independent methods converging on the same biology.

### 3.5 Co-expression network construction

#### Quality control and sample outliers

All 15,763 genes and all 183 samples passed `goodSamplesGenes`. Average-linkage
clustering with a data-driven cut height of 69.32 (mean + 3 SD of merge heights)
identified **ten outlying arrays**, which were removed: GSM2449655, GSM2449671,
GSM2449732, GSM2449748, GSM2449876, GSM2981062, GSM2981092, GSM2981125, GSM2981129 and
GSM2981299.

This removal is independently corroborated. **Six of the ten** are also among the ten
lowest-weighted arrays under the limma REML array quality weights of Section 3.3 — an
overlap far beyond chance for a 10-of-183 selection — even though the two criteria are
computed from entirely different quantities (between-sample Euclidean distance versus
gene-wise residual variance). The removal is therefore a defensible quality decision
rather than an arbitrary trim.

The network was accordingly built on **173 samples × 15,763 genes**: 95 RA and 78 HC;
138 female and 35 male; 80 RA-female, 15 RA-male, 58 HC-female and 20 HC-male.

#### Gene filtering

No variance filter was applied and all **15,763 genes entered the network**. This is a
change from the earlier configuration, which discarded the bottom 40 % of genes by
variance. That filter removed 6,305 genes and, with them, **14.9 % of the female DEGs
and 20.9 % of the male DEGs** — genes that could never subsequently appear in a disease
module and therefore could never become candidates. Disabling the filter both removes an
arbitrary threshold and recovers those candidates.

#### Soft-thresholding power

The scale-free topology criterion was satisfied: `pickSoftThreshold` returned an
estimate of **β = 3**, the first power reaching R² ≥ 0.85 (R² = 0.873). That power was
rejected on connectivity grounds — mean connectivity at β = 3 is **2,240**, which
produces a network far too dense for meaningful module resolution. The signed R² curve
plateaus at approximately 0.90 from β = 9 onwards while mean connectivity falls steeply
(Table 3.3). **β = 12 was used**, where R² = 0.901, slope = −1.78 and mean connectivity
= 47.4. This value also coincides with the WGCNA FAQ's recommendation for a signed
network with more than 40 samples, so the connectivity-based and the authors'
sample-size criteria agree; the FAQ table is reported as corroboration only, since it is
formally conditional on the scale-free fit failing, which it did not.

**Table 3.3 — Scale-free fit and connectivity across candidate powers (abridged).**

| Power | Signed R² | Slope | Mean k | Median k | Max k |
|---|---|---|---|---|---|
| 3 | 0.873 | −10.26 | 2,240.3 | 2,167.0 | 3,059 |
| 6 | 0.861 | −3.00 | 446.6 | 375.1 | 1,111 |
| 9 | 0.897 | −2.11 | 125.0 | 78.0 | 645 |
| 10 | 0.903 | −1.96 | 87.8 | 47.8 | 558 |
| **12** | **0.901** | **−1.78** | **47.4** | **18.7** | **432** |
| 14 | 0.908 | −1.69 | 28.2 | 7.8 | 344 |
| 16 | 0.909 | −1.65 | 18.0 | 3.4 | 280 |
| 20 | 0.914 | −1.60 | 8.5 | 0.7 | 193 |

#### Module detection

Single-block construction at β = 12 detected **12 modules**, with 6,019 genes remaining
unassigned (grey). Module sizes were: turquoise 2,584; blue 2,554; brown 1,878; yellow
708; green 639; red 525; black 295; pink 172; magenta 127; purple 120; greenyellow 108;
tan 34.

#### Module–trait relationships and disease-module selection

Correlating each module eigengene with the binary RA indicator across all 173 samples
gave the associations in Table 3.4.

**Table 3.4 — Module–trait correlations (all samples, n = 173).**

| Module | Size | cor(ME, RA) | p (RA) | cor(ME, Male) | p (Male) | Direction |
|---|---|---|---|---|---|---|
| **yellow** | **708** | **+0.556** | **1.9 × 10⁻¹⁵** | −0.176 | 0.021 | **UP in RA** |
| turquoise | 2,584 | +0.369 | 6.1 × 10⁻⁷ | −0.254 | 7.3 × 10⁻⁴ | up |
| red | 525 | +0.289 | 1.2 × 10⁻⁴ | +0.047 | 0.538 | up |
| grey | 6,019 | +0.161 | 0.034 | −0.099 | 0.197 | (unassigned) |
| pink | 172 | +0.073 | 0.337 | −0.020 | 0.796 | up |
| tan | 34 | −0.097 | 0.204 | +0.090 | 0.241 | down |
| magenta | 127 | −0.156 | 0.041 | +0.013 | 0.866 | down |
| black | 295 | −0.212 | 0.005 | +0.140 | 0.067 | down |
| greenyellow | 108 | −0.216 | 0.004 | +0.087 | 0.257 | down |
| blue | 2,554 | −0.286 | 1.3 × 10⁻⁴ | +0.265 | 4.3 × 10⁻⁴ | down |
| purple | 120 | −0.317 | 2.1 × 10⁻⁵ | −0.039 | 0.611 | down |
| green | 639 | −0.362 | 9.6 × 10⁻⁷ | −0.093 | 0.226 | down |
| **brown** | **1,878** | **−0.590** | **1.3 × 10⁻¹⁷** | +0.190 | 0.012 | **DOWN in RA** |

Two modules met the pre-specified rule (|cor| ≥ 0.5, p < 1 × 10⁻⁸):

- the **yellow module** (708 genes, r = +0.556, p = 1.9 × 10⁻¹⁵), up-regulated in RA;
- the **brown module** (1,878 genes, r = −0.590, p = 1.3 × 10⁻¹⁷), down-regulated in RA.

Pooling them gives a disease-module background of **2,586 genes (16.4 % of the
transcriptome)**. The gap to the next-strongest module is substantial — turquoise at
+0.369 and green at −0.362 fall well short of both thresholds — so the selection is not
sensitive to small changes in the cut-off.

Neither disease module is a sex module: their correlations with the male indicator are
−0.176 (yellow) and +0.190 (brown), an order of magnitude weaker than their disease
correlations. The modules index disease status, not sex, which is a precondition for
using the same module background in both strata.

**A note on module identity.** The two disease modules are reported by their WGCNA
colour names, but those names are assigned by module size rank and are not stable
identifiers: an earlier configuration of the same pipeline, differing only in the
variance filter, produced the same two disease programmes under the names *green* and
*brown*. Modules are therefore defined throughout by their correlation with RA, and no
analysis step in this chapter refers to a colour literal.

#### Sex-stratified module–trait associations

Recomputing the RA-versus-control correlation within each sex confirms that both
disease modules are associated with RA inside each stratum
(`module_trait_RAvsControl_{all,female,male}.csv`; Figure
`fig_module_trait_disease_selection`, Figure `fig_disease_module_selection_composite`).
Consistent with the analysis plan, the female and male columns are not compared with one
another: the male stratum contains 35 samples and its correlations are correspondingly
unstable, so a larger coefficient in men is not evidence of a stronger male effect.

#### Hub genes

Applying the hub definition (|kME| > 0.8 and |GS| > 0.2) identified **217 hub genes in
the yellow module** and **144 in the brown module**.

The yellow (up-in-RA) hubs are dominated by signalling and trafficking regulators:
*ZNF267* (kME 0.946, GS +0.534, connectivity 117.9), *ZFYVE16* (0.953, +0.490, 114.8),
*SP3* (0.946, +0.520, 113.6), *FAR1* (0.951, +0.519, 111.9), *SNX13* (0.937, +0.526,
111.8), *SPOPL* (0.955, +0.505, 110.7), *PPP1R12A*, *RRM2B*, *ROCK1*, *UBE2W*, *CPEB4*
and *CPEB2*.

The brown (down-in-RA) hubs are dominated by RNA-processing and spliceosomal machinery:
*KHSRP* (kME 0.911, GS −0.564, connectivity 130.1), *GPI* (0.907, −0.510, 124.5),
*SF3B3* (0.894, −0.567, 122.0), *CLSTN1* (0.896, −0.556, 119.1), *XRCC6*, *EXOSC10*,
*PTBP1*, *SARS*, *HCFC1*, *DDX24*, *NUMA1* and *SF1*.

#### Functional enrichment of the disease modules

Against the 15,763-gene network background, the pooled disease-module gene set was
enriched for **218 GO terms** (105 biological process, 86 cellular component, 27
molecular function) and **16 KEGG pathways** at FDR < 0.05.

The dominant theme is post-transcriptional RNA metabolism. The leading biological
processes were mRNA processing (118 genes, FDR 1.5 × 10⁻⁷), ribonucleoprotein complex
biogenesis (111), RNA splicing via transesterification reactions (83), mRNA splicing via
spliceosome (82), RNA splicing (104), rRNA metabolic process (70) and protein
localisation to nucleus (85). Cellular components implicated both the spliceosome
(spliceosomal complex, 57 genes; U2-type spliceosomal complex, 33) and mitochondria
(organelle inner membrane, 119; mitochondrial inner membrane, 105; mitochondrial matrix,
109). Molecular functions were led by catalytic activity acting on nucleic acid (132)
and on RNA (87), transcription coregulator activity (117) and helicase activity (43).

The 16 enriched KEGG pathways were headed by **spliceosome** (43 genes, FDR
3.7 × 10⁻⁵) and **citrate cycle (TCA cycle)** (17, 1.9 × 10⁻⁴), followed by Polycomb
repressive complex (28), ATP-dependent chromatin remodelling (29), 2-oxocarboxylic acid
metabolism (15), nucleocytoplasmic transport (31), carbon metabolism (34), mRNA
surveillance (27), amyotrophic lateral sclerosis (75), steroid biosynthesis (9),
biosynthesis of amino acids (22), ubiquitin-mediated proteolysis (37), human T-cell
leukaemia virus 1 infection (54), Parkinson disease (55), **Th17 cell differentiation**
(30) and proteasome (15).

Three of these — spliceosome, ATP-dependent chromatin remodelling and carbon
metabolism/TCA cycle — are the same programmes that gene-set enrichment analysis found
*suppressed* in RA at the transcriptome-wide level (Section 3.4), and they map onto the
brown module, which is the module *negatively* correlated with RA. The agreement between
a thresholded co-expression analysis and a ranked whole-transcriptome test, computed by
different methods on different summaries of the data, is mutually corroborating. The
appearance of Th17 cell differentiation is notable given the established role of Th17
biology in RA pathogenesis.

#### Sex-stratified networks and module preservation

Independent networks were constructed on the 138 female and 35 male samples at the same
soft-threshold power. The male stratum is at the lower bound of what WGCNA's authors
consider analysable (a minimum of 15 samples, with 20 or more preferred), and any
difference in module count between the two networks may reflect sample size rather than
biology — which is precisely why the formal, permutation-based preservation statistic
rather than a module-count comparison is used to answer the question.

**[PENDING]** — `modulePreservation` (female reference → male test, 200 permutations)
was still executing at the time of writing. To be reported: the Z-summary and
preservation status of each module, with particular reference to the two disease
modules, together with the number of genes passing quality control in both sexes.
*Expectation from the earlier configuration of this pipeline: all substantial modules
were strongly preserved (Z-summary > 10), i.e. co-expression structure was not found to
differ between the sexes; this is to be re-confirmed on the unfiltered network.*

### 3.6 Candidate gene identification

**[PENDING]** — the intersection depends on the network step above.

The structure of the result and the inputs are fixed: the disease-module background is
**2,586 genes** (yellow 708 + brown 1,878) and the sex-stratified DEG lists are
**5,131 female** and **5,820 male** genes. To be reported once the run completes:

- the number of female candidates (disease modules ∩ female DEGs) and male candidates
  (disease modules ∩ male DEGs), each split by module and by direction of change;
- the number of candidates shared by both sexes and the union carried into Mendelian
  randomisation;
- the directional-consistency table — for each module and sex, the overlap with the
  up-in-RA and down-in-RA DEG lists, with the biologically consistent direction flagged
  (up for yellow, down for brown);
- the count of DEGs lost before the network to the variance filter, which is now **zero
  by construction**, the filter having been disabled (it was 765 female and 1,216 male
  genes, 14.9 % and 20.9 % of the respective lists, under the previous configuration).

Figures: `fig_venn_female_disease_candidates`, `fig_venn_male_disease_candidates`,
`fig_disease_module_selection_composite`, `fig_wgcna_16_candidates`, and the per-module
directional Venns `fig_diseasemod_venn_{module}_{sex}_{up,down}`.

---

# PART 3 — WHAT CHANGED FROM THE PREVIOUS DRAFT, AND WHY

These are the substantive corrections to the methodology text as previously written.
Each reflects a change in the code, not a change of wording.

| # | Previous text | Corrected text | Reason |
|---|---|---|---|
| 1 | "the bottom 40 % of genes by variance were eliminated" | No variance filter; all 15,763 genes enter the network | Any gene filtered out before the network can never join a disease module and so can never become a candidate. The filter cost 14.9 % of female and 20.9 % of male DEGs. |
| 2 | "sample clustering was examined for outliers" | Ten outlying arrays were **identified and removed**; network built on 173 samples | The previous script printed a suggested cut height but never applied it. Outliers now removed, and the removal is corroborated against limma array weights. |
| 3 | "the soft-thresholding power was … set to β = 12 in accordance with the WGCNA authors' recommendation for a signed network with more than 40 samples" | β = 12 selected on scale-free fit **and** connectivity; the FAQ table is corroboration only | The FAQ's sample-size table is formally conditional on the scale-free fit *failing*. Here the fit succeeds (estimate β = 3), so the table cannot be the primary justification. The real argument is that β = 3 gives mean connectivity 2,240 whereas β = 12 gives 47.4 at equivalent R². |
| 4 | "maxBlockSize = 16000" | `maxBlockSize = 20000` | All 15,763 genes must fit in one block; the parameter was raised with the filter removed. |
| 5 | "modulePreservation … 100 permutations" | 200 permutations | 200 is the publication standard and is what the code now runs. |
| 6 | (not stated) | The preservation reference uses the **combined-network** module assignment, not the female-network modules | Otherwise the statistic tests female-network modules, which are a different gene set from the disease modules even where colour names coincide. |
| 7 | "The strongest 2 % of intramodular co-expression edges were used to create complementary betweenness-centrality hubs. When network access was permitted, a STRING protein–protein interaction network … was rebuilt" | **Removed** | Neither analysis exists in the current pipeline. Hub genes are defined by connectivity, kME and GS only. Re-add the text only if the STRING/betweenness steps are reinstated. |
| 8 | "the identical contrast was additionally refitted within each source series separately, and the per-series results are reported alongside the pooled analysis" | **Removed** | No per-series refit exists in `05_dge.R`. Either delete the claim or implement the refit — it cannot stand as written. |
| 9 | "McCarthy & Smyth … evaluated as sensitivity check" (no detail) | Explicit `treat` sensitivity analysis at τ = 0.1, log₂(1.2), 0.5 | The `treat` block now genuinely exists in `05_dge.R`, making the citation truthful. **It has not yet been run** — see Part 4. |
| 10 | Disease modules referred to as "green (up in RA) + brown (down in RA)" | Disease modules are **yellow (up) + brown (down)**, selected by rule, with an explicit statement that colour names are unstable | `labels2colors` assigns colours by size rank. The same biology was named green+brown under the filtered configuration and yellow+brown under the unfiltered one. |
| 11 | (not stated) | DEGs were fitted on 183 samples, the network on 173 | An examiner will notice the discrepancy. Stated openly with its justification. |
| 12 | (not stated) | Pearson used rather than `bicor`, with the memory constraint stated and sample-level outlier removal named as the substitute | The FAQ recommends `bicor`; departing from it silently would be a weakness. Stating the computational reason and the compensating measure is stronger than omitting it. |

---

# PART 4 — OUTSTANDING ITEMS BEFORE THIS SECTION IS FINAL

*Updated 2026-07-28. Items 1–3 of the previous list are RESOLVED; the tables were
regenerated on 28 July at 09:43–09:49 against the unfiltered network and the current
`05_dge.R`. `DEG_treat_sensitivity.csv` now exists. The disease modules are **yellow +
brown**, not green + brown.*

## Resolved since the previous draft

1. ~~`06_WGCNA.R` must finish~~ — complete; `WGCNA_10_module_preservation.csv`,
   `WGCNA_11_candidates_{female,male}.csv` and `WGCNA_12_candidate_summary.csv` are
   regenerated for the unfiltered network.
2. ~~Re-run 08, 09, 09b~~ — done. Disease modules are yellow (708 genes, cor 0.556,
   p = 1.9 × 10⁻¹⁵) and brown (1,878 genes, cor −0.590, p = 1.3 × 10⁻¹⁷); 2,045 female
   and 2,079 male candidates.
3. ~~Re-run `05_dge.R`~~ — done; the `treat` sensitivity is reportable.

## Added since the previous draft

4. **Sections 3.3 (interaction) and 3.3a (composition) are new** and are backed by
   `05d_interaction_report.R` and `05c_deconvolution.R`. The consolidated write-up of
   all five robustness analyses is `results/RESULTS_ROBUSTNESS.md`, which also covers
   the MHC sensitivity (`10c`), colocalisation (`10d`) and the panel-versus-composition
   benchmark (`13b`) — those belong to the MR and panel sections, which this sub-chapter
   does not yet reach.

## Still outstanding

5. **The causal language throughout the MR and panel sections must change.**
   Colocalisation (`10d`) shows that **no panel gene shares a causal variant with the RA
   association**, and six of nine show positive evidence of *distinct* causal variants
   (PP.H3 ≥ 0.8). "MR-causal" must become "MR-prioritised" wherever it appears. The MR
   step remains defensible as a genetically-informed filter on the candidate space; it
   is no longer evidence of causality.

6. **The male arm must be described as exploratory** in every table, figure and
   sentence. Train n = 38, internal n = 13, external n = 9. Its internal AUC of
   1.000 (1.000–1.000) is a perfect-separation artefact and is now flagged `SEPARATION`
   in `mr_final_panel_summary.csv`.

7. **Decide on the per-series heterogeneity analysis.** The previous draft claimed it;
   the code does not perform it. Either implement a per-series refit of the RA-vs-HC
   contrast in `05_dge.R` (worth doing — it directly addresses the differing origin
   normalisations noted in Section 2.1) or leave the claim deleted.

8. **Reconsider the |log₂FC| > 0.1 threshold.** `DEG_treat_sensitivity.csv` shows that
   formal TREAT testing at 1.2× reduces the female DEG count from 5,131 to 86, and at
   1.41× to 3. With 6,422 of 15,763 genes significant in the pooled contrast, the
   subsequent disease-module ∩ DEG intersection retains 2,045 of 2,586 module genes
   (79 %) and is therefore not doing the filtering work the Methods implies it does.
   This should be stated explicitly or the threshold raised.

9. **Figure numbering.** The Results text refers to figures by filename. These should be
   mapped to thesis figure numbers once the section order is settled. The mapping was
   maintained in `results/FIGURE_PROVENANCE.md`, **which does not exist in this tree** —
   regenerate it or remove every reference to it.
