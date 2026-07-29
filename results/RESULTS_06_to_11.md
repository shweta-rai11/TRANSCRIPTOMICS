# Results: co-expression network analysis, candidate selection and causal screening

Analysis run 27 July 2026. WGCNA cache key `00048f2e`. All figures and tables
referenced below were regenerated from this single parameter set; earlier
outputs keyed to `00018d30` have been removed.

---

## 1. Network input and quality control

Batch-corrected expression for the 70% training split comprised **15,763 genes
across 183 samples** (RA = 103, HC = 80; female = 145, male = 38). The held-out
30% was not seen at any point in this analysis.

`goodSamplesGenes` flagged no failing genes or samples. Hierarchical clustering
of samples on Euclidean distance identified **10 outlying arrays** at a cut
height of 69.32 (mean + 3 SD of merge heights), all of which were removed:

> GSM2449671, GSM2449748, GSM2449655, GSM2449732, GSM2449876, GSM2981062,
> GSM2981092, GSM2981125, GSM2981129, GSM2981299

The analysed matrix was therefore **173 samples × 15,763 genes** (RA = 95,
HC = 78; female = 138, male = 35).

**No variance filter was applied.** This is a deliberate departure from common
practice. Any gene removed before network construction can never be assigned to
a disease module and therefore can never enter the `disease-module ∩ sex-DEG`
intersection that defines the candidate set. A 40% variance filter, evaluated
previously, discarded 6,305 genes and with them 14.9% of female and 20.9% of
male DEGs. Retaining all genes eliminates an arbitrary threshold and recovers
those candidates; the cost is computational, and 38% of genes remain unassigned
(grey), which is reported transparently below.

*Figures: `fig_wgcna_01_gene_variance.png`, `fig_wgcna_02_sample_clustering.png`*

---

## 2. Soft-thresholding power

A signed network was constructed with Pearson correlation. Scale-free topology
fit was assessed across powers 1–20.

| Power | Signed R² | Mean connectivity |
|---|---|---|
| 3 (`pickSoftThreshold` estimate) | 0.873 | 2240.3 |
| 6 | 0.861 | 446.6 |
| 10 | 0.903 | 87.8 |
| **12 (used)** | **0.901** | **47.4** |
| 14 | 0.908 | 28.2 |
| 20 | 0.914 | 8.5 |

`pickSoftThreshold` returns the *first* power clearing the R² ≥ 0.85 criterion,
which here was 3. At that power the network is far too dense — mean connectivity
2,240 — producing modules with little biological coherence. The scale-free fit
plateaus from power 9 onward at R² ≈ 0.90, and **power 12** sits on that plateau
with a well-conditioned mean connectivity of 47.4. Power 12 is also the value
recommended by the WGCNA FAQ for signed networks with n > 40. Both criteria
agree, and power 12 was used.

*Figure: `fig_wgcna_03_soft_threshold.png`; Table: `WGCNA_01_soft_threshold.csv`*

---

## 3. Module detection

Blockwise network construction (single block, `maxBlockSize` = 20,000; minimum
module size 30; merge cut height 0.25) yielded **12 modules** plus 6,019
unassigned (grey) genes.

| Module | Size | Module | Size |
|---|---|---|---|
| turquoise | 2,584 | red | 525 |
| blue | 2,554 | black | 295 |
| brown | 1,878 | pink | 172 |
| yellow | 708 | magenta | 127 |
| green | 639 | purple | 120 |
| greenyellow | 108 | tan | 34 |
| *(grey, unassigned)* | *6,019* | | |

The high grey fraction (38%) is the direct consequence of retaining all genes;
it is stated as a limitation rather than concealed by a filter.

*Figure: `fig_wgcna_04_dendrogram.png`*

---

## 4. Module–trait relationships and disease-module selection

Module eigengenes were correlated with RA status across all 173 samples.

| Module | Size | cor(ME, RA) | p | Direction |
|---|---|---|---|---|
| **yellow** | 708 | **+0.556** | 1.91e-15 | UP in RA |
| turquoise | 2,584 | +0.369 | 6.09e-07 | UP in RA |
| red | 525 | +0.289 | 1.18e-04 | UP in RA |
| grey | 6,019 | +0.161 | 3.42e-02 | — |
| pink | 172 | +0.073 | 0.337 | ns |
| tan | 34 | −0.097 | 0.204 | ns |
| magenta | 127 | −0.156 | 4.06e-02 | DOWN in RA |
| black | 295 | −0.212 | 5.12e-03 | DOWN in RA |
| greenyellow | 108 | −0.216 | 4.23e-03 | DOWN in RA |
| blue | 2,554 | −0.286 | 1.34e-04 | DOWN in RA |
| purple | 120 | −0.317 | 2.10e-05 | DOWN in RA |
| green | 639 | −0.362 | 9.59e-07 | DOWN in RA |
| **brown** | 1,878 | **−0.590** | 1.29e-17 | DOWN in RA |

Disease modules were selected by a **pre-specified, data-driven rule**
(|cor(ME, RA)| ≥ 0.50 and p < 1e-8), never by colour name. WGCNA assigns colour
labels by module size rank, so colours are not stable identifiers across runs;
they are recorded here as output, not used as input anywhere in the pipeline.

Two modules met the rule:

- **yellow** — 708 genes, r = +0.556, p = 1.91e-15, up in RA
- **brown** — 1,878 genes, r = −0.590, p = 1.29e-17, down in RA

The next-strongest module (turquoise, r = +0.369) falls well below the
threshold, so the selection is unambiguous. Together the two modules contribute
**2,586 disease-module genes**.

*Figures: `fig_wgcna_05_module_trait.png`, `fig_module_trait_disease_selection.png`;
Tables: `WGCNA_02_module_trait.csv`, `WGCNA_03_disease_modules.csv`*

### 4.1 Within-sex module–trait association

The same contrast computed separately within each sex:

| Module | All (n = 173) | Female (n = 138) | Male (n = 35) |
|---|---|---|---|
| yellow | +0.56*** | +0.48*** | +0.77*** |
| brown | −0.59*** | −0.52*** | −0.81*** |

Both modules retain their direction and significance in each stratum. **The male
and female columns must not be compared directly.** The male stratum is roughly
one quarter the size of the female stratum, so its correlations are estimated
with far greater uncertainty and its point estimates are correspondingly more
extreme. Concluding a sex difference from the gap between +0.77 and +0.48 would
be the error described by Gelman & Stern (2006); a formal interaction test would
be required, and is not performed here.

*Figures: `fig_wgcna_module_trait_RAvsControl_{all,female,male}.png`;
Table: `module_trait_RAvsControl_ALLSTRATA.csv`*

---

## 5. Hub genes

Hub genes were defined within each disease module as |kME| > 0.80 and
|GS_RA| > 0.20.

| Module | Genes | Hubs |
|---|---|---|
| yellow | 708 | 217 |
| brown | 1,878 | 144 |

Highest-connectivity hubs:

- **yellow** — SPOPL, ZFYVE16, FAR1, ZNF267, SP3, RRM2B, CPEB4, SNX13, MTMR6,
  UBE2W, EXOC8, SLMAP
- **brown** — KHSRP, GPI, EXOSC10, PTBP1, CLSTN1, SCAMP3, SF3B3, RNPS1, ANAPC5,
  SARS, NUMA1, NUDCD3

*Figures: `fig_wgcna_08_kME_vs_GS.png`, `fig_wgcna_09_hub_genes.png`,
`fig_wgcna_10_disease_heatmap.png`; Tables: `WGCNA_06_disease_module_hubs.csv`,
`WGCNA_07_hub_genes_only.csv`*

---

## 6. Functional enrichment

Over-representation analysis of the 2,586 disease-module genes against a
background of all 15,763 network genes (the correct universe: those genes that
*could* have been assigned to a module) returned **218 GO terms** and
**16 KEGG pathways** at BH-adjusted p < 0.05.

**Top GO terms**

| Term | Gene ratio | p.adj |
|---|---|---|
| mRNA processing | 118/2342 | 1.52e-07 |
| ribonucleoprotein complex biogenesis | 111/2342 | 1.52e-07 |
| RNA splicing, via transesterification | 83/2342 | 1.81e-07 |
| mRNA splicing, via spliceosome | 82/2342 | 1.81e-07 |
| spliceosomal complex | 57/2418 | 6.38e-06 |
| mitochondrial inner membrane | 105/2418 | 1.93e-05 |

**Top KEGG pathways**

| Pathway | Gene ratio | p.adj |
|---|---|---|
| Spliceosome | 43/1224 | 3.72e-05 |
| Citrate cycle (TCA cycle) | 17/1224 | 1.93e-04 |
| Polycomb repressive complex | 28/1224 | 5.91e-03 |
| ATP-dependent chromatin remodeling | 29/1224 | 6.06e-03 |
| Nucleocytoplasmic transport | 31/1224 | 6.06e-03 |
| Ubiquitin mediated proteolysis | 37/1224 | 3.36e-02 |

The enrichment signal is dominated by RNA processing and spliceosome biology,
driven by the larger brown module. The yellow module's gene content is by
contrast strongly myeloid — ARG1, CLEC4D, BCL2A1, LY96 (the TLR4 co-receptor
MD-2), MS4A4A, MS4A3, ANXA3 — consistent with the innate immune activation
expected in RA whole blood. The two disease modules therefore capture
biologically distinct axes: an up-regulated innate/myeloid programme and a
down-regulated RNA-processing programme.

*Figure: `fig_wgcna_12_enrichment.png`; Tables: `WGCNA_08_disease_GO.csv`,
`WGCNA_09_disease_KEGG.csv`*

---

## 7. Sex-specific networks and module preservation

Independent networks were constructed within each sex at the same soft power,
yielding **15 modules in females** (n = 138) and **11 in males** (n = 35). All
15,763 genes passed QC in both strata.

Preservation of the combined-network modules was assessed from female
(reference) to male (test) with 200 permutations. Reference module labels were
taken from the **combined** network, not the female-specific network, so that
the statistic answers the question actually being asked — whether the yellow and
brown *disease* modules hold in both sexes.

| Module | Zsummary | Status |
|---|---|---|
| turquoise | 51.60 | Strong |
| **yellow** | **48.14** | Strong |
| blue | 44.43 | Strong |
| black | 38.38 | Strong |
| **brown** | **36.95** | Strong |
| *gold (random benchmark)* | *34.19* | — |
| red | 27.24 | Strong |
| green | 26.00 | Strong |
| magenta | 23.03 | Strong |
| purple | 21.02 | Strong |
| pink | 16.42 | Strong |
| greenyellow | 11.84 | Strong |
| tan | 8.63 | Moderate |
| grey | 4.61 | Moderate |

**Interpretation requires care on two points.**

First, `gold` is WGCNA's synthetic benchmark module composed of randomly
selected genes; it indexes the preservation expected in the absence of genuine
module structure. At Z = 34.19 the background is itself highly preserved.
Against that reference, yellow (48.14) and turquoise (51.60) are clearly
preserved, whereas **brown (36.95) exceeds the random benchmark only
marginally**. The conventional Z > 10 "strong" label overstates the evidence for
brown, and the appropriate comparison is against gold rather than against the
fixed cut-off.

Second, the `moduleSize` reported by `modulePreservation` is capped at 1,000
genes; large modules are subsampled. Brown's statistic is therefore computed on
a 1,000-gene subset of its 1,878 genes, and the tabulated sizes for brown, blue,
turquoise and grey should not be quoted as module sizes.

*Figures: `fig_wgcna_13_dendro_female.png`, `fig_wgcna_14_dendro_male.png`,
`fig_wgcna_15_preservation.png`; Table: `WGCNA_10_module_preservation.csv`*

---

## 8. Candidate genes: disease module ∩ sex-stratified DEG

Candidates were defined as genes both (i) assigned to a disease module and
(ii) differentially expressed between RA and HC *within* the sex in question.

| Sex | Significant DEGs | Candidates | Up (yellow) | Down (brown) |
|---|---|---|---|---|
| Female | 5,131 | **2,045** | 595 | 1,450 |
| Male | 5,820 | **2,079** | 585 | 1,494 |

- Shared by both sexes: **1,773**
- Female-specific: **272**
- Male-specific: **306**
- Union: **2,351**

Because no variance filter was applied, **zero DEGs were lost before network
construction** in either sex — the design objective stated in §1 is met exactly.

### 8.1 Directional consistency

Every candidate's direction of differential expression agrees with the direction
of its module's association with RA: all 595/585 yellow candidates are up in RA,
all 1,450/1,494 brown candidates are down, and the inconsistent overlap is
**zero in both sexes**. This is a coherence check that the assignment passes,
but it is a weak one: in a *signed* network all members of a module correlate
positively with the module eigengene by construction, so near-perfect
directional agreement is close to guaranteed and should not be presented as
independent validation.

### 8.2 Leading candidates

Ranked by |log2FC|, the strongest candidates in both sexes are predominantly
yellow-module myeloid genes:

| Gene | Module | logFC (F) | logFC (M) |
|---|---|---|---|
| BCL2A1 | yellow | +1.16 | +1.84 |
| LY96 | yellow | +0.92 | +1.55 |
| ARG1 | yellow | +1.36 | +1.43 |
| CLEC4D | yellow | +1.17 | +1.35 |
| ZNF267 | yellow | +0.79 | +1.42 |
| MS4A3 | yellow | +0.82 | +1.44 |

The male fold-changes are uniformly larger than the female ones. **This is not
evidence of a stronger effect in men.** With only 15 male RA cases and 20 male
controls, effect sizes among genes passing a significance threshold are subject
to the well-known upward bias of small samples ("winner's curse"). The ratio of
male to female log-fold-change carries no interpretable meaning here.

*Figures: `fig_wgcna_16_candidates.png`, `fig_venn_{female,male}_disease_candidates.png`,
`fig_diseasemod_venn_{yellow,brown}_{female,male}_{up,down}.png`;
Tables: `WGCNA_11_candidates_{female,male}.csv`, `candidates_{female,male}_disease.csv`,
`diseasemod_DEG_direction_summary.csv`*

### 8.3 Limitation of the candidate step

The DEG threshold inherited from the differential expression analysis is
|log2FC| > 0.10 with FDR < 0.05 — a 7% fold change, which calls 33–37% of the
transcriptome differentially expressed. As a consequence the intersection
retains **2,045 of 2,586 disease-module genes (79%) in females** and a similar
proportion in males. The step therefore performs little prioritisation, and the
large candidate set imposes a heavy multiple-testing burden on the causal
screen that follows. A more stringent effect-size threshold would materially
sharpen this stage.

---

## 9. Causal screening by two-sample Mendelian randomisation

### 9.1 Design

The candidate genes above are defined by *association*. Two-sample MR was used
to ask which of them show evidence of a **causal** effect on RA risk, using
germline genetic variants as instrumental variables. Because alleles are
allocated at conception and precede disease onset, they are not subject to the
reverse causation or confounding (by medication, smoking or inflammation itself)
that limits observational transcriptomics.

- **Exposure** — *cis*-eQTLs from eQTLGen (Võsa et al. 2021; ~31,000 whole-blood
  samples), accessed as `eqtl-a-<ENSG>`, p < 5e-8, LD-clumped. Cis instruments
  were used in preference to trans because their mechanism is direct, reducing
  violation of the exclusion-restriction assumption.
- **Outcome** — rheumatoid arthritis, Okada et al. (2014) European ancestry
  GWAS (`ieu-a-832`; 14,361 cases / 43,923 controls), fetched in 250-SNP
  batches with LD proxies.
- **Instrument strength** — SNPs retained only at F = (β/SE)² ≥ 10.
- **Harmonisation** — effect alleles aligned with `action = 2`, inferring strand
  for palindromic variants from allele frequency and discarding ambiguous ones.
- **Estimators** — chosen by instrument count: ≥3 SNPs, inverse-variance
  weighted (IVW) with MR-Egger and weighted median as sensitivity analyses;
  2 SNPs, IVW; 1 SNP, Wald ratio. Primary estimate ordered
  IVW > Wald > weighted median > MR-Egger.
- **Pleiotropy and heterogeneity** — Cochran's Q and the MR-Egger intercept
  where ≥3 instruments permitted.
- **Multiple testing** — Benjamini–Hochberg FDR computed **within each sex
  stratum**, the denominator matching the stratified design. A pooled correction
  across the union of both sexes is also retained for backward compatibility but
  is not the column reported.

### 9.2 Instruments

Of the 2,351 union candidates, **1,980 (84%) had at least one genome-wide
significant cis-eQTL surviving clumping and the F ≥ 10 filter**, yielding
**4,932 instrument SNPs**. The remaining 371 genes are not testable by this
design, as no adequately strong cis instrument exists for them in eQTLGen.

### 9.3 Results

After harmonisation against the Okada outcome (3,363 of 4,932 instrument SNPs
retrieved), MR estimates were obtained for **1,490 genes in the female candidate
set and 1,493 in the male set**.

| | Female | Male |
|---|---|---|
| Candidates | 2,045 | 2,079 |
| MR-tested | 1,490 | 1,493 |
| Nominal p < 0.05 | 113 | 115 |
| *Expected false positives at p < 0.05* | *74.5* | *74.7* |
| **Surviving within-stratum BH-FDR < 0.05** | **33** | **27** |

Of the FDR-surviving genes, **25 are shared**, 8 are female-only
(FCGR2B, GNL1, HLA-DPA1, IKZF4, LAX1, LSM2, RETSAT, ZBTB9) and 2 are male-only
(CISH, TAB1).

**Strongest associations (female stratum, ordered by p):**

| Gene | OR [95% CI] | p | FDR | nSNP | Method |
|---|---|---|---|---|---|
| HLA-DRB1 | 0.34 [0.32–0.36] | 1.0e-250 | 1.5e-247 | 1 | Wald ratio |
| HNRNPM | 287.3 [195.9–421.4] | 1.4e-184 | 1.0e-181 | 1 | Wald ratio |
| FOXP3 | 263.9 [173.4–401.6] | 3.1e-149 | 1.5e-146 | 1 | Wald ratio |
| WDR46 | 4.53 [3.64–5.64] | 5.2e-42 | 1.9e-39 | 1 | Wald ratio |
| AIF1 | 0.56 [0.50–0.62] | 2.5e-24 | 7.5e-22 | 2 | IVW |
| VPS52 | 0.41 [0.35–0.49] | 4.2e-24 | 1.0e-21 | 1 | Wald ratio |
| GNL1 | 4.63 [3.26–6.57] | 9.0e-18 | 1.9e-15 | 1 | Wald ratio |
| VARS2 | 1.20 [1.15–1.26] | 1.1e-16 | 2.0e-14 | 3 | IVW |
| AP4B1 | 1.42 [1.29–1.57] | 1.3e-12 | 2.2e-10 | 1 | Wald ratio |
| CSNK2B | 0.55 [0.46–0.65] | 1.3e-11 | 1.9e-09 | 1 | Wald ratio |

*Tables: `MR_{female,male}_primary_okada.csv`,
`MR_causal_FDR_{female,male}.csv`, `FS_input_{female,male}.csv`,
`MR_{female,male}_TABLE1-4.csv`, `MR_{sex}_all_tables.xlsx`*

### 9.4 Critical limitations of the MR results

Four features of these results materially constrain what may be concluded from
them, and each should be stated in the thesis rather than left for a reader to
discover.

**(a) The nominal screen is close to uninformative.** At p < 0.05 the female
stratum returns 113 genes against an expectation of **74.5 by chance alone**;
the male stratum returns 115 against 74.7. Roughly two thirds of the nominally
screened set is therefore noise. Because `FS_input_{female,male}.csv` — the
input to downstream feature selection — is populated from the *nominal* screen
rather than the FDR-surviving set, the feature-selection stage inherits a list
that is predominantly false positives. Only the 33 female and 27 male
FDR-surviving genes should be described as causal candidates.

**(b) The signal is dominated by the MHC, and the genes there are almost
certainly not independent.** A large fraction of the FDR-surviving set lies in
the chromosome 6p21 MHC region: HLA-DRB1, HLA-DMA, HLA-DPA1, BTN3A3, C6orf136,
CSNK2B, ZBTB9, GNL1, LSM2, VARS2, WDR46, AIF1. HLA-DRB1 is the canonical RA
susceptibility locus (the shared-epitope alleles), and its recovery as the
single strongest signal (p = 1e-250) is a reassuring positive control. However,
linkage disequilibrium across the MHC is exceptionally long-range and complex,
so a *cis*-eQTL for any gene in the region will tag the causal HLA haplotype.
These twelve genes are best interpreted as **one signal observed twelve times**,
not twelve independent causal genes. Any claim about a specific MHC gene
requires conditional or colocalisation analysis that has not been performed
here.

**(c) Several odds ratios are not biologically interpretable.** HNRNPM
(OR = 287) and FOXP3 (OR = 264) are implausible as causal effect sizes. Both are
single-instrument Wald-ratio estimates, which are the ratio β_outcome /
β_exposure and therefore inflate without bound as the denominator approaches
zero — a weak-instrument artefact, notwithstanding that both passed the F ≥ 10
filter. These estimates should be reported as direction-of-effect only, or
excluded, and must not be quoted as effect magnitudes.

**(d) Most estimates rest on a single instrument.** Eight of the ten strongest
female associations derive from one SNP, which permits neither the Cochran's Q
heterogeneity test nor the MR-Egger intercept test. For those genes no
pleiotropy assessment is possible, and the core MR assumption of exclusion
restriction is untested rather than satisfied. Only AIF1 (2 SNPs), HLA-DMA
(2 SNPs) and VARS2 (3 SNPs) among the leading hits carry any sensitivity
analysis at all.

Taken together, the defensible statement is that the causal screen recovers the
established HLA-DRB1 association and nominates a small number of non-MHC
candidates — among them FOXP3, IKZF3, CD74, SH2B1, SREBF2, KLF2 and CDC37 —
whose individual effect estimates remain provisional pending multi-instrument
replication.

### 9.5 Interpretation limit: the MR step is not sex-stratified

This must be stated explicitly. Both MR inputs — the eQTLGen exposure and the
Okada outcome — are **sex-combined**. A survey of all 37 rheumatoid arthritis
datasets available in OpenGWAS confirmed that **no sex-stratified RA GWAS
exists**. Consequently the MR estimate for any given gene is identical whether
that gene is reported in the female or the male table; the two tables differ
only in *which* genes appear (272 female-specific, 306 male-specific) and in the
FDR denominator applied within each stratum.

The design is therefore accurately described as **sex-stratified candidate
discovery followed by sex-combined causal validation**. Sex-specificity in this
work derives entirely from the within-sex differential expression analysis, not
from the causal screen. No claim of a sex-specific causal effect is made or
could be supported by these data. Establishing one would require sex-stratified
outcome summary statistics and a formal test of the difference between strata,
z = (β_F − β_M) / √(SE_F² + SE_M²), rather than the comparison of separate
significance in each sex.

---

## References

- Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network
  analysis. *BMC Bioinformatics* 2008;9:559.
- Zhang B, Horvath S. A general framework for weighted gene co-expression
  network analysis. *Stat Appl Genet Mol Biol* 2005;4:Article 17.
- Langfelder P, Luo R, Oldham MC, Horvath S. Is my network module preserved and
  reproducible? *PLoS Comput Biol* 2011;7(1):e1001057.
- Langfelder P, Zhang B, Horvath S. Defining clusters from a hierarchical
  cluster tree: the Dynamic Tree Cut library. *Bioinformatics* 2008;24:719–720.
- Yu G, et al. clusterProfiler 4.0. *Innovation (Camb)* 2021;2(3):100141.
- Võsa U, et al. Large-scale cis- and trans-eQTL analyses identify thousands of
  genetic loci and polygenic scores that regulate blood gene expression.
  *Nat Genet* 2021;53:1300–1310.
- Okada Y, et al. Genetics of rheumatoid arthritis contributes to biology and
  drug discovery. *Nature* 2014;506:376–381.
- Gelman A, Stern H. The difference between "significant" and "not significant"
  is not itself statistically significant. *Am Stat* 2006;60(4):328–331.
- Burgess S, Butterworth A, Thompson SG. Mendelian randomization analysis with
  multiple genetic variants using summarized data. *Genet Epidemiol*
  2013;37:658–665.
- Bowden J, Davey Smith G, Burgess S. Mendelian randomization with invalid
  instruments: effect estimation and bias detection through Egger regression.
  *Int J Epidemiol* 2015;44:512–525.
