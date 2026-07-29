# References for the pipeline (methods, tools, datasets)

Citations for every method, R/Bioconductor package, and dataset used across the
scripts, organised by pipeline step. Package versions are printed at the end of
individual scripts (`sessionInfo`). **⚠️ = verify the primary publication /
accession against GEO before citing** (dataset papers not fabricated here).

---

## Data acquisition & preprocessing (01–04, 27b)
- **GEOquery** — Davis S, Meltzer PS. GEOquery: a bridge between the Gene
  Expression Omnibus (GEO) and BioConductor. *Bioinformatics* 2007;23(14):1846–7.
- **Quantile normalisation (preprocessCore)** — Bolstad BM, Irizarry RA,
  Åstrand M, Speed TP. A comparison of normalization methods for high density
  oligonucleotide array data. *Bioinformatics* 2003;19(2):185–93.
- **ComBat batch correction (sva)** — Johnson WE, Li C, Rabinovic A. Adjusting
  batch effects in microarray expression data using empirical Bayes methods.
  *Biostatistics* 2007;8(1):118–27. Leek JT, et al. The sva package.
  *Bioinformatics* 2012;28(6):882–3.
- **edgeR / TMM normalisation (RNA-seq synovium)** — Robinson MD, McCarthy DJ,
  Smyth GK. edgeR. *Bioinformatics* 2010;26(1):139–40.
- **voom** — Law CW, Chen Y, Shi W, Smyth GK. voom: precision weights unlock
  linear model analysis tools for RNA-seq read counts. *Genome Biol* 2014;15:R29.

## Differential expression (05)
- **limma** — Ritchie ME, et al. limma powers differential expression analyses
  for RNA-sequencing and microarray studies. *Nucleic Acids Res* 2015;43(7):e47.
- **Empirical array quality weights** — Ritchie ME, Diyagama D, Neilson J, et al.
  Empirical array quality weights in the analysis of microarray data.
  *BMC Bioinformatics* 2006;7:261.

## Co-expression network (06, 08)
- **Benchmark protocol (WGCNA + module preservation + enrichment)** — Nguyen P,
  Zeng E. A protocol for weighted gene co-expression network analysis with module
  preservation and functional enrichment analysis for tumor and normal
  transcriptomic data. *Bio Protoc* 2025;15(18):e5447.
  doi:10.21769/BioProtoc.5447. PMID: 41000162; PMCID: PMC12457846.
  *(Cited in Section 2.3 as the published protocol against which parameter choices
  are benchmarked. Agrees: signed network, deepSplit = 2, minModuleSize = 30,
  mergeCutHeight = 0.25, joint R²/connectivity power criterion, Z-summary bands.
  Diverges: gene filtering, permutation count, preservation direction — each
  divergence justified in text.)*
- **WGCNA framework** — Zhang B, Horvath S. A general framework for weighted gene
  co-expression network analysis. *Stat Appl Genet Mol Biol* 2005;4:Article 17.
- **WGCNA** — Langfelder P, Horvath S. WGCNA: an R package for weighted
  correlation network analysis. *BMC Bioinformatics* 2008;9:559.
  *(Also the source cited in text for the package FAQ and tutorials; if the
  department's style requires the online resources to be listed separately, add:
  Langfelder P, Horvath S. WGCNA package FAQ / Tutorials for the WGCNA package.
  University of California, Los Angeles. Available at:
  https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
  [accessed DD Month YYYY].)*
- Module–trait correlation / eigengenes — Langfelder P, Horvath S. Eigengene
  networks for studying the relationships between co-expression modules.
  *BMC Syst Biol* 2007;1:54.
- **Dynamic Tree Cut (module detection)** — Langfelder P, Zhang B, Horvath S.
  Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package
  for R. *Bioinformatics* 2008;24(5):719–20.
- **Module preservation** — Langfelder P, Luo R, Oldham MC, Horvath S. Is my
  network module preserved and reproducible? *PLoS Comput Biol* 2011;7(1):e1001057.
- **Hierarchical clustering (Ward's method, cited as the rejected alternative)** —
  Ward JH Jr. Hierarchical grouping to optimize an objective function.
  *J Am Stat Assoc* 1963;58(301):236–44.
- **Standardised-connectivity sample outlier detection (Z.k)** — Oldham MC,
  Langfelder P, Horvath S. Network methods for describing sample relationships in
  genomic datasets: application to Huntington's disease. *BMC Syst Biol* 2012;6:63.
- **arrayQualityMetrics (rejected alternative)** — Kauffmann A, Gentleman R,
  Huber W. arrayQualityMetrics — a bioconductor package for quality assessment of
  microarray data. *Bioinformatics* 2009;25(3):415–6.
- **RLE / NUSE array quality statistics (rejected alternative)** —
  Brettschneider J, Collin F, Bolstad BM, Speed TP. Quality assessment for short
  oligonucleotide microarray data. *Technometrics* 2008;50(3):241–64.

## Candidate integration (09)
- 3-set intersection (DEG × sex-DEG × module) design follows the workflow of
  the reference biomarker paper. ⚠️ **Wang et al.** (the paper you are mirroring;
  insert full citation — Cancer Inform / the SLE platelet ML paper).

## Mendelian randomisation (10, 11)
- **TwoSampleMR** — Hemani G, et al. The MR-Base platform supports systematic
  causal inference across the human phenome. *eLife* 2018;7:e34408.
- **OpenGWAS / ieugwasr** — Elsworth B, et al. The MRC IEU OpenGWAS data
  infrastructure. *bioRxiv* 2020. doi:10.1101/2020.08.10.244293.
- **eQTLGen (cis-eQTL exposures)** — Võsa U, et al. Large-scale cis- and
  trans-eQTL analyses identify thousands of genetic loci and polygenic scores
  that regulate blood gene expression. *Nat Genet* 2021;53(9):1300–10.
- **RA GWAS outcome (Okada, ieu-a-832)** — Okada Y, et al. Genetics of
  rheumatoid arthritis contributes to biology and drug discovery. *Nature*
  2014;506(7488):376–81.
- **RA GWAS replication (Stahl, ieu-a-833/834)** — Stahl EA, et al. Genome-wide
  association study meta-analysis identifies seven new rheumatoid arthritis risk
  loci. *Nat Genet* 2010;42(6):508–14.
- **IVW estimator** — Burgess S, Butterworth A, Thompson SG. Mendelian
  randomization analysis with multiple genetic variants using summarized data.
  *Genet Epidemiol* 2013;37(7):658–65.
- **MR-Egger** — Bowden J, Davey Smith G, Burgess S. Mendelian randomization
  with invalid instruments: effect estimation and bias detection through Egger
  regression. *Int J Epidemiol* 2015;44(2):512–25.
- **Weighted median** — Bowden J, et al. Consistent estimation in Mendelian
  randomization with some invalid instruments using a weighted median estimator.
  *Genet Epidemiol* 2016;40(4):304–14.
- **Weak-instrument / F-statistic** — Burgess S, Thompson SG. Avoiding bias from
  weak instruments in Mendelian randomization studies. *Int J Epidemiol*
  2011;40(3):755–64.
- **LD clumping (1000 Genomes ref)** — 1000 Genomes Project Consortium. A global
  reference for human genetic variation. *Nature* 2015;526:68–74.

## Feature selection & classification (12, 15, 16)
- **LASSO / elastic-net (glmnet)** — Friedman J, Hastie T, Tibshirani R.
  Regularization paths for generalized linear models via coordinate descent.
  *J Stat Softw* 2010;33(1):1–22.
- **Random forest** — Breiman L. Random forests. *Machine Learning* 2001;45:5–32.
  Liaw A, Wiener M. Classification and regression by randomForest. *R News* 2002.
- **SVM-RFE** — Guyon I, Weston J, Barnhill S, Vapnik V. Gene selection for
  cancer classification using support vector machines. *Machine Learning*
  2002;46:389–422. e1071: Meyer D, et al. (libsvm — Chang CC, Lin CJ. ACM TIST 2011).
- **caret (hyperparameter tuning)** — Kuhn M. Building predictive models in R
  using the caret package. *J Stat Softw* 2008;28(5):1–26.
- **XGBoost** — Chen T, Guestrin C. XGBoost: a scalable tree boosting system.
  *KDD* 2016:785–94.

## Validation & honest evaluation (12b, 13, 14, 24)
- **Feature-selection bias / nested CV** — Ambroise C, McLachlan GJ. Selection
  bias in gene extraction on the basis of microarray gene-expression data.
  *PNAS* 2002;99(10):6562–6. Simon R, et al. Pitfalls in the use of DNA
  microarray data for diagnostic and prognostic classification. *JNCI*
  2003;95(1):14–8.
- **ROC / AUC (pROC)** — Robin X, et al. pROC: an open-source package for R and
  S+ to analyze and compare ROC curves. *BMC Bioinformatics* 2011;12:77.
- **DeLong CI** — DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas
  under two or more correlated ROC curves: a nonparametric approach.
  *Biometrics* 1988;44(3):837–45.
- **Diagnosis × sex interaction (limma)** — as limma above; interaction contrast
  per Smyth GK. *Stat Appl Genet Mol Biol* 2004;3:Article3.

## Clinical prediction model (24)
- **Nomogram / calibration (rms)** — Harrell FE Jr. *Regression Modeling
  Strategies*, 2nd ed. Springer; 2015. (rms R package.)
- **Decision curve analysis** — Vickers AJ, Elkin EB. Decision curve analysis:
  a novel method for evaluating prediction models. *Med Decis Making*
  2006;26(6):565–74.
- **Clinical impact curve** — Kerr KF, Brown MD, Zhu K, Janes H. Assessing the
  clinical impact of risk prediction models with decision curves. *J Clin Oncol*
  2016;34(21):2534–40.

## Pathway enrichment (25)
- **clusterProfiler** — Wu T, et al. clusterProfiler 4.0: A universal enrichment
  tool for interpreting omics data. *Innovation* 2021;2(3):100141. (orig: Yu G,
  et al. *OMICS* 2012;16(5):284–7.)
- **Gene Ontology** — Ashburner M, et al. Gene Ontology. *Nat Genet*
  2000;25(1):25–9. Gene Ontology Consortium. *Nucleic Acids Res* 2021;49:D325–34.
- **KEGG** — Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes.
  *Nucleic Acids Res* 2000;28(1):27–30.
- **org.Hs.eg.db** — Carlson M. org.Hs.eg.db: Genome wide annotation for Human.
  Bioconductor.

## Immune-cell deconvolution and composition adjustment (05c, 13b)
- **CIBERSORT** — Newman AM, et al. Robust enumeration of cell subsets from
  tissue expression profiles. *Nat Methods* 2015;12(5):453–7.
- **IOBR (wrapper used)** — Zeng D, et al. IOBR: multi-omics immuno-oncology
  biological research to decode tumor microenvironment. *Front Immunol*
  2021;12:687975.
- **MCP-counter (independent marker-based cross-check)** — Becht E, et al.
  Estimating the population abundance of tissue-infiltrating immune and stromal
  cell populations using gene expression. *Genome Biol* 2016;17:218.
- **Compositional data / centred log-ratio** — Aitchison J. The statistical
  analysis of compositional data. *J R Stat Soc B* 1982;44(2):139–77.
- **Deconvolution in immunology, general** — Shen-Orr SS, Gaujoux R. Computational
  deconvolution: extracting cell type-specific information from heterogeneous
  samples. *Curr Opin Immunol* 2013;25(5):571–8.
- **Bootstrap confidence intervals (small-n AUC)** — Carpenter J, Bithell J.
  Bootstrap confidence intervals: when, which, what? *Stat Med*
  2000;19(9):1141–64.

## MHC sensitivity and colocalisation (10c, 10d)
- **coloc / approximate Bayes factor colocalisation** — Giambartolomei C, et al.
  Bayesian test for colocalisation between pairs of genetic association studies
  using summary statistics. *PLoS Genet* 2014;10(5):e1004383.
- **Colocalisation priors, p12** — Wallace C. Eliciting priors and relaxing the
  single causal variant assumption in colocalisation analyses. *PLoS Genet*
  2020;16(4):e1008720.
- **SMR / HEIDI — rationale for distinguishing pleiotropy from linkage** —
  Zhu Z, et al. Integration of summary data from GWAS and eQTL studies predicts
  complex trait gene targets. *Nat Genet* 2016;48(5):481–7.
- **MHC fine-mapping in rheumatoid arthritis** — Raychaudhuri S, et al. Five
  amino acids in three HLA proteins explain most of the association between MHC
  and seropositive rheumatoid arthritis. *Nat Genet* 2012;44(3):291–6.
- **Long-range LD across the MHC** — Trynka G, et al. Dense genotyping identifies
  and localizes multiple common and rare variant association signals in celiac
  disease. *Nat Genet* 2011;43(12):1193–201.

## Sex-differential expression (05d)
- **Sex effects on gene expression across human tissues** — Oliva M, et al. The
  impact of sex on gene expression across human tissues. *Science*
  2020;369(6509):eaba3066.
- **Difference between "significant" and "significantly different"** — Gelman A,
  Stern H. The difference between "significant" and "not significant" is not
  itself statistically significant. *Am Stat* 2006;60(4):328–31.

## Figures
- **ggplot2** — Wickham H. *ggplot2: Elegant Graphics for Data Analysis*.
  Springer; 2016.
- **ComplexHeatmap** — Gu Z, Eils R, Schlesner M. Complex heatmaps reveal
  patterns and correlations in multidimensional genomic data. *Bioinformatics*
  2016;32(18):2847–9.
- **patchwork / cowplot / ggrepel / ggVennDiagram** — Pedersen TL (patchwork);
  Wilke CO (cowplot); Slowikowski K (ggrepel); Gao CH, et al. (ggVennDiagram,
  *Front Genet* 2021).

## Datasets (GEO) — ⚠️ verify primary publications
- **GSE93272** (whole blood, RA) — Tasaki S, et al. Multi-omic disease and drug
  signatures of rheumatoid arthritis. *Nat Commun* 2018;9:2755. ⚠️ confirm.
- **GSE110169** (whole blood) — ⚠️ insert primary publication.
- **GSE15573** (PBMC, external blood) — ⚠️ insert primary publication.
- **GSE89408** (synovial tissue RNA-seq, cross-tissue) — ⚠️ insert primary
  publication (GEO accession GSE89408).
- **East Asian RA GWAS (BBJ, bbj-a-151)** — Ishigaki K, et al. Large-scale
  genome-wide association study in a Japanese population identifies novel
  susceptibility loci across multiple diseases. *Nat Genet* 2020;52:669–79. ⚠️ confirm.

## R / Bioconductor
- R Core Team. R: A Language and Environment for Statistical Computing. R
  Foundation for Statistical Computing; cite the version from `sessionInfo()`.
- Gentleman RC, et al. Bioconductor: open software development for computational
  biology and bioinformatics. *Genome Biol* 2004;5:R80.
