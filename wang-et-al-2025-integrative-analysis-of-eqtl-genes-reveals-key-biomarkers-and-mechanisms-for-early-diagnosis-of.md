```
https://doi.org/10.1177/11769351251400465
```
Cancer Informatics
Volume 24: 1 –14
```
© The Author(s) 2025
```
Article reuse guidelines:
sagepub.com/journals-permissions
```
DOI: 10.1177/11769351251400465
```
journals.sagepub.com/home/cix
Creative Commons Non Commercial CC BY-NC: This article is distributed under the terms of the Creative Commons Attribution-NonCommercial
```
4.0 License (https://creativecommons.org/licenses/by-nc/4.0/) which permits non-commercial use, reproduction and distribution of the work without
```
```
further permission provided the original work is attributed as specified on the SAGE and Open Access pages (https://us.sagepub.com/en-us/nam/open-access-at-sage).
```
Integrative Analysis of eQTL
Genes Reveals Key Biomarkers and
Mechanisms for Early Diagnosis of
Pancreatic Ductal Adenocarcinoma
Xuebo Wang1* , Xusheng Zhang1*, Shicai Liang1, Jialong Wang1,
Yannan Xie 1, Jiawei Wang1 and Bendong Chen2
Abstract
```
Background: Pancreatic ductal adenocarcinoma (PDAC) is a highly lethal malignancy with a dismal 5-year survival rate,
```
largely due to the absence of reliable biomarkers for early detection. The molecular mechanisms underpinning PDAC
pathogenesis remain incompletely understood, highlighting the urgent need for novel diagnostic strategies.
```
Objective: This study aimed to integrate eQTL-driven Mendelian randomization (MR) with transcriptomic and genome-
```
wide association data to identify causal PDAC-associated genes and construct a diagnostic nomogram based on 5 hub
```
genes (CTSC, SMYD3, MFGE8, IGFBP7, POC1B) for early detection of pancreatic ductal adenocarcinoma (PDAC).
```
```
Methods: Transcriptomic data from GSE62165 and GSE25471 were retrieved from the Gene Expression Omnibus
```
```
(GEO) and processed for differential expression using LIMMA and GEO2R, followed by batch correction and weighted
```
```
gene co-expression network analysis (WGCNA). Summary-level eQTL statistics were obtained from OpenGWAS, and
```
```
GWAS data included over 5000 PDAC cases. MR analysis was performed using inverse variance weighted (IVW) as
```
the primary approach, supplemented with MR-Egger, weighted median, weighted mode, and MR-PRESSO. Instrument
strength, pleiotropy, and heterogeneity were assessed via F-statistics, Egger intercept, and Cochran’s Q test. Candidate
```
genes were filtered using a consensus approach combining random forest (RF), support vector machine-recursive feature
```
```
elimination (SVM-RFE), and Lasso regression. Diagnostic performance was evaluated via ROC curves, C-index, calibration
```
plots, and decision curve analysis. Mechanistic insights were derived from KEGG and GO enrichment analyses, as well
```
as protein-protein interaction (PPI) network analyses.
```
```
Results: Five eQTL-associated hub genes—CTSC, SMYD3, MFGE8, IGFBP7, and POC1B—were identified as
```
causally linked to PDAC via robust MR analysis with minimal evidence of pleiotropy or heterogeneity. These genes
```
demonstrated high diagnostic potential (AUC > 0.85, P < .001). A diagnostic nomogram incorporating these genes
```
```
achieved strong predictive performance (C-index = 0.92) with favorable clinical decision curve results. Functional
```
enrichment and PPI analyses implicated these genes, particularly CTSC, in modulating the ITGAV/ITGB3–PI3K–Akt
signaling axis, contributing to PDAC cell cycle regulation and apoptosis resistance.
```
Conclusions: This study presents a multi-omics, MR-informed framework for identifying eQTL-regulated biomarkers
```
of PDAC. The identified hub genes offer promising avenues for early detection, while the mechanistic mapping of the
PI3K–Akt pathway provides translational insights. These findings warrant further validation in clinical and experimental
settings and hold potential to reshape PDAC diagnostic strategies.
```
Pancreatic ductal adenocarcinoma (PDAC) remains a formidable clinical challenge due to its aggressive nature and lack
```
of effective early diagnostic biomarkers. To address this, we integrated transcriptomic data, genome-wide association
```
studies (GWAS), and expression quantitative trait loci (eQTL) information using Mendelian randomization (MR) to
```
identify genes causally associated with PDAC risk. Differentially expressed genes were identified across 2 GEO datasets
```
(GSE62165, GSE25471) and prioritized using weighted gene co-expression network analysis (WGCNA). MR analysis
```
employing IVW, MR-Egger, weighted median, and MR-PRESSO identified 5 hub genes—CTSC, SMYD3, MFGE8, IGFBP7,
1Ningxia Medical University, Yinchuan, China
2Department of Hepatobiliary Surgery, General Hospital of Ningxia
Medical University, Yinchuan, China
*These authors contributed equally to this work and share the first
authorship.
1400465CIX Cancer InformaticsWang et al
Original Research
Corresponding author:
Bendong Chen, Department of Hepatobiliary Surgery, General Hospital
of Ningxia Medical University, Yinchuan 750004, China.
```
Email: chenbendong@nxmu.edu.cn
```
2 Cancer Informatics
and POC1B—as significant causal drivers of PDAC. These genes were incorporated into a diagnostic model constructed
```
using machine learning approaches (random forest, SVM-RFE, Lasso), which achieved strong classification performance
```
```
(AUC > 0.85) and excellent calibration (C-index = 0.92). Functional enrichment and protein-protein interaction analyses
```
revealed that CTSC regulates the ECM-integrin–PI3K–Akt signaling pathway, contributing to tumor cell proliferation
and survival. The findings establish a multi-omics-based biomarker panel with strong diagnostic utility and mechanistic
relevance, suggesting a potential framework for future translational validation in clinical cohorts.
Keywords
pancreatic ductal adenocarcinoma, eQTL, Mendelian randomization, diagnostic biomarkers, PI3K/AKT pathway,
bioinformatics integration
```
Received: 23 July 2025; accepted: 7 November 2025
```
Highlights
• Integration of GWAS, eQTL, and transcriptomic
data for PDAC risk gene discovery.
• MR identifies CTSC, SMYD3, MFGE8, IGFBP7,
POC1B as causal hub genes.
• High-accuracy diagnostic nomogram model devel-
```
oped (C-index = 0.92).
```
• PI3K/AKT signaling implicated in PDAC
progression.
• Translational potential for early detection and thera-
peutic targeting of PDAC.
Introduction
```
Pancreatic ductal adenocarcinoma (PDAC) ranks as the
```
seventh most common malignancy globally and carries an
exceptionally poor 5-year survival rate of approximately
11%.1 As the fourth leading cause of cancer-related mortality
in Western populations,2,3 PDAC is projected to become the
second most fatal cancer by 2030.4 Its incidence varies
widely across geographic and ethnic groups, reflecting a
complex interplay of genetic predisposition, lifestyle factors,
and environmental exposures.2 Although established risk
factors—such as smoking, chronic pancreatitis, obesity, type
2 diabetes, and family history—have been identified,5-8 the
molecular mechanisms initiating PDAC remain incom-
pletely characterized. Therefore, identifying biomarkers
capable of enabling early detection is a critical unmet clinical
need. Current biomarkers like CA19-9 exhibit limited sensi-
```
tivity (50%-70%) for early-stage PDAC,53 underscoring the
```
urgency for novel molecular tools. Genome-wide association
```
studies (GWASs) have revealed both rare high-penetrance
```
mutations and common low-penetrance variants that contrib-
ute to PDAC susceptibility.9-17 Meta-analyses and large-scale
multicenter studies have identified additional risk loci,18-27
but their functional relevance—particularly their regulatory
effects on gene expression—remains poorly understood,
limiting clinical translation.
To date, GWASs have identified over 30 genomic loci
associated with PDAC risk.9 However, the precise genes or
regulatory elements through which these variants exert
functional effects remain largely undefined. Expression
```
quantitative trait loci (eQTLs)—genomic regions that
```
modulate gene expression levels—represent a valuable tool
to bridge this gap. eQTLs are categorized as cis-eQTLs,
located within ± 1 Mb of their target genes, or trans-eQTLs,
which influence distal gene expression via long-range regu-
latory interactions.28,29 By linking non-coding genetic varia-
tion to transcriptomic changes, eQTLs provide mechanistic
insight into cancer susceptibility and tumor biology.30-32
Several studies have demonstrated that eQTLs contribute to
oncogenesis in various cancers,33-36 such as NTN4 eQTLs in
```
breast cancer37 and splicing QTLs (sQTLs) driving tran-
```
scriptomic dysregulation in non-small cell lung cancer.38
When integrated with GWAS summary statistics, Mendelian
```
randomization (MR) enables causal inference between gene
```
expression traits and disease outcomes, providing a robust
framework for identifying functionally relevant genes.39
Unlike GWAS, which identifies disease-associated
```
genomic loci, expression quantitative trait loci (eQTLs)
```
map genetic variants that directly regulate gene expression
```
levels. Cis-eQTLs localize near target genes (±1 Mb),
```
while trans-eQTLs exert distal effects. 28,29 Integrating
eQTLs with GWAS via MR enables causal inference
between gene expression and disease.
Beyond risk association, eQTLs influence key onco-
genic processes, including tumor progression, treatment
response, and outcomes of immunotherapy. As such, they
are increasingly recognized as promising biomarkers for
both early diagnosis and therapeutic targeting. We hypoth-
esize that specific pancreatic eQTLs mediate PDAC risk by
driving transcriptional dysregulation of pathogenic genes.
By integrating eQTL and GWAS data using MR, we aim to
identify causal PDAC-associated genes and elucidate the
regulatory pathways underpinning early tumor develop-
ment. This approach holds potential to uncover clinically
actionable biomarkers that enable timely intervention in
high-risk populations.
Despite growing understanding of PDAC pathogenesis,
no standardized screening protocols currently exist for
asymptomatic individuals.40 Comprehensive identification
of both heritable and non-heritable risk factors is essential for
stratifying population-level risk and informing targeted
screening strategies. Systematic characterization of regula-
```
tory single-nucleotide polymorphisms (SNPs) that influence
```
gene expression will advance mechanistic insight into PDAC
and potentially transform diagnostic paradigms.
Wang et al 3
To our knowledge, this study represents the first integra-
tive investigation linking eQTLs to PDAC pathogenesis
through a comprehensive Mendelian randomization frame-
work. Leveraging publicly available GWAS and transcrip-
tomic datasets, we performed 2-sample MR analyses to
identify causal eQTL-gene-disease relationships, supported
by sensitivity and pleiotropy testing. Weighted gene co-
```
expression network analysis (WGCNA) and pathway
```
enrichment analysis were employed to characterize func-
tional modules relevant to PDAC. To translate these find-
ings into clinical practice, we constructed a diagnostic
model using machine learning approaches—random forest
```
(RF), support vector machine-recursive feature elimination
```
```
(SVM-RFE), and Lasso regression—and evaluated its pre-
```
dictive accuracy using receiver operating characteristic
```
(ROC) curves, calibration plots, and decision curve analy-
```
```
sis (DCA). A final nomogram integrating key diagnostic
```
genes was developed for individualized risk assessment.
By integrating genomic, transcriptomic, and causal
inference methods, this study advances both mechanistic
understanding and diagnostic precision in PDAC. The
identification of eQTL-mediated hub genes and their
involvement in key oncogenic pathways, particularly the
ITGAV/ITGB3–PI3K–Akt signaling axis, provides a foun-
dation for translational biomarker development. Our find-
ings establish a framework for future validation using
patient-derived models, ultimately supporting the develop-
ment of precision diagnostics and targeted therapeutics in
pancreatic cancer.
Materials and Methods
Data Sources and Preprocessing
Transcriptomic data for pancreatic ductal adenocarcinoma
```
(PDAC) were obtained from the Gene Expression Omnibus
```
```
(GEO) under accession numbers GSE62165 and GSE25471.
```
Raw expression profiles were preprocessed using standard
R packages, and batch effects were corrected using ComBat
in the sva package, with empirical Bayes adjustment for
site-specific biases in GSE62165 and GSE25471. Genome-
wide association summary statistics for PDAC were derived
from a large-scale meta-analysis encompassing 5430 histo-
logically confirmed cases of PDAC. eQTL summary-level
```
data were retrieved from the OpenGWAS Project (https://
```
```
gwas.mrcieu.ac.uk/; Supplemental Table 1), comprising
```
19 942 individuals of European ancestry. Only cis-eQTLs—
defined as SNPs located within ±2 Mb of a gene’s tran-
```
scription start site (TSS)—were included for downstream
```
```
Mendelian randomization (MR) analysis. SNPs reaching
```
```
genome-wide significance (P < 5 × 10 −8) were retained for
```
further instrument selection.
Instrument Selection and MR Design
Genetic instruments were selected based on genome-wide
```
significance (P < 5 × 10 −8) and independence, as deter-
```
```
mined by linkage disequilibrium (LD) clumping (r2 < 0.001,
```
```
clumping window = 10 000 kb). The final instrumental
```
```
variables (IVs) were restricted to cis-eQTLs to minimize
```
horizontal pleiotropy. Analyses used R v4.4.0 with pack-
```
ages: TwoSampleMR (v0.5.6), MendelianRandomization
```
```
(v0.9.0), sva (v3.48.0) for batch correction. The inverse
```
```
variance weighted (IVW) method served as the primary
```
estimator. Supplementary MR estimators included
MR-Egger, weighted median, weighted mode, and simple
mode approaches to assess the consistency and robustness
of causal effects. Instrument strength was evaluated using
```
F-statistics, and weak instruments were excluded (F < 10).
```
```
Effect alleles for exposure (eQTL) and outcome (PDAC
```
```
GWAS) were harmonized to ensure consistent directional-
```
ity. Palindromic SNPs with intermediate allele frequencies
were excluded.
Assessment of Heterogeneity and Pleiotropy
To evaluate heterogeneity and horizontal pleiotropy among
instruments, we performed Cochran’s Q test and MR-Egger
intercept analysis, respectively. The HEIDI-outlier test was
conducted to exclude variants with potential single-SNP-
driven effects. The sensitivity of causal estimates was
assessed through leave-one-out analyses, where each SNP
was iteratively excluded. Additionally, MR-PRESSO
```
(Mendelian Randomization Pleiotropy RESidual Sum and
```
```
Outlier) was used to detect and correct for global horizontal
```
```
pleiotropy. A false discovery rate (FDR) adjustment was
```
applied using the Benjamini-Hochberg method to control
for multiple testing across gene-trait associations.
Weighted Gene Co-expression Network
```
Analysis (WGCNA) and Functional
```
Enrichment
```
A total of 2090 differentially expressed genes (DEGs) were
```
identified via LIMMA-based analysis of the GEO datasets,
followed by batch correction. WGCNA was then applied to
construct co-expression modules associated with PDAC
status. Hub modules were selected based on module-trait
correlation coefficients. Functional annotation of prior-
itized gene modules was performed using the clusterPro-
filer, enrichplot, and ggplot2 R packages. Gene Ontology
```
(GO) terms (biological process, cellular component, molec-
```
```
ular function) and Kyoto Encyclopedia of Genes and
```
```
Genomes (KEGG) pathway enrichment were assessed,
```
with statistical significance defined as FDR < 0.05.
Machine Learning-Based Diagnostic Model
Development
To optimize diagnostic gene selection, we applied 3 super-
```
vised machine learning algorithms: Random Forest (RF),
```
Support Vector Machine Recursive Feature Elimination
```
(SVM-RFE), and Least Absolute Shrinkage and Selection
```
```
Operator (LASSO) regression. Model residuals were com-
```
pared via reverse cumulative residual boxplots and distri-
bution plots. Genes identified by all 3 algorithms were
retained for model construction. Predictive performance
was evaluated using receiver operating characteristic
4 Cancer Informatics
```
(ROC) curves and area under the curve (AUC) values.
```
Model overfitting was minimized via cross-validation.
Nomogram Construction and Clinical
Evaluation
A diagnostic nomogram was developed to visualize indi-
vidual-level PDAC risk based on the selected feature genes.
Gene-specific coefficients were transformed into risk
points, and total risk scores were computed to stratify
patients. Model calibration was evaluated using calibration
plots, while clinical utility was assessed via decision curve
```
analysis (DCA) and clinical impact curves (CICs). External
```
validation was performed using an independent dataset
from GEO to assess generalizability.
Mechanistic Exploration and Pathway
Mapping
To further investigate underlying mechanisms, hub genes
identified via MR and WGCNA were subjected to protein-
```
protein interaction (PPI) network reconstruction using
```
```
GeneMANIA (https://genemania.org/). Candidate genes,
```
including CTSC and ITGB3, were annotated with GO and
KEGG databases via clusterProfiler. Special emphasis was
placed on the ECM-integrin–PI3K–Akt signaling cascade to
explore oncogenic processes, including apoptosis resistance,
extracellular matrix remodeling, and cell cycle progression.
The overall study design is illustrated in Figure 1.
Results
```
MR analysis prioritized eight causal genes (ASNS, CTSC,
```
```
GLIPR2, IGFBP2, IGFBP7, MFGE8, POC1B, SMYD3),
```
```
of which 5 (CTSC, SMYD3, MFGE8, IGFBP7, POC1B)
```
```
were consolidated into a diagnostic model (Figure 3).
```
Transcriptomic Data Processing and Co-
expression Network Construction
To characterize PDAC-related transcriptional signatures, 2
```
GEO datasets (GSE62165 and GSE25471) were analyzed
```
using LIMMA and GEO2R, identifying 1898 overlapping
```
differentially expressed genes (DEGs; Figure 2A-C). After
```
```
batch effect correction via the sva package (Figure 2D and E),
```
```
a merged expression matrix yielded 2090 DEGs (Figure 2F).
```
```
Weighted gene co-expression network analysis (WGCNA)
```
```
identified 9 gene modules (Figure 2G and H), among which
```
Figure 1. Overall study design. Schematic flowchart illustrating the integrated multi-omics workflow of this study. Transcriptomic
```
data (GEO), eQTL mapping (OpenGWAS), Mendelian Randomization (MR), WGCNA, functional enrichment, and machine learning
```
algorithms were combined to identify diagnostic biomarkers and causal pathways in PDAC.
Wang et al 5
```
the blue module (n = 412 genes) exhibited the highest correla-
```
```
tion with PDAC status (r = .82, P < 1 × 10−6). Within this
```
```
module, gene significance (GS) and module membership
```
```
(MM) analysis confirmed strong intra-module connectivity
```
```
for PDAC-associated hub genes (Figure 2I and J).
```
Instrument Selection and Strength
Evaluation
Cis-eQTLs reaching genome-wide significance
```
(P < 5 × 10−8) were selected as candidate instrumental vari-
```
```
ables (IVs). SNPs in high linkage disequilibrium (r2 > 0.001,
```
```
window size = 10 000 kb) were excluded via clumping. All
```
```
retained instruments showed strong strength (F-statistics:
```
```
12.5–89.3), surpassing the conventional threshold (F > 10),
```
thus mitigating weak instrument bias. These criteria ensured
the validity of assumptions regarding relevance, independ-
ence, and exclusion restriction for MR analysis.
Primary Mendelian Randomization Findings
```
Using a 2-sample Mendelian Randomization (MR) frame-
```
```
work with the inverse-variance weighted (IVW) method, we
```
identified 125 gene-trait pairs significantly associated with
Figure 3. Integration of DEGs and MR results to identify
```
signature genes. (A) Venn diagram showing the intersection of
```
1898 DEGs, 2090 limma-derived DEGs, and 412 WGCNA blue
module genes, identifying 1721 PDAC-associated candidates.
```
(B) The overlap of upregulated DEGs and MR-derived risk loci
```
```
(OR > 1) identified 6 disease-associated genes. (C) Overlap of
```
```
downregulated DEGs and MR-derived protective loci (OR < 1)
```
yielded 2 protective genes.
```
Figure 2. Differential expression and co-expression network analysis. (A–B) Volcano plots depicting differentially expressed genes
```
```
(DEGs) between PDAC and normal tissues. (C) Venn diagram showing the overlap of key gene sets. (D–E) PCA plots of GSE62165
```
```
and GSE25471 datasets before and after batch correction. (F) Heatmap of upregulated and downregulated DEGs. (G) Model fit
```
```
R² curve. (H) Dendrogram and module color assignment via WGCNA. (I) Heatmap of module–trait correlations. (J) Correlation
```
```
between module membership (MM) and gene significance (GS) in the blue module (r = 0.8, P = 5.4e−93).
```
6 Cancer Informatics
```
pancreatic ductal adenocarcinoma (PDAC) (P < .05;
```
```
Supplemental Table 2). To prioritize biologically relevant can-
```
didates, we performed integrative filtering across 3 data lay-
```
ers: 1898 differentially expressed genes (DEGs) from
```
cross-dataset analysis, 2090 DEGs identified via LIMMA, and
412 WGCNA-derived blue module genes. This yielded a
```
refined set of 1721 PDAC-associated candidate genes (Figure
```
```
3A). Among these, 8 genes demonstrated statistically signifi-
```
cant causal associations with PDAC based on MR analysis,
```
including ASNS (odds ratio [OR] = 0.708, P = .011), CTSC
```
```
(OR = 1.273, P = .041), GLIPR2 (OR = 1.850, P = .049),
```
```
IGFBP2 (OR = 0.704, P = .038), IGFBP7 (OR = 1.378,
```
```
P = .034), MFGE8 (OR = 1.387, P = .011), POC1B (OR = 1.824,
```
```
P = .014), and SMYD3 (OR = 1.399, P = .041; Figure 4A and
```
```
B). These findings highlight both risk-promoting and protec-
```
tive gene candidates, supporting their potential roles as eQTL-
associated modulators of PDAC susceptibility.
Sensitivity Analyses and Pleiotropy
Assessment
Sensitivity analyses confirmed the robustness of the
causal estimates. Leave-one-out analysis showed that no
```
single SNP drove the association signal (Figure 5B).
```
MR-Egger regression yielded intercepts not significantly
```
different from zero (P > .05 for all genes), indicating no
```
evidence of directional pleiotropy. The weighted median
estimator confirmed significant associations for all 8
```
genes (P < .05), supporting consistent effect directions
```
```
across analytical methods (Figure 6). MR-PRESSO
```
detected no outlier SNPs that influenced causal infer-
ence, further affirming the validity. Cochran’s Q test
revealed no significant heterogeneity among instruments
```
(IVW Q P = .15; MR-Egger Q P = .22), indicating homo-
```
geneity in variant-level effects.
```
Figure 4. Mendelian randomization identifies causal genes for PDAC. (A) Forest plot showing causal effects (odds ratios with 95%
```
```
CI) of 8 eQTL-associated genes on PDAC risk. (B) Scatter plot visualizing SNP-level effects used in MR analysis.
```
Wang et al 7
Functional Annotation of Causal Genes
```
Gene Ontology (GO) and KEGG pathway enrichment anal-
```
yses revealed that the 8 MR-identified genes were signifi-
cantly enriched in glucocorticoid signaling, growth factor
binding, extracellular matrix organization, and insulin-like
```
growth factor activity (Figure 7A). Cellular localization
```
analysis highlighted enrichment in the endoplasmic reticu-
lum lumen and extracellular matrix. KEGG analysis indi-
cated involvement in amino acid biosynthesis, lysine
degradation, and carbon metabolism pathways, suggesting
functional convergence on metabolic reprograming and
immune modulation during PDAC progression.
Diagnostic Model Construction and
Validation
```
The dataset (n = 320 samples) was split into 70% training
```
```
(n = 224) and 30% validation (n = 96) sets for all machine
```
```
learning algorithms. Random Forest (RF) and Support
```
```
Vector Machine (SVM) algorithms were compared via
```
```
residual analysis (Figure 7B and C), with both showing
```
```
comparable diagnostic performance (AUC differ-
```
```
ence < 0.02; Figure 7D). Feature selection using RF, SVM-
```
RFE, and LASSO identified 5 consensus diagnostic genes
```
(Figure 7E-I). A nomogram based on these genes achieved
```
excellent calibration and high diagnostic accuracy
```
(C-index = 0.91; Figure 8A and B). Decision curve analysis
```
```
(DCA) and clinical impact curves (CICs) confirmed the
```
model’s clinical benefit, particularly at risk thresholds of
```
20% to 60% (Figure 8C and D), supporting its feasibility
```
for population-based screening in high-risk cohorts.
Mechanistic Exploration via Network and
Pathway Analysis
```
Protein–protein interaction (PPI) network analysis identi-
```
fied 25 highly interconnected PDAC hub genes, including
```
CTSC, ITGAV, and ITGB3 (Figure 9A). Functional annota-
```
tion revealed that these genes are involved in insulin-like
```
Figure 5. Sensitivity analysis of MR results. (A) Funnel plot evaluating heterogeneity among SNP instruments. (B) Leave-one-out
```
analysis assessing the influence of individual SNPs on MR estimates.
8 Cancer Informatics
growth factor binding and lipid metabolic regulation.
Pearson correlation analysis demonstrated significant co-
```
expression of CTSC with ITGB3 (r = .477, P < .05) and
```
```
ITGAV (r = .568, P < .05; Figure 9B). KEGG pathway
```
mapping showed convergence on the ECM–integrin–
```
PI3K–Akt signaling axis (Figure 9C), suggesting CTSC as
```
a central modulator of apoptosis resistance, ECM remode-
ling, and PDAC progression.
Discussion
```
Pancreatic ductal adenocarcinoma (PDAC) is a polygenic
```
and multifactorial malignancy with complex and incom-
pletely understood genetic architecture. Unlike other com-
mon cancers, the genetic susceptibility of PDAC remains
underexplored, particularly regarding regulatory variants.
```
While genome-wide association studies (GWAS) have
```
identified risk loci, many of these lie in non-coding regions,
suggesting that transcriptional regulation rather than pro-
tein-coding alterations may be involved. Expression quan-
```
titative trait loci (eQTLs), which modulate gene expression,
```
serve as functional links between non-coding variants and
disease-relevant transcriptional activity. In this study, we
integrated transcriptomic data from GEO with GWAS-
derived eQTL information using Mendelian randomization
```
(MR) and weighted gene co-expression network analysis
```
```
(WGCNA) to identify eQTL-driven genes and signaling
```
pathways potentially causal in PDAC pathogenesis.
Our findings nominate CTSC, SMYD3, and MFGE8 as
causal mediators of PDAC pathogenesis through distinct
oncogenic mechanisms. Whereas ABO and NR5A2 are
established PDAC risk genes,48 our study implicates novel
```
candidates: POC1B (centrosome integrity49) and ASNS
```
```
(asparagine synthetase). Unlike IGFBP7 (a known tumor
```
```
suppressor50), POC1B has not been previously linked to
```
PDAC, suggesting its potential as a unique early biomarker.
CTSC, a lysosomal protease, has previously been implicated
in metastasis through the formation of neutrophil extracellu-
```
lar traps (NETs).41 Its potential role in PDAC may involve
```
the recruitment of immune cells and the promotion of tumor-
promoting inflammation. SMYD3, an epigenetic regulator,
has been shown to activate Ras/ERK signaling and drive
transcriptional amplification in multiple cancers, including
PDAC.42 MFGE8, known for its immune-modulatory prop-
erties, may promote an immunosuppressive microenviron-
ment by facilitating PD-L1 trafficking on extracellular
vesicles, thereby contributing to immune escape and resist-
ance to anti-PD-1 therapy.43 Collectively, these genes high-
light the interplay between immune modulation, epigenetic
regulation, and tumor progression in PDAC.
Mechanistically, CTSC activates the ECM-integrin-
```
PI3K/AKT axis (Figure 9C). This pathway drives PDAC
```
```
progression by: (i) Enhancing cell survival via apoptosis
```
```
resistance,44,45 (ii) Facilitating ECM remodeling to promote
```
metastasis.60,61 Notably, CTSC was strongly associated
with integrin family members ITGAV and ITGB3, and
pathway enrichment confirmed convergence on the PI3K/
AKT cascade.44,45 This axis is known to drive resistance to
apoptosis and enhance oncogenic signaling in aggressive
PDAC subtypes. Inhibitors of PI3K have shown efficacy in
preclinical PDAC models, supporting the translational
potential of targeting this network.46,47
Figure 6. Summary of Mendelian randomization findings. Forest plot displaying causal associations of candidate genes with PDAC
from MR analysis, highlighting effect size and confidence intervals.
Wang et al 9
```
In addition to well-established PDAC risk genes (eg,
```
```
ABO, NR5A2), 48 our study identifies novel candidates,
```
including POC1B and IGFBP7. POC1B is crucial for
centrosome integrity and mitotic stability, 49 while
IGFBP7, a well-established tumor suppressor, regulates
TGF-β signaling and cell proliferation. 50 Their inclusion
in our diagnostic model suggests potential utility as early
biomarkers or therapeutic targets. The additive influence
of these genes reflects the polygenic nature of PDAC
and underscores the importance of integrating eQTL-
based functional annotation into biomarker discovery
frameworks.
Our machine learning-derived diagnostic model, inte-
grating RF, SVM-RFE, and LASSO feature selection, out-
```
performed conventional CA19-9 assays (AUC = 0.92). 51
```
Compared with multi-omics diagnostic tools in colorectal
cancer,52 our approach uniquely leverages eQTL-informed
genetic drivers of PDAC, enhancing diagnostic sensitivity
for early-stage disease. Future studies will extend valida-
tion to ctDNA-based platforms and integrate single-cell
```
RNA sequencing (scRNA-seq) to account for tumor micro-
```
environment heterogeneity.53,54
To our knowledge, 55 this is the first study to integrate
GEO transcriptomic data, GWAS summary statistics, and
```
Figure 7. GO/KEGG enrichment score and establishment and of diagnostic model through machine learning. (A) GO/KEGG-
```
```
enriched bar plot showing the top 3 pathways with the highest significance for Biological process (BP), Cellular component (CC),
```
```
and Molecular function (MF). (B) Residual boxplots of Random Forest (RF) and Support Vector Machine (SVM). (C) Reverse
```
```
cumulative residual distribution plots of Random Forest (RF) and Support Vector Machine (SVM). (D) ROC curves of Random
```
```
Forest (RF) and Support Vector Machine (SVM) models. (E and F) Random Forest Plot and Importance Ranking Plot. Importance
```
ranking of gene features in the random forest model based on the “Mean Decreased Gini” index. Gene importance ranking, with
the x-axis showing gene names and the y-axis indicating importance scores, reflecting the discriminative power of genes in the
```
model. (G) SVM Gene Selection Line Chart. Gene selection results using Support Vector Machine Recursive Feature Elimination
```
```
(SVM-RFE). The x-axis represents the number or rank of features, and the y-axis shows classification performance metrics (eg,
```
```
accuracy or error), optimizing the selection of critical gene subsets. (H) Confirmation intervals under each lambda in the LASSO
```
```
regression. (I) Venn diagram identifying five5 consensus signature genes through intersection analysis.
```
10 Cancer Informatics
eQTL information via MR to identify causally relevant
genes in PDAC. Using TwoSampleMR, we established
```
robust associations between eQTL-regulated genes (eg,
```
```
CTSC, SMYD3, MFGE8) and susceptibility to PDAC.
```
These findings provide a foundation for functional investi-
gation of candidate genes such as CTSC, potentially link-
ing neutrophil activation and NETosis to metastatic
progression.56 Similarly, SMYD3′s epigenetic function
may intersect with FAK pathway activity,57 and MFGE8
may contribute to immune evasion by reshaping extracel-
lular vesicle signaling. 58,59
Our results also support the ECM-integrin-PI3K/AKT
axis as a core signaling pathway in PDAC susceptibility.
Prior studies demonstrate that ECM remodeling drives
PDAC aggressiveness through mechanotransduction and
resistance to therapy.60,61 Integrins ITGAV and ITGB3,
enriched in basal-like PDAC subtypes, interact with focal
```
adhesion kinase (FAK) to promote chemoresistance. 62,63
```
CTSC may modulate this axis through regulation of integ-
rin expression or activity, contributing to cell adhesion and
invasion. These findings warrant in vitro and in
vivo validation to confirm the mechanistic link between
eQTL-mediated gene regulation and kinase pathway
activation.64-66
The nomogram developed here achieved diagnostic
accuracy comparable to cfDNA-based methylation pan-
```
els (eg, EpiPanGI Dx, AUC = 0.88), 67 offering a poten-
```
tially cost-effective screening alternative. Nevertheless,
validation in neoadjuvant-treated cohorts is essential, as
emerging evidence suggests that therapy-induced tran-
scriptional plasticity—such as the reprograming of neu-
roendocrine progenitor-like cells—can confound
diagnostic models. Future model iterations will incorpo-
rate spatial transcriptomics and proteomics to dissect
intratumoral heterogeneity and improve robustness
across clinical contexts. 68
This study has limitations. First, blood-derived eQTLs
may not fully represent pancreas-specific expression, par-
ticularly within the diverse tumor microenvironment. 69,70
Second, regulatory mechanisms such as histone modifica-
tion or RNA methylation were not integrated, possibly
missing key post-transcriptional effects.71 Third, the
```
Figure 8. Validation of the diagnostic model. (A) Nomogram illustrating the predicted probability of PDAC based on selected
```
```
gene expression levels. (B) Calibration plot comparing predicted versus observed probabilities. (C) Decision curve analysis (DCA)
```
```
assessing clinical net benefit. (D) Clinical impact curves (CICs) evaluating population-level diagnostic performance.
```
Wang et al 11
European ancestry of the reference population may limit
the generalizability of the results. Future studies should pri-
oritize eQTL datasets from diverse populations and explore
multi-omic QTLs—including mQTLs, sQTLs, and
pQTLs—to comprehensively characterize genetic regula-
tion in PDAC.72
Conclusion
In summary, we identified 5 eQTL-associated hub genes—
CTSC, SMYD3, MFGE8, IGFBP7, and POC1B—as puta-
tive causal mediators of PDAC development. Mechanistic
analyses support a model in which CTSC promotes apop-
totic resistance and cell cycle deregulation via the integrin-
PI3K/AKT axis. These findings provide a foundation for
biomarker development and therapeutic targeting, with
future validation in organoid and animal models expected
to bridge the gap between genetic epidemiology and trans-
lational oncology.73
ORCID iDs
Xuebo Wang https://orcid.org/0009-0004-1095-397X
Bendong Chen https://orcid.org/0009-0001-5618-6962
Author Contributions
All authors have read and approved the final manuscript. The spe-
cific contributions of each author are as follows:
Xuebo Wang※: Conceptualization, Formal analysis, Investigation,
Writing – original draft.
Ningxia Medical University, Yinchuan 750004, China,
```
Email: byron_wang89@163.com
```
Xusheng Zhang※: Data curation, Software, Validation,
Visualization, Writing – original draft.
Ningxia Medical University, Yinchuan 750004, China,
```
Email: zxswle1995@163.com
```
Shicai Liang: Methodology, Resources, Data curation.
Ningxia Medical University, Yinchuan 750004, China,
```
Email: liang199810@yeah.net
```
```
Figure 9. Protein–protein interaction (PPI) network and mechanistic validation. (A) PPI network of top 25 hub genes based on
```
```
interaction degree. (B) Correlation plots showing associations between CTSC expression and integrins ITGB3 and ITGAV (P < .05).
```
```
Pearson coefficients and linear regression lines are shown. (C) Schematic of ECM–integrin–PI3K–AKT pathway inferred from
```
eQTL–gene interactions.
12 Cancer Informatics
Jialong Wang: Investigation, Validation.
Ningxia Medical University, Yinchuan 750004, China,
```
Email: 18809535935@163.com
```
Yannan Xie: Investigation, Validation.
Ningxia Medical University, Yinchuan 750004, China,
```
Email: xieyanan0301@163.com
```
Jiawei Wang: Investigation, Validation.
Ningxia Medical University, Yinchuan 750004, China,
```
Email: Wjw9964@126.com
```
Bendong Chen*: Supervision, Project administration, Funding
acquisition, Writing – review & editing.
Department of Hepatobiliary Suegery, General Hospital of
Ningxia Medical University, Yinchuan 750004, China,
```
Email: chenbendong@nxmu.edu.cn
```
Funding
The authors received no financial support for the research, author-
ship, and/or publication of this article.
Declaration of Conflicting Interests
The authors declared no potential conflicts of interest with respect
to the research, authorship, and/or publication of this article.
Data Availability Statement
```
The data used in this study were obtained from the GEO (https://
```
```
www.ncbi.nlm.nih.gov/geo/) database and the genome-wide asso-
```
```
ciation study (GWAS) summary dataset. All of these databases are
```
publicly available. This study adheres to their respective data
usage and publication rules.
Clinical Trial Number
Not applicable.
Supplemental Material
Supplemental material for this article is available online.
References
1. Siegel RL, Miller KD, Fuchs HE, Jemal A. Cancer statistics,
2022. CA Cancer J Clin. 2022;72(1):7-33.
2. Ferlay J, Colombet M, Soerjomataram I, et al. Estimating the
global cancer incidence and mortality in 2018: GLOBOCAN
```
sources and methods. Int J Cancer. 2019;144:1941-1953.
```
3. Bray F, Ferlay J, Soerjomataram I, Siegel RL, Torre LA,
Jemal A. Global cancer statistics 2018: GLOBOCAN esti-
mates of incidence and mortality worldwide for 36 cancers
```
in 185 countries. CA Cancer J Clin. 2018;68:394-424.
```
4. Rahib L, Smith BD, Aizenberg R, Rosenzweig AB, Fleshman
JM, Matrisian LM. Projecting cancer incidence and deaths to
2030: the unexpected burden of thyroid, liver, and pancreas
```
cancers in the United States. Cancer Res. 2014;74:2913-
```
2921.
5. Afghani E, Klein AP. Pancreatic adenocarcinoma: trends in
epidemiology, risk factors, and outcomes. Hematol Oncol
```
Clin N Am. 2022;36(5):879-895.
```
6. Klein AP. Pancreatic cancer epidemiology: understand-
ing the role of lifestyle and inherited risk factors. Nat Rev
```
Gastroenterol Hepatol. 2021;18(7):493-502.
```
7. Maisonneuve P, Lowenfels AB. Risk factors for pancreatic
```
cancer: a summary review of meta-analytical studies. Int J
```
```
Epidemiol. 2015;44:186-198.
```
8. Lu Y, Gentiluomo M, Lorenzo-Bermejo J, et al. Mendelian
randomisation study of the effects of known and putative risk
```
factors on pancreatic cancer. J Med Genet. 2020;57:820-828.
```
9. Gentiluomo M, Canzian F, Nicolini A, Gemignani F, Landi
S, Campa D. Germline genetic variability in pancreatic can-
```
cer risk and prognosis. Semin Cancer Biol. 2022;79(20):105-
```
131.
10. Amundadottir L, Kraft P, Stolzenberg-Solomon RZ, et al.
Genome-wide association study identifies variants in the
ABO locus associated with susceptibility to pancreatic can-
```
cer. Nat Genet. 2009;41:986-990.
```
11. Petersen GM, Amundadottir L, Fuchs CS, et al. A genome-
wide association study identifies pancreatic cancer suscep-
tibility loci on chromosomes 13q22.1, 1q32.1 and 5p15.33.
```
Nat Genet. 2010;42:224-228.
```
12. Wolpin BM, Rizzato C, Kraft P, et al. Genome-wide associa-
tion study identifies multiple susceptibility loci for pancre-
```
atic cancer. Nat Genet. 2014;46:994-1000.
```
13. Childs EJ, Mocci E, Campa D, et al. Common variation at
2p13.3, 3q29, 7p13 and 17q25.1 associated with susceptibil-
```
ity to pancreatic cancer. Nat Genet. 2015;47:911-916.
```
14. Zhang M, Wang Z, Obazee O, et al. Three new pancreatic can-
cer susceptibility signals identified on chromosomes 1q32.1,
```
5p15.33 and 8q24.21. Oncotarget. 2016;7:66328-66343.
```
15. Klein AP, Wolpin BM, Risch HA, et al. Genome-wide meta-
analysis identifies five new susceptibility loci for pancreatic
```
cancer. Nat Commun. 2018;9:556.
```
16. Galeotti AA, Gentiluomo M, Rizzato C, et al. Polygenic and
multifactorial scores for pancreatic ductal adenocarcinoma
```
risk prediction. J Med Genet. 2021;58:369-377.
```
17. Lin Y, Nakatochi M, Hosono Y, et al. Genome-wide asso-
ciation meta-analysis identifies GP2 gene risk variants for
```
pancreatic cancer. Nat Commun. 2020;11:3175.
```
18. Campa D, Pastore M, Gentiluomo M, et al. Functional sin-
gle nucleotide polymorphisms within the cyclin-dependent
kinase inhibitor 2A/2B region affect pancreatic cancer risk.
```
Oncotarget. 2016;7:57011-57020.
```
19. Campa D, Rizzato C, Stolzenberg-Solomon R, et al. TERT
gene harbors multiple variants associated with pancreatic
```
cancer susceptibility. Int J Cancer. 2015;137:2175-2183.
```
20. Campa D, Rizzato C, Bauer AS, et al. Lack of replication
of seven pancreatic cancer susceptibility loci identified in
two Asian populations. Cancer Epidemiol Biomarkers Prev.
```
2013;22:320-323.
```
21. Gentiluomo M, Peduzzi G, Lu Y, Campa D, Canzian F.
Genetic polymorphisms in inflammatory genes and pancre-
atic cancer risk: a two-phase study on more than 14 000 indi-
```
viduals. Mutagenesis. 2019;34:395-401.
```
22. Gentiluomo M, Lu Y, Canzian F, Campa D. Genetic vari-
ants in taste-related genes and risk of pancreatic cancer.
```
Mutagenesis. 2019;34:391-394.
```
23. Xu X, Qian D, Liu H, et al. Genetic variants in the liver
kinase B1-AMP-activated protein kinase pathway genes and
```
pancreatic cancer risk. Mol Carcinog. 2019;58:1338-1348.
```
24. Feng Y, Liu H, Duan B, et al. Potential functional variants in
SMC2 and TP53 in the AURORA pathway genes and risk of
```
pancreatic cancer. Carcinogenesis. 2019;40:521-528.
```
25. Yang W, Liu H, Duan B, et al. Three novel genetic variants in
NRF2 signaling pathway genes are associated with pancre-
```
atic cancer risk. Cancer Sci. 2019;110:2022-2032.
```
26. Campa D, et al. Genome-wide association study identifies
an early onset pancreatic cancer risk locus. Int J Cancer.
```
2020;147:2065-2074.
```
27. Corradi C, Gentiluomo M, Gajdán L, et al. Genome-wide scan
of long noncoding RNA single nucleotide polymorphisms and
Wang et al 13
```
pancreatic cancer susceptibility. Int J Cancer. 2021;148:2779-
```
2788.
28. Han H, Su H, Lv Z, Zhu C, Huang J. Identifying MTHFD1
and LGALS4 as potential therapeutic targets in prostate can-
cer through multi-omics Mendelian randomization analysis.
```
Biomedicines. 2025;13(1):185.
```
29. Jain L, Fadason T, Schierding W, Vickers MH, O’Sullivan
JM, Perry JK. 3D interactions with the growth hormone
locus in cellular signalling and cancer-related pathways. J
```
Mol Endocrinol. 2020;64(4):209-222.
```
30. Strunz T, Grassmann F, Gayán J, et al. A mega-analysis of
```
expression quantitative trait loci (eQTL) provides insight
```
into the regulatory architecture of gene expression variation
```
in liver. Sci Rep. 2018;8:5865.
```
31. Gilad Y, Rifkin SA, Pritchard JK. Revealing the architec-
ture of gene regulation: the promise of eQTL studies. Trends
```
Genet. 2008;24:408-415.
```
32. Nica AC, Dermitzakis ET. Expression quantitative trait loci:
present and future. Philos. Trans. R. Soc. Lond. B Biol. Sci.
```
2013;368:20120362.
```
33. Dai JY, Wang X, Wang B, et al. DNA methylation and cis-
regulation of gene expression by prostate cancer risk SNPs.
```
PLoS Genet. 2020;16:e1008667.
```
34. Ferreira MA, et al. Genome-wide association and transcrip-
tome studies identify target genes and risk loci for breast can-
```
cer. Nat Commun. 2019;10:1741.
```
35. Fung JN, Mortlock S, Girling JE, et al. Genetic regulation
of disease risk and endometrial gene expression highlights
potential target genes for endometriosis and polycystic ovar-
```
ian syndrome. Sci Rep. 2018;8:11424.
```
36. Loo LWM, Lemire M, Le Marchand L. In silico pathway
analysis and tissue specific ciseQTL for colorectal cancer
```
GWAS risk variants. BMC Genomics. 2017;18:381.
```
37. Luan Y, Xian D, Zhao C, et al. Therapeutic targets for lung
```
cancer: genome-wide Mendelian randomization and colo-
```
```
calization analyses. Front Pharmacol. 2024;15:1441233.
```
38. Chen J, Wang Y, Jiang R, Qu Y, Li Y, Zhang Y. Application of
Mendelian randomized research method in oncology research:
```
bibliometric analysis. Front Oncol. 2024;14:1424812.
```
39. Zhu Z, Zhang F, Hu H, et al. Integration of summary data
from GWAS and eQTL studies predicts complex trait gene
```
targets. Nat Genet. 2016;48:481-487.
```
40. Neoptolemos JP, et al. Pancreatic cancer. Pancreat. Cancer.
```
2018;2:1-1661.
```
41. Xiao Y, Cong M, Li J, et al. Cathepsin C promotes breast
cancer lung metastasis by modulating neutrophil infiltration
and neutrophil extracellular trap formation. Cancer Cell.
```
2021;39(3):423-437.e7.
```
42. Giakountis A, Moulos P, Sarris ME, Hatzis P, Talianidis I.
Smyd3-associated regulatory pathways in cancer. Semin
```
Cancer Biol. 2017;42:70-80.
```
43. Wang W, Chen J, Wang S, et al. MFGE8 induces anti-PD-1
therapy resistance by promoting extracellular vesicle sorting
```
of PD-L1. Cell Rep Med. 2025;6(2):101922.
```
44. Yamamoto S, Tomita Y, Hoshida Y, et al. Prognostic signifi-
cance of activated Akt expression in pancreatic ductal adeno-
```
carcinoma. Clin Cancer Res. 2004;10:2846-2850.
```
45. Mehra S, Deshpande N, Nagathihalli N. Targeting PI3K
pathway in pancreatic ductal adenocarcinoma: rationale and
```
progress. Cancers. 2021;13(17):4434.
```
46. Bondar VM, Sweeney-Gotsch B, Andreeff M, Mills GB,
McConkey DJ. Inhibition of the phosphatidylinositol 3’-kinase-
AKT pathway induces apoptosis in pancreatic carcinoma cells
```
in vitro and in vivo. Mol Cancer Ther. 2002;1:989-997.
```
47. Alam H, Gu B, Lee MG. Histone methylation modi-
fiers in cellular signaling pathways. Cell Mol Life Sci.
```
2015;72(23):4577-4592.
```
48. Rizzato C, Campa D, Giese N, et al. Pancreatic cancer
susceptibility loci and their role in survival. PLoS One.
```
2011;6(11):e27921.
```
49. Venoux M, Tait X, Hames RS, Straatman KR, Woodland HR,
Fry AM. Poc1A and Poc1B act together in human cells to
```
ensure centriole integrity. J Cell Sci. 2013;126(1):163-175.
```
50. An W, Ben QW, Chen HT, et al. Low expression of IGFBP7
is associated with poor outcome of pancreatic ductal adeno-
```
carcinoma. Ann Surg Oncol. 2012;19(12):3971-3978.
```
51. Boyd LN, Ali M, Kam L, et al. The diagnostic value of the
CA19-9 and bilirubin ratio in patients with pancreatic cancer,
distal bile duct cancer and benign periampullary diseases, a
```
novel approach. Cancers. 2022;14(2):344.
```
52. Thomas AM, Manghi P, Asnicar F, et al. Author Correction:
Metagenomic analysis of colorectal cancer datasets identi-
fies cross-cohort microbial diagnostic signatures and a link
```
with choline degradation. Nat Med. 2019;25(12):1948.
```
53. Lapitz A, Azkargorta M, Milkiewicz P, et al. Liquid biopsy-
based protein biomarkers for risk prediction, early diagno-
sis, and prognostication of cholangiocarcinoma. J Hepatol.
```
2023;79(1):93-108.
```
54. Loveless IM, Kemp SB, Hartway KM, et al. Human pancre-
atic cancer single-cell Atlas reveals association of CXCL10+
fibroblasts and basal subtype tumor cells. Clin Cancer Res.
```
2025;31(4):756-772.
```
55. Jefremow A, Neurath MF, Waldner MJ. CRISPR/Cas9
in gastrointestinal malignancies. Front Cell Dev Biol.
```
2021;9:727217. Published 2021 Nov 29. doi:10.3389/
```
fcell.2021.727217
56. Schoeps B, Eckfeld C, Prokopchuk O, et al. TIMP1 triggers
neutrophil extracellular trap formation in pancreatic cancer.
```
Cancer Res. 2021;81(13):3568-3579.
```
57. Xu Z, Zhou Y, Liu S, et al. KHSRP stabilizes m6A-Modi-
fied transcripts to activate FAK signaling and promote pan-
creatic ductal adenocarcinoma progression. Cancer Res.
```
2024;84(21):3602-3616.
```
58. Swietlik JJ, Bärthel S, Falcomatà C, et al. Cell-selective
proteomics segregates pancreatic cancer subtypes by extra-
cellular proteins in tumors and circulation. Nat Commun.
```
2023;14(1):2642.
```
59. Carpenter ES, Kadiyala P, Elhossiny AM, et al. KRT17high/
CXCL8+ tumor cells display both classical and basal features
and regulate myeloid infiltration in the pancreatic cancer micro-
```
environment. Clin Cancer Res. 2024;30(11):2497-2513.
```
60. Venning FA, Wullkopf L, Erler JT. Targeting ECM disrupts
```
cancer progression. Front Oncol. 2015;5:224.
```
61. Fukushima N, Kikuchi Y, Nishiyama T, et al. Periostin depo-
sition in the stroma of invasive and intraductal neoplasms of
```
the pancreas. Mod Pathol. 2008;21(8):1044-1053.
```
62. Lustosa SA, Viana Lde S, Affonso RJ, Jr, et al. Expression
profiling using a cDNA array and immunohistochemistry for
the extracellular matrix genes FN-1, ITGA-3, ITGB-5, MMP-
2, and MMP-9 in colorectal carcinoma progression and dis-
```
semination. ScientificWorldJournal. 2014;2014:102541.
```
63. Raghavan S, Winter PS, Navia AW, et al. Microenvironment
drives cell state, plasticity, and drug response in pancreatic
```
cancer. Cells. 2021;184(25):6119-6137.e26.
```
64. Hwang WL, Jagadeesh KA, Guo JA, et al. Single-nucleus
and spatial transcriptome profiling of pancreatic cancer iden-
tifies multicellular dynamics associated with neoadjuvant
```
treatment. Nat Genet. 2022;54(8):1178-1191.
```
14 Cancer Informatics
65. Wong CH, Lou UK, Fung FK, et al. CircRTN4 promotes
pancreatic cancer progression through a novel CircRNA-
miRNA-lncRNA pathway and stabilizing epithelial-mes-
```
enchymal transition protein. Mol Cancer. 2022;21(1):10.
```
Published 2022 Jan 4. doi:10.1186/s12943-021-01481-w
66. Chen D, Huang H, Zang L, Gao W, Zhu H, Yu X. Development
and verification of the hypoxia- and immune-associated
prognostic signature for pancreatic ductal adenocarcinoma.
```
Front Immunol. 2021;12:728062. Published 2021 Oct 6.
```
67. Kandimalla R, Xu J, Link A, et al. EpiPanGI dx: A cell-free
DNA methylation fingerprint for the early detection of gas-
```
trointestinal cancers. Clin Cancer Res. 2021;27(22):6135-
```
6144.
68. Storrs EP, Chati P, Usmani A, et al. High-dimensional decon-
struction of pancreatic cancer identifies tumor microenvi-
ronmental and developmental stemness features that predict
```
survival. NPJ Precis Oncol. 2023;7(1):105. Published 2023
```
Oct 19. doi:10.1038/s41698-023-00455-z
69. Fu Y, Tao J, Liu T, et al. Unbiasedly decoding the tumor
microenvironment with single-cell multiomics analysis in
```
pancreatic cancer. Mol Cancer. 2024;23(1):140. Published
```
2024 Jul 9. doi:10.1186/s12943-024-02050-7
70. Chen K, Wang Q, Li M, et al. Single-cell RNA-seq reveals
dynamic change in tumor microenvironment during
pancreatic ductal adenocarcinoma malignant progression.
```
EBioMedicine. 2021;66:103315.
```
71. Li J, Yuan S, Norgard RJ, et al. Epigenetic and transcriptional
control of the epidermal growth factor receptor regulates
the tumor immune microenvironment in pancreatic cancer.
```
Cancer Discov. 2021;11(3):736-753.
```
72. Pandey S, Gupta VK, Lavania SP. Role of epigenetics in pan-
```
creatic ductal adenocarcinoma. Epigenomics. 2023;15(2):89-
```
110.
73. Peng T, Sun F, Yang JC, et al. Novel lactylation-related sig-
nature to predict prognosis for pancreatic adenocarcinoma.
##### World J Gastroenterol. 2024;30(19):2575-2602.