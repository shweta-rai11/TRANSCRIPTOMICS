#!/usr/bin/env Rscript
# =============================================================================
# 23_methodology_flowchart.R  —  publication-ready methodology flowchart for the
# sex-stratified RA biomarker pipeline (base graphics; corrected & complete).
# =============================================================================
fig <- "results/figures"; dir.create(fig, showWarnings = FALSE, recursive = TRUE)

png(file.path(fig, "fig_methodology_flowchart.png"),
    width = 2000, height = 2500, res = 200)
par(mar = c(0, 0, 0, 0)); plot.new(); plot.window(xlim = c(0, 100), ylim = c(0, 100))

# ---- helpers ---------------------------------------------------------------
box <- function(x, y, w, h, label, fill = "white", border = "#1F3B5C",
                cex = 0.72, font = 1, tcol = "black", lwd = 1.6) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2, col = fill, border = border, lwd = lwd)
  text(x, y, label, cex = cex, font = font, col = tcol)
}
arr <- function(x0, y0, x1, y1, col = "#2E5E8C", lwd = 2.4)
  arrows(x0, y0, x1, y1, length = 0.09, angle = 22, lwd = lwd, col = col)

# palette by stage
c_data <- "#D6E4F0"; c_pre <- "#EDEDED"; c_deg <- "#D9EAD3"
c_wg <- "#E6DAF0"; c_ml <- "#FCE5CD"; c_val <- "#FFF2CC"
c_down <- "#F4CCCC"; c_final <- "#FFD54A"

# ================= MAIN SPINE (centre x = 38) ==============================
cx <- 38
box(cx, 95, 34, 6.5, "Discovery cohorts (NCBI GEO)\nGSE93272 + GSE110169  — whole blood, RA vs HC", c_data, cex=0.68, font=2)
box(cx, 87, 30, 4.5, "Merge (common genes)", c_pre)
box(cx, 80.5, 30, 4.5, "Quantile normalization", c_pre)
box(cx, 74, 30, 5, "ComBat batch correction\n(protect group + sex)", c_pre, cex=0.68)
box(cx, 66, 34, 5.5, "Differential expression (DEG)\nAll  •  Female  •  Male", c_deg, cex=0.72, font=2)
arr(cx,91.8,cx,89.3); arr(cx,84.8,cx,82.8); arr(cx,78.3,cx,76.5); arr(cx,71.5,cx,68.8)

box(cx, 57.5, 34, 5.5, "Feature selection (sex-stratified)\nLASSO + RF + SVM-RFE + Boruta", c_ml, cex=0.68)
box(cx, 50.5, 30, 4.8, "Female panel (15)  •  Male panel (8)", c_ml, cex=0.7, font=2)
box(cx, 43.5, 30, 5, "Diagnostic models\nLR / RF / SVM / KNN / ANN", c_ml, cex=0.68)
box(cx, 36.5, 30, 4.5, "Nested cross-validation", c_ml)
box(cx, 29, 34, 5.5, "External validation\nBlood GSE15573  •  Synovium GSE89408", c_val, cex=0.66)
box(cx, 22, 30, 4.5, "Cross-tissue transferability", c_val)
arr(cx,63.2,cx,60.3); arr(cx,54.8,cx,52.9); arr(cx,48.1,cx,46);
arr(cx,41,cx,38.8); arr(cx,34.2,cx,31.8); arr(cx,26.2,cx,24.3)

# ================= WGCNA BRANCH (right x = 80) =============================
wx <- 80
box(wx, 57.5, 30, 4.8, "WGCNA co-expression network", c_wg, cex=0.7, font=2)
box(wx, 50.5, 30, 5, "Modules (blue = down, yellow = up)\n+ module–trait", c_wg, cex=0.66)
# Hub genes are defined by kME, GS and intramodular connectivity ONLY.
# No STRING PPI and no betweenness-centrality analysis exists in this
# pipeline; the previous label advertised an analysis that was never run.
box(wx, 43.5, 30, 4.5, "Hub genes (kME, GS, connectivity)", c_wg, cex=0.62)
box(wx, 36.5, 30, 4.5, "DEG ∩ WGCNA candidates", c_wg)
# DEG -> WGCNA branch
arr(cx+17, 66, wx-15, 57.5)
arr(wx,55.1,wx,52.9); arr(wx,48,wx,45.8); arr(wx,41.2,wx,38.8)

# ================= DOWNSTREAM (left x = 12) ================================
lx <- 12
box(lx, 45, 20, 5, "CIBERSORT\nimmune deconvolution", c_down, cex=0.66)
box(lx, 33, 20, 5.5, "Mendelian randomization\n(eQTL -> RA, causal): GGA2", c_down, cex=0.62)
box(lx, 21, 21, 6, "Comparative ancestry analysis\nEUR (Okada/Stahl) vs EAS (Sakaue/BBJ)\nGGA2 direction-heterogeneous",
    "#F9CB9C", cex=0.58, font=2)                       # ancestry comparison box
box(wx, 29, 26, 4.8, "DEG n WGCNA n MR\nintegration", c_wg, cex=0.66)
# inputs: signatures -> left rail -> CIBERSORT -> MR -> ancestry  (orthogonal)
seg <- function(x0,y0,x1,y1) segments(x0,y0,x1,y1, col="#2E5E8C", lwd=2.4)
seg(cx-15, 50.5, lx, 50.5); arr(lx, 50.5, lx, 47.6)   # signatures -> CIBERSORT
arr(lx, 42.5, lx, 35.9)                               # CIBERSORT -> MR
arr(lx, 30.2, lx, 24.1)                               # MR -> ancestry comparison
arr(wx, 34.2, wx, 31.4)                                # candidates -> integration

# ---- PLANNED multi-omics extension (dashed = future work, not yet done) ----
c_plan <- "#ECECEC"
box_d <- function(x, y, w, h, label, cex = 0.56) {
  rect(x-w/2, y-h/2, x+w/2, y+h/2, col=c_plan, border="#7A7A7A", lwd=1.8, lty=2)
  text(x, y, label, cex=cex, font=3, col="grey20") }
darr <- function(x0,y0,x1,y1)
  arrows(x0,y0,x1,y1, length=0.08, angle=22, lwd=1.8, col="#7A7A7A", lty=2)
box_d(26, 14, 28, 6.5, "Single-cell RNA-seq  (PLANNED)\ncell-type resolution of the\nsex-stratified biomarkers")
c_meth <- "#CDE7DF"
box(64, 14.2, 40, 7.6,
    "DNA-METHYLATION VALIDATION  (designed - see methylomics figure)\nGSE42861 blood  |  GSE111942 blood-F  |  GSE80071 synovium-F\nsex-stratified EWAS  ->  intersect 23 biomarker loci  ->  concordance\nQ1 methylation?   Q2 cross-tissue (blood<->synovium)?   Q3 consistent in both sexes?",
    c_meth, cex=0.46, font=1, border="#2E7D6E", lwd=2)

# ================= FINAL: two biomarker categories =========================
box(27, 4.3, 40, 7,
    "TISSUE-SPECIFIC BIOMARKERS  (non-cross-tissue)\nMale blood: GGA2*, AP1B1, UBAP2\nMale synovium: STARD3NL, POM121, TTC7A",
    "#F6C34D", cex=0.54, font=2, border="#B8860B", lwd=2.4)
box(72, 4.3, 40, 7,
    "CROSS-TISSUE BIOMARKERS  (blood <-> synovium)\nFemale: CNIH4, EXOC6, HMGB2,\nBMX, BCAT1, C5orf30",
    c_final, cex=0.54, font=2, border="#B8860B", lwd=2.4)

# cross-tissue transferability -> multi-omics boxes (split rail)
seg(cx, 19.7, cx, 18); seg(26, 18, 64, 18)
arr(26, 18, 26, 17.3); arr(64, 18, 64, 18.1)
# scRNA (dashed, planned) + methylomics (solid, designed) -> final biomarkers
darr(26, 10.75, 27, 8);  arr(64, 10.4, 72, 8)
# solid causal / integration -> final
arr(lx, 18, lx, 7.9)                                  # ancestry/MR (GGA2) -> tissue-specific
arr(wx, 26.6, wx, 7.9)                                # integration -> cross-tissue
text(50, 0.4, "* GGA2 = causally validated (MR).   Methylomics = designed (3 datasets verified, not yet run).   Single-cell RNA-seq = planned.",
     cex=0.48, font=3, col="grey30")

# title + legend
text(50, 99.4, "Sex-stratified rheumatoid arthritis biomarker discovery pipeline",
     cex=0.95, font=2)
legend("bottomright", inset=c(0.01,0.01), bty="n", cex=0.56,
  legend=c("Data","Preprocessing","DEG","WGCNA branch","ML biomarker","Validation","Causal / immune","Final","Methylomics (designed)","Planned (future)"),
  fill=c(c_data,c_pre,c_deg,c_wg,c_ml,c_val,c_down,c_final,c_meth,c_plan), border="grey40")
dev.off()
cat("Wrote fig_methodology_flowchart.png\n")
