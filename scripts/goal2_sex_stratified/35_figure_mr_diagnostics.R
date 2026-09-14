#!/usr/bin/env Rscript
# Rebuild MR diagnostic figures (forest, SNP support, funnel, leave-one-out) from the current MR objects, per sex.
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})
proc <- "data/processed"; tab <- "results/tables"
figd <- "results/figures/current"; dir.create(figd, showWarnings=FALSE, recursive=TRUE)

o    <- readRDS(file.path(proc,"new","MR_primary_objects.rds"))
prim <- as.data.table(o$primary); dat <- as.data.table(o$dat)
fs   <- list(female=fread(file.path(tab,"FS_input_female.csv"))$gene,
             male  =fread(file.path(tab,"FS_input_male.csv"))$gene)
th <- theme_bw(base_size=12) + theme(
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.border=element_rect(colour="black", linewidth=.4, fill=NA),
        axis.text=element_text(colour="black"),
        axis.title=element_text(face="plain"),
        legend.background=element_blank(),
        legend.key=element_blank(),
        strip.background=element_rect(fill="grey92", colour="black", linewidth=.4))

# instrument availability (drives what is drawable)
nsnp_tab <- dat[mr_keep==TRUE, .N, by=gene]
nsnp_v   <- setNames(as.integer(nsnp_tab$N), nsnp_tab$gene)
lookup_n <- function(g){ v <- as.integer(nsnp_v[g]); v[is.na(v)] <- 0L; v }
avail <- rbindlist(lapply(names(fs), function(s){
  d <- data.table(sex=s, gene=fs[[s]]); d[, n_snp := lookup_n(gene)]; d }))
avail[, diagnostics_possible := fifelse(n_snp>=3, "funnel/LOO/Egger",
                                fifelse(n_snp==2, "IVW only", "Wald ratio only"))]
fwrite(avail, file.path(tab,"MR_diagnostics_availability.csv"))
cat("=== instrument availability among prioritised genes ===\n")
print(avail[, .N, by=.(sex, diagnostics_possible)][order(sex, -N)])

for (s in names(fs)){
  P <- prim[gene %in% fs[[s]]]
  P[, n_snp := lookup_n(gene)]
  P[, dirn := fifelse(OR>1, "risk (OR>1)", "protective (OR<1)")]

  # FIG 1: forest of all prioritised genes
  Pf <- P[order(OR)]; Pf[, gene := factor(gene, levels=gene)]
  g <- ggplot(Pf, aes(OR, gene, colour=dirn)) +
    geom_vline(xintercept=1, linetype=2, colour="grey40", linewidth=.4) +
    geom_errorbarh(aes(xmin=OR_lo, xmax=OR_hi), height=0, linewidth=.5) +
    geom_point(size=2) + scale_x_log10() +
    scale_colour_manual(values=c("risk (OR>1)"="#B2182B","protective (OR<1)"="#2166AC"), name=NULL) +
    labs(x="Odds ratio for RA (log scale)", y=NULL) +
    th + theme(axis.text.y=element_text(size=7.5), legend.position="top",
               legend.text=element_text(size=10))
  ggsave(file.path(figd, sprintf("FIG_MR_01_forest_%s.png", s)), g,
         width=7.2, height=max(5, nrow(Pf)*0.16), dpi=600, limitsize=FALSE)
  ggsave(file.path(figd, sprintf("FIG_MR_01_forest_%s.pdf", s)), g,
         width=7.2, height=max(5, nrow(Pf)*0.16), limitsize=FALSE)

  # FIG 2: instrument support
  A <- avail[sex==s]
  g2 <- ggplot(A, aes(factor(n_snp))) + geom_bar(fill=NA, colour="#2166AC", linewidth=1.1, width=.3) +
    geom_text(stat="count", aes(label=after_stat(count)), vjust=-.5, size=3.5) +
    scale_y_continuous(expand=expansion(mult=c(0, .08))) +
    labs(x="Number of cis-eQTL instruments retained after harmonisation",
         y="Number of genes") +
    th + theme(panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey90", linewidth=.3))
  ggsave(file.path(figd, sprintf("FIG_MR_02_snp_support_%s.png", s)), g2,
         width=6.4, height=4.4, dpi=600)
  ggsave(file.path(figd, sprintf("FIG_MR_02_snp_support_%s.pdf", s)), g2,
         width=6.4, height=4.4)

  # FIG 3/4: funnel + leave-one-out for genes with >=3 SNPs
  multi <- A[n_snp>=3]$gene
  D <- dat[gene %in% multi & mr_keep==TRUE]
  if (nrow(D) > 0){
    D[, wald := beta.outcome/beta.exposure]
    D[, prec := abs(beta.exposure)/se.outcome]      # instrument strength / precision
    gf <- ggplot(D, aes(wald, prec, colour=gene)) +
      geom_point(size=2, alpha=.85) +
      geom_vline(xintercept=0, linetype=2, colour="grey45", linewidth=.4) +
      labs(x="SNP-specific Wald ratio (beta outcome / beta exposure)",
           y="Instrument precision |beta exposure| / SE outcome") +
      th + theme(legend.position=if(length(multi)<=12) "right" else "none",
                 panel.grid.major.x=element_line(colour="grey90", linewidth=.3))
    ggsave(file.path(figd, sprintf("FIG_MR_03_funnel_%s.png", s)), gf,
           width=7.6, height=5.4, dpi=600)
    ggsave(file.path(figd, sprintf("FIG_MR_03_funnel_%s.pdf", s)), gf,
           width=7.6, height=5.4)

    # leave-one-out: recompute IVW omitting each SNP, per gene
    loo <- rbindlist(lapply(multi, function(gn){
      d <- D[gene==gn]
      rbindlist(lapply(seq_len(nrow(d)), function(i){
        dd <- d[-i]
        w  <- 1/(dd$se.outcome^2)
        b  <- sum(w*dd$beta.outcome*dd$beta.exposure)/sum(w*dd$beta.exposure^2)
        data.table(gene=gn, omitted=d$SNP[i], OR=exp(b)) }))
    }), fill=TRUE)
    allp <- rbindlist(lapply(multi, function(gn){
      d <- D[gene==gn]; w <- 1/(d$se.outcome^2)
      data.table(gene=gn, omitted="(all SNPs)",
                 OR=exp(sum(w*d$beta.outcome*d$beta.exposure)/sum(w*d$beta.exposure^2))) }))
    L <- rbind(loo, allp)
    gl <- ggplot(L, aes(OR, omitted)) +
      geom_vline(xintercept=1, linetype=2, colour="grey45", linewidth=.4) +
      geom_point(aes(colour=omitted=="(all SNPs)"), size=2.1) +
      scale_colour_manual(values=c(`TRUE`="#B2182B",`FALSE`="grey30"), guide="none") +
      facet_wrap(~gene, scales="free_y") + scale_x_log10() +
      labs(x="IVW odds ratio omitting each SNP (log scale)", y=NULL) +
      th + theme(axis.text.y=element_text(size=7), strip.text=element_text(size=8))
    ggsave(file.path(figd, sprintf("FIG_MR_04_leaveoneout_%s.png", s)), gl,
           width=9, height=max(4.5, ceiling(length(multi)/3)*2.1), dpi=600, limitsize=FALSE)
    ggsave(file.path(figd, sprintf("FIG_MR_04_leaveoneout_%s.pdf", s)), gl,
           width=9, height=max(4.5, ceiling(length(multi)/3)*2.1), limitsize=FALSE)
    cat(sprintf("  %-6s : forest(%d genes), funnel+LOO(%d genes >=3 SNPs)\n",
                s, nrow(Pf), length(multi)))
  } else cat(sprintf("  %-6s : forest(%d genes); NO gene has >=3 SNPs -> funnel/LOO not drawable\n",
                     s, nrow(Pf)))
}
cat("\nWrote MR diagnostic figures to", figd, "\n")
