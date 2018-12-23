#!/applications/R/R-3.3.2/bin/Rscript

# Plot bar chart of log2(observed:expected) peaks overlapping other features

# Usage:
# /applications/R/R-3.3.2/bin/Rscript other_features_vs_MSH4_Rep1_genome_wide_peaks.R "MSH4 peaks" "wt MSH4 Rep1" "wt_MSH4_Rep1" 10000

library(ggplot2)
library(ggthemes)

dataName <- "MSH4 peaks"
pt1MSH4Name <- "wt MSH4 Rep1"
pt1MSH4LibName <- "wt_MSH4_Rep1"
# Number of permutations (randomisations) performed
perms <- 10000

args <- commandArgs(trailingOnly = T)
dataName <- as.character(args[1])
pt1MSH4Name <- as.character(args[2])
pt1MSH4LibName <- as.character(args[3])
# Number of permutations (randomisations) performed
perms <- as.numeric(args[4])

MSH4Dir1 <- "/home/ajt200/analysis/160902_Sasha_ChIP_MSH4_Rep1/peaks/PeakRanger1.18/ranger/p0.001_q0.01/MSH4_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
plotDir <- "/home/ajt200/analysis/160902_Sasha_ChIP_MSH4_Rep1/peaks/PeakRanger1.18/ranger/p0.001_q0.01/MSH4_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/plots/"

otherNames <- c("REC8_HA_Rep1GR", "REC8_HA_Rep2GR",
                "REC8_MYC_Rep1GR", "kss_REC8_HA_Rep1GR",
                "nucleRnucsGR", "rangernucsGR",
                "SPO11GR", "SPO11_ChIP4GR", "SPO11_ChIP13GR",
                "H3K4me3GR", "H3K9me2GR", "H3K9me2GRbcp", "COsGR",
                "genesGR", "promotersGR", "terminatorsGR",
                "TSSdownstream500GR", "TTSupstream500GR",
                "exonsGR", "intronsGR", "TEsGR")
otherNamesPlot <- c("REC8-HA Rep1 peaks", "REC8-HA Rep2 peaks",
                    "REC8-MYC Rep1 peaks", "kss REC8-HA Rep1 peaks",
                    "Nucleosomes", "Nucleosomes (ranger)",
                    "SPO11-1-oligo hotspots", "SPO11-1 ChIP peaks", "SPO11-1 ChIP13 peaks",
                    "H3K4me3 peaks", "H3K9me2 peaks", "H3K9me2 peaks (bcp)", "Crossovers",
                    "Genes", "Gene promoters", "Gene terminators",
                    "Gene 5' ends", "Gene 3' ends",
                    "Gene exons", "Gene introns", "Transposons")

otherNames <- otherNames[c(-2:-4, -6, -9, -12)]
otherNamesPlot <- otherNamesPlot[c(-2:-4, -6, -9, -12)]

load(paste0(MSH4Dir1, "permTest_MSH4_Rep1_rangerPeaks_vs_others.RData"))
pt1_MSH4 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL
pt1_MSH4 <- pt1_MSH4[c(-2:-4, -6, -9, -12)]

# pt1_MSH4
pt1_MSH4_Pval <- lapply(seq_along(pt1_MSH4), function(x) {
  pt1_MSH4[[x]]$numOverlaps$pval
})
pt1_MSH4_Obs <- lapply(seq_along(pt1_MSH4), function(x) {
  pt1_MSH4[[x]]$numOverlaps$observed
})
pt1_MSH4_Perm <- lapply(seq_along(pt1_MSH4), function(x) {
  pt1_MSH4[[x]]$numOverlaps$permuted
})
pt1_MSH4_Exp <- lapply(seq_along(pt1_MSH4), function(x) {
  mean(pt1_MSH4[[x]]$numOverlaps$permuted)
})
pt1_MSH4_log2ObsExp <- lapply(seq_along(pt1_MSH4_Obs), function(x) {
  log2(pt1_MSH4_Obs[[x]]/pt1_MSH4_Exp[[x]])
})
pt1_MSH4_Zscore <- lapply(seq_along(pt1_MSH4), function(x) {
  pt1_MSH4[[x]]$numOverlaps$zscore
})
pt1_MSH4_AltHyp <- lapply(seq_along(pt1_MSH4), function(x) {
  pt1_MSH4[[x]]$numOverlaps$alternative
})
pt1_MSH4_alpha0.05 <- lapply(seq_along(pt1_MSH4_Perm), function(x) {
  if(pt1_MSH4_AltHyp[[x]] == "greater") {
    quantile(pt1_MSH4_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1_MSH4_Perm[[x]], 0.05)[[1]]
  }
})
pt1_MSH4_log2alpha0.05 <- lapply(seq_along(pt1_MSH4_alpha0.05), function(x) {
  log2(pt1_MSH4_alpha0.05[[x]]/pt1_MSH4_Exp[[x]])
})

pt1_MSH4_log2ObsExp_sorted <- unlist(pt1_MSH4_log2ObsExp[sort.int(unlist(pt1_MSH4_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_MSH4_log2alpha0.05_sorted <- unlist(pt1_MSH4_log2alpha0.05[sort.int(unlist(pt1_MSH4_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_MSH4_otherNames_sorted <- otherNames[sort.int(unlist(pt1_MSH4_log2ObsExp), decreasing = T, index.return = T)$ix]
pt1_MSH4_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(pt1_MSH4_log2ObsExp), decreasing = T, index.return = T)$ix]

df <- data.frame(Sample = rep(c(pt1MSH4Name),
                              each = length(pt1_MSH4_log2ObsExp_sorted)),
                 Annotation_feature = pt1_MSH4_otherNamesPlot_sorted,
                 log2ObsExp = c(pt1_MSH4_log2ObsExp_sorted),
                 log2alpha0.05 = c(pt1_MSH4_log2alpha0.05_sorted))

df$Annotation_feature <- factor(df$Annotation_feature,
                                levels = c(pt1_MSH4_otherNamesPlot_sorted))
df$Sample <- factor(df$Sample,
                    levels = c(pt1MSH4Name))

bp <- ggplot(data = df,
             mapping = aes(x = Annotation_feature,
                           y = log2ObsExp,
                           fill = Sample)) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(name = "Sample",
                        values = c("black"),
                        labels = c(pt1MSH4Name)) +
      geom_point(mapping = aes(Annotation_feature, log2alpha0.05),
                 position = position_dodge(0.9),
                 shape = "-", colour  = "red", size = 7) +
      labs(x = "Annotation feature",
           y = expression("Log"[2]*"(observed:expected) peak overlap")) +
      theme_bw() +
      theme(axis.line.y = element_line(size = 0.5, colour = "black"),
            axis.ticks.y = element_line(size = 0.25, colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 8),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0(dataName, " (", as.character(perms), " permutations)"))
ggsave(paste0(plotDir, "barplot_other_features_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_",
              pt1MSH4LibName, "_peaks.pdf"),
       plot = bp,
       height = 4.5, width = 7)
save(bp,
     file = paste0(plotDir, "barplot_other_features_permTestResults_",
                   as.character(perms), "perms_",
                   "log2_Observed_Expected_",
                   pt1MSH4LibName, "_peaks.RData"))
