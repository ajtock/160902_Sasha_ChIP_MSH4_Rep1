#!/applications/R/R-3.3.2/bin/Rscript

# Plot bar chart of log2(observed:expected) peaks overlapping other features

# Usage:
# /applications/R/R-3.3.2/bin/Rscript TE_family_vs_MSH4_Rep1_vs_SPO11oligo_genome_wide_peaks.R "MSH4 peaks and SPO11-1-oligo hotspots" "wt MSH4 Rep1" "wt_MSH4_Rep1" "wt SPO11-1-oligos" "wt_SPO11oligo" 10000

library(ggplot2)
library(ggthemes)

#dataName <- "MSH4 peaks and SPO11-1-oligo hotspots"
#pt1MSH4Name <- "wt MSH4 Rep1"
#pt1MSH4LibName <- "wt_MSH4_Rep1"
#pt1Name <- "wt SPO11-1-oligos"
#pt1LibName <- "wt_SPO11oligo"
## Number of permutations (randomisations) performed
#perms <- 10000

args <- commandArgs(trailingOnly = T)
dataName <- as.character(args[1])
pt1MSH4Name <- as.character(args[2])
pt1MSH4LibName <- as.character(args[3])
pt1Name <- as.character(args[4])
pt1LibName <- as.character(args[5])
# Number of permutations (randomisations) performed
perms <- as.numeric(args[6])

MSH4Dir1 <- "/home/ajt200/analysis/160902_Sasha_ChIP_MSH4_Rep1/peaks/PeakRanger1.18/ranger/p0.001_q0.01/MSH4_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
inDir1 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/"
plotDir <- "/home/ajt200/analysis/160902_Sasha_ChIP_MSH4_Rep1/peaks/PeakRanger1.18/ranger/p0.001_q0.01/MSH4_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/plots/"

famNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger",
              "rna", "gypsy", "copia", "linel1", "sine")
famNamesPlot <- c("DNA", "Helitron", "Pogo/Tc1/mariner", "MuDR", "EnSpm", "hAT", "Harbinger",
                  "RNA", "Gypsy LTR", "Copia LTR", "LINE-1", "SINE")

load(paste0(MSH4Dir1, "permTest_MSH4_Rep1_rangerPeaks_vs_TEsDNA.RData"))
pt1_MSH4_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(MSH4Dir1, "permTest_MSH4_Rep1_rangerPeaks_vs_TEsRNA.RData"))
pt1_MSH4_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt1_MSH4 <- c(pt1_MSH4_DNA, pt1_MSH4_RNA)

load(paste0(inDir1, "permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_TEsDNA.RData"))
pt1_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(inDir1, "permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_TEsRNA.RData"))
pt1_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt1 <- c(pt1_DNA, pt1_RNA)

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
pt1_MSH4_famNames_sorted <- famNames[sort.int(unlist(pt1_MSH4_log2ObsExp), decreasing = T, index.return = T)$ix]
pt1_MSH4_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(pt1_MSH4_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt1
pt1_Pval <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$pval
})
pt1_Obs <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$observed
})
pt1_Perm <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$permuted
})
pt1_Exp <- lapply(seq_along(pt1), function(x) {
  mean(pt1[[x]]$numOverlaps$permuted)
})
pt1_log2ObsExp <- lapply(seq_along(pt1_Obs), function(x) {
  log2(pt1_Obs[[x]]/pt1_Exp[[x]])
})
pt1_Zscore <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$zscore
})
pt1_AltHyp <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$alternative
})
pt1_alpha0.05 <- lapply(seq_along(pt1_Perm), function(x) {
  if(pt1_AltHyp[[x]] == "greater") {
    quantile(pt1_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1_Perm[[x]], 0.05)[[1]]
  }
})
pt1_log2alpha0.05 <- lapply(seq_along(pt1_alpha0.05), function(x) {
  log2(pt1_alpha0.05[[x]]/pt1_Exp[[x]])
})

pt1_log2ObsExp_sorted <- unlist(pt1_log2ObsExp[sort.int(unlist(pt1_MSH4_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_log2alpha0.05_sorted <- unlist(pt1_log2alpha0.05[sort.int(unlist(pt1_MSH4_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_famNames_sorted <- famNames[sort.int(unlist(pt1_MSH4_log2ObsExp), decreasing = T, index.return = T)$ix]
pt1_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(pt1_MSH4_log2ObsExp), decreasing = T, index.return = T)$ix]


df <- data.frame(Sample = rep(c(pt1MSH4Name, pt1Name),
                              each = length(pt1_MSH4_log2ObsExp_sorted)),
                 Transposon_family = rep(pt1_MSH4_famNamesPlot_sorted, 2),
                 log2ObsExp = c(pt1_MSH4_log2ObsExp_sorted,
                                pt1_log2ObsExp_sorted),
                 log2alpha0.05 = c(pt1_MSH4_log2alpha0.05_sorted,
                                   pt1_log2alpha0.05_sorted))

df$Transposon_family <- factor(df$Transposon_family,
                               levels = c(pt1_MSH4_famNamesPlot_sorted))
df$Sample <- factor(df$Sample,
                    levels = c(pt1MSH4Name, pt1Name))

bp <- ggplot(data = df,
             mapping = aes(x = Transposon_family,
                           y = log2ObsExp,
                           fill = Sample)) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(name = "Sample",
                        values = c("red",
                                   "blue"),
                        labels = c(pt1MSH4Name,
                                   pt1Name)) +
      geom_point(mapping = aes(Transposon_family, log2alpha0.05),
                 position = position_dodge(0.9),
                 shape = "-", colour  = "grey70", size = 7) +
      labs(x = "Transposon superfamily",
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
ggsave(paste0(plotDir, "barplot_TE_family_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_",
              pt1MSH4LibName, "_peaks_",
              pt1LibName, "_peaks.pdf"),
       plot = bp,
       height = 4.5, width = 7)
save(bp,
     file = paste0(plotDir, "barplot_TE_family_permTestResults_",
                   as.character(perms), "perms_",
                   "log2_Observed_Expected_",
                   pt1MSH4LibName, "_peaks_",
                   pt1LibName, "_peaks.RData"))
