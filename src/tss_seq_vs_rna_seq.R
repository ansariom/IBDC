library(tidyr)
all_peaks_root <- read.csv("~/Downloads/ibdc/Aug2018/aligned.peaks.annotated_root.capped")
all_peaks_leaf <- read.csv("~/Downloads/ibdc/Aug2018/aligned.peaks.annotated_leaf.capped")

peaks <- read.csv("~/Downloads/ibdc/Aug2018/aligned.peaks.annotated.capped.filtered")
rnaseq <- read.table("~/Downloads/ibdc/Aug2018/ath_root_leaf_rsem_deseq_diff_expr_results.txt", header = T)
gff <- read.table("~/Downloads/ibdc/July2018/TAIR10_GFF3_genes_geneIds_geneType.txt", col.names = c("Chr", "V2", "V3", "Start", "End", "dot1", "Strand", "dot2", "type", "geneID"))

root <- peaks[peaks$tissue == "root",]
leaf <- peaks[peaks$tissue == "leaf",]


m1 <- merge(root[,c("ReadCount", "TranscriptID")], rnaseq, by.x = "TranscriptID", by.y = "Accession", all.y = T)
m1[is.na(m1)] <- 0
m1$tissue.y <- "root"
nonZero <- merge(leaf[,c("ReadCount", "TranscriptID")], m1, by = "TranscriptID", all.y = T, suffixes = c(".leaf", ".root"))
nonZero[is.na(nonZero)] <- 0
nonZero$tissue.x = "leaf"
all_zer0 <- nonZero[nonZero$ReadCount.root == 0 & nonZero$ReadCount.leaf == 0 & nonZero$baseMean == 0,]
nonZero <- nonZero[!nonZero$TranscriptID %in% all_zer0$TranscriptID,]
nonZero <- extract(nonZero, TranscriptID, into = c("geneID"), regex = "([^.]+).\\d", remove = F)
nonZero <- merge(nonZero, gff[,c("geneID", "type")], by = "geneID")
aggregate(nonZero$type, list(nonZero$type), length)


noTSS_yRNA <- rnaseq[!(rnaseq$Accession %in% root$TranscriptID | rnaseq$Accession %in% leaf$TranscriptID), ]
noTSS_yRNA <- extract(noTSS_yRNA, Accession, into = c("geneID"), regex = "([^.]+).\\d", remove = F)
noTSS_yRNA <- noTSS_yRNA[noTSS_yRNA$baseMean > 0,]
noTSS_yRNA <- merge(noTSS_yRNA, gff[,c("geneID", "type")], by = "geneID")
aggregate(noTSS_yRNA$type, list(noTSS_yRNA$type), length)
pcg <- noTSS_yRNA[noTSS_yRNA$type == "protein_coding_gene",]
pcg <- as.data.frame(table(cut(noTSS_yRNA$baseMean, breaks = seq(0, max(pcg$baseMean), by = 10))))
pcg <- pcg[1:80,]

no_rnaseq <- nonZero[nonZero$baseMean == 0 & (nonZero$ReadCount.root > 0 | nonZero$ReadCount.leaf > 0),]
no_rnaseq$tss_fc <- log10(no_rnaseq$ReadCount.root / no_rnaseq$ReadCount.leaf)

# Focus on the following data
all_pcg_noTSS <- noTSS_yRNA[noTSS_yRNA$type == "protein_coding_gene" & (noTSS_yRNA$mean_leaf_norm > 300 | noTSS_yRNA$mean_root_norm > 300),]
all_pcg_noTSS <- merge(all_peaks_leaf[,c("ReadCount", "TranscriptID", "TranscriptLocation", "X..Capped")], all_pcg_noTSS,  by.x = "TranscriptID", by.y = "Accession")
all_pcg_noTSS <- merge(all_peaks_root[,c("ReadCount", "TranscriptID", "TranscriptLocation", "X..Capped")], all_pcg_noTSS,  by.x = "TranscriptID", by.y = "Accession")

all_pcg_noTSS$eff_reads <- all_pcg_noTSS$ReadCount * (all_pcg_noTSS$X..Capped/100)
all_pcg_noTSS <- all_pcg_noTSS[!all_pcg_noTSS$TranscriptLocation %in% c('coding', "3'utr", "intron"),]
all_pcg_noTSS[all_pcg_noTSS$eff_reads > 50,]
###------

all_main <- nonZero[!nonZero$TranscriptID %in% noTSS_yRNA$Accession,]
all_main <- all_main[!all_main$TranscriptID %in% no_rnaseq$TranscriptID,]
all_main$ReadCount.leaf <- ifelse(all_main$ReadCount.leaf == 0, 1, all_main$ReadCount.leaf)
all_main$ReadCount.root <- ifelse(all_main$ReadCount.root == 0, 1, all_main$ReadCount.root)
all_main$tss_fc <- all_main$ReadCount.root / all_main$ReadCount.leaf
all_main$tss_diff <- log2(all_main$tss_fc)
aggregate(all_main$type, list(all_main$type), length)

inf <- all_main[all_main$b %in% c("Inf", "-Inf"),]
all_main$mean_leaf_norm <- all_main$mean_leaf_norm + 1
all_main$mean_root_norm <- all_main$mean_root_norm + 1
all_main$foldChange <- all_main$mean_root_norm / all_main$mean_leaf_norm
all_main$b <- log2(all_main$foldChange)

final_output <- unique(all_main[,c("TranscriptID", "baseMean", "mean_leaf_norm", "mean_root_norm", "foldChange", "b", "pval", "qval")])
colnames(final_output) <- c("Accession", "baseMean", "mean_leaf_norm", "mean_root_norm", "foldChange", "b", "pval", "qval")
write.table(final_output, file = "ath_root_leaf_rsem_deseq_diff_expr_results_filtered.txt", quote = F, row.names = F, sep = "\t")

#plot(all_main$b, all_main$b1)
plot(all_main$tss_diff, all_main$b, xlab = "TSS foldChange (log2)", ylab = "RNA-seq foldChange (log2)", main = "Correlation plot between RNA-seq and TSS-seq expression")

agree <- all_main[(all_main$b > 0 & all_main$tss_diff > 0) | (all_main$b < 0 & all_main$tss_diff < 0),]
disagree <- all_main[(all_main$b > 0 & all_main$tss_diff < 0) | (all_main$b < 0 & all_main$tss_diff >0),]

disagree <- all_main[abs(all_main$b - all_main$tss_diff) > 10 & all_main$baseMean > 100,]

# Disagreements can be for two reasons: 1- two tramscripts that have same peak assigned but realy diff expressoin 2- peak pile up count error
# Let's check the first issue
disagree_fun <- function(i, disagree, all_main) {
  tID = disagree[i, "TranscriptID"]
  gID = disagree[i, "geneID"]
  print(tID)
  a <- all_main[all_main$geneID == gID & all_main$TranscriptID != tID,]
  print(a)
  count <-  nrow(a[abs(a$b1 - a$tss_diff) > 5,])
  print(count)
  return(count)
}

lapply(seq(1,nrow(disagree)), disagree_fun, disagree, all_main)
do.call(rbind, lapply(seq(1,nrow(disagree)), disagree_fun, all_main))

plot(all_main$ReadCount.root, all_main$mean_root_norm)
cor.test(all_main$tss_diff, all_main$b)

library(ggplot2)
ggplot(pcg) + geom_col(aes(pcg$Var1, pcg$Freq)) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + xlab("mean expression")+ ylab("Frequency") +
  ggtitle("Expression values for protein coding genes in which there is no TSS-peak associated with them")


