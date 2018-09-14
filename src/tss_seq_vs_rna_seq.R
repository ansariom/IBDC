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
m <- merge(leaf[,c("ReadCount", "TranscriptID")], m1, by = "TranscriptID", all.y = T, suffixes = c(".leaf", ".root"))
m[is.na(m)] <- 0
m$tissue.x = "leaf"
all_zer0 <- m[m$ReadCount.root == 0 & m$ReadCount.leaf == 0 & m$baseMean == 0,]
m <- m[!m$TranscriptID %in% all_zer0$TranscriptID,]
m <- extract(m, TranscriptID, into = c("geneID"), regex = "([^.]+).\\d", remove = F)
m <- merge(m, gff[,c("geneID", "type")], by = "geneID")
aggregate(m$type, list(m$type), length)


noTSS_yRNA <- rnaseq[!(rnaseq$Accession %in% root$TranscriptID | rnaseq$Accession %in% leaf$TranscriptID), ]
noTSS_yRNA <- extract(noTSS_yRNA, Accession, into = c("geneID"), regex = "([^.]+).\\d", remove = F)
noTSS_yRNA <- noTSS_yRNA[noTSS_yRNA$baseMean > 0,]
noTSS_yRNA <- merge(noTSS_yRNA, gff[,c("geneID", "type")], by = "geneID")
aggregate(noTSS_yRNA$type, list(noTSS_yRNA$type), length)
pcg <- noTSS_yRNA[noTSS_yRNA$type == "protein_coding_gene",]
pcg <- as.data.frame(table(cut(noTSS_yRNA$baseMean, breaks = seq(0, max(pcg$baseMean), by = 10))))
pcg <- pcg[1:80,]

no_rnaseq <- m[m$baseMean == 0 & (m$ReadCount.root > 0 | m$ReadCount.leaf > 0),]
no_rnaseq$tss_fc <- log10(no_rnaseq$ReadCount.root / no_rnaseq$ReadCount.leaf)

# Focus on the following data
all_pcg_noTSS <- noTSS_yRNA[noTSS_yRNA$type == "protein_coding_gene" & (noTSS_yRNA$mean_leaf_norm > 300 | noTSS_yRNA$mean_root_norm > 300),]
all_pcg_noTSS <- merge(all_peaks_leaf[,c("ReadCount", "TranscriptID", "TranscriptLocation", "X..Capped")], all_pcg_noTSS,  by.x = "TranscriptID", by.y = "Accession")
all_pcg_noTSS <- merge(all_peaks_root[,c("ReadCount", "TranscriptID", "TranscriptLocation", "X..Capped")], all_pcg_noTSS,  by.x = "TranscriptID", by.y = "Accession")

all_pcg_noTSS$eff_reads <- all_pcg_noTSS$ReadCount * (all_pcg_noTSS$X..Capped/100)
all_pcg_noTSS <- all_pcg_noTSS[!all_pcg_noTSS$TranscriptLocation %in% c('coding', "3'utr", "intron"),]
all_pcg_noTSS[all_pcg_noTSS$eff_reads > 50,]
###------

all_main <- m[!m$TranscriptID %in% noTSS_yRNA$Accession,]
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
all_main$b1 <- log2(all_main$foldChange)
plot(all_main$b, all_main$b1)
plot(all_main$tss_diff, all_main$b1)

disagree <- all_main[abs(all_main$b1 - all_main$tss_diff) > 10 & all_main$baseMean > 300,]

# Disagreements can be for two reasons: 1- two tramscripts that have same peak assigned but realy diff expressoin 2- peak pile up count error
# Let's check the first issue
disagree_fun <- function(i, all_main) {
  print(i)
  tID = all_main[i, "TranscriptID"]
  print(tID)
  gID = all_main[i, "GeneID"]
  a <- all_main[all_main$geneID == gID & all_main$TranscriptID != tID,]
  return(nrow(a[abs(a$b1 - a$tss_diff) > 5,]))
}
do.call(rbind, lapply(seq(1,nrow(disagree)), disagree_fun, all_main))

plot(all_main$ReadCount.root, all_main$mean_root_norm)
cor.test(all_main$tss_diff, all_main$b)

library(ggplot2)
ggplot(pcg) + geom_col(aes(pcg$Var1, pcg$Freq)) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + xlab("mean expression")+ ylab("Frequency") +
  ggtitle("Expression values for protein coding genes in which there is no TSS-peak associated with them")


