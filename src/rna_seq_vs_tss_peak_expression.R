library(ggplot2)
library(gridExtra)
library(tidyr)

root_dir = "~/Downloads/ibdc/"
expr_file = paste(root_dir, "mean_norm_leaf_root.txt", sep = "")
peaks_file = paste(root_dir, "July2018/aligned.peaks.annotated.capped.filtered", sep = "")
gff_file = paste(root_dir, "July2018/TAIR10_GFF3_genes_geneIds_geneType.txt", sep = "")

leaf_ann_file <- "~/Downloads/ibdc/July2018/aligned.peaks.annotated_leaf.capped"
root_ann_file <- "~/Downloads/ibdc/July2018/aligned.peaks.annotated_root.capped"

expr <- read.table(expr_file, sep = "\t", header = T)
peaks <- read.delim(peaks_file, sep = ",", header = T)
gff <- read.delim(gff_file, sep = "\t", header = F)

leaf_ann <- read.delim(leaf_ann_file, sep = ",", header = T)[,c("TranscriptID", "TranscriptLocation","X..Capped")]
root_ann <- read.delim(root_ann_file, sep = ",", header = T)[,c("TranscriptID", "TranscriptLocation","X..Capped")]

leaf_ann$tissue = "leaf"
root_ann$tissue = "root"

expr <- expr[!grepl("ATC|ATM", expr$target_id),]
noTss <- expr[!expr$target_id %in% peaks$TranscriptID & expr$mean_leaf_norm > 1 & expr$mean_root_norm > 1,]
noTss <- extract(noTss, target_id, "gene_id", regex = "([^.]+).\\d", remove = F)
noTss <- merge(noTss, gff, by.x = "gene_id", by.y = "V10")
gene_groups <- aggregate(noTss$V9, list(noTss$V9), length)
noTss <- noTss[noTss$V9 == "protein_coding_gene",]

noTss$avg_expr <- (noTss$mean_leaf_norm + noTss$mean_root_norm )/ 2
noTss <- noTss[,c(1,2,3,4,13,14)]
m <- table(cut(noTss$avg_expr, breaks = c(1, 5, 50, 1000)))
m

high_expr <- noTss[noTss$avg_expr > 100,]
h <- merge(high_expr, leaf_ann, by.x = "target_id", by.y = "TranscriptID", all.x = T)
h <- merge(h, root_ann, by.x = "target_id", by.y = "TranscriptID", all.x = T)
hist(h$X..Capped.x)
hist(h$X..Capped.y)

aggregate(h$TranscriptLocation.x, list(h$TranscriptLocation.x), length)

plot(noTss$mean_leaf_norm, noTss$mean_root_norm)
absent_tss <- noTss[noTss$mean_leaf_norm > 1000 | noTss$mean_root_norm > 1000,]


df <- merge(peaks, expr, by.x = "TranscriptID", by.y = "target_id")
plot(df$mean_leaf_norm, df[df$tissue == "leaf",]$ReadCount)

R = cor(leaf_tss_rna_seq$ReadCount, leaf_tss_rna_seq$mean_leaf_norm)
R = format(round(R, 2), nsmall = 2)
leaf <- ggplot(leaf_tss_rna_seq, aes(ReadCount, mean_leaf_norm)) + 
  geom_point() + geom_smooth(method = lm, se = FALSE) + 
  ggtitle("TSS-seq expression vs. RNA-seq Expression in LEAF") +
  xlab("TSS peak read count") + ylab("RNA-seq mean-normalized count") + 
  annotate(geom = "text", label = paste("pearson-cor = ", R, sep = ""), color = "blue", x = 1500000, y = 5000)

R = cor(root_tss_rna_seq$ReadCount, root_tss_rna_seq$mean_leaf_norm)
R = format(round(R, 2), nsmall = 2)
root <- ggplot(root_tss_rna_seq, aes(ReadCount, mean_leaf_norm)) + 
  geom_point() + geom_smooth(method = lm, se = FALSE) + 
  ggtitle("TSS-seq expression vs. RNA-seq Expression in ROOT") +
  xlab("TSS peak read count") + ylab("RNA-seq mean-normalized count") + 
  annotate(geom = "text", label = paste("pearson-cor = ", R, sep = ""), color = "blue", x = 1000000, y = 1000)
grid.arrange(leaf, root, nrow = 1)

#------
leaf_tss_rna_seq <- leaf_tss_rna_seq[leaf_tss_rna_seq$ReadCount > 0 & leaf_tss_rna_seq$mean_leaf_norm > 0,]
root_tss_rna_seq <- root_tss_rna_seq[root_tss_rna_seq$ReadCount > 0 & root_tss_rna_seq$mean_root_norm > 0,]

a <- leaf_tss_rna_seq[leaf_tss_rna_seq$ReadCount > 10 & leaf_tss_rna_seq$mean_leaf_norm < 10,]
a <- root_tss_rna_seq[root_tss_rna_seq$ReadCount > 10 & root_tss_rna_seq$mean_root_norm > 10,]

par(mfrow=c(2,2))
hist(leaf_tss_rna_seq[leaf_tss_rna_seq$ReadCount < 200,]$ReadCount, breaks = 100, xlab = "TSS reads", main = "TSS-peaks leaf")
hist(leaf_tss_rna_seq[leaf_tss_rna_seq$mean_leaf_norm < 200,]$mean_leaf_norm, breaks = 100, xlab = "rna-seq reads", main = "RNA-seq leaf")
hist(root_tss_rna_seq[root_tss_rna_seq$ReadCount < 200,]$ReadCount, breaks = 100, xlab = "TSS reads", main = "RNA-seq root")
hist(root_tss_rna_seq[root_tss_rna_seq$mean_root_norm < 200,]$mean_root_norm, breaks = 100, xlab = "rna-seq reads", main = "RNA-seq root")

plot(leaf_tss_rna_seq[leaf_tss_rna_seq$mean_leaf_norm < 200,]$ReadCount, leaf_tss_rna_seq[leaf_tss_rna_seq$mean_leaf_norm < 200,]$mean_leaf_norm)


leaf_tss_rna_seq$tissue = "leaf"
leaf_tss_rna_seq$rnaSeq_reads <- leaf_tss_rna_seq$mean_leaf_norm
root_tss_rna_seq$tissue = "root"
root_tss_rna_seq$rnaSeq_reads <- root_tss_rna_seq$mean_root_norm
all <- rbind(leaf_tss_rna_seq, root_tss_rna_seq)

tbl <- aggregate(all$ReadCount, list(all$TranscriptID, all$tissue), max)
colnames(tbl) <- c("TranscriptID", "tissue", "ReadCount")
library(tidyr)
a <- spread(tbl, tissue, ReadCount, fill = 0)
leaf_g <- ggplot(a, aes(a$leaf, a$root)) +  geom_point() + xlab("TSS Read Count (Leaf)") + ylab("TSS Read Count (Root)") +
  xlim(1,max(a$root))
leaf_g

tbl <- aggregate(all$rnaSeq_reads, list(all$TranscriptID, all$tissue), max)
colnames(tbl) <- c("TranscriptID", "tissue", "ReadCount")
a <- spread(tbl, tissue, ReadCount, fill = 0)
root_g <- ggplot(a, aes(a$leaf, a$root)) +  geom_point() + xlab("RNA-seq Read Count (Leaf)") + ylab("RNA-seq Read Count (Root)") 
root_g
grid.arrange(leaf_g, root_g, nrow=1)

######## diff expression of TSSs in the absence of RNA-seq

