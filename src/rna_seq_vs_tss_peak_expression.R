library(ggplot2)
library(gridExtra)

load("~/Downloads/leaf_tss_rna_seq.rdat")
load("~/Downloads/root_tss_rna_seq.rdat")

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

