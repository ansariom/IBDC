root_rsem <- read.delim("~/Downloads/ibdc/July2018/root.isoforms.results", header = T, sep = "\t")
leaf_rsem <- read.delim("~/Downloads/ibdc/July2018/leaf.isoforms.results", header = T, sep = "\t")
kallisto <- read.delim("~/Downloads/ibdc/mean_norm_leaf_root.txt", sep = "\t", header = T)

rsem <- merge(leaf_rsem[,c("transcript_id", "TPM")], root_rsem[,c("transcript_id", "TPM")], by = "transcript_id", suffixes = c(".leaf",".root"))
library(tidyr)
rsem <- extract(rsem, transcript_id, c("transcript_id"), regex = "([^_]+)[^_]*")

all <- merge(rsem, kallisto, by.x = "transcript_id", by.y = "target_id")

a = cor.test(all$TPM.leaf, all$mean_leaf_norm)
b = cor.test(all$TPM.root, all$mean_root_norm)

par(mfrow = c(1,2))
plot(all$TPM.leaf, all$mean_leaf_norm, xlab = "RSEM gene expr", ylab = "Kallisto gene expression", main = paste("Leaf pearson cor-coef = ", a$estimate[[1]], sep = ""))
plot(all$TPM.root, all$mean_root_norm, xlab = "RSEM gene expr", ylab = "Kallisto gene expression", main = paste("Root pearson cor-coef = ", b$estimate[[1]], sep = ""))



