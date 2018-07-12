#!/usr/bin/Rscript

args = commandArgs(trailingOnly = T)

infile = args[1]
name = args[2]
outfile_auc = paste(infile, ".plots.auroc.png", sep = "")
outfile_prc = paste(infile, ".plots.auprc.png", sep = "")

df <- read.table(infile, header = F, col.names = c("auROC", "auPRC", "Expr_Level", "fold_change"))

library(ggplot2)

p <- ggplot(df, aes(fold_change, auROC)) + geom_boxplot() + facet_wrap(~Expr_Level) + ggtitle(paste(name, " AUROC", sep = ""))
ggsave(p, file = outfile_auc, dpi = 320)

p <- ggplot(df, aes(fold_change, auPRC)) + geom_boxplot() + facet_wrap(~Expr_Level) + ggtitle(paste(name, " AUPRC", sep = ""))
ggsave(p, file = outfile_prc, dpi = 320)
