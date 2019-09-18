#!/usr/bin/Rscript
library(tidyr)
args = commandArgs(trailingOnly = T)

indir = args[1]
outbase = args[2]

fs <- list.files(indir)
all <- c(auROC=c(), auPRC=c(), group = c())
for (f in fs) {
  infile = paste(indir, "/", f, sep = "")
  df <- read.table(infile, header = F, col.names = c("auROC", "auPRC"))
  if (nrow(df) == 0) next
  df$group <- basename(infile)
  all <- rbind(df, all)
}
print(head(all))
all <- extract(all, group, into = "group", regex = "performances_top(\\d+).txt")
all$group <- as.factor(as.numeric(all$group))
library(ggplot2)
outfile = paste(outbase, "_performance_plot_auROC.png", sep = "")
g = ggplot(all, aes(group, auROC)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(g, file = outfile)

outfile = paste(outbase, "_performance_plot_auPRC.png", sep = "")
g = ggplot(all, aes(group, auPRC)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(g, file = outfile)


#-------
