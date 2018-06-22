#!/usr/bin/Rscript
args = commandArgs(trailingOnly = T)

indir = args[1]
outbase = args[2]

fs <- list.files(indir)
all <- c(auROC=c(), auPRC=c(), group = c())
for (f in fs) {
  infile = paste(indir, "/", f, sep = "")
  df <- read.table(infile, header = F, col.names = c("auROC", "auPRC"))
  df$group <- basename(infile)
  all <- rbind(df, all)
}

library(ggplot2)
outfile = paste(outbase, "_performance_plot_auROC.png", sep = "")
png(outfile)
boxplot(all$group, all$auROC)
dev.off()



#-------
