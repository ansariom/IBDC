#!/usr/bin/Rscript

args = commandArgs(trailingOnly = T)

infile1 = "~/Downloads/ibdc/jan2018/tile_performances.txt"
infile2 = "~/Downloads/ibdc/jan2018/roe_performances.txt"
name1 = "Tile"
name2 = "ROE"


infile1 = args[1]
infile2 = args[2]
name1 = args[3]
name2 = args[4]
model_type= args[5]
outfile_auc = paste(model_type, "performance.plots.ROE-Tile.png", sep = "")
#outfile_prc = paste(infile, ".plots.auprc.png", sep = "")

df1 <- read.table(infile1, header = F, col.names = c("auROC", "auPRC", "Expr_Level", "fold_change"))
df2 <- read.table(infile2, header = F, col.names = c("auROC", "auPRC", "Expr_Level", "fold_change"))

library(ggplot2)
library(tidyr)

df1$fold_change <- as.factor(df1$fold_change)
df1 <- gather(df1, perf_metric, value, -fold_change, -Expr_Level)
p1 = ggplot(df1, aes(fold_change, value, color = fold_change)) + geom_boxplot() + facet_grid(perf_metric ~ Expr_Level) + ggtitle(paste(name1, " Performance", sep = "")) + ylim(0.8,0.96)

df2$fold_change <- as.factor(df2$fold_change)
df2 <- gather(df2, perf_metric, value, -fold_change, -Expr_Level)
p2 = ggplot(df2, aes(fold_change, value, color = fold_change)) + geom_boxplot() + facet_grid(perf_metric ~ Expr_Level) + ggtitle(paste(name2, " Performance", sep = "")) + ylim(0.8,0.96)

library(grid)
library(gridExtra)

g = grid.arrange(p1, p2, nrow = 1)
g
ggsave(g, file = outfile_auc, dpi = 320, width = 12)

