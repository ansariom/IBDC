#############################################################################
### Testing root only 3peat model on leaf  only test sets and vice versa ####
#############################################################################
library(ggplot2)
library(plotly)
# ----- Test1: What are misclassified genes in each model? -----#
coef_root <- read.table("Downloads/ALL100_Root-not-Leaf_arab_allPWMs_Train_-5000_5000.20_negfrom_-2000_-200.TAIR10_cds_227.simple.SumScoreVars.0.00215.model.UnModCoefs.plot.table.fac.summary")
coef_leaf <- read.table("Downloads/ALL100_Leaf-not-Root_arab_allPWMs_Train_-5000_5000.20_negfrom_-2000_-200.TAIR10_cds_270.simple.SumScoreVars.0.00144.model.UnModCoefs.plot.table.fac.summary")

coef_merged <- merge(coef_root, coef_leaf, by.x = "V1", by.y = "V1", all = T)

coef_merged$V2.x <- NULL
coef_merged$V2.y <- NULL

coef_merged[is.na(coef_merged)] <- 0

colnames(coef_merged) <- c("factor", "root_coef", "leaf_coef")
ggplot(coef_merged, aes(root_coef, leaf_coef)) + geom_point(aes(colour = factor)) 
ggplotly()

root_facs <- coef_merged[coef_merged[,"leaf_coef"] == 0,]
root_facs <- root_facs[order(-root_facs$root_coef),]
root_facs$leaf_coef <- NULL
write.table(root_facs, file="Downloads/root_only_top_factors.txt", sep = "\t", quote = F, row.names = F)

leaf_facs <- coef_merged[coef_merged[,"root_coef"] == 0,]
leaf_facs <- leaf_facs[order(-leaf_facs$leaf_coef),]
leaf_facs$root_coef <- NULL
write.table(leaf_facs, file="Downloads/leaf_only_top_factors.txt", sep = "\t", quote = F, row.names = F)
