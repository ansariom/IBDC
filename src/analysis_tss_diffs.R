
load("all_features_diffs_wide.rdat")
all_labeled <- subset(all_features_diffs_wide, (all_features_diffs_wide[4] > 4 | all_features_diffs_wide[4] < -4))[1:14]
all_labeled <- all_labeled[,c(1,4,12,13,14)]
all_labeled$modLoc <- as.numeric(all_labeled$loc) + 2000

root_data <- read.csv("aligned.peaks.annotated.capped_root.filtered")
leaf_data <- read.csv("aligned.peaks.annotated.capped_leaf.filtered")

root_data <- root_data[, c(3,4,5,6,7,8,10)]
leaf_data <- leaf_data[, c(3,4,5,6,7,8,10)]

all_root <- merge(all_labeled, root_data, by.x = c("gene_id", "modLoc"), by.y = c("TranscriptID", "ModeLocation"), all.x = T)
all_root_leaf <- merge(all_root, leaf_data, by.x = c("gene_id", "modLoc"), by.y = c("TranscriptID", "ModeLocation"), all.x = T)

########## Feature extraction ########

all_labeled <- subset(all_features_diffs_wide, (all_features_diffs_wide[4] > 4 | all_features_diffs_wide[4] < -4))
feature_all <- all_labeled[16:ncol(all_labeled)]
feature_all$gene_id <- all_labeled[,1]
feature_all$foldchange <- all_labeled[,4]

root_spec <- subset(feature_all, feature_all$foldchange > 4)
leaf_spec <- subset(feature_all, feature_all$foldchange < -4)


root_spec$gene_id <- NULL
feature_all$foldchange <- all_labeled[,4] 
root_mean <- apply(root_spec, 2, mean)
root_mean <- sort(root_mean, decreasing=T)

dr <- as.data.frame(root_mean)
dr$fac_win <- row.names(dr)
row.names(dr) <- NULL
dr$factor <- lapply(dr$fac_win, function(x) { gsub("_FWD_\\d+", "", x)} )
dr$factor <- lapply(dr$factor, function(x) { gsub("_REV_\\d+", "", x)} )
#dr$factor <- lapply(dr$factor, function(x) { gsub("_OC_P_", "_OC", x)} )

colnames(dr) = c("loglik", "factor_win", "factor")
m = cbind(dr$loglik, dr$factor)
write.table(m, file="root_mean.txt", sep = "\t", quote=F, row.names=F)
root_mean_data <- read.table("root_mean.txt", sep="\t", header=T)
colnames(root_mean_data) <- c("loglik_score", "fac")
rmax_g <- aggregate(root_mean_data[,1], by=list(root_mean_data$fac), max)
colnames(rmax_g) <- c("fac_root", "max_score")
ordered_rmax <- rmax_g[order(-rmax_g$max_score),]
write.table(ordered_rmax, "top_facs_root.txt", sep = "\t", quote=F, row.names=F)

#---- Leaf
leaf_mean <- apply(leaf_spec, 2, mean)
leaf_mean <- sort(leaf_mean, decreasing=T)

dl <- as.data.frame(leaf_mean)
dl$fac_win <- row.names(dl)
row.names(dl) <- NULL
dl$factor <- lapply(dl$fac_win, function(x) { gsub("_FWD_\\d+", "", x)} )
dl$factor <- lapply(dl$factor, function(x) { gsub("_REV_\\d+", "", x)} )
#dr$factor <- lapply(dr$factor, function(x) { gsub("_OC_P_", "_OC", x)} )

colnames(dl) = c("loglik", "factor_win", "factor")
m = cbind(dl$loglik, dl$factor)
write.table(m, file="leaf_mean.txt", sep = "\t", quote=F, row.names=F)
leaf_mean_data <- read.table("leaf_mean.txt", sep="\t", header=T)
colnames(leaf_mean_data) <- c("loglik_score", "fac")
lmax_g <- aggregate(leaf_mean_data[,1], by=list(leaf_mean_data$fac), max)
colnames(lmax_g) <- c("fac_leaf", "max_score")
ordered_lmax <- lmax_g[order(-lmax_g$max_score),]
write.table(ordered_lmax, "top_facs_leaf.txt", sep = "\t", quote=F, row.names=F)
