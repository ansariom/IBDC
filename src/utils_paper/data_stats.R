load("~/Downloads/all_features_diffs_wide.rdat")

which(is.na(all_features_diffs_wide), arr.ind = T)

diffs_colnames <- c("gene_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs",
                    "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq",
                    "tss_name", "chr", "loc", "strand", "offset?", "")
d <- all_features_diffs_wide[, !colnames(all_features_diffs_wide) %in% diffs_colnames]
mean_d <- as.data.frame(apply(d, MARGIN = 2, FUN = mean))
mean_d$feature <- rownames(mean_d)
colnames(mean_d)[1] <- "mean_value"
# oc only
oc_mean_leaf <- mean_d[grepl("_P_LEAF", mean_d$feature),]
oc_mean_leaf$tissue <- "LEAF"
oc_mean_root <- mean_d[grepl("_P_ROOT", mean_d$feature),]
oc_mean_root$tissue <- "ROOT"
all_oc <- rbind(oc_mean_leaf, oc_mean_root)

library(ggplot2)
ggplot(all_oc, aes(x = mean_value)) + geom_histogram() + facet_wrap(~tissue) +
  ggtitle("Avg %openness in all diff expressed promoters")

library(tidyr)
all_oc <- all_oc[!grepl("OVERALL", all_oc$feature),]
all_oc <- extract(all_oc, feature, into = c("strand", "win"), regex = "\\D+(FWD|REV)_(.)_\\D+", remove = F)
ggplot(all_oc, aes(x = mean_value)) + geom_histogram() + facet_wrap(~tissue + win) +
  ggtitle("Avg %openness in all diff expressed promoters")
