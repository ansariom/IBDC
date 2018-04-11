load("~/Downloads/all_features_diffs_wide.rdat")

which(is.na(all_features_diffs_wide), arr.ind = T)

diffs_colnames <- c("gene_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs",
                    "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq",
                    "tss_name", "chr", "loc", "strand", "offset?", "")
d <- all_features_diffs_wide[, !colnames(all_features_diffs_wide) %in% diffs_colnames]
mean_d <- as.data.frame(apply(d, MARGIN = 2, FUN = mean))
mean_d$feature <- rownames(mean_d)
colnames(mean_d)[1] <- "mean_value"

library(tidyr)
# Extract feature types
mean_d <- extract(mean_d, feature, into = c("pwm" ,"strand", "win", "type"), regex = "(.+)_(FWD|REV)_(.)(\\D*)", remove = F)
mean_d$type[mean_d$type == ""] <- "TFBS"
mean_d$type[mean_d$type == "_OC_P_LEAF"] <- "OC_LEAF"
mean_d$type[mean_d$type == "_OC_P_ROOT"] <- "OC_ROOT"

tfbs <- mean_d[mean_d$type == "TFBS",]
hist(tfbs$mean_value, main = "Histogram of TFBS scores", xlab = "mean log-lik score (normalized)")

oc_closed_leaf <- mean_d[mean_d$mean_value < 0.3 & mean_d$type == "OC_LEAF",]
oc_open_leaf <- mean_d[mean_d$mean_value > 0.5 & mean_d$type == "OC_LEAF",]
oc_closed_root <- mean_d[mean_d$mean_value < 0.3 & mean_d$type == "OC_ROOT",]

# close OC and TFBS score
m_leaf <- merge(tfbs, oc_closed_leaf, by = c("pwm", "strand", "win"))
hist(m_leaf$mean_value.x)
high_tfbs_in_closed_leaf <- m_leaf[m_leaf$mean_value.x > 0.7,]

m_root <- merge(tfbs, oc_closed_root, by = c("pwm", "strand", "win"))
hist(m_root$mean_value.x)
high_tfbs_in_closed_root <- m_root[m_root$mean_value.x > 0.7,]
write.table(high_tfbs_in_closed_leaf, file = "~/Downloads/closed_strong_tfbs_sites_leaf.txt", quote = F, row.names = F, sep = "\t")
write.table(high_tfbs_in_closed_root, file = "~/Downloads/closed_strong_tfbs_sites_root.txt", quote = F, row.names = F, sep = "\t")

ggplot(high_tfbs_in_closed_leaf, aes(feature.x, mean_value.y)) + geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Top weighted features and Avg of % openness (LEAF)") + xlab("feature") + ylab("avg %open")
ggplot(high_tfbs_in_closed_root, aes(feature.x, mean_value.y)) + geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Top weighted features and Avg of % openness (ROOT)") + xlab("feature") + ylab("avg %open")

# OC stats for top 50 features
coef_table <- read.table("~/Downloads/ibdc/coefs/top0_coefs.txt", header = F, col.names = c("feature", "coef"))
topX_coef <- coef_table[1:50,]
topX_coef <- extract(topX_coef, feature, into = c("pwm" ,"strand", "win"), regex = "(.+)_(FWD|REV)_(.)\\D*", remove = F)

oc_leaf <- mean_d[mean_d$type == "OC_LEAF",]
oc_root <- mean_d[mean_d$type == "OC_ROOT",]
top_mean_oc_leaf <- merge(topX_coef, oc_leaf, by = c("pwm", "strand", "win"))
top_mean_oc_root <- merge(topX_coef, oc_root, by = c("pwm", "strand", "win"))

hist(top_mean_oc_leaf$mean_value, main = "OC openness for top 50 features (LEAF OC)", xlab = "mean %OC Openness")
hist(top_mean_oc_root$mean_value, main = "OC openness for top 50 features (ROOT OC)", xlab = "mean %OC Openness")

closed_features_leaf <- top_mean_oc_leaf[top_mean_oc_leaf$mean_value < 0.3,]
closed_features_root <- top_mean_oc_root[top_mean_oc_root$mean_value < 0.3,]

# oc only
oc_mean_leaf <- mean_d[grepl("_P_LEAF", mean_d$feature),]
oc_mean_leaf$tissue <- "LEAF"
oc_mean_root <- mean_d[grepl("_P_ROOT", mean_d$feature),]
oc_mean_root$tissue <- "ROOT"
all_oc <- rbind(oc_mean_leaf, oc_mean_root)

library(ggplot2)
library(tidyr)
all_oc <- all_oc[!grepl("OVERALL", all_oc$feature),]
all_oc <- extract(all_oc, feature, into = c("strand", "win"), regex = "\\D+(FWD|REV)_(.)_\\D+", remove = F)
ggplot(all_oc, aes(x = mean_value)) + geom_histogram() + facet_wrap(~tissue + win) +
  ggtitle("Avg %openness in all diff expressed promoters")
