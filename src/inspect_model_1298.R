library(tidyr)
library(ggplot2)
library(NeatMap)

# given two dataframes with similar rownames, merge (join) them using those rownames
merge_by_rownames <- function(df1, df2, ...) {
  df1$asdfasdf <- rownames(df1)
  df2$asdfasdf <- rownames(df2)
  res <- merge(df1, df2, by = "asdfasdf", ...)
  rownames(res) <- res$asdfasdf
  res$asdfasdf <- NULL
  return(res)
}

# load the data
if(TRUE) {
  setwd("~/Documents/cgrb/pis/Megraw/mitra_output/")
  load("featureinfo_features_diffsclasses.Rdata")
}

# sort the feature info by abs(coefficient)
feature_info_sorted <- feature_info[rev(order(abs(feature_info$coefficient))), ]

# grab the top 50 features, 
features_top_50 <- features[, feature_info_sorted$feature[1:50]]

# make a dataframe of just calls & true class information
prob1 <- diffs_classes[, "prob1", drop = F]
true_class <- diffs_classes[diffs_classes$class %in% c(0,1), "class", drop = F]
colnames(true_class) <- "true_class"
info <- merge_by_rownames(prob1, true_class)

# normize the data feature cols
features_top_50_scaled <- as.data.frame(scale(features_top_50))
# multiply them by the coefficients
features_top_50_scaled_by_coeffs <- as.data.frame(scale(features_top_50_scaled, scale = feature_info_sorted$coefficient[1:50]))
#features_top_50_rescaled <- as.data.frame(scale(features_top_50_scaled_by_coeffs))

# merge that data along with the class info above
features_top_50_info <- merge_by_rownames(info, features_top_50_scaled_by_coeffs)
# scale the true class and class calls into a larger range for plotting
features_top_50_info$true_class <- (features_top_50_info$true_class - 0.5) * 600
features_top_50_info$prob1 <- (features_top_50_info$prob1 - 0.5) * 600

# make a heatmap
mat <- as.matrix(features_top_50_info)
f <- make.heatmap1(profiles = mat,
              row.method = "average.linkage",
              column.method = "average.linkage",
              row.normalize = FALSE,
              column.labels = colnames(mat))
f <- f + expand_limits(y = c(-500, 2250))

plot(f)



## for molly: look for pairs of differentially expressed TSSs that are similar in profile
# the strategy here is to loop over all rows that are root-specific, and compare them to all rows that are 
# leaf-specific. That's a lot of comparisons, and we don't want to store all of them; rather we just store
# top N. This is done by keeping a dataframe of results with N rows, and if the current comparison is larger (in similarity, abs(correlation coefficient))
# than any existing row in the table, replace that row with the current comparison.
n <- 50
res_df_top_n <- data.frame(root_spec_rownum = seq(1, n), leaf_spec_rownum = seq(1, n), correlation = rep(0, n))
features_top_50 <- features_top_50_info[, !colnames(features_top_50_info) %in% c("prob1", "true_class")]
for(root_spec_rownum in which(features_top_50_info$true_class < 0)) {
  for(leaf_spec_rownum in which(features_top_50_info$true_class > 0)) {
    root_spec_values <- as.numeric(features_top_50[root_spec_rownum,])
    leaf_spec_values <- as.numeric(features_top_50[leaf_spec_rownum,])
    correlation <- cor.test(root_spec_values, leaf_spec_values)$estimate
    if(abs(correlation) > min(abs(res_df_top_n$correlation))) {
      res_df_top_n[which.min(abs(res_df_top_n$correlation)), ] <- c(root_spec_rownum, leaf_spec_rownum, correlation)
    }
  }
  print(paste(root_spec_rownum/nrow(features_top_50), "%"))
}


print(res_df_top_n)
# order the above by the correlations
res_df_top_n <- res_df_top_n[order(abs(res_df_top_n$correlation)), ]

# Assign to the results various useful stuff for plotting - especially their location along the y axis (numid)
# so they can be plotted in groups
root_patterns <- features_top_50[res_df_top_n$root_spec_rownum, ]
root_patterns$numid <- seq(1, nrow(root_patterns))
root_patterns$id <- paste(rownames(features_top_50[res_df_top_n$root_spec_rownum, ]), "root_specific")
root_patterns$correlation <- res_df_top_n$correlation
root_patterns$pattern = "root_specific"

leaf_patterns <- features_top_50[res_df_top_n$leaf_spec_rownum, ]
leaf_patterns$numid <- seq(1.5, nrow(leaf_patterns) + 0.5)
leaf_patterns$id <- paste(rownames(features_top_50[res_df_top_n$leaf_spec_rownum, ]), "leaf_specific")
leaf_patterns$correlation <- res_df_top_n$correlation
leaf_patterns$pattern = "leaf_specific"

all_patterns <- rbind(root_patterns, leaf_patterns)

all_patterns_long <- tidyr::gather(all_patterns, feature, value, -id, -correlation, -pattern, -numid)
head(all_patterns_long)

# order the features by the coefficients for plotting (sorting info gathered from feature_info_sorted above)
all_patterns_long$feature <- factor(all_patterns_long$feature, 
                                    labels = feature_info_sorted$feature, 
                                    levels = feature_info_sorted$feature,
                                    ordered = T)
  
p <- ggplot(all_patterns_long) +
  geom_tile(aes(x = feature, y = numid, fill = value)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(
    data = unique(all_patterns_long[, c("id", "numid")]),
    aes(x = -1, y = numid, label = id), hjust = 1) +
    expand_limits(x = c(-20, 50)) +
  geom_hline(data = data.frame(y = seq(0.75, nrow(leaf_patterns)+1, 1)),
             aes(yintercept = y)) +
  geom_text(data = feature_info_sorted[1:50,],
            aes(x = feature, y = 0.1, label = round(coefficient, 3), angle = 90, hjust = 1)) +
  scale_fill_gradient2() +
  scale_y_continuous(name = "Top 50 pairs of leaf & root specific TSSs that are correlated by top 50 weighted features")
plot(p)

ggsave("correlated_root_vs_leaf_pairs.pdf", p, width = 15, height = 18)



