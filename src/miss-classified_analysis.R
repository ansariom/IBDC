library(ggplot2)

load("~/Downloads/featureInfo_hardCodedSoftCoded.rdat")
low_diffs_classes <- diffs_classes
low_diffs_classes <- low_diffs_classes[low_diffs_classes$mean_root_norm >=30 |low_diffs_classes$mean_leaf_norm >=30, ]
load("~/Downloads/ibdc/Aug2018/med_high/2.5/147842/featureInfo_hardCodedSoftCoded.rdat")
high_diffs_classes <- diffs_classes
high_diffs_classes <- high_diffs_classes[high_diffs_classes$mean_root_norm >=300 |high_diffs_classes$mean_leaf_norm >=300, ]

l_classified_data = low_diffs_classes[(low_diffs_classes["class"] != -1000),]
l_classified_data$data_type <- "med-low"
h_classified_data = high_diffs_classes[(high_diffs_classes["class"] != -1000),]
h_classified_data$data_type <- "med-high"
classified_data <- rbind(l_classified_data, h_classified_data)

ggplot(classified_data, aes(mean_root_norm, mean_leaf_norm)) + geom_point(shape = 1) +
  ggtitle("mean normalized expression levele for diff expressed genes") + facet_wrap(~data_type)

ex_low <- l_classified_data
ex_high <- h_classified_data
ex_low$mean_root_leaf_expr <- (ex_low$mean_leaf_norm + ex_low$mean_root_norm)/2
ex_high$mean_root_leaf_expr <- (ex_high$mean_leaf_norm + ex_high$mean_root_norm)/2
ex_low <- data.frame(table(cut(ex_low$mean_root_leaf_expr, breaks = seq(1, 500, 10))))
ex_low$bin <- seq(1, nrow(ex_low))
ex_high <- data.frame(table(cut(ex_high$mean_root_leaf_expr, breaks = seq(1, 500, 10))))
ex_high$bin <- seq(1, nrow(ex_high))
ex_low$data_type <- "med-low"
ex_high$data_type <- "med-high"
expr_levels <- rbind(ex_high, ex_low)
ggplot(expr_levels, aes(bin, Freq)) + geom_bar(stat = "identity") + xlab("mean normalized expression")+ facet_wrap(~data_type)


missclassified_1 = diffs_classes[(diffs_classes["prob0"]>0.5) & (diffs_classes["class"]==1),]
missclassified_0 = diffs_classes[(diffs_classes["prob1"]>0.5) & (diffs_classes["class"]==0),]
missclassified = rbind(missclassified_0, missclassified_1)
expr_miss <- missclassified
expr_miss$expr_mean <- (expr_miss$mean_root_norm+expr_miss$mean_leaf_norm)/2
ggplot(expr_miss, aes(mean_root_norm, mean_leaf_norm)) + geom_point(shape = 1) +
  ggtitle("mean normalized expression level for misclassified genes")+ facet_wrap(~data_type)

highest_diff_misclassed = expr_miss[expr_miss$mean_leaf_norm > 500 | expr_miss$mean_root_norm > 500,]

correct_class_rowids = setdiff(row.names(classified_data),row.names(missclassified))
correctly_classified = classified_data[correct_class_rowids,]

#### Plot the distribution of probabolities that leads to miss-classification
library(ggplot2)
missclassified <- within(missclassified, class_diff <- abs(prob1 - prob0))
ggplot(missclassified, aes(prob0, prob1)) + geom_point(shape = 1, aes(colour = class_diff)) + 
  ggtitle("Prob. of root (0) vs leaf(1) predictions for miss-classified data")  + facet_wrap(~data_type)

ggplot(missclassified) + geom_histogram(aes(prob0, fill = data_type, alpha = 0.6), position = "identity") + ggtitle("Histpgram of prob0 in misclassified")
ggplot(missclassified) + geom_histogram(aes(prob1, fill = data_type, alpha = 0.6), position = "identity") + ggtitle("Histpgram of prob1 in misclassified")
ggplot(missclassified) + geom_density(aes(prob0, fill = data_type, alpha = 0.6), position = "identity", stat = "density") + ggtitle("Density of prob0 in misclassified")
ggplot(missclassified) + geom_density(aes(prob1, fill = data_type, alpha = 0.6), position = "identity", stat = "density") + ggtitle("Density of prob1 in misclassified")


extremely_missclassified <- missclassified[missclassified["class_diff"] > 0.7, ]
extremely_missclassified <- highest_diff_misclassed
miss_ids = row.names(extremely_missclassified)

m1 <- missclassified[missclassified$mean_root_norm > 300 & missclassified$mean_leaf_norm > 300,]
aggregate(m1$data_type, list(m1$data_type), length)
# ---------------------- ---------------------- ---------------------- ----------------------
### Where the features fall for missclassified as compared to correctly-classified feature dist.
# ---------------------- ---------------------- ---------------------- ----------------------
library(reshape)
correct_features <- features[correct_class_rowids,]
correct_features = cbind(correct_features, class = correctly_classified[,"class"] [match(row.names(correct_features), row.names(correctly_classified))])

extremely_missclassified <- missclassified[missclassified["class_diff"] > 0.7, ]
extremely_missclassified <- highest_diff_misclassed
miss_ids = row.names(extremely_missclassified)

feature_means <- aggregate(correct_features, by = list(correct_features[,"class"]), mean)
feature_means <- as.data.frame(t(feature_means))
colnames(feature_means) <- c("true_mean_root", "true_mean_leaf")
#colnames(feature_means) <- feature_means[1,]
feature_means = feature_means[-1,]
feature_means$feature_id <- row.names(feature_means)
#feature_means = melt(feature_means, id=c("feature_id"))

miss_features <- features[miss_ids,]
miss_features <- cbind(miss_features, class = extremely_missclassified[,"class"] [match(row.names(extremely_missclassified), row.names(miss_features))])
miss_features$gene_id <- row.names(miss_features)
miss_features = melt(miss_features, id=c("gene_id"))
colnames(miss_features) <- c("gene_id", "feature_id", "feature_value")
d <- extremely_missclassified[, c("tss_name", "class")]
miss_features <- merge(miss_features, d, by.x = "gene_id", by.y = "tss_name")

# Here you go! All we need for analysis is here!
miss_feature_merged1 <- merge(feature_means, miss_features, by.x = "feature_id", by.y = "feature_id")
miss_feature_merged <- merge(miss_feature_merged1, feature_info, by.x = "feature_id", by.y = "feature")

miss_feature_merged$mean_root_p <- miss_feature_merged$true_mean_root * miss_feature_merged$coefficient
miss_feature_merged$mean_leaf_p <- miss_feature_merged$true_mean_leaf * miss_feature_merged$coefficient
miss_feature_merged$featureval_p <- miss_feature_merged$feature_value * miss_feature_merged$coefficient

miss_feature_merged$featurevalp_diff_root = abs(miss_feature_merged$featureval_p - miss_feature_merged$mean_root_p)
miss_feature_merged$featurevalp_diff_leaf = abs(miss_feature_merged$featureval_p - miss_feature_merged$mean_leaf_p)
miss_feature_merged$featurevalp_diff_both = miss_feature_merged$featurevalp_diff_leaf + miss_feature_merged$featurevalp_diff_root

# Now separate OC from TFBS features and threshold them separately
OC_subset <- miss_feature_merged[miss_feature_merged$type == "OC",]
TFBS_subset <- miss_feature_merged[miss_feature_merged$type == "SLL",]

ordered_oc = OC_subset[order(-OC_subset[,"featurevalp_diff_both"]),]
ordered_tfbs = TFBS_subset[order(-TFBS_subset[,"featurevalp_diff_both"]),]

# Does OC in specific window matter?
selected_ordered_oc <- ordered_oc[1:100,]
ggplot(selected_ordered_oc, aes(as.numeric(selected_ordered_oc$window), fill = factor(tissue))) + geom_histogram(bins = 14) +
  ggtitle("Freq of observation of highly different OC values between \n correctly classified vs. miss-classified genes ") + 
  xlab("window_number") + ylab("frequency")

ggplot(selected_ordered_oc, aes(selected_ordered_oc$pwm, fill = factor(tissue))) + geom_bar() +
  ggtitle("Freq of observation of highly different OC values between \n correctly classified vs. miss-classified genes ") + 
  xlab("PWM") + ylab("frequency") + facet_wrap(~gene_id) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


selected_ordered_tfbs <- ordered_tfbs[1:100,]
ggplot(selected_ordered_tfbs, aes(as.numeric(selected_ordered_tfbs$window), fill = factor(tissue))) + geom_histogram(bins = 14) +
  ggtitle("Freq of observation of very different TFBS scores between \n correctly classified vs. miss-classified genes ") + 
  xlab("window_number") + ylab("frequency")

# Merge   selected OC and tfbs and look at gene level
selected_both <- rbind(selected_ordered_oc, selected_ordered_tfbs)
selected_both = selected_both[order(selected_both[,"gene_id"], selected_both[,"feature_id"]),]
selected_both[is.na(selected_both[,"tissue"]), "tissue"] = "tfbs"
selected_both[selected_both[,"tissue"] == "ROOT", "tissue"] = "OC_root"
selected_both[selected_both[,"tissue"] == "LEAF", "tissue"] = "OC_leaf"

selected_root <- selected_both[selected_both[,"class"] == 0,]
selected_leaf <- selected_both[selected_both[,"class"] == 1,]

ggplot(selected_root, aes(factor(pwm) ,fill = tissue)) + geom_bar() +
  facet_wrap(~gene_id, ncol=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("root-specific miss_classified genes")

ggplot(selected_leaf, aes(factor(pwm) ,fill = tissue)) + geom_bar() +
  facet_wrap(~gene_id, ncol=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("leaf-specific miss_classified genes")

selected_both$label <- paste(selected_both$gene_id, selected_both$class, sep = "*")
ggplot(selected_both, aes(factor(pwm) ,fill = tissue)) + geom_bar() +
  facet_wrap(~label, ncol=2) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("miss_classified genes")

unique(selected_both$gene_id)


# ---------------------- ---------------------- ---------------------- ----------------------
# TSS analysis for miss-classified data
# ---------------------- ---------------------- ---------------------- ----------------------
leaf_peaks = read.csv("~/Downloads/aligned.peaks.annotated.capped_leaf.filtered", header = T)
root_peaks = read.csv("~/Downloads/aligned.peaks.annotated.capped_root.filtered", header = T)

root_peaks = root_peaks[, c(3,4,5,6,7,8,10)]
leaf_peaks = leaf_peaks[, c(3,4,5,6,7,8,10)]

root_related = merge(missclassified, root_peaks, by.x = c("gene_id"), by.y = c("TranscriptID"))

#----s1=features["AT2G01830.3_Chr2_366165_0",]
s1=features["AT2G01830.3_Chr2_366165_0",]
s1$id = row.names(s1)
library(reshape)
sm = melt(s1, id=(c("id")))
colnames(sm) <- c("id", "fac", "coef")
