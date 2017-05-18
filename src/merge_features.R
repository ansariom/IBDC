#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 10) {
	print("This script is meant to join 4 very specifically-formatted 'feature' files into a single file for modeling and analysis.")
	print("It uses maybe like... 25G of RAM.")
	print("")
	print("Usage: features.rdat tiled_features.rdat features_long_oc_p_leaf.txt features_long_oc_p_root.txt tiled_features_long_oc_p_leaf.txt tiled_features_long_oc_p_root.txt  features_long_oc_p_overall_leaf.txt features_long_oc_p_overall_root.txt  diff_exp_results.txt output.txt")
	print("")
	quit()
}

feat_wide_file <- args[1]
tiled_feat_wide_file <- args[2]
oc_leaf_long_file <- args[3]
oc_root_long_file <- args[4]
tiled_oc_leaf_long_file <- args[5]
tiled_oc_root_long_file <- args[6]
oc_leaf_overall_file <- args[7]
oc_root_overall_file <- args[8]
diff_file <- args[9]
output_file <- args[10]




############# 
## Read and merge features
#############

features <- read.table(feat_wide_file, header = TRUE, stringsAsFactors = FALSE)
features$tss_name <- rownames(features)
features_long <- gather(features, feature, value, -tss_name)

tiled_features <- read.table(tiled_feat_wide_file, header = TRUE, stringsAsFactors = FALSE)
tiled_features$tss_name <- rownames(tiled_features)
tiled_features_long <- gather(tiled_features, feature, value, -tss_name)

oc_features_leaf <- read.table(oc_leaf_long_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)
oc_features_root <- read.table(oc_root_long_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)

tiled_oc_features_leaf <- read.table(tiled_oc_leaf_long_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)
tiled_oc_features_root <- read.table(tiled_oc_root_long_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)

oc_features_overall_leaf <- read.table(oc_leaf_overall_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)
oc_features_overall_root <- read.table(oc_root_overall_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)

all_features_long <- rbind(features_long, tiled_features_long, oc_features_leaf, oc_features_root, tiled_oc_features_leaf, tiled_oc_features_root, oc_features_overall_leaf, oc_features_overall_root)

# to here: ~ 15 mins, 15 gigs (before tiled features addition...)
print("Step1")
all_features_wide <- spread(all_features_long, feature, value)
print("all_features_wide is ready" )

## cleanup a bit to save RAM
rm(features)
rm(features_long)
rm(oc_features_leaf)
rm(oc_features_root)
rm(all_features_long)
gc()

#####################################
### Make the diff table nice
#####################################

diff <- read.table(diff_file, header = TRUE, stringsAsFactors = FALSE)

# Mitra: lots of NAs in df > omit them
diff <- na.omit(diff)
diff <- diff[(abs(diff$b) > 4),]

# Mitra: Now we have less infomation to dig into => runtime cost will be reduced dramtically
diff$gene_id <- diff$Accession
diff$Accession <- NULL

# this stuff is for the old version from Jason
#diff <- separate(diff, Accession, c("gene_id", "gene_version"), sep = "\\.")

#diff$mean_root_norm <- diff[, grepl("root_norm", colnames(diff))] %>% apply(1, mean)
#diff$sd_root_norm <- diff[, grepl("root_norm", colnames(diff))] %>% apply(1, sd)
#diff$mean_leaf_norm <- diff[, grepl("leaf_norm", colnames(diff))] %>% apply(1, mean)
#diff$sd_leaf_norm <- diff[, grepl("leaf_norm", colnames(diff))] %>% apply(1, sd)
#diff$root_fold_leaf <- log2(diff$mean_root_norm/diff$mean_leaf_norm)
#diff$bonf_pval <- p.adjust(diff$P.val, method = "bonferroni")

# delete raw counts, is diff, gene_version, individual norm counts
#diff <- diff[,!grepl("_raw", colnames(diff))]
#diff <- diff[,!grepl("^X.", colnames(diff))]
#diff$Is_DIFF <- NULL
#diff$gene_version <- NULL


#####################################
### merge to big-ass table
#####################################

all_features_wide <- extract(all_features_wide, remove = FALSE,
                             tss_name, c("gene_id", "chr", "loc", "offset?"), 
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
gc()
all_features_diffs_wide <- merge(diff, all_features_wide, by = "gene_id")  # about 4500 genes had tss peaks but no diff information; lots (18056) also had diffs info but no peaks. This is an inner join.
print("Done! We want to save them now!")
rm(all_features_wide)
rm(diff)
gc()

save(all_features_diffs_wide, file = output_file)




