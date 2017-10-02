#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 11) {
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
all_TSS_diff_outfile <- args[11]

#####################################
### Make the diff table nice
#####################################

diff <- read.table(diff_file, header = TRUE, stringsAsFactors = FALSE)

# A copy of complete features
diff_all <- diff

# Mitra: lots of NAs in df > omit them
diff <- na.omit(diff)
diff <- diff[(abs(diff$b) > 4),]

# Mitra: Now we have less infomation to dig into => runtime cost will be reduced dramtically
diff$gene_id <- diff$Accession
diff$Accession <- NULL


#############
## Read and merge features
#############
library(data.table)
drops <- c("gene_id", "chr", "loc", "offset?")

features <- as.data.frame(fread(feat_wide_file, header = TRUE, stringsAsFactors = FALSE))
features$tss_name <- features$V1
features <- extract(features, remove = FALSE,
                             tss_name, drops,
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
features <- subset(features, features$gene_id %in% diff$gene_id)
features$V1 <- NULL
rownames(features) <- features$tss_name
features <- features[, !(names(features) %in% drops)]
features_long <- gather(features, feature, value, -tss_name)
print(paste("features.rdat   ", nrow(features), sep = ""))
#-------------------
tiled_features <- as.data.frame(fread(tiled_feat_wide_file, header = TRUE, stringsAsFactors = FALSE))
tiled_features$tss_name <- tiled_features$V1
tiled_features <- extract(tiled_features, remove = FALSE,
                             tss_name, drops,
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
tiled_features <- subset(tiled_features, tiled_features$gene_id %in% diff$gene_id)
tiled_features$V1 <- NULL
rownames(tiled_features) <- tiled_features$tss_name
tiled_features <- tiled_features[, !(names(tiled_features) %in% drops)]
tiled_features_long <- gather(tiled_features, feature, value, -tss_name)
print(paste("tiled_features.rdat   ", nrow(tiled_features), sep = ""))
#-------------------
oc_features_leaf <- read.table(oc_leaf_long_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)
oc_features_leaf <- extract(oc_features_leaf, remove = FALSE,
                             tss_name, drops,
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
oc_features_leaf <- subset(oc_features_leaf, oc_features_leaf$gene_id %in% diff$gene_id)
oc_features_leaf <- oc_features_leaf[, !(names(oc_features_leaf) %in% drops)]
print(paste("c_features_leaf   ", nrow(oc_features_leaf), sep = ""))
#-------------------
oc_features_root <- read.table(oc_root_long_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)
oc_features_root <- extract(oc_features_root, remove = FALSE,
                             tss_name, drops,
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
oc_features_root <- subset(oc_features_root, oc_features_root$gene_id %in% diff$gene_id)
oc_features_root <- oc_features_root[, !(names(oc_features_root) %in% drops)]
print(paste("c_features_root   ", nrow(oc_features_root), sep = ""))

#-------------------

oc_features_overall_leaf <- read.table(oc_leaf_overall_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)
oc_features_overall_leaf <- extract(oc_features_overall_leaf, remove = FALSE,
                             tss_name, drops,
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
oc_features_overall_leaf <- subset(oc_features_overall_leaf, oc_features_overall_leaf$gene_id %in% diff$gene_id)
oc_features_overall_leaf <- oc_features_overall_leaf[, !(names(oc_features_overall_leaf) %in% drops)]
print(paste("oc_features_overall_leaf   ", nrow(oc_features_overall_leaf), sep = ""))
#-------------------

oc_features_overall_root <- read.table(oc_root_overall_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)
oc_features_overall_root <- extract(oc_features_overall_root, remove = FALSE,
                             tss_name, drops,
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
oc_features_overall_root <- subset(oc_features_overall_root, oc_features_overall_root$gene_id %in% diff$gene_id)
oc_features_overall_root <- oc_features_overall_root[, !(names(oc_features_overall_root) %in% drops)]
print(paste("c_features_overall_root   ", nrow(oc_features_overall_root), sep = ""))

#-------------------
tiled_oc_features_leaf <- as.data.frame(fread(tiled_oc_leaf_long_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE))

# Break large file into multiple smaller files and merge them later 
all <- data.frame(tss_name=c(), feature=c(), value=c())
s = 10000000
nparts <- floor(nrow(tiled_oc_features_leaf)/s)
print(nparts)
for (i in 1:nparts) {
	start = ((i-1) * s ) + 1
	end = i * s
	print(paste(start, "-", end, sep = ""))
	parti <- tiled_oc_features_leaf[start:end,]
	parti <- extract(parti, remove = FALSE,
                             tss_name, drops,
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
	parti <- subset(parti, parti$gene_id %in% diff$gene_id)
	parti <- parti[, !(names(parti) %in% drops)]
	all <- rbind(all, parti)
	rm(parti)
}
#remain <- nrow(tiled_oc_features_leaf) - (nparts * s)
lastPart <- tiled_oc_features_leaf[((nparts * s) + 1):nrow(tiled_oc_features_leaf),]
lastPart <- extract(lastPart, remove = FALSE,
                             tss_name, drops,
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
lastPart <- subset(lastPart, lastPart$gene_id %in% diff$gene_id)
lastPart <- lastPart[, !(names(lastPart) %in% drops)]
all <- rbind(all, lastPart)

#tiled_oc_features_leaf <- extract(tiled_oc_features_leaf, remove = FALSE,
#                             tss_name, drops,
#                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
#tiled_oc_features_leaf <- subset(tiled_oc_features_leaf, tiled_oc_features_leaf$gene_id %in% diff$gene_id)
#tiled_oc_features_leaf <- tiled_oc_features_leaf[, !(names(tiled_oc_features_leaf) %in% drops)]
tiled_oc_features_leaf <- all
rm(all)
print(paste("tiled_oc_features_leaf   ", nrow(tiled_oc_features_leaf), sep = ""))

#-------------------
tiled_oc_features_root <- as.data.frame(fread(tiled_oc_root_long_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE))
all <- data.frame(tss_name=c(), feature=c(), value=c())
s = 10000000
nparts <- floor(nrow(tiled_oc_features_root)/s)
print(nparts)
for (i in 1:nparts) {
        start = ((i-1) * s ) + 1
        end = i * s
        print(paste(start, "-", end, sep = ""))
        parti <- tiled_oc_features_root[start:end,]
        parti <- extract(parti, remove = FALSE,
                             tss_name, drops,
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
        parti <- subset(parti, parti$gene_id %in% diff$gene_id)
        parti <- parti[, !(names(parti) %in% drops)]
        all <- rbind(all, parti)
        rm(parti)
}
#remain <- nrow(tiled_oc_features_leaf) - (nparts * s)
lastPart <- tiled_oc_features_root[((nparts * s) + 1):nrow(tiled_oc_features_root),]
lastPart <- extract(lastPart, remove = FALSE,
                             tss_name, drops,
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
lastPart <- subset(lastPart, lastPart$gene_id %in% diff$gene_id)
lastPart <- lastPart[, !(names(lastPart) %in% drops)]
all <- rbind(all, lastPart)

#tiled_oc_features_leaf <- extract(tiled_oc_features_leaf, remove = FALSE,
#                             tss_name, drops,
#                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
#tiled_oc_features_leaf <- subset(tiled_oc_features_leaf, tiled_oc_features_leaf$gene_id %in% diff$gene_id)
#tiled_oc_features_leaf <- tiled_oc_features_leaf[, !(names(tiled_oc_features_leaf) %in% drops)]
tiled_oc_features_root <- all
rm(all)
print(paste("tiled_oc_features_root   ", nrow(tiled_oc_features_root), sep = ""))

#--------------------
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
rm(tiled_oc_features_root)
rm(tiled_oc_features_leaf)
rm(all_features_long)
gc()

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

diff_all$gene_id <- diff_all$Accession
diff_all$Accession <- NULL
all_tss_diffs_wide <- merge(diff_all, all_features_wide, by = "gene_id")
save(all_tss_diffs_wide, file = all_TSS_diff_outfile)

diff_all <- NULL
all_tss_diffs_wide <- NULL
gc()












