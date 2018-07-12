#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

feature_rdat <- args[1]
filter_pwms_file <- args[2]
outdir <-  args[3]

outfile <- paste(outdir, "/all_features_diffs_wide.rdat", sep = "")

load(feature_rdat)
filter_pwms <- read.table(filter_pwms_file, col.names=c("pwm"))

diffs_colnames <- c("gene_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs",
                    "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq",
                    "tss_name", "chr", "loc", "strand", "offset?")
  missing_cols <- all_features_diffs_wide[, colnames(all_features_diffs_wide) %in% diffs_colnames]

cols_to_exclude <- all_features_diffs_wide[, grep(paste(filter_pwms$pwm, collapse='|'), colnames(all_features_diffs_wide)),]

all_features_diffs_wide <- all_features_diffs_wide[, colnames(all_features_diffs_wide) %in% colnames(cols_to_exclude)]
all_features_diffs_wide <- cbind(all_features_diffs_wide, missing_cols)

save(all_features_diffs_wide, file = outfile)
