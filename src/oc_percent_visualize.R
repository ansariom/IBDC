#!/usr/bin/Rscript

load("ibdc_roe-only_fixedROE/features_long_oc_p_root_wide.rdat")
oc_root <- all_features_wide

load("ibdc_roe-only_fixedROE/features_long_oc_p_leaf_wide.rdat")
oc_leaf <- all_features_wide

load("ibdc_roe-only_fixedROE/all_features_diffs_wide.rdat")

r <- all_features_diffs_wide[all_features_diffs_wide$b > 4,]
l <- all_features_diffs_wide[all_features_diffs_wide$b < -4,]

oc_root_root <- oc_root[oc_root$tss_name %in% r$tss_name,]
oc_root_leaf <- oc_root[oc_root$tss_name %in% l$tss_name,]
oc_leaf_leaf <- oc_leaf[oc_leaf$tss_name %in% l$tss_name,]
oc_leaf_root <- oc_leaf[oc_leaf$tss_name %in% r$tss_name,]


s1 <- oc_root_root[, grepl("SORLREP5", colnames(oc_root_root))]
s2 <- oc_leaf_leaf[, grepl("SORLREP5", colnames(oc_leaf_leaf))]

apply(s1, MARGIN= 2, FUN = mean)

apply(s2, MARGIN= 2, FUN = mean)
