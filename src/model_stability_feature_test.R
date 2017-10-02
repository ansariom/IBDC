

#coefs <- c("Y_Patch_FWD_3" ,"ATHB2_binding_site_motif_FWD_3", "ATHB2_binding_site_motif_REV_4", "CCA1_binding_site_motif_REV_5", "CCA1_motif1_BS_in_CAB1_REV_2", "Ibox_promoter_motif_REV_3")
load("~/Downloads/roe_collapsed_oc/featureinfo_features_diffsclasses.Rdata")
diffs_classes$tss_name <- row.names(diffs_classes)

base_dir <- "~/Downloads/crossVals_roe//"
seeds1 <- c(16023)
#seeds1 <- c(1298,16023,1899418)
#seeds2 <- c(12421, 5203)
seeds2 <- c(12421)

coefs <- c("ATHB2_binding_site_motif_FWD_3", "ATHB2_binding_site_motif_REV_4", "CCA1_binding_site_motif_REV_5", "CCA1_motif1_BS_in_CAB1_REV_2", "Ibox_promoter_motif_REV_3")
#result <- data.frame(seed1=c(), seed2=c(), nleaf1=c(), nleaf2=c(), nroot1=c(), nroot2=c())

rootness <- data.frame(feature=c(), group1_rootsum=c(),group1_leafsum=c(), group2_leaf_sum=c(), group2_rootsum=c(), diff_root=c(), diff_leaf=c())

feature_mean <- data.frame()

for (seed1 in seeds1) {
  for (seed2 in seeds2) {
    print(paste(seed1, "-", seed2, sep = ""))
    p1 <- paste(base_dir, "/", seed1, "/", "roe_only_x_train_set_scaled.rdat", sep = "")
    p2 <- paste(base_dir, "/", seed2, "/", "roe_only_x_train_set_scaled.rdat", sep = "")
    
    coef_file <- paste(base_dir, "/", seed1, "/", "roe_only_coef-table.sorted.txt", sep = "")
    coef_table <- read.table(coef_file, header = T)
    coef_table <- coef_table[100:1000,]
    
    selected <- features[,colnames(features) %in% coef_table$Feature]
    selected$tss_name <- row.names(features)
    sf <- merge(selected, diffs_classes, by = "tss_name")
    sf_root <- sf[sf$class==0,]
    sf_leaf <- sf[sf$class==1,]
    
    load(p1)
    a <- as.data.frame(x_train_set_scaled)
    
    load(p2)
    b <- as.data.frame(x_train_set_scaled)
    
#    a_c <- a[, colnames(a) %in% coefs]
#    b_c <- b[, colnames(b) %in% coefs]
    
    a_c <- a[, colnames(a) %in% coef_table$Feature]
    b_c <- b[, colnames(b) %in% coef_table$Feature]
    
    
   # print(apply(a_c, 2, mean))
    #print(apply(b_c, 2, mean))
    
    a_c$tss_name <- rownames(a_c)
    b_c$tss_name <- rownames(b_c)
    
    a_c <- extract(a_c, tss_name, into=c("fold", "tss_name"), regex ="(.)\\.(.+)")
    b_c <- extract(b_c, tss_name, into=c("fold", "tss_name"), regex ="(.)\\.(.+)")
    
    c1 <- merge(a_c, diffs_classes, by = "tss_name")
    c1$g <- "g1"
    c2 <- merge(b_c, diffs_classes, by = "tss_name")
    c2$g <- "g2"
    
    #------ Number of samples in each set root vs. leaf -------
    #a_leaf_no = nrow(c1[c1$class == 1,])
    #b_leaf_no = nrow(c2[c2$class == 1,])
    
    #a_root_no = nrow(c1[c1$class == 0,])
    #b_root_no = nrow(c2[c2$class == 0,])
    
    #d <- data.frame(seed1=seed1, seed2=seed2, nleaf1=a_leaf_no, nleaf2=b_leaf_no, nroot1=a_root_no, nroot2=b_root_no)
    #result <- rbind(d, result)
    #---------
    #--------- Plot coefficients ----------
    leaf_a <- c1[c1$class == 1,]
    root_a <- c1[c1$class == 0,]
    
    leaf_b <- c2[c2$class == 1,]
    root_b <- c2[c2$class == 0,]
    
    leaf_both <- rbind(leaf_a, leaf_b)
    root_both <- rbind(root_a, root_b)
    
    #all <- rbind(root_both, leaf_both)
    
    for(coef_name in coef_table$Feature) {
      fname <- paste("~/Downloads/coef_plots/", coef_name, ".png", sep = "")
      #p = ggplot(all, aes(all[,coef_name], fill = factor(class))) + geom_histogram(alpha = 0.5, position = 'identity') + facet_wrap(~g, ncol = 2)
      #ggsave(fname, p)
      
      #print(coef_name)
      leaf_a_hist <- hist(leaf_a[,coef_name])
      sa_leaf <- sum(leaf_a_hist$mids * leaf_a_hist$counts)
      root_a_hist <- hist(root_a[,coef_name])
      sa_root <- sum(root_a_hist$mids * root_a_hist$counts)
      
      leaf_b_hist <- hist(leaf_b[,coef_name])
      sb_leaf <- sum(leaf_b_hist$mids * leaf_b_hist$counts)
      root_b_hist <- hist(root_b[,coef_name])
      sb_root <- sum(root_b_hist$mids * root_b_hist$counts)
      
      
      d <- data.frame(feature=coef_name, sa_root, sa_leaf, sb_leaf, sb_root, diff_root = abs(sa_root - sb_root), diff_leaf=abs(sa_leaf - sb_leaf))
      #print(d)
      rootness <- rbind(rootness, d)
    }
    
 #   for(coef_name in coef_table$Feature) {
#      fname <- paste("~/Downloads/coef_plots/group1.", coef_name, ".png", sep = "")
#      p <- ggplot(leaf_both, aes(leaf_both[,coef_name], fill = factor(class))) + geom_histogram(alpha = 0.5, position = 'identity')
#      ggsave(fname, p)
#      fname <- paste("~/Downloads/coef_plots/group2.", coef_name, ".png", sep = "")
#      p <- ggplot(root_both, aes(root_both[,coef_name], fill =  factor(class))) + geom_histogram(alpha = 0.5, position = 'identity')
#      ggsave(fname, p)
      #fname <- paste("~/Downloads/coef_plots/both.", coef_name, ".png", sep = "")
      #p <- ggplot(all, aes(all[,coef_name], fill = g)) + geom_histogram(alpha = 0.5, position = 'identity')
      #ggsave(fname, p)
  #  }
  }
}

diff_features <- rootness[rootness$diff_root > 50 | rootness$diff_leaf > 50,]

seed1 <- 16023
seed2 <- 12421

p1 <- paste(base_dir, "/", seed1, "/", "roe_only_x_train_set_scaled.rdat", sep = "")
p2 <- paste(base_dir, "/", seed2, "/", "roe_only_x_train_set_scaled.rdat", sep = "")

selected <- features[,colnames(features) %in% diff_features$feature]
selected$tss_name <- row.names(features)
sf <- merge(selected, diffs_classes, by = "tss_name")
sf_root <- sf[sf$class==0,]
sf_leaf <- sf[sf$class==1,]

load(p1)
a <- as.data.frame(x_train_set_scaled)

load(p2)
b <- as.data.frame(x_train_set_scaled)
a_c <- a[, colnames(a) %in% coef_table$Feature]
b_c <- b[, colnames(b) %in% coef_table$Feature]
a_c$tss_name <- rownames(a_c)
b_c$tss_name <- rownames(b_c)
a_c <- extract(a_c, tss_name, into=c("fold", "tss_name"), regex ="(.)\\.(.+)")
b_c <- extract(b_c, tss_name, into=c("fold", "tss_name"), regex ="(.)\\.(.+)")
c1 <- merge(a_c, diffs_classes, by = "tss_name")
c1$g <- "g1"
c2 <- merge(b_c, diffs_classes, by = "tss_name")
c2$g <- "g2"

#--------- Plot coefficients ----------
leaf_a <- c1[c1$class == 1,]
root_a <- c1[c1$class == 0,]

leaf_b <- c2[c2$class == 1,]
root_b <- c2[c2$class == 0,]

leaf_both <- rbind(leaf_a, leaf_b)
root_both <- rbind(root_a, root_b)

all <- rbind(root_both, leaf_both)

for(coef_name in diff_features$feature) {
  fname <- paste("~/Downloads/coef_plots/", coef_name, ".png", sep = "")
  p = ggplot(all, aes(all[,coef_name], fill = factor(class))) + geom_histogram(alpha = 0.5, position = 'identity') + facet_wrap(~g, ncol = 2)
  ggsave(fname, p)
}


