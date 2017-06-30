library(tidyr)
library(ggplot2)
#library(NeatMap)

# given two dataframes with similar rownames, merge (join) them using those rownames
merge_by_rownames <- function(df1, df2, ...) {
  df1$asdfasdf <- rownames(df1)
  df2$asdfasdf <- rownames(df2)
  res <- merge(df1, df2, by = "asdfasdf", ...)
  rownames(res) <- res$asdfasdf
  res$asdfasdf <- NULL
  return(res)
}

normalize_vector_01 <- function(x) {
  (x-min(x))/(max(x)-min(x))
}

# load the data
if(TRUE) {
  setwd("~/Downloads/")
  load("featureinfo_features_diffsclasses.Rdata")
}

# make a dataframe of just calls & true class information
prob1 <- diffs_classes[, "prob1", drop = F]
true_class <- diffs_classes[diffs_classes$class %in% c(0,1), "class", drop = F]
colnames(true_class) <- "true_class"
info <- merge_by_rownames(prob1, true_class)

# normalize the data feature cols. Use scaling relative to the trained model
load("roe_only_x_train_set_scaled.rdat")
features_scaled <- as.data.frame(scale(features, attr(x_train_set_scaled, "scaled:center"), attr(x_train_set_scaled, "scaled:scale")))

# multiply them by the coefficients
#features_scaled_by_coeffs <- as.data.frame(scale(features_scaled, scale = feature_info$coefficient))
features_scaled_by_coeffs <- features_scaled * feature_info$coefficient

# normalize data such that each element is between 0 and 1
features_list_scaled_by_coeffs_norm <- lapply(features_scaled_by_coeffs, normalize_vector_01)
features_scaled_by_coeffs_norm <- as.data.frame(features_list_scaled_by_coeffs_norm)
row.names(features_scaled_by_coeffs_norm) <- row.names(features_scaled_by_coeffs)


#features_scaled_by_coeffs <- scale(features_scaled_by_coeffs, scale = F)

# merge that data along with the class info above
features_product_info <- merge_by_rownames(info, features_scaled_by_coeffs_norm)
# scale the true class and class calls into a larger range for plotting
features_product_info$true_class <- (features_product_info$true_class - 0.5) * 600
features_product_info$prob1 <- (features_product_info$prob1 - 0.5) * 600

## for molly: look for pairs of differentially expressed TSSs that are similar in profile
# the strategy here is to loop over all rows that are root-specific, and compare them to all rows that are 
# leaf-specific. That's a lot of comparisons, and we don't want to store all of them; rather we just store
# top N. This is done by keeping a dataframe of results with N rows, and if the current comparison is larger (in similarity, abs(correlation coefficient))
# than any existing row in the table, replace that row with the current comparison.
root_specific <- features_product_info[which(features_product_info$true_class < 0), !colnames(features_product_info) %in% c("prob1", "true_class")]
leaf_specific <- features_product_info[which(features_product_info$true_class > 0), !colnames(features_product_info) %in% c("prob1", "true_class")]

#r <- as.numeric(root_specific[1,])
#l <- as.numeric(leaf_specific[1,])
#diff <- abs(r-l)
#qplot(diff, data=data.frame(diff), geom="histogram", alpha=I(.3), col=I("red"),fill=I("blue")) + 
#  ggtitle(paste("Difference between ", 1,1,sep = ","))
#hist(diff)
#cutoff <- 0.2
#d <- diff[which(diff > cutoff)]
#rn <- r[which(diff > cutoff)]
#ln <- l[which(diff > cutoff)]

outdir <- "~/Downloads/plots"

compute_diff <- function(list_of_pairs, root_specific, leaf_specific) {
  #print(list_of_pairs)
  cutoffs <- c(0.2, 0.3, 0.4)
  i = as.numeric(list_of_pairs$i)
  j = as.numeric(list_of_pairs$j)
  #print(i)
  #print(j)
  root_spec_values <- as.numeric(root_specific[i,])
  leaf_spec_values <- as.numeric(leaf_specific[j,])
  diff <- abs(root_spec_values - leaf_spec_values)
  
  if (sd(diff) == 0) { 
    #print(paste(i,j, sep = ","))
    next
  }
  
  results <- list()
  #results <- data.frame(root_index=c(), leaf_index=c(), cutoff=c(), num_diffs=c(), cor=c())
  for (cutoff in cutoffs) {
    d <- diff[which(diff > cutoff)]
    
    if (length(d) < 10 ) {
      q <- qplot(diff, data=data.frame(diff), geom="histogram", alpha=I(.3), col=I("red"),fill=I("blue")) + 
        ggtitle(paste("Difference between ", i, j,sep = ","))
      outfile <- paste(outdir, "/", i, "_", j, ".png", sep = "");
      ggsave(filename = outfile, q)
    }
    root_items <- root_spec_values[which(diff > cutoff)]
    leaf_items <- leaf_spec_values[which(diff > cutoff)]
  
    pres <- list(root_index=i, leaf_index=j, cutoff=cutoff, num_diffs=length(d), cor=as.numeric(cor.test(root_items, leaf_items)$estimate))
    #pres <- data.frame(roo_index=i, leaf_index=j, cutoff=cutoff, num_diffs=length(d), cor=as.numeric(cor.test(root_items, leaf_items)$estimate))
    #print(pres)
    results <- c(results, pres)
    #print(pres)
    #results <- rbind(results, pres)
    
  }
  return(results)
}


library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

root_index <- nrow(root_specific)
leaf_index <- nrow(leaf_specific)
list_of_pairs = list()
for(i in 1:3) {
  for (j in 1:2) {
    l <- list(i=i, j=j)
    list_of_pairs <- c(list_of_pairs,list(l))
  }
}
result <- parLapply(cl, list_of_pairs, compute_diff, root_specific, leaf_specific)
all <- as.data.frame(matrix(unlist(result), byrow = T, ncol = 5))
colnames(all) <- c("root_index", "leaf_index", "cutoff", "num_diffs", "cor_coef")
write.table(all, file = table_outfile, row.names = F, quote = F, sep = "\t")
#print(result)
#all_results <- rbind(all_results, unlist(result))
#all_results <- c(all_results, result)

stopCluster(cl)
