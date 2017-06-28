#!/usr/bin/env Rscript

library(LiblineaR)
library(ggplot2)
library(dplyr)
library(tidyr)
#library(e1071) # for svm with prob output
library(PRROC)
library(rstackdeque)
library(purrr)
library(parallel)


rm(list = ls())
#options(scipen=999)
geometric_mean <- function(x) exp(mean(log(x)))


# adds a "class" column of either "leaf_specific" or "root_specific" based
# on the threshold, assumes the present of column "b" (fold change) and "qval" (Q value)
# this also removes all unclassed examples
add_class <- function(input, qval_thresh, fold_thresh) {
  input$class <- NA
  input$class[input$qval <= qval_thresh & input$b < -1 * fold_thresh] <- 1#"leaf_specific"
  input$class[input$qval <= qval_thresh & input$b > fold_thresh] <- 0#"root_specific"
  input <- input[!is.na(input$class), ]
  return(input)
}



# split dataframe into train (and this into folds), final test
# works on any data frame
# rows will be randomized first
# returns a list; first el being a list of folds (data frames), second being the final test data frame
split_data <- function(input, percent_train = 0.8, folds = 5) {
  randomized_features_diffs_wide <- input[order(runif(nrow(input))), ]
  train_index <- as.integer(percent_train * nrow(randomized_features_diffs_wide))
  train_features_diffs_wide <- randomized_features_diffs_wide[seq(1, train_index), ]
  final_test_features_diffs_wide <- randomized_features_diffs_wide[seq(train_index + 1, nrow(randomized_features_diffs_wide)), ]
  
  folds_list <- split(train_features_diffs_wide, seq(1,nrow(train_features_diffs_wide)) %% folds)
  ret_list <- list(train_folds = folds_list, final_test = final_test_features_diffs_wide)
  return(ret_list)
}



# given a name for a fold, and a (named) list of all the folds, extracts
# that name as the test, the next name as the validation,
# and the rest as training (collapsed with rbind)
folds_to_train_validate_test <- function(test_fold_name, all_folds) {
  validation_fold_index <- (which(names(train_folds) == test_fold_name) + 1) %% length(all_folds) + 1
  validation_fold_name <- names(all_folds)[validation_fold_index]
  test_fold <- all_folds[[test_fold_name]]
  validation_fold <- all_folds[[validation_fold_name]]
  other_folds <- all_folds[!names(all_folds) %in% c(test_fold_name, validation_fold_name)]
  other_folds_bound <- do.call(rbind, other_folds)
  ret_list <- list(train_set = other_folds_bound, test_set = test_fold, validation_set = validation_fold)
  return(ret_list)
}

# given a name for a fold, and a (named) list of all the folds, extracts
# that name as the test, the rest as training (collapsed with rbind)
folds_to_train_test <- function(test_fold_name, all_folds) {
  test_fold <- all_folds[[test_fold_name]]
  other_folds <- all_folds[!names(all_folds) %in% c(test_fold_name)]
  other_folds_bound <- do.call(rbind, other_folds)
  ret_list <- list(train_set = other_folds_bound, test_set = test_fold)
  return(ret_list)
}

# runs the model on the given train_validate list (2 data frames), using regularization parameter param.
# note that this ASSUMES that train_validate_test_list is a list of 2 data frames,
# and each of those data frames has only numeric columns and a "class" column (factor) to be
# predicted
# returns: list of 6: param, confusion_matrix, coeffs_df, model, auroc, auprc
run_and_validate_model <- function(param, train_validate) {
  #print(paste("run_and_validate_model : model_type = ", model_type, "param = ", param, sep = ""))
  train_set <- train_validate[[1]]
  test_set <- train_validate[[2]]
  # hm, we gotta get rid of the cols that are all identical if there are any (0s sometimes)
  # but will this F it up? Might make comparisons after the fact tricky...
  # also we gotta not do this determination based on the "class" column
  train_keep_logical <- !unlist(lapply(train_set[,colnames(train_set) != "class"], function(col){sd(col) == 0}))
  train_keep_logical <- c(train_keep_logical, TRUE)
  #print(length(train_keep_logical))
  #print(dim(train_set))
  train_set <- train_set[, train_keep_logical]
  test_set <- test_set[, train_keep_logical]
  
  x_train_set <- train_set[, !colnames(train_set) %in% "class"]
  class_train_set <- train_set$class
  
  x_test_set <- test_set[, !colnames(test_set) %in% "class"]
  class_test_set <- test_set$class
  
  # scale the data
  x_train_set_scaled <- scale(x_train_set, center = TRUE, scale = TRUE)
  # also scale the test set by the same scale factor
  x_test_set_scaled <- scale(x_test_set, attr(x_train_set_scaled, "scaled:center"), attr(x_train_set_scaled, "scaled:scale"))
  
  model <- LiblineaR::LiblineaR(data = x_train_set_scaled, 
                                target = class_train_set,
                                type = 7,   # L2 regularized logistic regression (dual) = 7
                                cost = param,
                                bias = TRUE,  # ?? (recommended by vignette)
                                # cross = 10, # built-in cross validation; probably better to do it ourselves
                                verbose = TRUE)
  
  if (save_model == TRUE) {
    print("Saving Model ..")
    model_outfile <- paste(outdir, "/", model_type, "_model.rdat", sep = "")
    save(model, file = model_outfile)
    print(attr(x_train_set_scaled, "scaled:center"))
    scale_center_outfile <- paste(outdir, "/", model_type, "_x_train_set_scaled.rdat", sep = "")
    save(x_train_set_scaled, file = scale_center_outfile)
}
  coefficients <- model$W
  # drop bias coefficient
  coefficients <- coefficients[1:(length(coefficients) - 1)]
  
  p <- predict(model, x_test_set_scaled, proba = TRUE, decisionValues = TRUE)
  # produce a confusion matrix
  confusion_matrix <- table(predictions = p$predictions, actuals = class_test_set)
  
  probabilities <- as.data.frame(p$probabilities)
  rownames(p$probabilities) <- rownames(test_set)
  
  coeffs_df <- data.frame(coefficients, Feature = colnames(x_train_set_scaled), stringsAsFactors = FALSE)
  auroc <- PRROC::roc.curve(probabilities[,"1"], weights.class0 = class_test_set, curve = TRUE)$auc
  auprc <- PRROC::pr.curve(probabilities[,"1"], weights.class0 = class_test_set, curve = TRUE)$auc.davis.goadrich
  
  retlist <- list(param = param, confusion_matrix = confusion_matrix, coeffs_df = coeffs_df, model = p, auroc = auroc, auprc = auprc)
  
  return(retlist)
}


# given a fold name (length-1 character vector), uses it to extract train, validate, test sets,
# also for each param in the params list, tries that param.
# folds list is a list of data frames
# returns a list with lots of goodies
find_pstar <- function(fold_name, params_list, folds_list) {
  train_valid_test <- folds_to_train_validate_test(fold_name, folds_list)
  train_validate <- train_valid_test[c("train_set", "validation_set")]
  
  param_results <- lapply(params_list, run_and_validate_model, train_validate)
  
  best_result <- param_results[[1]]
  best_auroc <- param_results[[1]]$auroc
  for(result in param_results) {
    auroc <- result$auroc
    auprc <- result$auprc
    if(auroc > best_auroc) {
      best_result <- result
      best_auroc <- auroc
    }
  }
  
  all_params <- unlist(params_list)
  all_fold_names <- rep(fold_name, length(all_params))
  all_aurocs <- unlist(lapply(param_results, function(x) { return(x$auroc) } ))
  all_auprcs <- unlist(lapply(param_results, function(x) { return(x$auprc) } ))
  within_params_list <- list(fold_name = all_fold_names, param = all_params, auroc = all_aurocs, auprcs = all_auprcs)
  
  best_param <- best_result$param
  test_result <- run_and_validate_model(best_param, train_valid_test[c("train_set", "test_set")])
  ret_list <- list(fold_name = fold_name, 
                   train_set = train_valid_test$train_set,
                   test_set = train_valid_test$test_set,
                   best_model = best_result$model,
                   best_param = best_param, 
                   best_auroc = best_result$auroc, 
                   best_auprc = best_result$auprc, 
                   test_auroc = test_result$auroc, 
                   test_auprc = test_result$auprc, 
                   test_coeffs_df = test_result$coeffs_df,
                   mean_auroc = mean(all_aurocs), 
                   sd_auroc = sd(all_aurocs), 
                   within_params_list = within_params_list)
  return(ret_list)
  
}


# train_folds is a named list of data frames, possible_params is a named list of param values
# breaks out computation by fold 
# returns a list (of lists) for each fold with lots of goodies
n_fold_cross <- function(train_folds, possible_params) {
  train_folds_names <- as.list(names(train_folds))
  #print(train_folds_names)
  bests_by_fold <- lapply(train_folds_names, find_pstar, possible_params, train_folds)
  return(bests_by_fold)
}



####################################
##  End functions, begin script
####################################
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 5) {
  print("./cross_validation.R [all_features.rdat] [nfolds] [outdir] [ncpus] [tile_only|roe_only|tile_roe]")
  quit()
}

input_features <- args[1]
nfolds <- as.numeric(args[2])
outdir <- args[3]
ncpu <- as.numeric(args[4])
model_type = args[5]
save_model <- FALSE

# chosen by fair die roll, gauranteed to be random
set.seed(1298)
# load the features and differential expression data
# into all_features_diffs_wide
load(input_features)

if (model_type == "tile_only") {
  ### : remove all but tiled features
  all_features_diffs_wide <- all_features_diffs_wide[, !(grepl("(FWD|REV)", colnames(all_features_diffs_wide)) & !grepl("tile", colnames(all_features_diffs_wide))) ]
} else if (model_type == "roe_only") {
  ### or, remove tiled features
  all_features_diffs_wide <- all_features_diffs_wide[, !grepl("tile", colnames(all_features_diffs_wide)) ]
} else if (model_type == "oc_only") {
  missing_cols <- all_features_diffs_wide[,1:15]
  all_features_diffs_wide <- all_features_diffs_wide[, (grepl("_OC_", colnames(all_features_diffs_wide)) & !grepl("tile", colnames(all_features_diffs_wide))) ]
  all_features_diffs_wide <- cbind(missing_cols, all_features_diffs_wide)
} else if (model_type == "tfbs_only" ) {
  all_features_diffs_wide <- all_features_diffs_wide[, (!grepl("_OC_", colnames(all_features_diffs_wide)) & !grepl("tile", colnames(all_features_diffs_wide))) ]
} 

# set rownames to tss names
rownames(all_features_diffs_wide) <- all_features_diffs_wide$tss_name

# define classes, get rid of unclassed columns
classed_features_diffs_wide <- add_class(all_features_diffs_wide, qval_thresh = 0.05, fold_thres = 4)

# NAs were introduced because many TSSs have overall OC features but not others
## todo: why is this again?
classed_features_diffs_wide <- classed_features_diffs_wide[complete.cases(classed_features_diffs_wide), ]
print("Overall class sizes:")
print(table(classed_features_diffs_wide$class))

# strip out the differential expression stuff
diffs_colnames <- c("gene_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs", 
                    "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq", 
                    "tss_name", "chr", "loc", "offset?")

# differential expression data
classed_diffs_info <- classed_features_diffs_wide[, diffs_colnames]
# features and class only
classed_features_class <- classed_features_diffs_wide[, !colnames(classed_features_diffs_wide) %in% diffs_colnames]


# split into 80% 8-fold set, and 20% final test
folds_final_test <- split_data(classed_features_class, percent_train = 0.8, folds = nfolds)
train_folds <- folds_final_test$train_folds

# # make it all parallel...
library(parallel)
cl <- makeCluster(ncpu)
clusterExport(cl, list("folds_to_train_validate_test", "train_folds",
                       "run_and_validate_model", "model_type", "save_model"))
# replace lapply with parLapply
lapply <- function(...) {parLapply(cl, ...)}


# we'll try a bunch of different params
#possible_params <- as.list(10^seq(-6,-1,0.2))
possible_params <- as.list(seq(0.0001, 0.003, 0.0002))

print("trying params:")
print(unlist(possible_params))
names(possible_params) <- as.character(possible_params)

print("Start running cross validation ..")
bests_by_fold <- n_fold_cross(train_folds, possible_params)
#str(bests_by_fold[1])

print("cross-validation done!")
################################
####### Cross val done.
################################

######### Line plots start
print("Create Plots ....")

# make it into a table
# grab everything but the "within_params_list" entries and build a table
# map_df -> turns some parts of a list into a dataframe - from purrr library
bestfolds_table <- paste(outdir, "/", model_type ,"_bests_by_folds_param_vs_aurocs.txt", sep = "")
bests_by_fold_table <- map_df(bests_by_fold, .f = function(x) {return(x[!names(x) %in% c("within_params_list", "test_coeffs_df", "train_set", "test_set", "best_model")])} )
write.table(bests_by_fold_table, file = bestfolds_table, quote = F, sep = "\t", row.names = F)
print(bests_by_fold_table)

pstar_avg <- mean(bests_by_fold_table$best_param)
print(paste("mean_param = ", pstar_avg, sep = ""))
param_file <- paste(outdir, "/", model_type, "_param_mean.txt", sep = "")
write.table(pstar_avg, param_file, quote = F, row.names = F, col.names = F)

# grab the "within_params_list" entries and build a table
foldout_table <- paste(outdir, "/", model_type ,"_within_folds_param_vs_aurocs.txt", sep = "")
within_folds_table <- map_df(map(bests_by_fold, "within_params_list"), I)
print(as.data.frame(within_folds_table), row.names = FALSE)
write.table(within_folds_table, file = foldout_table, quote = F, sep = "\t", row.names = F)

foldout_plot <- paste(outdir, "/", model_type ,"_within_folds_param_vs_aurocs.png", sep = "")
df <- within_folds_table
df$fold_name <- as.character(df$fold_name)
plot_title <- paste(model_type, " Features", sep = "")
g <- ggplot(df) + geom_line(aes(x = param, y = auroc, color = fold_name)) +
  expand_limits(y = c(0.80, 1.0)) +
  ggtitle(plot_title) 
ggsave(g, filename = foldout_plot)

######### Line plots end

######## folds_final_test is a list of 2; first is a list of data frames (folds), second is the held out df
print("Start running heldout test...")
all_train <- do.call(rbind, folds_final_test$train_folds)
print(paste("train set size : ", dim(all_train), sep = ""))
final_test <- folds_final_test$final_test
print(paste("test set size: ", dim(final_test), sep = ""))
save_model <- TRUE
final_res <- run_and_validate_model(pstar_avg, list(all_train, final_test))

perf_out <- paste(outdir, "/", model_type, "_heldoutTest_performance.txt", sep = "")
perf_df <- data.frame(auroc = final_res$auroc, auprc = final_res$auprc)
write.table(perf_df, file = perf_out, sep = "\t", row.names = F, quote = F)

model_file <- paste(outdir, "/", model_type, "_model.rdat", sep = "")
save(final_res, file = model_file)
#print("model saved at: ", model_file, sep = "")

### Test calls
test_calls_features <- cbind(final_res$model$probabilities, final_res$model$predictions)
colnames(test_calls_features)[1] <- "prob_class_0"
colnames(test_calls_features)[2] <- "prob_class_1"
colnames(test_calls_features)[3] <- "class_call"
all_input_test_rows <- classed_features_diffs_wide[rownames(classed_features_diffs_wide) %in% rownames(final_test), ]

test_outfile <- paste(outdir, "/", model_type, "_testout.table.txt", sep = "")
write.table(all_input_test_rows, file = test_outfile, quote = F, sep = "\t")
print("test is done!")

###### Coefficients table #######
print("Writing coef tables..")
coef_df <- final_res$coeffs_df
coef_df <- coef_df[order(-abs(coef_df$coefficients)),]
coef_table <- paste(outdir, "/", model_type, "_coef-table.sorted.txt", sep = "")
write.table(coef_df, file = coef_table, quote = F, row.names = F, sep = "\t")





