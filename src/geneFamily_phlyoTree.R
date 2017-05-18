
add_class <- function(input, qval_thresh, fold_thresh, strip_unclassed = TRUE) {
  input$class <- "none"
  input$class[input$qval <= qval_thresh & input$b < -1 * fold_thresh] <- "leaf" #"leaf_specific"
  input$class[input$qval <= qval_thresh & input$b > fold_thresh] <- "root" #"root_specific"
  if(strip_unclassed) {
    input <- input[!is.na(input$class), ]
  }
  return(input)
}

load("~/Downloads/all_features_diffs_wide.rdat")
families <- read.table("~/Downloads/gene_families_sep_29_09_update.txt", sep = "\t", header = T, fill = T)
families$Genomic_Locus_Tag <- toupper(families$Genomic_Locus_Tag)

all <- add_class(all_features_diffs_wide, 0.05, 4)
geneInfo <- all[,c("gene_id", "class")]
geneInfo <- geneInfo[order(geneInfo$gene_id),]
geneInfo$gene_id <- unlist(lapply(geneInfo$gene_id, function(x){ gsub("\\.\\d+", "", x)  }))
geneInfo <- as.data.frame(unique(geneInfo))

genes <- merge(geneInfo, families, by.x = "gene_id", by.y = "Genomic_Locus_Tag")

