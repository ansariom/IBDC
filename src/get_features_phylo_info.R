#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 3) {
	print("./get_features_phylo_info.R [tair10_gene_family_file] [features_info.RDat] [outfile] ")
	quit()
}

gene_family_file <- args[1]
features_info <- args[2]
outfile <- args[3]

families <- read.table(gene_family_file, sep = "\t", header = T, fill = T)
families <- families[,c("Gene_Family", "Genomic_Locus_Tag")]
families$Genomic_Locus_Tag <- toupper(families$Genomic_Locus_Tag)

load(features_info)
classified_data = diffs_classes[(diffs_classes["class"] != -1000),]
geneInfo <- classified_data[,c("gene_id","class", "prob1", "prob0")]

geneInfo <- geneInfo[order(geneInfo$gene_id),]
geneInfo$tair_id <- unlist(lapply(geneInfo$gene_id, function(x){ gsub("\\.\\d+", "", x)  }))
geneInfo <- as.data.frame(unique(geneInfo))

genes <- merge(geneInfo, families, by.x = "tair_id", by.y = "Genomic_Locus_Tag")
genes$prob <- with(genes, ifelse(class == "1", prob1, prob0))
genes$class <- with(genes, ifelse(class == "1", "leaf", "root"))


genes <- genes[order(genes$Gene_Family),]
genes$ID <- paste(genes$gene_id, "_", genes$class, "_", genes$prob, sep = "")

write.table(genes, file = outfile, quote = F, row.names = F, sep = "\t")
