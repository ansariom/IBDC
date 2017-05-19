
families <- read.table("~/Downloads/gene_families_sep_29_09_update.txt", sep = "\t", header = T, fill = T)
families <- families[,c("Gene_Family", "Genomic_Locus_Tag")]
families$Genomic_Locus_Tag <- toupper(families$Genomic_Locus_Tag)

load("~/Downloads/featureinfo_features_diffsclasses.Rdata")
classified_data = diffs_classes[(diffs_classes["class"] != -1000),]
geneInfo <- classified_data[,c("gene_id","class", "prob1", "prob0")]

geneInfo <- geneInfo[order(geneInfo$gene_id),]
geneInfo$gene_id <- unlist(lapply(geneInfo$gene_id, function(x){ gsub("\\.\\d+", "", x)  }))
geneInfo <- as.data.frame(unique(geneInfo))

genes <- merge(geneInfo, families, by.x = "gene_id", by.y = "Genomic_Locus_Tag")
genes$prob <- with(genes, ifelse(class == "1", prob1, prob0))
genes$class <- with(genes, ifelse(class == "1", "leaf", "root"))


genes <- genes[order(genes$Gene_Family),]
genes$ID <- paste(genes$gene_id, "_", genes$class, "_", genes$prob, sep = "")

write.table(genes, file = "~/Downloads/genes.txt", quote = F, row.names = F, sep = "\t")
