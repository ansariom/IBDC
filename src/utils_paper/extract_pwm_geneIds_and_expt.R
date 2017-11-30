huges_genes <- read.delim("~/Downloads/huges_genes_mats.txt", fill = T, col.names = c("pwm", "gene_id"), sep = "\t")
huges_genes <- extract(huges_genes, pwm, c("matrix_id", "suffix"), regex = "(M\\d+)_(.*)", remove = F)
transcfac_genes <- read.delim("~/Downloads/ath.mats_genes.txt", fill = T, sep = "\t", col.names = c("pwm", "gene_id"))
coef_table <- read.table("~/Downloads/roe_only_coef_table.txt", col.names = c("feature", "coef"))
diff <- read.table("~/Downloads/diff_exp_results.txt", header = T)


huges_genes <- na.omit(huges_genes)

coef_table <- extract(coef_table, feature, c("pwm", "strand", "window"), regex = "(.+?)_(FWD|REV)_(\\d+)", remove = F)
coef_table <- extract(coef_table, pwm, c("matrix_id", "suffix"), regex = "(M\\d+)_(.*)", remove = F)
c1 <- merge(coef_table, transcfac_genes, by.x = "matrix_id", by.y = "pwm", all.x = T)

c1 <- merge(c1, huges_genes, by = "matrix_id", all.x = T)

c1$pwm.y <- NULL
c1$strand <- NULL
c1$suffix.x <- NULL
c1$suffix.y <- NULL
c1$matrix_id <- NULL

c1$gene_ids <- paste(c1$gene_id.x, c1$gene_id.y, sep = "")
c1$gene_id.x <- NULL
c1$gene_id.y <- NULL
c1 <- c1[order(-abs(c1$coef)),]

