library(DESeq)
base = "~/Downloads/"
appx = ".isoforms.results"
m = ""
read_counts <- function(i){
  root = paste(base, "root", i, appx, sep = "")
  leaf = paste(base, "leaf", i, appx, sep = "")
  r = read.delim(root, sep = "\t")[,c("transcript_id", "expected_count")]
  rownames(r) <- r$transcript_id
  r$transcript_id <- NULL
  colnames(r) <- paste("root", i, sep = "")
  l = read.delim(leaf, sep = "\t")[,c("transcript_id", "expected_count")]
  rownames(l) <- l$transcript_id
  l$transcript_id <- NULL
  colnames(l) <- paste("leaf", i, sep = "")
  cbind(r,l)
}
counts <- do.call(cbind, lapply(seq(1,3), read_counts))[,c( "root1", "root2", "root3",
                                                            "leaf1", "leaf2", "leaf3" )]
counts = (counts * 100) %/% 100

meta_data = data.frame(
  row.names = colnames( counts ),
  condition = c( "root", "root", "root",
                 "leaf", "leaf", "leaf" ),
  libType = c( "single-end", "single-end", "single-end",
               "single-end", "single-end", "single-end") )

head(meta_data)
condition = factor(meta_data$condition)
cds = newCountDataSet(counts, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )
cds = estimateDispersions( cds )
plotDispEsts( cds )
res = nbinomTest( cds, "leaf", "root" )
colnames(res) <- c("Accession", "baseMean", "leaf_tpm", "root_tpm", "foldChange", "b", "pval", "qval")
head(res)
write.table(res, file = "ath_root_leaf_rsem_deseq_diff_expr_results.txt", quote = F, row.names = F, sep = "\t")

plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
resSig = res[ res$padj < 0.1, ]
head(resSig)

k_diff = read.delim("~/Downloads/ibdc/July2018/kallisto_diff_exp_results.txt", sep ="\t")
a <- merge(k_diff, res, by = "Accession")[,c("Accession", "foldChange", "b.x",  "b.y", "leaf_tpm", "root_tpm", "pval.x", "pval.y", "qval.x", "qval.y")]

diff <- a[abs(a$b.x) > 2 & abs(a$b.y) > 2 & a$qval.y < 0.1 & a$qval.x < 0.05,]
plot(a$b, a$log2FoldChange, xlab = "Kallisto_logFoldChange", ylab = "rsem_logFoldChange")
a = na.omit(a)
cor(a$b, a$log2FoldChange)
