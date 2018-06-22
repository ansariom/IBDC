leaf_peaks_file <- "~/Downloads/ibdc/peaks/aligned.peaks.annotated.capped_leaf.filtered"
root_peaks_file <- "~/Downloads/ibdc/peaks/aligned.peaks.annotated.capped_root.filtered"
diff_file <- "~/Downloads/ibdc/peaks/diff_exp_results.txt"
gff_file <- "~/Downloads/ibdc/peaks/TAIR10_GFF3_genes.gff"
min_reads <- 50
mod_count <- 10

leaf_peaks <- read.delim(leaf_peaks_file, sep = ",", header = T)
root_peaks <- read.delim(root_peaks_file, sep = ",", header = T)
diff_info <- read.delim(diff_file, sep = "\t", header = T)
gff <- read.delim(gff_file, sep = "\t", header = T)

gff <- gff[, c(1,3,4,5,7,9)]
colnames(gff) <- c("chr", "type", "start", "end", "strand", "ID")
gff <- gff[gff$type == "mRNA",]
library(tidyr)
gff <- extract(gff, ID, c("gene_id", "transcript_id"), regex = "([^;]+);([^;]+);\\D+")
gff <- extract(gff, gene_id, "TranscriptID", regex = "ID=([^;]+)")
gff <- extract(gff, transcript_id, "GeneName", regex = "Parent=([^;]+)")

leaf_peaks <- leaf_peaks[leaf_peaks$ReadCount > min_reads & leaf_peaks$ModeReadCount > mod_count,]
root_peaks <- root_peaks[root_peaks$ReadCount > min_reads & root_peaks$ModeReadCount > mod_count,]
a <- aggregate(leaf_peaks$GeneName, list(leaf_peaks$GeneName), length)
a <- leaf_peaks[leaf_peaks$GeneName == "AT"]

leaf_peaks$tissue <- "leaf"
root_peaks$tissue <- "root"
peaks <- rbind(leaf_peaks, root_peaks)

