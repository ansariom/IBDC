library(ggplot2)

desktop_mode = TRUE

args = commandArgs(trailingOnly = T)
leaf_open_promoter = args[1]
root_open_promoter = args[2]
all_tss <- args[3]
diff_expr_file = args[4]

if (desktop_mode) {
  indir = "./"
  leaf_open_promoter <- paste(indir, "leaf_promoter_open_regions_3000-3000.txt", sep = "")
  root_open_promoter <- paste(indir, "root_promoter_open_regions_3000-3000.txt", sep = "")
  all_tss <- paste(indir, "aligned.peaks.annotated.capped.filtered", sep = "")
  diff_expr_file = paste(indir,"ath_root_leaf_rsem_deseq_diff_expr_results.txt", sep = "")
}

col_names = c("tss_id", "oc_id", "chr", "start", "end", "rel_start", "rel_end")
leaf_oc <- read.table(leaf_open_promoter, col.names = col_names)
root_oc <- read.table(root_open_promoter, col.names = col_names)

de_trx <- read.table(diff_expr_file, header = T)
up_down_trx <- de_trx[abs(de_trx$foldChange) > 3 & (de_trx$mean_leaf_norm > 300 | de_trx$mean_root_norm > 300),]
df_tss <- read.delim(all_tss, sep = ",", header = T)
head(df_tss)
#######

df_tss$tss_id <- paste(df_tss$Chromosome, df_tss$TranscriptID, df_tss$Strand, df_tss$ModeLocation, df_tss$tissue, sep = "_")
aggregate(df_tss$tss_id, list(df_tss$tissue), length)
tss_per_trx <- aggregate(df_tss$tss_id, list(df_tss$tissue, df_tss$TranscriptID), length)
head(tss_per_trx)
a = as.data.frame(table(tss_per_trx[tss_per_trx$Group.1 == "root",]$x))
a$tissue <- "root"
b = as.data.frame(table(tss_per_trx[tss_per_trx$Group.1 == "leaf",]$x))
b$tissue = "shoot"
counts = rbind(a,b)
colnames(counts) <- c("num_transcripts", "tss_count", "tissue")
g <- ggplot(counts, aes(num_transcripts, tss_count, col = tissue, shape = tissue)) + geom_point() + theme_bw() + 
  ggtitle("Number of mapped TSS peaks  and transcript counts")
g
outfile = paste(indir, "tss_ntranscripts.png", sep = "")
ggsave(g, file = outfile)

mapped_loc <- aggregate(df_tss$tss_id, list(df_tss$tissue, df_tss$TranscriptLocation), length)
colnames(mapped_loc) <- c("Tissue", "Location", "No_TSSs")
g <- ggplot(mapped_loc, aes(x = Location, y = No_TSSs, fill = Tissue)) +
  geom_col() + ylab("Number of TSS peaks") + xlab("Transcript Location") +
  theme_bw() + ggtitle("Promoter location for mapped TSS peaks relative to Transcripts")

outfile = paste(indir, "tss_mapped_locs_promoter.png", sep = "")
ggsave(g, file = outfile)

length(unique(df_tss$TranscriptID))

aggregate(df_tss$TranscriptID, list(df_tss$tissue, df_tss$TranscriptID), length)

root_tss <- df_tss[df_tss$tissue == "root",]
leaf_tss <- df_tss[df_tss$tissue == "leaf",]

length(unique(root_tss$TranscriptID))
length(unique(leaf_tss$TranscriptID))

both <- merge(root_tss, leaf_tss, by = "TranscriptID", suffixes = c(".root", ".leaf"))
both$left <- ifelse(both$Start.root > both$Start.leaf, both$Start.root, both$Start.leaf) 
both$right <- ifelse(both$End.root < both$End.leaf, both$End.root, both$End.leaf)
both$overlap <- ifelse(both$Start.root > both$End.leaf,0, ifelse(both$End.root < both$Start.leaf, 0, abs(both$left - both$right)))
both$mod_dist <- abs(both$ModeLocation.root - both$ModeLocation.leaf)
both$tss_dist_groups <- cut(both$mod_dist, breaks = c(seq(0,4, by=5), seq(5,100, by=96), seq(101,max(both$mod_dist), by=(max(both$mod_dist)-101))), 
                            include.lowest = T)
#ggplot(both, aes(factor(tss_dist_groups))) + geom_bar(stat = "count", fill="steelblue") + xlab("TSS mode distance (base pair)") + 
#  geom_text(aes(label=count), vjust=1.6, color="white", size=3.5) + theme_minimal()

g <- ggplot(both, aes(factor(tss_dist_groups), fill = tss_dist_groups)) + geom_bar(stat = "count", col = "black") + xlab("TSS mode distance (base pair)") + 
  theme_minimal() + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#A69F00")) + geom_text(stat='count', aes(label=..count..), vjust=-1) +
  scale_x_discrete(labels = c("[0,5]" = "0-5","(5,101]"  ="5-100", "(101,3.42e+03]" = ">100"))

outfile = paste(indir, "tss_mod_distance_barplot.png", sep = "")
ggsave(g, file = outfile)

de_tss <- both[both$TranscriptID %in% up_down_trx$Accession,]

outfile = paste(indir, "DE_transcripts_TSS_mod_distance.txt", sep = "")
write.table(de_tss, file = outfile, quote = F, row.names = F, col.names = T, sep = "\t")

g <- ggplot(de_tss, aes(factor(tss_dist_groups), fill = tss_dist_groups)) + geom_bar(stat = "count", col = "black") + xlab("TSS mode distance (base pair)") + 
  theme_minimal() + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#A69F00")) + geom_text(stat='count', aes(label=..count..), vjust=-1) +
  scale_x_discrete(labels = c("[0,5]" = "0-5","(5,101]"  ="5-100", "(101,3.42e+03]" = ">100"))
outfile = paste(indir, "diffexpr_tss_mod_distance_barplot.png")
ggsave(g, file = outfile)

g <- ggplot(de_tss, aes(mod_dist)) + geom_histogram(binwidth = 2) + xlab("TSS mode distance (base pair)") + 
  theme_minimal() 

outfile = paste(indir, "diffexpr_tss_mod_distance_histogram.png")
ggsave(g, file = outfile)


both$tss <- "tss"
ggplot(both, aes(tss, overlap)) + geom_violin() + ylab( "Size of overlap between TSS peaks (basepair)") +
  xlab("") + ggtitle("Size of overlap between root and shoot TSS peaks in each promoter")
mean(both$overlap)
median(both$overlap)
hist(both$overlap)
both$range <- cut(both$overlap, c(seq(-1,50, by = 10)))
table(both$range)

#------------ OC stats on not-DE transcripts not having diff TSS ---------
nde_tss <- de_tss[de_tss$mod_dist < 30,]
length(unique(nde_tss$tss_id.root))

nde_tss$tss_id.root <- paste(nde_tss$TranscriptID, nde_tss$Chromosome.root, nde_tss$ModeLocation.root, nde_tss$Strand.root, sep = "_")
nde_tss$tss_id.leaf <- paste(nde_tss$TranscriptID, nde_tss$Chromosome.leaf, nde_tss$ModeLocation.leaf, nde_tss$Strand.leaf, sep = "_")
nde_tss <- nde_tss[,c("tss_id.root", "tss_id.leaf", "mod_dist")]
#colnames(nde_tss) <- c("tss_id.root", "tss_id.leaf", "mod_dist")

orig_nde_tss <- nde_tss

nde_tss <- merge(nde_tss, root_oc[,c("tss_id", "start", "end")], by.x = "tss_id.root", by.y = "tss_id")
length(unique(nde_tss$tss_id.root))

nde_tss <- merge(nde_tss, leaf_oc[,c("tss_id", "start", "end")], by.x = "tss_id.leaf", by.y = "tss_id", suffixes = c(".root", ".leaf"))
length(unique(nde_tss$tss_id.root))

not_open_in_root <- orig_nde_tss[!orig_nde_tss$tss_id.root %in% root_oc$tss_id & orig_nde_tss$tss_id.leaf %in% leaf_oc$tss_id,]
not_open_in_leaf <- orig_nde_tss[orig_nde_tss$tss_id.root %in% root_oc$tss_id & !orig_nde_tss$tss_id.leaf %in% leaf_oc$tss_id,]
nrow(orig_nde_tss[!orig_nde_tss$tss_id.leaf %in% nde_tss$tss_id.leaf,])

nde_tss$left <- ifelse(nde_tss$start.root > nde_tss$start.leaf, nde_tss$start.root, nde_tss$start.leaf) 
nde_tss$right <- ifelse(nde_tss$end.root < nde_tss$end.leaf, nde_tss$end.root, nde_tss$end.leaf)
nde_tss$l <- ifelse(nde_tss$start.root > nde_tss$start.leaf, nde_tss$start.leaf, nde_tss$start.root) 
nde_tss$r <- ifelse(nde_tss$end.root < nde_tss$end.leaf, nde_tss$end.leaf, nde_tss$end.root)

nde_tss$overlap <- ifelse(nde_tss$start.root > nde_tss$end.leaf,0, ifelse(nde_tss$end.root < nde_tss$start.leaf, 0, 
                                                                          abs(nde_tss$start.root - nde_tss$end.leaf)))
nde_tss$percent_overlap <- nde_tss$overlap / (nde_tss$r - nde_tss$l)
all <- aggregate(nde_tss$percent_overlap, list(nde_tss$tss_id.root, nde_tss$tss_id.leaf, nde_tss$start.root, nde_tss$end.root), max)
head(all)

hist(all$x, xlab = "% OC overlap between root and shoot promoters", main = "Histogram of open regions comparision between DE transcripts \n having differenet TSS locations in root and shoot")

all_avg <- aggregate(all$x, list(all$Group.1, all$Group.2), mean)
hist(all_avg$x, xlab = "% OC overlap between root and shoot promoters", main = "Histogram of open regions comparision between DE transcripts \n having similar TSS locations in root and shoot")


a = all[all$x < 0.1,]
nrow(all[all$x < 0.1,])
length(unique(a$Group.2))

length(unique(nde_tss$TranscriptID))
length(unique(de_tss$TranscriptID))
nde_tss <- nde_tss[, c("tss_id.leaf", "tss_id.root", "mod_dist", "overlap")]
