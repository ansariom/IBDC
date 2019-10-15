indir = "~/Downloads/ibdc/aug2019/"
leaf_open_promoter <- paste(indir, "leaf_promoter_open_regions_3000-3000.txt", sep = "")
root_open_promoter <- paste(indir, "root_promoter_open_regions_3000-3000.txt", sep = "")
all_tss <- paste(indir, "/aligned.peaks.annotated.capped.filtered", sep = "")

#######
df_tss <- read.delim(all_tss, sep = ",", header = T)
head(df_tss)

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
ggplot(counts, aes(num_transcripts, tss_count, col = tissue, shape = tissue)) + geom_point() + theme_bw() + 
  ggtitle("Number of mapped TSS peaks  and transcript counts")

mapped_loc <- aggregate(df_tss$tss_id, list(df_tss$tissue, df_tss$TranscriptLocation), length)
colnames(mapped_loc) <- c("Tissue", "Location", "No_TSSs")
ggplot(mapped_loc, aes(x = Location, y = No_TSSs, fill = Tissue)) +
  geom_col() + ylab("Number of TSS peaks") + xlab("Transcript Location") +
  theme_bw() + ggtitle("Promoter location for mapped TSS peaks relative to Transcripts")

length(unique(df_tss$TranscriptID))

aggregate(df_tss$TranscriptID, list(df_tss$tissue, df_tss$TranscriptID), length)

root_tss <- df_tss[df_tss$tissue == "root",]
leaf_tss <- df_tss[df_tss$tissue == "leaf",]

length(unique(root_tss$TranscriptID))
length(unique(leaf_tss$TranscriptID))

both <- merge(root_tss, leaf_tss, by = "TranscriptID")
both$left <- ifelse(both$Start.x > both$Start.y, both$Start.x, both$Start.y) 
both$right <- ifelse(both$End.x < both$End.y, both$End.x, both$End.y)
both$overlap <- ifelse(both$Start.x > both$End.y,0, ifelse(both$End.x < both$Start.y, 0, abs(both$left - both$right)))
both$mod_dist <- abs(both$ModeLocation.x - both$ModeLocation.y)
both$tss_dist_groups <- cut(both$mod_dist, breaks = c(seq(0,4, by=5), seq(5,30, by=26), seq(31,100, by=70), seq(101,max(both$mod_dist), by=(max(both$mod_dist)-101))), 
                            include.lowest = T)
ggplot(both, aes(factor(tss_dist_groups))) + geom_bar(stat = "count", fill="steelblue") + xlab("TSS mode distance (base pair)") + 
  geom_text(aes(label=count), vjust=1.6, color="white", size=3.5) + theme_minimal()

ggplot(both, aes(factor(tss_dist_groups), fill = tss_dist_groups)) + geom_bar(stat = "count", col = "black") + xlab("TSS mode distance (base pair)") + 
  theme_minimal() + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#A69F00")) + geom_text(stat='count', aes(label=..count..), vjust=-1)

both$tss <- "tss"
ggplot(both, aes(tss, overlap)) + geom_violin() + ylab( "Size of overlap between TSS peaks (basepair)") +
  xlab("") + ggtitle("Size of overlap between root and shoot TSS peaks in each promoter")
mean(both$overlap)
median(both$overlap)
hist(both$overlap)
both$range <- cut(both$overlap, c(seq(-1,50, by = 10)))
table(both$range)
#######


all_peaks <- read.delim(all_tss, sep = ",", header = T)
all_peaks$pstart <- all_peaks$Start - 3000
all_peaks$tss_name <- paste(all_peaks$TranscriptID, all_peaks$Chromosome, all_peaks$Start - 3000, all_peaks$Strand, sep = "_")

col_names = c("tss_id", "oc_id", "chr", "start", "end", "rel_start", "rel_end")
leaf_oc <- read.table(leaf_open_promoter, col.names = col_names)
root_oc <- read.table(root_open_promoter, col.names = col_names)

m = merge(all_peaks, root_oc, by.x = "tss_name", by.y = "tss_id")
close_root <- all_peaks[all_peaks$tss_name %in% root_oc$tss_id,]
length(unique(root_oc$tss_id))
length(unique(leaf_oc$tss_id))

#root_oc <- leaf_oc

#all
left = -1000
right = 1000
root_oc$allL <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$allR <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)


# -100 to +100
left = -200
right = 200

root_oc$tssL <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$tssR <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)

# 100 to 500
left = 200
right = 600
root_oc$down500L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$down500R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)
  
# 1000 to 3000
left = 600
right = 1000
root_oc$down1000L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$down1000R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)

left = 1000
right = 3000
root_oc$down3000L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$down3000R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)

# -500 to -1000
left = -200
right = -600
root_oc$up500L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$up500R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)
  
# 1000 to 3000
left = -600
right = -1000
root_oc$up1000L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$up1000R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)

left = -1000
right = -3000
root_oc$up3000L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$up3000R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)

###########
df <- data.frame(location = c(), bp=c())
#all
left = -1000
right = 1000
root_oc$all <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$allL - root_oc$allR) ))
r = root_oc[root_oc$all > 0, "all"]
d = data.frame(location = "All(-1000_+1000)", bp = r)
df <- rbind(df, d)

left = -1000
right = -3000
root_oc$up3000 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$up3000L - root_oc$up3000R) ))
r = root_oc[root_oc$up3000 > 0, "up3000"]
d = data.frame(location = "-1000_-3000", bp = r)
#df <- rbind(df, d)

left = -600
right = -1000
root_oc$up1000 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$up1000L - root_oc$up1000R) ))
r = root_oc[root_oc$up1000 > 0, "up1000"]
d = data.frame(location = "-600_-1000", bp = r)
df <- rbind(df, d)

left = -200
right = -600
root_oc$up500 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$up500L - root_oc$up500R) ))
r = root_oc[root_oc$up500 > 0, "up500"]
d = data.frame(location = "-200_-600", bp = r)
df <- rbind(df, d)

left = -200
right = 200
root_oc$TSS <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$tssR - root_oc$tssL) ))
r = root_oc[root_oc$TSS > 0, "TSS"]
d = data.frame(location = "TSS", bp = r)
df <- rbind(df, d)

left = 200
right = 600
root_oc$down500 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$down500L - root_oc$down500R) ))
r = root_oc[root_oc$down500 > 0, "down500"]
d = data.frame(location = "200-600", bp = r)
df <- rbind(df, d)

left = 600
right = 1000
root_oc$down1000 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$down1000L - root_oc$down1000R) ))
r = root_oc[root_oc$down1000 > 0, "down1000"]
d = data.frame(location = "600-1000", bp = r)
df <- rbind(df, d)

left = 1000
right = 3000
root_oc$down3000 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$down3000L - root_oc$down3000R) ))
r = root_oc[root_oc$down3000 > 0, "down3000"]
d = data.frame(location = "1000-3000", bp = r)
#df <- rbind(df, d)

ggplot(df, aes(x = location, y = bp, fill = location)) +
  geom_violin() + ylab("Size of Open region (basepair)") + xlab("Relative Location to TSS") +
  theme_bw() + ggtitle("Coverage Regions of Open Chromatin (Shoot)")

table(df$location)

v = root_oc[root_oc$down3000 >0,]

