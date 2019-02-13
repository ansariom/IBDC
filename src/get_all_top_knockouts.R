
library(tidyr)


get_all_data <- function(d) {
  print(basename(d))
  f = list.files(d, pattern = "n_top*")
  a <- read.table(paste(d, f, sep = "/"))
  a <- a[,c(1,2,5)]
  colnames(a) <- c("expr_change", "feature", "promoter")
  a$type = basename(d)
  a <- extract(a, type, "type", regex = "\\d_\\d+_([^_]+)_roe")
  return(a)
}

dirs <- dir(".", pattern="3_")
all <- do.call(rbind, lapply(dirs, get_all_data))

get_groups <- function(df) {
  r = data.frame(feature=df$feature[1], promoter=df$promoter[1], nobserved=nrow(df), models=paste(df$type, collapse = ","))
  r
}

g <- do.call(rbind, by(all, all[,c("feature", "promoter")], get_groups))
g <- g[order(-g$nobserved),]
write.table(g, file = "all_top_candidates_all_models_ranked_noBoundary.txt", sep = "\t", quote=F, row.names=F, col.names=T)

a <- aggregate(g$models, list(g$models), length)
colnames(a) <- c("models", "count")
a <- a[order(a$count),]



#####
df <- read.csv("seqContentTest/ATContent_classed.csv",header = T)
r = df[df$class == 0,]
l = df[df$class == 1,]
l$class <- NULL
r$class <- NULL
r = data.frame(apply(r, 2, mean))
l = data.frame(apply(l, 2, mean))
r$feature <- rownames(r)
l$feature <- rownames(l)
colnames(r) <- c("mean_ATContent", "feature")
colnames(l) <- c("mean_ATContent", "feature")
r <- extract(r, feature, into = c("seq", "tile_no"), regex = "(.+Content)(\\d+)")
l <- extract(l, feature, into = c("seq", "tile_no"), regex = "(.+Content)(\\d+)")
r$type <- "root"
l$type <- "leaf"
df <- rbind(r,l)
#g <- ggplot(df, aes(factor(tile_no,levels = tile_no), mean_ATContent, color = type)) + geom_line()
g <- ggplot(df, aes(as.numeric(tile_no), mean_ATContent, color = type)) + geom_line() + scale_x_discrete(name ="tile no", breaks=seq(1,28,1)) + geom_point()
ggsave(g, file = "mean_AT_content.png")

df <- read.csv("seqContentTest/GCContent_classed.csv",header = T)
r = df[df$class == 0,]
l = df[df$class == 1,]
l$class <- NULL
r$class <- NULL
r = data.frame(apply(r, 2, mean))
l = data.frame(apply(l, 2, mean))
r$feature <- rownames(r)
l$feature <- rownames(l)
colnames(r) <- c("mean_GCContent", "feature")
colnames(l) <- c("mean_GCContent", "feature")
r <- extract(r, feature, into = c("seq", "tile_no"), regex = "(.+Content)(\\d+)")
l <- extract(l, feature, into = c("seq", "tile_no"), regex = "(.+Content)(\\d+)")
r$type <- "root"
l$type <- "leaf"
df <- rbind(r,l)
#g <- ggplot(df, aes(factor(tile_no,levels = tile_no), mean_ATContent, color = type)) + geom_line()
g <- ggplot(df, aes(as.numeric(tile_no), mean_GCContent, color = type)) + geom_line() + scale_x_discrete(name ="tile no", breaks=seq(1,28,1)) + geom_point()
ggsave(g, file = "mean_GC_content.png")



