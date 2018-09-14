library(tidyr)
library(plotly)

read_files <- function(indir, top_count, regx) {
  all_models <- data.frame(feature=c(), coef=c())
  all_coefs_all <- data.frame(feature=c(), coef=c())
  a <- c()
  fs = list.files(indir)
  count = length(fs)
  for (f in fs) {
    print(f)
    all_coefs = read.table(paste(indir, "/",  f, sep = ""), col.names = c("feature", "coef"))
    top_tfs <- all_coefs[1:top_count,]
    alist <- extract(top_tfs, feature, c("pwm", "strand", "win", "type"), regex = regx, remove = FALSE)
    a <- c(a,unique(alist$pwm))
    all_models <- rbind(all_models, top_tfs)
    all_coefs_all <- rbind(all_coefs_all, all_coefs)
  }
  p = data.frame(pwm=a)
  pwm_repeats = aggregate(p, by = list(p$pwm), FUN = length)
  colnames(pwm_repeats) <- c("pwm", "nobserved")
  sel_tfs <- pwm_repeats[pwm_repeats$nobserved > (count/2), "pwm"]
  
  all_models_tmp <- extract(all_models, feature, c("pwm", "strand", "win", "type"), regex = regx, remove = FALSE)
  all_models_tmp <- all_models_tmp[all_models_tmp$pwm %in% sel_tfs,]
  all_models <- all_models[all_models$feature %in% all_models_tmp$feature,]
  
  all_coefs_tmp <- extract(all_coefs_all, feature, c("pwm", "strand", "win", "type"), regex = regx, remove = FALSE)
  all_coefs_tmp <- all_coefs_tmp[all_coefs_tmp$pwm %in% sel_tfs,]
  all_coefs_all <- all_coefs_all[all_coefs_all$feature %in% all_coefs_tmp$feature,]
  
  return(list("all_models" = all_models, "all_coefs" = all_coefs_all))
}

avg_coef_computation <- function(all_models, all_coefs, regx) {
  all_models$fullpwm <- paste(all_models$pwm, all_models$type, sep = "")
  tfs_top50 <- unique(all_models$fullpwm)
  all_coefs$fullpwm <- paste(all_coefs$pwm, all_coefs$type, sep = "")
  coefs_top50_all <- all_coefs[all_coefs$fullpwm %in% tfs_top50,]
  selected_tfs <- coefs_top50_all[, colnames(coefs_top50_all) %in% c("feature", "coef")]
  avg_coefs_top50 <- aggregate(selected_tfs, by = list(selected_tfs$feature), FUN = mean)
  avg_coefs_top50$feature <- NULL
  colnames(avg_coefs_top50) <- c("feature", "coef_avg")
  avg_coefs_top50 <- extract(avg_coefs_top50, feature, c("pwm", "strand", "win", "type"), regex = regx, remove = FALSE)
  return(avg_coefs_top50)
}

prepare_plot_data <- function(avg_coefs_top50) {
  #avg_coefs_top50 <- avg_coefs_top50_roe
  avg_coefs_top50$pwm <- paste(avg_coefs_top50$pwm, avg_coefs_top50$type, sep = "")
  avg_coefs_top50_fwd <- avg_coefs_top50[avg_coefs_top50$strand == "FWD", colnames(avg_coefs_top50) %in% c("pwm", "coef_avg", "left")]
  avg_coefs_top50_rev <- avg_coefs_top50[avg_coefs_top50$strand == "REV", colnames(avg_coefs_top50) %in% c("pwm", "coef_avg", "left")]
  
  d_fwd <- spread(avg_coefs_top50_fwd, key = left, value = coef_avg, fill = 0)
  d_rev <- spread(avg_coefs_top50_rev, key = left, value = coef_avg, fill = 0)
  
  mat_fwd <- d_fwd[,2:length(d_fwd)]
  mat_fwd <- as.matrix(mat_fwd, rownames.force = F)
  y_fwd <- d_fwd[,1]
  x_fwd <- colnames(d_fwd)[2:length(d_fwd)]
  mat_rev <- d_rev[,2:length(d_rev)]
  mat_rev <- as.matrix(mat_rev, rownames.force = F)
  y_rev <- d_rev[,1]
  x_rev <- colnames(d_rev)[2:length(d_rev)]
  return(list("x_fwd" = x_fwd, "y_fwd" = y_fwd, "mat_fwd" = mat_fwd, "x_rev" = x_rev, "y_rev" = y_rev, "mat_rev" = mat_rev))
}

  f <- list(
    family = "Courier New, monospace",
    size = 1,
    color = "#7f7f7f"
  )
  a <- list(
    title = "PWM",
    showticklabels = TRUE,
    tickangle = 0,
    exponentformat = "E",
    pad = 1
  )
  xaxis <- list(
    title = "Promoter region",
    range = c(-1000,2000)
  )
  m <- list(
    l = 200,
    r = 50,
    b = 50,
    t = 50,
    pad = 1
  )

tile_dir <- "~/Downloads/ibdc/Aug2018/tile_vs_roe/tiles/"
roe_dir = "~/Downloads/ibdc/Aug2018/tile_vs_roe/roes/"
top_count = 100
tile_regx <- "(.+?)_(FWD|REV)_(.+)_tile100(.*)"
roe_regx <- "(.+?)_(FWD|REV)_(.)(.*)"
roe_win_coords = read.table("~/Downloads/ibdc/Aug2018/tile_vs_roe/roe_win_coords.txt", col.names = c("feature", "left", "right")) # coordinates of each window in ROE table

# ----- Tiled Models ------
l = read_files(tile_dir, top_count, tile_regx)
all_models_tile <-  extract(l$all_models, feature, c("pwm", "strand", "win", "type"), regex = tile_regx, remove = FALSE)
all_coefs_tile <-  extract(l$all_coefs, feature, c("pwm", "strand", "win", "type"), regex = tile_regx, remove = FALSE)
avg_coefs_top50_tile = avg_coef_computation(all_models_tile, all_coefs_tile, tile_regx)
avg_coefs_top50_tile$left <- -1000 + (as.numeric(avg_coefs_top50_tile$win) - 1) * 100
avg_coefs_top50_tile$right <- avg_coefs_top50_tile$left + 100

p_title = paste("Average model weights for top", top_count, " TFs in 10 Tiled models (FWD)", sep = "")
p_tile <- prepare_plot_data(avg_coefs_top50_tile)
p <- p_tile
plot_ly( z = p$mat_fwd, type = "heatmap", y = p$y_fwd, x = as.numeric(p$x_fwd), colors = c( "#641E16","white", "white", "#1D8348")) %>%
  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = p_title)
#plot_ly( z = p$mat_fwd, type = "heatmap", y = p$y_fwd, x = as.numeric(p$x_fwd), colors = c( "#641E16","white")) %>%
#  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = p_title)
p_title = paste("Average model weights for top", top_count, " TFs in 10 Tiled models (REV)", sep = "")
plot_ly( z = p$mat_rev, type = "heatmap", y = p$y_rev, x = as.numeric(p$x_rev), colors = c( "#641E16","white", "white", "#1D8348")) %>%
  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = p_title)
#plot_ly( z = p$mat_fwd, type = "heatmap", y = p$y_fwd, x = as.numeric(p$x_fwd), colors = c( "#641E16","white")) %>%
#  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = p_title)

# ----- ROE Models ----
l = read_files(roe_dir, top_count, roe_regx)
all_models_roe <-  extract(l$all_models, feature, c("pwm", "strand", "win", "type"), regex = roe_regx, remove = FALSE)
all_coefs_roe <-  extract(l$all_coefs, feature, c("pwm", "strand", "win", "type"), regex = roe_regx, remove = FALSE)
avg_coefs_top50_roe = avg_coef_computation(all_models_roe, all_coefs_roe, roe_regx)
avg_coefs_top50_roe$fullPwm <- paste(avg_coefs_top50_roe$pwm, avg_coefs_top50_roe$strand, avg_coefs_top50_roe$win, sep = "_")
avg_coefs_top50_roe = merge(avg_coefs_top50_roe, roe_win_coords, by.x = "fullPwm", by.y = "feature")
avg_coefs_top50_roe$fullPwm <- NULL

p_roe <- prepare_plot_data(avg_coefs_top50_roe)
p <- p_roe
p_title = paste("Average model weights for top", top_count, " TFs in 10 ROE models (FWD)", sep = "")
plot_ly( z = p$mat_fwd, type = "heatmap", y = p$y_fwd, x = as.numeric(p$x_fwd), colors = c( "#641E16","white", "white","white", "white" ,"#1D8348")) %>%
  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = p_title)

p_title = paste("Average model weights for top", top_count, " TFs in 10 ROE models (REV)", sep = "")
plot_ly( z = p$mat_rev, type = "heatmap", y = p$y_rev, x = as.numeric(p$x_rev), colors = c( "#641E16","white", "white", "#1D8348")) %>%
  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = p_title)

#------ Copmare Tiled and ROE
avg_coefs_top50_roe$coef_avg <- avg_coefs_top50_roe$coef_avg/max(abs(avg_coefs_top50_roe$coef_avg))
write.table(avg_coefs_top50_tile[order(-abs(avg_coefs_top50_tile$coef_avg)),], file = "~/Downloads/tile_top_coefs_sorted.txt", sep = "\t", quote = F, row.names = F)
avg_coefs_top50_roe$model_type = "ROE"
avg_coefs_top50_tile$model_type = "Tiled"
avg_coefs_top50_tile$coef_avg <- avg_coefs_top50_tile$coef_avg / max(abs(avg_coefs_top50_tile$coef_avg))
both <- rbind(avg_coefs_top50_roe, avg_coefs_top50_tile)
both$pwm <- paste(both$pwm, both$type, sep = "")
# ggplot option for zoomed in view
data_fwd = both[both$strand == "FWD",]
data_rev = both[both$strand == "REV",]

## common vs different
length(unique(both$pwm))
common_80run <- merge(avg_coefs_top50_roe, avg_coefs_top50_tile, by = c("pwm", "strand"))
length(unique(common_80run$pwm))
thresh <- 100
common_80run$in_roe <- ifelse(common_80run$right.x - common_80run$left.y < (-1 * thresh), 0, ifelse(common_80run$left.x - common_80run$right.y > thresh, 0, 1))
not_in_roe <- common_80run[common_80run$in_roe == 0,]

## plots
data <- data_fwd
ggplot(data, aes(y = pwm)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = pwm, xend = data$right, yend = pwm, color = model_type, alpha = abs(coef_avg)), size = 4)  +
#  geom_point(aes(x = data$left, color = model_type, alpha = abs(coef_avg)), size = 5, shape = 'I') +
#  geom_point(aes(x = data$right, color = model_type, alpha = abs(coef_avg)), size = 5, shape = 'I') + 
  scale_size(guide = "none") +
  scale_alpha_continuous(guide = FALSE) +
  #scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
  xlim(-1000,1500) +
  scale_color_manual(values = c("red", "blue")) + 
  theme_bw() +
  theme(axis.text = element_text(size=7, color="black",face="bold"), axis.text.x = element_text(angle = 45)) + 
  theme(legend.text=element_text(size=9,face="bold")) +
  theme(strip.text = element_text(size=9,face="bold")) +
  ggtitle("Comparison between Top 50 Highly Weighted Features in ROE vs. Tile Models (FWD)")+  facet_wrap(~model_type) 

data <- data_rev
ggplot(data, aes(y = pwm)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = pwm, xend = data$right, yend = pwm, color = model_type, alpha = abs(coef_avg)), size = 4)  +
#  geom_point(aes(x = data$left, color = model_type, alpha = abs(coef_avg)), size = 5, shape = 'I') +
#  geom_point(aes(x = data$right, color = model_type, alpha = abs(coef_avg)), size = 5, shape = 'I') + 
  scale_size(guide = "none") +
  scale_alpha_continuous(guide = FALSE) +
  #scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
  xlim(-1000,1500) +
  scale_color_manual(values = c("red", "blue")) + 
  theme_bw() +
  theme(axis.text = element_text(size=7, color="black",face="bold"), axis.text.x = element_text(angle = 45)) + 
  theme(legend.text=element_text(size=9,face="bold")) +
  theme(strip.text = element_text(size=9,face="bold")) +
  ggtitle("Comparison between Top 50 Highly Weighted Features in ROE vs. Tile Models (REV)")+
  facet_wrap(~model_type) 
# ----- numerical stats -----

