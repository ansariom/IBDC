## Peak finding
#=============================================
argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

## closest index finding
#=============================================
find_closest_to_max <- function(max_idx, mod_idxes) {
  distance_to_max <- 10000 
  closest <- -1
  for (idx in mod_idxes) { 
    distance = abs(idx - max_idx) 
    if (distance < distance_to_max) { 
      distance_to_max <- distance 
      closest <- idx 
    } 
  } 
  return(closest)
}

## closes bin finding
#=============================================
find_closest_bin <- function(search_idx, bin_indxes, win = 50) {
  distance_to_max <- 10000 
  closest <- ""
  for (idx in bin_indxes) { 
    distance = abs(search_idx - idx) 
    if (distance < distance_to_max && distance < win) { 
      distance_to_max <- distance 
      closest <- idx 
    } 
  } 
  return(closest)
}

# Write outputs to the table
#=============================================
printout_NA_table <- function(outfname) {
  outTable <- matrix(nrow=1, ncol=4)
  outTable[1,1] <- NA
  outTable[1,2] <- NA
  outTable[1,3] <- NA
  outTable[1,4] <- NA
  dimnames(outTable)[[2]] <- c("MaxPeakLoc", "HalfWidth", "Left", "Right")
  pwmName <- unlist(strsplit(f, "[.]"))[1]
  dimnames(outTable)[[1]] <- pwmName
  write.table(outTable, file=outfname, row.names=T, col.names=T, quote=F, sep="	")
}


# Plot smoothed peaks
#=============================================
plot_modes <- function() {
  smooth_plotfile <- paste(outdir, f, ".smoothed.jpg", sep = "")
  jpeg(file=smooth_plotfile) 
  plot(locs, scores, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep="")) 
  lines(locs, peaks$y.hat,  lwd=2) 
  scores.min <- min(scores) 
  sapply(peaks$i, function(i) lines(c(locs[i],locs[i]), c(scores.min, peaks$y.hat[i]),col="Red", lty=2)) 
  points(locs[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, cex=1.25) 
  dev.off()
}

# Plot ROE peaks
#=============================================
plot_roe <- function() {
  left <- locs[leftind]
  right <- locs[rightind]
  peak_mode_loc <- locs[peak_mod_index]
  halfWidth <- (right - left)/2.0;
  if (halfWidth > 300) {
    halfWidth <- 250 
    left <- peak_mode_loc - halfWidth 
    right <- peak_mode_loc + halfWidth 
  } 
  
  outTable <- matrix(nrow=1, ncol=4)
  outTable[1,1] <- format(peak_mode_loc, digits=2) 
  outTable[1,2] <- format(halfWidth, digits=2) 
  outTable[1,3] <- format(left, digits=2) 
  outTable[1,4] <- format(right, digits=2) 
  
  dimnames(outTable)[[2]] <- c("MaxPeakLoc", "HalfWidth", "Left", "Right")
  pwmName <- unlist(strsplit(f, "[.]"))[1]
  dimnames(outTable)[[1]] <- pwmName
  write.table(outTable, file=outfname, row.names=T, col.names=T, quote=F, sep="	")
  # Plot ROEs
  subT <- paste("Interval = [", toString(format(left, digits=2)), ",", toString(format(right, digits=2)), "]")
  mainT <- paste(pwmName, subT, sep="\n")
  plot_filename <- paste(outdir, f, ".jpg", sep="")
  jpeg(file=plot_filename)
  par(mfrow=c(2,1))
  plot(locs, scores, main=mainT, xlab="locs", ylab="scores")
  points(locs[locs >= left & locs <= right], scores[locs >= left & locs <= right], col='red', pch=1)
  plot(locs[(locs >= left - 20) & (locs <= right + 20)], scores[(locs >= left - 20) & (locs <= right + 20)], main=mainT, xlab="locs", ylab="scores")
  points(locs[locs >= left & locs <= right], scores[locs >= left & locs <= right], col='red', pch=1)
  points(peak_mode_loc, maxbgval, pch=20, col='green')
  arrows(left, maxbgval, right, maxbgval, col='green', length=0.1, code=3)
  results <- dev.off()
}

# Initial screeing for existance of ROE
#=============================================
is_roe_exist <- function(mod_idxes) {
  all_bins <- data.frame(idx = c(), count = c(), sum = c())
  for (peak_candidate_idx in mod_idxes ) {
    item <- data.frame(idx = peak_candidate_idx, count = 1, sum = 0)
    all_bins <- rbind(all_bins, item)
  }
  for (j in 1:length(scores)) {
    closest_peak_idx <- find_closest_bin(j, all_bins$idx)
    if (closest_peak_idx > 0) {
      bin_idx <- which(all_bins$idx == closest_peak_idx)
      all_bins[bin_idx,2] <- all_bins[bin_idx, "count"] + 1
      all_bins[bin_idx,3] <- scores[j] + all_bins[bin_idx, "sum"]
    }
  }
  # normalize values
  sum_sd <- sd(all_bins$sum / max(all_bins$sum))
  if (sum_sd > 0.04)
    return(TRUE)
  return(FALSE)
}

# get the peak that is more enriched by high scores
#=============================================
get_max_peak <- function(above_avg_mod_idxes, extreme_max_idxs) {
  # initialize the bins for each peak mod
  bins <- data.frame(idx = c(), count = c(), sum = c())
  for (peak_candidate_idx in above_avg_mod_idxes ) {
    item <- data.frame(idx = peak_candidate_idx, count = 1, sum = 0)
    bins <- rbind(bins, item)
  }
  
  # add scores to each bin that are within 50 bp of mod
  for (j in extreme_max_idxs) {
    closest_peak_idx <- find_closest_bin(j, bins$idx)
    if (closest_peak_idx > 0) {
      bin_idx <- which(bins$idx == closest_peak_idx)
      bins[bin_idx,2] <- bins[bin_idx, "count"] + 1
      bins[bin_idx,3] <- scores[j] + bins[bin_idx, "sum"]
    }
  }
  
  # get bins that have more max in them
  max_bins <- bins[bins$count == max(bins$count),]
  
  if (nrow(max_bins) > 1) {
    max_idxs <- max_bins[max_bins$sum == max(max_bins$sum), "idx"]
    top_idx <-  which(scores == max(scores))
    closest_to_max_idx <- find_closest_to_max(top_idx, max_idxs)
    max_idx <- closest_to_max_idx
    # if even the sum is the same then choose first one
  } else {
    max_idx <- max_bins$idx
  }
  
  return(max_idx)
}

datadir <- "~/Downloads/test/dists_all///"
outdir <- "~/Downloads/test/plots/"

flist <- list.files(datadir)

# Init default parameters 
w <- 100 # win length for signal mode finder
span <- 0.01 # span size [0,1] (small values = more sensitive)
buffer_left <- 50
buffer_right <- 50
peak_mode_loc <- -6000
max_neighbor_dist <- 500

#------------ Start Main Loop Over all dists --------------

for (f in flist) {
  #f <- flist[1]
  print(f)
  dist_tbl <- read.table(paste(datadir, f, sep = ""), header = F)
  locs <- dist_tbl[,1] 
  scores <- dist_tbl[,2] 
  maxdist <- length(locs)/2 
  outfname <- paste(outdir, f, ".table", sep = "")

  
  # find the modes and plot them
  peaks <- argmax(locs, scores, w=w, span=span)
  plot_modes()
  
  ## process found peaks
  ##########################
  mod_idxes <- peaks$i 
  
  if (length(mod_idxes) < 1 ) {
    print("no peaks are found!")
    printout_NA_table()
    next
  }
  
  smooth_scores <- peaks$y.hat
  
  # Don't consider peaks beyond 1.5kb frim TSS
  mod_idxes = mod_idxes[which(locs[mod_idxes] > -1500 & locs[mod_idxes] < 1500)]
  
  # Pre filter peaks that are almost equal to each other
  #if (!is_roe_exist(mod_idxes)) {
  #  print("No ROE can be found for this distribution of modes")
  #  printout_NA_table(outfname)
  #  next
  #}
  
  # get top two mode locations
  if (length(mod_idxes > 2)) {
    last_idx = 2
  } else {
    last_idx = length(mod_idxes)
  }
  cutoff <- sort(smooth_scores[mod_idxes], decreasing = T)[last_idx]
  above_avg_mod_idxes <- mod_idxes[which((smooth_scores[mod_idxes] >= cutoff))]
  if (length(above_avg_mod_idxes) < 1) {
    print("No good peak were found!")
    printout_NA_table(outfname)
    next
  }
  
  # which region have more strog binding site accumulation ( scores above Q3)
  thr <- quantile(density(scores)$x)[[4]] + abs(quantile(density(scores)$x)[[4]] - max(density(scores)$x)) / 2
  extreme_max_idxs <- which(scores > thr)
  max_idx <- get_max_peak(above_avg_mod_idxes, extreme_max_idxs)
  
  # set initial lef and right indexes for peak finding
  peak_mod_index <- max_idx
  leftind <- peak_mod_index
  rightind <- peak_mod_index
  peaks_scores <- smooth_scores[above_avg_mod_idxes]
  # if the next highest peak is close, then set the right index at the second peak
  if (length(above_avg_mod_idxes) == 2) {
    if (abs(above_avg_mod_idxes[1] - above_avg_mod_idxes[2]) < max_neighbor_dist) {
      if (max_idx == above_avg_mod_idxes[1]) {
        peak_mod_index <- above_avg_mod_idxes[1]
        leftind <- peak_mod_index
        rightind <- above_avg_mod_idxes[2]
      } else {
        peak_mod_index <- above_avg_mod_idxes[2]
        leftind <- above_avg_mod_idxes[1]
        rightin <- peak_mod_index
      }
    }
  } 
  
  #### Use old method to find region but using smoothed values
  # Find Left 
  maxbgval <- mean(density(smooth_scores[1:leftind])$x)
  nBelow <- 0 
  no_roe <- FALSE
  while (nBelow < buffer_left) { 
    if (leftind > length(scores) | leftind < 1) { 
      print("No ROE! left index is out of bound!")
      printout_NA_table(outfname)
      no_roe <- TRUE
    } 
    if (smooth_scores[leftind] < maxbgval) { 
      nBelow <- nBelow + 1 
    } 
    leftind <- leftind - 1 
  } 
  
  if (no_roe) 
    next
  
  # Find Right
  maxbgval <- mean(density(smooth_scores[rightind:length(smooth_scores)])$x)
  print(maxbgval)
  nBelow <- 0
  while (nBelow < buffer_right) {
    if (rightind > length(scores)) {
      print("No ROE! right index is out of bound!")
      printout_NA_table(outfname)
      no_roe <- TRUE
      break
    }
    if (smooth_scores[rightind] < maxbgval) { 
      nBelow <- nBelow + 1 
    } 
    rightind <- rightind + 1 
  } 
  if (no_roe)
    next
  
  if ((leftind == peak_mod_index) || (rightind == peak_mod_index)) { 
    print("No ROE found! right or left index got stuck on the mode!")
    printout_NA_table(outfname)
    next
  } 
  
  ## print final ROE
  plot_roe()
}

#------------ End of the loop #------------------------ 






