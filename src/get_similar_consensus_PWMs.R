library(seqLogo)

indir = "~/Downloads/sim_pwms/"
files <- list.files(indir)
df <- data.frame(name=c(), s=c())
for (f in files) {
  n = basename(f)
  f = paste(indir, f, sep = "")
  m <- t(read.table(f))
  pwm <- makePWM(m)
  
  #
  a <- unlist(strsplit(pwm@consensus, ""))
  idxs <- which(pwm@ic > 0.4)
  idxs_p <- which(pwm@ic <= 0.4)
  a[idxs] <- a[idxs]
  a[idxs_p] <- 'N'
  print(a)
  #
  df <- rbind(df, data.frame(name=n, s=str_flatten(a)))
}
outfile = "~/Downloads/consensus_pwms_withIC.txt"
write.table(na.omit(df), file = outfile, sep = "\t", quote = F, row.names = F, col.names = F)

#0------
cons <- read.table(outfile, col.names = c("pwm", "consensus"))
head(cons)

df <- data.frame(m1=c(), m2=c(), c1=c(), c2=c())
for (i in seq(1,nrow(cons))) {
  for (j in seq(1, nrow(cons))) {
    x <- cons[i,2]
    #print(x)
    y <- cons[j,2]
    a <- agrep(x,y,max.distance = 0)
    if (!is.na(a[1])) {
      df = rbind(df, data.frame(m1=cons[i,1], m2=cons[j,1], c1=cons[i,2], c2=cons[j,2]))
    }
  }
}
outfile = "~/Downloads/consensus_pwms_result_mismatch0.txt"
df <- df[df$m1 != df$m2,]
write.table(df, file = outfile, sep = "\t", quote = F, row.names = F, col.names = F)


####
f = "M0119_1.02.txt"
f = "M0118_1.02.txt"
f <- paste(indir, f, sep = "")
m <- t(read.table(f))
pwm <- makePWM(m)
pwm@consensus
pwm@ic
#pwm

a <- unlist(strsplit(pwm@consensus, ""))
idxs <- which(pwm@ic > 0.4)
idxs_p <- which(pwm@ic <= 0.4)
a[idxs] <- a[idxs]
a[idxs_p] <- 'N'
a
str_flatten(a)
