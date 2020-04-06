indir="~/Downloads/ibdc/aug2019/3peat/coefs/"
froot_coefs = paste(indir, "ALL_50_ROOT_IBDC_Train_-5000_5000.20_negfrom_-2000_-200.TAIR10_cds_16189.simple.SumScoreVars.0.00022.model.UnModCoefs.plot.table.fac.summary", sep = "")
fshoot_coefs = paste(indir, "ALL_50_SHOOT_IBDC_Train_-5000_5000.20_negfrom_-2000_-200.TAIR10_cds_14839.simple.SumScoreVars.0.00026.model.UnModCoefs.plot.table.fac.summary" , sep = "")

shoot = read.table(fshoot_coefs, header = F)
root = read.table(froot_coefs, header = F)
colnames(shoot) <- c("PWM", "N", "ModelWeight")
colnames(root) <- c("PWM", "N", "ModelWeight")

root$tissue <- "Root"
shoot$tissue <- "Shoot"


library(ggplot2)
both <- rbind(root[1:20,], shoot[1:20,])

ggplot(both, aes(PWM, ModelWeight, color = tissue)) + geom_point() + theme_bw()
