classes1 <- read.table("~/Downloads/1298.txt", header = F)
load("~/Downloads/roe_only_x_train_set_scaled_1298.rdat")
load("~/Downloads/crossVals_roe/1298/roe_only_x_test_set_scaled.rdat")
x_test_set_scaled_1298 = as.data.frame(x_test_set_scaled)
x_train_set_scaled_1298 <- x_train_set_scaled
test_1298_classes <- read.table("~/Downloads/crossVals_roe/1298/test_set_classes.txt")

model <- LiblineaR::LiblineaR(data = x_train_set_scaled_1298,
                              target = classes1$V2,
                              type = 7,   # L2 regularized logistic regression (dual) = 7
                              cost = 0.00224,
                              bias = T,  # ?? (recommended by vignette)
                              # cross = 10, # built-in cross validation; probably better to do it ourselves
                              verbose = TRUE)

p <- predict(model, x_test_set_scaled, proba = TRUE, decisionValues = TRUE)
confusion_matrix <- table(predictions = p$predictions, actuals = test_1298_classes$x)

auroc <- PRROC::roc.curve(p$probabilities[,"1"], weights.class0 = test_1298_classes$x, curve = TRUE)$auc
auprc <- PRROC::pr.curve(p$probabilities[,"1"], weights.class0 = test_1298_classes$x, curve = TRUE)$auc.davis.goadrich


d <- as.data.frame(colnames(x_train_set_scaled_1298))
write.table(d, "labels.txt", col.names = F, quote = F)

x_train_set_scaled_1298 <- as.data.frame(x_train_set_scaled_1298)
x_train_set_scaled_1298 <- cbind(x_train_set_scaled_1298, class=classes1$V2)
x_test_set_scaled_1298 <- cbind(x_test_set_scaled_1298, test_1298_classes$x)
write.table(x_test_set_scaled_1298, file = "test_1298.csv", quote = F, sep = ",", col.names = F, row.names = T)

classes2 <- read.table("~/Downloads/12421.txt", header = F)
load("~/Downloads/roe_only_x_train_set_scaled_12421.rdat")
x_train_set_scaled_12421 <- x_train_set_scaled
load("~/Downloads/crossVals_roe/12421//roe_only_x_test_set_scaled.rdat")
x_test_set_scaled_12421 = as.data.frame(x_test_set_scaled)
x_train_set_scaled_12421 <- x_train_set_scaled
test_12421_classes <- read.table("~/Downloads/crossVals_roe/12421/test_set_classes.txt")

model2 <- LiblineaR::LiblineaR(data = x_train_set_scaled_12421,
                              target = classes2$V2,
                              type = 7,   # L2 regularized logistic regression (dual) = 7
                              cost = 0.004,
                              bias = T,  # ?? (recommended by vignette)
                              # cross = 10, # built-in cross validation; probably better to do it ourselves
                              verbose = TRUE)
p <- predict(model, x_test_set_scaled_1298, proba = TRUE, decisionValues = TRUE)
confusion_matrix <- table(predictions = p$predictions, actuals = test_12421_classes$x)
confusion_matrix

x_train_set_scaled_12421 <- as.data.frame(x_train_set_scaled_12421)
x_train_set_scaled_12421 <- cbind(x_train_set_scaled_12421, class=classes2$V2)
write.table(x_train_set_scaled_12421, "~/data_12421.csv", quote = F,col.names = F, sep = ",")
#----------- 
df1 <- as.data.frame(model1$W)
b1 <- df1$Bias
df1$Bias <- NULL
m1 <- gather(df1, key = feature, value = value)

df2 <- as.data.frame(model2$W)
m2 <- gather(df2, key = feature, value = value)

#-----------
#---
m <- merge(classes1, classes2, by = "V1")
#c1 <- classes1[!(classes1$V1 %in% m$V1),]
c1 <- classes1[(classes1$V1 %in% m$V1),]
x_train_set_scaled_1298 <- as.data.frame(x_train_set_scaled_1298)
x_train_set_scaled_1298$tss_name <- rownames(x_train_set_scaled_1298)
x1 <- extract(x_train_set_scaled_1298, tss_name, into=c("fold", "tss_name"), regex ="(.)\\.(.+)")
x1$fold <- NULL

q <- "AT1G50732.1_Chr1_18799651_+_0" #1
t1 <- x1[x1$tss_name == q,]
t1$tss_name <- NULL
t1 <- gather(t1, key = feature, value = value)
multi1 <- merge(m1, t1, by= "feature")
sum(multi1$value.x * multi1$value.y)

x1 <- merge(x1, c1, by.x = "tss_name", by.y = "V1")
x1$tss_name <- NULL
clas <- x1$V2
x1$V2 <- NULL
model3 <- LiblineaR::LiblineaR(data = x1,
                              target = clas,
                              type = 7,   # L2 regularized logistic regression (dual) = 7
                              cost = 0.004,
                              bias = T,  # ?? (recommended by vignette)
                              # cross = 10, # built-in cross validation; probably better to do it ourselves
                              verbose = TRUE)
p <- predict(model2, proba = T, x1)
confusion_matrix <- table(predictions = p$predictions, actuals = clas)
#---
c2 <- classes2[!(classes2$V1 %in% m$V1),]
x_train_set_scaled_12421 <- as.data.frame(x_train_set_scaled_12421)
x_train_set_scaled_12421$tss_name <- rownames(x_train_set_scaled_12421)
x2 <- extract(x_train_set_scaled_12421, tss_name, into=c("fold", "tss_name"), regex ="(.)\\.(.+)")
x2$fold <- NULL
x2 <- merge(x2, c2, by.x = "tss_name", by.y = "V1")
x2$tss_name <- NULL
clas <- x2$V2
x2$V2 <- NULL
model1 <- LiblineaR::LiblineaR(data = x2,
                               target = clas,
                               type = 7,   # L2 regularized logistic regression (dual) = 7
                               cost = 0.004,
                               bias = T,  # ?? (recommended by vignette)
                               # cross = 10, # built-in cross validation; probably better to do it ourselves
                               verbose = TRUE)
#---
tss_id <- "0.AT1G62800.1_Chr1_23254523_-_0"
a1 <- as.data.frame(x_train_set_scaled_1298[rownames(x_train_set_scaled_1298) == tss_id,])
colnames(a1) <- "coef"
a1$feature <- rownames(a1)

tss_id <- "4.AT1G62800.1_Chr1_23254523_-_0"
a2 <- as.data.frame(x_train_set_scaled_12421[rownames(x_train_set_scaled_12421) == tss_id,])
colnames(a2) <- "coef"
a2$feature <- rownames(a2)

sum(a1 * model$W)
sum(a2 * model2$W)

plot(a1$coef, a2$coef)


m <- merge(a1, a2, by = "feature")
