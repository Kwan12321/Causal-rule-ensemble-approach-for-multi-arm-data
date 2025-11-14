library(speff2trial)
library(doParallel)
library(ggplot2)
#########################

dir <- getwd()
source(paste0(dir, "/mcre_agenet.R"))
source(paste0(dir, "/imp.R"))

###############
### Dataset ###
###############

data(ACTG175)
dat <- ACTG175

seed <-123
set.seed(seed)

Y <- as.numeric(dat$cd420 - dat$cd40)
W <- dat$arm
X <- dat[, c("age", "wtkg", "hemo", "homo", "drugs", "karnof", "race", "gender", "str2", "symptom",  "cd40", "cd80")]
dat <- data.frame(Y = Y, W = W, X)
dat <- dat[dat$cd40>= 200 & dat$cd40 <= 500,]
Y <- dat$Y
W <- dat$W
X <- dat[,-c(1:2)]

######################
### Model Building ###
######################

### Estimate the propensity score ###

PS.RF <- nnet::multinom(W~., data = data.frame(W = as.factor(W), X))
W.hat <- predict(PS.RF, data = data.frame(X), type = "prob")

### Fit the model ### 
fit <- mcre_agenet(X = X, Y = Y, W = W, W.hat = W.hat, gen.tree = "gbm", tune.init = "CV", tune.fin = "CV", seed = seed, nfolds = 5, learnrate = 0.01, ntrees = 1000, depth = 3)
pred <- predict.mcre(fit, newdata = data.frame(Y = Y, W = as.factor(W), X))
tau.beta <- pred$tau

###############
### Table 2 ###
###############

EST <- tau.beta[c(9, 16, 333),]
EST <- data.frame(EST)
EST$c <- EST[,3] - EST[,2]
colnames(EST) <- c("ZDV vs ZDV + DID", "ZDV vs ZDV + ZAL", "ZDV vs DID", "ZDV + ZAL vs DID")
EST <- cbind(ID = ACTG175$pidnum[c(9, 16, 333)], round(EST,1))
write.csv(EST, paste0(dir,"/Table_Figure/Section5/Table2.csv"))

################
### Figure 6 ###
################

opt.id1 <- (tau.beta[,1] > 0 & dat$W == 1)|(tau.beta[,1] < 0 & dat$W == 0)
opt.id2 <- (tau.beta[,2] > 0 & dat$W == 2)|(tau.beta[,2] < 0 & dat$W == 0)
opt.id3 <- (tau.beta[,3] > 0 & dat$W == 3)|(tau.beta[,3] < 0 & dat$W == 0)

df1 <- data.frame(Y = dat$Y, W = dat$W, id = opt.id1)[(dat$W == 0|dat$W == 1),]
df2 <- data.frame(Y = dat$Y, W = dat$W, id = opt.id2)[(dat$W == 0|dat$W == 2),]
df3 <- data.frame(Y = dat$Y, W = dat$W, id = opt.id3)[(dat$W == 0|dat$W == 3),]

df1$W <- ifelse(df1$W == 0, "ZDV (Control)", "ZDV + DID")
df2$W <- ifelse(df2$W == 0, "ZDV (Control)", "ZDV + ZAL")
df3$W <- ifelse(df3$W == 0, "ZDV (Control)", "DID")

df1$id <- ifelse(df1$id, "Optimal", "Non-optimal")
df2$id <- ifelse(df2$id, "Optimal", "Non-optimal")
df3$id <- ifelse(df3$id, "Optimal", "Non-optimal")

df <- rbind(df1,df2,df3)
df$W <- factor(df$W, levels = c("ZDV (Control)", "ZDV + DID", "ZDV + ZAL", "DID"))

g0 <- ggplot(df, aes(x = id, y = Y, fill = id)) + geom_boxplot() + facet_wrap(~W, ncol=4)
g0 <- g0 + ylab("Change in CD4 Count (cells/mm3)") + xlab("")
g0 <- g0 + theme_bw()  + theme(legend.position="none", text = element_text(size = 20))
ggsave(paste0(dir, "/Table_Figure/Section5/Fig7.eps"),g0 , width = 14, height = 7)

###############
### Table 3 ###
###############

### Base importance ------------------------------------------------------------
vip.can <- variable.importance(fit, data.frame(Y = Y, W = W, X))
r.id <- which(apply(vip.can$base.imp,1,sum)!= 0) 
vip.all <- vip.can$base.imp
supp <- apply(tran_rules(dat, rownames(vip.can$base.imp)),2, mean)
beta.all <- matrix(fit$beta, nrow = length(fit$beta)/4)

### ZDV vs ZDV + DID -----------------------------------------------------------
beta1 <- round((beta.all[,2] - beta.all[,1])[r.id],1)
vip1 <- vip.all[r.id,1]
rules1 <- rownames(vip.can$base.imp)[r.id]
supp1 <- round(supp[r.id],2)
sum1 <- data.frame(Rules = rules1, Importance = vip1, Coef = round(beta1,1), Support = round(supp1,2))
sum1$Importance <- round(100*(sum1$Importance)/max(sum1$Importance),0)
sum1 <- sum1[order(sum1$Importance, decreasing = TRUE),]
rownames(sum1) <- 1:nrow(sum1)
r.select1 <- which(sum1$Support >= 0.1)
sum1 <- sum1[r.select1,][1:5,]

### ZDV vs ZDV + ZAL -----------------------------------------------------------
beta2 <- round((beta.all[,3] - beta.all[,1])[r.id],1)
vip2 <- vip.all[r.id,2]
rules2 <- rownames(vip.can$base.imp)[r.id]
supp2 <- round(supp[r.id],2)
sum2 <- data.frame(Rules = rules2, Importance = vip2, Coef = round(beta2,1), Support = round(supp2,2))
sum2$Importance <- round(100*(sum2$Importance)/max(sum2$Importance),0)
sum2 <- sum2[order(sum2$Importance, decreasing = TRUE),]
rownames(sum2) <- 1:nrow(sum2)
r.select2 <- which(sum2$Support >= 0.1)
sum2 <- sum2[r.select2,][1:5,]

### ZDV vs DID -----------------------------------------------------------------
beta3 <- round((beta.all[,4] - beta.all[,1])[r.id],1)
vip3 <- vip.all[r.id,3]
rules3 <- rownames(vip.can$base.imp)[r.id]
supp3 <- round(supp[r.id],2)
sum3 <- data.frame(Rules = rules3, Importance = vip3, Coef = round(beta3,1), Support = round(supp3,2))
sum3$Importance <- round(100*(sum3$Importance)/max(sum3$Importance),0)
sum3 <- sum3[order(sum3$Importance, decreasing = TRUE),]
rownames(sum3) <- 1:nrow(sum3)
r.select3 <- which(sum3$Support >= 0.1)
sum3 <- sum3[r.select3,][1:5,]

### ZDV + ZAL vs DID -----------------------------------------------------------
beta4 <- round((beta.all[,4] - beta.all[,3])[r.id], 1)
vip4 <- vip.all[r.id,6]
rules4 <- rownames(vip.can$base.imp)[r.id]
supp4 <- round(supp[r.id],2)
sum4 <- data.frame(Rules = rules4, Importance = vip4, Coef = round(beta4,1), Support = round(supp4,2))
sum4$Importance <- round(100*(sum4$Importance)/max(sum4$Importance),0)
sum4 <- sum4[order(sum4$Importance, decreasing = TRUE),]
rownames(sum4) <- 1:nrow(sum4)
r.select4 <- which(sum4$Support >= 0.1)
sum4 <- sum4[r.select4,][1:5,]

colnames(sum1) <- colnames(sum2) <- colnames(sum3) <- colnames(sum4) <- c("Rules","Importance","Support","Coefficients", "Actual HTEs")[-5]
rownames(sum1) <- paste0("Rule", rownames(sum1))
rownames(sum2) <- paste0("Rule", rownames(sum2))
rownames(sum3) <- paste0("Rule", rownames(sum3))
rownames(sum4) <- paste0("Rule", rownames(sum4))

write.csv(sum1, paste0(dir,"/Table_Figure/Section5/Table3-A.csv"))
write.csv(sum2, paste0(dir,"/Table_Figure/Section5/Table3-B.csv"))
write.csv(sum3, paste0(dir,"/Table_Figure/Section5/Table3-C.csv"))
write.csv(sum4, paste0(dir,"/Table_Figure/Section5/Table3-D.csv"))

################
### Figure 7 ###
################

### Variable importance (Proposed method) ###

vip.imp.sep <- vip.can$var.imp[,c(1:3,6)]
for(r in 1:ncol(vip.imp.sep)){
  vip.imp.sep[,r] <- round(100*vip.imp.sep[,r]/max(vip.imp.sep[,r]),1)
}

d1 <- vip.imp.sep[,1]
d1 <- data.frame(names = names(d1), value = d1)
d1$names <- factor(d1$names, levels = d1$names[order(d1$value, decreasing = TRUE)])
g1 <- ggplot(d1, aes(x = names, y= value)) + geom_bar(stat = "identity", position = "dodge")
g1 <- g1 + xlab("") + ylab("Variable importance") + theme(text = element_text(size = 15))
g1 <- g1 + theme_minimal() + theme(text = element_text(size = 15), axis.text.x = element_text(face = ifelse(levels(d1$names) %in% c("cd40", "cd80", "wtkg", "age"), "bold", "plain"))
)

d2 <- vip.imp.sep[,2]
d2 <- data.frame(names = names(d2), value = d2)
d2$names <- factor(d2$names, levels = d2$names[order(d2$value, decreasing = TRUE)])
g2 <- ggplot(d2, aes(x = names, y= value)) + geom_bar(stat = "identity", position = "dodge")
g2 <- g2 + xlab("") + ylab("Variable importance") + theme(text = element_text(size = 15))
g2 <- g2 + theme_minimal() + theme(text = element_text(size = 15), axis.text.x = element_text(face = ifelse(levels(d2$names) %in% c("cd40", "cd80", "wtkg", "age"), "bold", "plain"))
)

d3 <- vip.imp.sep[,3]
d3 <- data.frame(names = names(d3), value = d3)
d3$names <- factor(d3$names, levels = d3$names[order(d3$value, decreasing = TRUE)])
g3 <- ggplot(d3, aes(x = names, y= value)) + geom_bar(stat = "identity", position = "dodge")
g3 <- g3 + xlab("") + ylab("Variable importance") + theme(text = element_text(size = 15))
g3 <- g3 + theme_minimal() + theme(text = element_text(size = 15), axis.text.x = element_text(face = ifelse(levels(d3$names) %in% c("cd40", "cd80", "wtkg", "age"), "bold", "plain"))
)

d4 <- vip.imp.sep[,4]
d4 <- data.frame(names = names(d4), value = d4)
d4$names <- factor(d4$names, levels = d4$names[order(d4$value, decreasing = TRUE)])
g4 <- ggplot(d4, aes(x = names, y= value)) + geom_bar(stat = "identity", position = "dodge")
g4 <- g4 + xlab("") + ylab("Variable importance")
g4 <- g4 + theme_minimal() + theme(text = element_text(size = 15), axis.text.x = element_text(face = ifelse(levels(d4$names) %in% c("cd40", "cd80", "wtkg", "age"), "bold", "plain"))
)

### Variable importance (Causal Forest) ###

dat1 <- dat[(dat$W == 1|dat$W==0),]
fit1 <- grf::causal_forest(X = dat1[,-c(1:2)] ,Y = dat1$Y, W = dat1$W, seed = seed)
vip_cf1 <- grf::variable_importance(fit1)
vip.cf1 <- data.frame(names = colnames(X), value = vip_cf1*100/max(vip_cf1))
vip.cf1$names <- factor(vip.cf1$names, levels = vip.cf1$names[order(vip.cf1$value, decreasing = TRUE)])
gc1 <- ggplot(vip.cf1, aes(x = names, y= value)) + geom_bar(stat = "identity", position = "dodge")
gc1 <- gc1 + xlab("") + ylab("Variable importance")
gc1 <- gc1 + theme_minimal() + theme(text = element_text(size = 15), axis.text.x = element_text(face = ifelse(levels(vip.cf1$names) %in% c("cd40", "cd80", "wtkg", "age"), "bold", "plain"))
)


dat2 <- dat[(dat$W == 2|dat$W==0),]
dat2$W[dat2$W == 2] <- 1
fit2 <- grf::causal_forest(X = dat2[,-c(1:2)] ,Y = dat2$Y, W = dat2$W, seed = seed)
vip_cf2 <- grf::variable_importance(fit2)
vip.cf2 <- data.frame(names = colnames(X), value = vip_cf2*100/max(vip_cf2))
vip.cf2$names <- factor(vip.cf2$names, levels = vip.cf2$names[order(vip.cf2$value, decreasing = TRUE)])
gc2 <- ggplot(vip.cf2, aes(x = names, y= value)) + geom_bar(stat = "identity", position = "dodge")
gc2 <- gc2 + xlab("") + ylab("Variable importance")
gc2 <- gc2 + theme_minimal() + theme(text = element_text(size = 15), axis.text.x = element_text(face = ifelse(levels(vip.cf2$names) %in% c("cd40", "cd80", "wtkg", "age"), "bold", "plain"))
)

dat3 <- dat[(dat$W == 3|dat$W==0),]
dat3$W[dat3$W == 3] <- 1
fit3 <- grf::causal_forest(X = dat3[,-c(1:2)] ,Y = dat3$Y, W = dat3$W, seed = seed)
vip_cf3 <- grf::variable_importance(fit3)
vip.cf3 <- data.frame(names = colnames(X), value = vip_cf3*100/max(vip_cf3))
vip.cf3$names <- factor(vip.cf3$names, levels = vip.cf3$names[order(vip.cf3$value, decreasing = TRUE)])
gc3 <- ggplot(vip.cf3, aes(x = names, y= value)) + geom_bar(stat = "identity", position = "dodge")
gc3 <- gc3 + xlab("") + ylab("Variable importance")
gc3 <- gc3 + theme_minimal() + theme(text = element_text(size = 15), axis.text.x = element_text(face = ifelse(levels(vip.cf3$names) %in% c("cd40", "cd80", "wtkg", "age"), "bold", "plain"))
)

dat4 <- dat[(dat$W == 2|dat$W==3),]
dat4$W[dat4$W == 2] <- 0
dat4$W[dat4$W == 2] <- 1
fit4 <- grf::causal_forest(X = dat4[,-c(1:2)] ,Y = dat4$Y, W = dat4$W, seed = seed)
vip_cf4 <- grf::variable_importance(fit4)
vip.cf4 <- data.frame(names = colnames(X), value = vip_cf4*100/max(vip_cf4))
vip.cf4$names <- factor(vip.cf4$names, levels = vip.cf4$names[order(vip.cf4$value, decreasing = TRUE)])
gc4 <- ggplot(vip.cf4, aes(x = names, y= value)) + geom_bar(stat = "identity", position = "dodge")
gc4 <- gc4 + xlab("") + ylab("Variable importance") 
gc4 <- gc4 + theme_minimal() + theme(text = element_text(size = 15), axis.text.x = element_text(face = ifelse(levels(vip.cf4$names) %in% c("cd40", "cd80", "wtkg", "age"), "bold", "plain"))
)


gg <-  grid.arrange(arrangeGrob(g1 + labs(title = "A1) ZDV vs ZDV + DID (Propposed method)"),
                                gc1 + labs(title = "A2) ZDV vs ZDV + DID (Causal forest)"),
                                g2 + labs(title = "B1) ZDV vs ZDV + ZAL (Propposed method)"),
                                gc2 + labs(title = "B2) ZDV vs ZDV + ZAL (Causal forest)"),
                                g3 + labs(title = "C1) ZDV vs DID (Propposed method)"),
                                gc3 + labs(title = "C2) ZDV vs DID (Causal forest)"),
                                g4 + labs(title = "D1) ZDV + ZAL vs DID (Propposed method)"),
                                gc4 + labs(title = "D2) ZDV + ZAL vs DID (Causal forest)"),nrow = 4), nrow = 1)

ggsave(paste0(dir,"/Table_Figure/Section5/Fig8.eps"),gg , width = 18, height = 15)





