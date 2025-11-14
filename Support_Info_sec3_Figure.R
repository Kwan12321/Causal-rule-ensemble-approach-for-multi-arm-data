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

#################
### Figure 18 ###
#################

### ZDV vs ZDV + DID -----------------------------------------------------------
dat1 <- data.frame(Y = Y, W = W, tau = tau.beta[,1])
dat1 <- dat1[which(dat1$W == 0 | dat1$W == 1),]
dat1 <- dat1[order(dat1$tau),]
id1 <- 1:nrow(dat1)
split_vec1 <- split(id1, cut(seq_along(id1), 5, labels = FALSE))
ate1 <- c()
ate1.est <- c()
for(l in 1:length(split_vec1)){
  dat.sub1 <- dat1[split_vec1[[l]],]
  ate1.t <- mean(dat.sub1$Y[dat.sub1$W == 1]) - mean(dat.sub1$Y[dat.sub1$W == 0])
  ate1.est.t <- mean(dat.sub1$tau)
  ate1 <- c(ate1, ate1.t)
  ate1.est <- c(ate1.est, ate1.est.t)
}
ate1.dat <- data.frame(senario = "ZDV vs ZDV + DID", 
                       ate = c(ate1, ate1.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p1 <- clinfun::jonckheere.test(ate1.dat$ate[1:5], ate1.dat$ate[6:10])$p.value
cor1 <- cor(ate1.dat$ate[1:5], ate1.dat$ate[6:10], method = "spearman")

### ZDV vs ZDV + ZAL -----------------------------------------------------------
dat2 <- data.frame(Y = Y, W = W, tau = tau.beta[,2])
dat2 <- dat2[which(dat2$W == 0 | dat2$W == 2),]
dat2 <- dat2[order(dat2$tau),]
id2 <- 1:nrow(dat2)
split_vec2 <- split(id2, cut(seq_along(id2), 5, labels = FALSE))
ate2 <- c()
ate2.est <- c()
for(l in 1:length(split_vec2)){
  dat.sub2 <- dat2[split_vec2[[l]],]
  ate2.t <- mean(dat.sub2$Y[dat.sub2$W == 2]) - mean(dat.sub2$Y[dat.sub2$W == 0])
  ate2.est.t <- mean(dat.sub2$tau)
  ate2 <- c(ate2, ate2.t)
  ate2.est <- c(ate2.est, ate2.est.t)
}
ate2.dat <- data.frame(senario = "ZDV vs ZDV + ZAL", 
                       ate = c(ate2, ate2.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p2 <- clinfun::jonckheere.test(ate2.dat$ate[1:5], ate2.dat$ate[6:10])$p.value
cor2 <- cor(ate2.dat$ate[1:5], ate2.dat$ate[6:10], method = "spearman")

### ZDV vs DID -----------------------------------------------------------
dat3 <- data.frame(Y = Y, W = W, tau = tau.beta[,3])
dat3 <- dat3[which(dat3$W == 0 | dat3$W == 3),]
dat3 <- dat3[order(dat3$tau),]
id3 <- 1:nrow(dat3)
split_vec3 <- split(id3, cut(seq_along(id3), 5, labels = FALSE))
ate3 <- c()
ate3.est <- c()
for(l in 1:length(split_vec3)){
  dat.sub3 <- dat3[split_vec3[[l]],]
  ate3.t <- mean(dat.sub3$Y[dat.sub3$W == 3]) - mean(dat.sub3$Y[dat.sub3$W == 0])
  ate3.est.t <- mean(dat.sub3$tau)
  ate3 <- c(ate3, ate3.t)
  ate3.est <- c(ate3.est, ate3.est.t)
}
ate3.dat <- data.frame(senario = "ZDV vs DID", 
                       ate = c(ate3, ate3.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p3 <- clinfun::jonckheere.test(ate3.dat$ate[1:5], ate3.dat$ate[6:10])$p.value
cor3 <- cor(ate3.dat$ate[1:5], ate3.dat$ate[6:10], method = "spearman")

ate.dat <- rbind(ate1.dat, ate2.dat, ate3.dat)
ate.dat$senario <- factor(ate.dat$senario, levels = unique(ate.dat$senario))


library(ggplot2)
library(gridExtra)
g1 <- ggplot(ate.dat[1:10,], aes(x = sub, y = ate, fill = sub.c))
g1 <- g1 + geom_bar(stat = "identity", position = "dodge")
g1 <- g1 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g1 <- g1 + coord_cartesian(ylim = c(-100, 250))
g1 <- g1 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g1 <- g1 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

g2 <- ggplot(ate.dat[11:20,], aes(x = sub, y = ate, fill = sub.c))
g2 <- g2 + geom_bar(stat = "identity", position = "dodge")
g2 <- g2 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g2 <- g2 + coord_cartesian(ylim = c(-100, 250))
g2 <- g2 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g2 <- g2 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

g3 <- ggplot(ate.dat[21:30,], aes(x = sub, y = ate, fill = sub.c))
g3 <- g3 + geom_bar(stat = "identity", position = "dodge")
g3 <- g3 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g3 <- g3 + coord_cartesian(ylim = c(-100, 250))
g3 <- g3 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g3 <- g3 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(g1)

g <- grid.arrange(arrangeGrob(g1 + labs(title = "A) ZDV vs ZDV + DID"), 
                              g2 + labs(title = "B) ZDV vs ZDV + ZAL"), 
                              g3 + labs(title = "C) ZDV vs DID"), 
                              nrow=1))

ggsave(paste0(dir,"/Table_Figure_support/sup_section3/s-Fig18.eps"),g , width = 16, height = 7)

###################################################################
### Figure 19-b  gbm + adaptive group lasso (Same to Figure 18) ###
###################################################################
ggsave(paste0(dir,"/Table_Figure_support/sup_section3/s-Fig19-b.eps"),g , width = 16, height = 6)


######################################
### Figure 19-a  gbm + group lasso ###
######################################

### Fit the model ### 
fit <- mcre_agenet(X = X, Y = Y, W = W, W.hat = W.hat, gen.tree = "gbm", tune.init = "CV", tune.fin = "CV", seed = seed, nfolds = 5, learnrate = 0.1, ntrees = 333, depth = 4)
pred <- predict.mcre(fit, newdata = data.frame(Y = Y, W = as.factor(W), X))
tau.beta <- pred$tau.first

### ZDV vs ZDV + DID -----------------------------------------------------------
dat1 <- data.frame(Y = Y, W = W, tau = tau.beta[,1])
dat1 <- dat1[which(dat1$W == 0 | dat1$W == 1),]
dat1 <- dat1[order(dat1$tau),]
id1 <- 1:nrow(dat1)
split_vec1 <- split(id1, cut(seq_along(id1), 5, labels = FALSE))
ate1 <- c()
ate1.est <- c()
for(l in 1:length(split_vec1)){
  dat.sub1 <- dat1[split_vec1[[l]],]
  ate1.t <- mean(dat.sub1$Y[dat.sub1$W == 1]) - mean(dat.sub1$Y[dat.sub1$W == 0])
  ate1.est.t <- mean(dat.sub1$tau)
  ate1 <- c(ate1, ate1.t)
  ate1.est <- c(ate1.est, ate1.est.t)
}
ate1.dat <- data.frame(senario = "ZDV vs ZDV + DID", 
                       ate = c(ate1, ate1.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p1 <- clinfun::jonckheere.test(ate1.dat$ate[1:5], ate1.dat$ate[6:10])$p.value
cor1 <- cor(ate1.dat$ate[1:5], ate1.dat$ate[6:10], method = "spearman")

### ZDV vs ZDV + ZAL -----------------------------------------------------------
dat2 <- data.frame(Y = Y, W = W, tau = tau.beta[,2])
dat2 <- dat2[which(dat2$W == 0 | dat2$W == 2),]
dat2 <- dat2[order(dat2$tau),]
id2 <- 1:nrow(dat2)
split_vec2 <- split(id2, cut(seq_along(id2), 5, labels = FALSE))
ate2 <- c()
ate2.est <- c()
for(l in 1:length(split_vec2)){
  dat.sub2 <- dat2[split_vec2[[l]],]
  ate2.t <- mean(dat.sub2$Y[dat.sub2$W == 2]) - mean(dat.sub2$Y[dat.sub2$W == 0])
  ate2.est.t <- mean(dat.sub2$tau)
  ate2 <- c(ate2, ate2.t)
  ate2.est <- c(ate2.est, ate2.est.t)
}
ate2.dat <- data.frame(senario = "ZDV vs ZDV + ZAL", 
                       ate = c(ate2, ate2.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p2 <- clinfun::jonckheere.test(ate2.dat$ate[1:5], ate2.dat$ate[6:10])$p.value
cor2 <- cor(ate2.dat$ate[1:5], ate2.dat$ate[6:10], method = "spearman")

### ZDV vs DID -----------------------------------------------------------
dat3 <- data.frame(Y = Y, W = W, tau = tau.beta[,3])
dat3 <- dat3[which(dat3$W == 0 | dat3$W == 3),]
dat3 <- dat3[order(dat3$tau),]
id3 <- 1:nrow(dat3)
split_vec3 <- split(id3, cut(seq_along(id3), 5, labels = FALSE))
ate3 <- c()
ate3.est <- c()
for(l in 1:length(split_vec3)){
  dat.sub3 <- dat3[split_vec3[[l]],]
  ate3.t <- mean(dat.sub3$Y[dat.sub3$W == 3]) - mean(dat.sub3$Y[dat.sub3$W == 0])
  ate3.est.t <- mean(dat.sub3$tau)
  ate3 <- c(ate3, ate3.t)
  ate3.est <- c(ate3.est, ate3.est.t)
}
ate3.dat <- data.frame(senario = "ZDV vs DID", 
                       ate = c(ate3, ate3.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p3 <- clinfun::jonckheere.test(ate3.dat$ate[1:5], ate3.dat$ate[6:10])$p.value
cor3 <- cor(ate3.dat$ate[1:5], ate3.dat$ate[6:10], method = "spearman")

ate.dat <- rbind(ate1.dat, ate2.dat, ate3.dat)
ate.dat$senario <- factor(ate.dat$senario, levels = unique(ate.dat$senario))


library(ggplot2)
library(gridExtra)
g1 <- ggplot(ate.dat[1:10,], aes(x = sub, y = ate, fill = sub.c))
g1 <- g1 + geom_bar(stat = "identity", position = "dodge")
g1 <- g1 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g1 <- g1 + coord_cartesian(ylim = c(-100, 250))
g1 <- g1 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g1 <- g1 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

g2 <- ggplot(ate.dat[11:20,], aes(x = sub, y = ate, fill = sub.c))
g2 <- g2 + geom_bar(stat = "identity", position = "dodge")
g2 <- g2 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g2 <- g2 + coord_cartesian(ylim = c(-100, 250))
g2 <- g2 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g2 <- g2 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

g3 <- ggplot(ate.dat[21:30,], aes(x = sub, y = ate, fill = sub.c))
g3 <- g3 + geom_bar(stat = "identity", position = "dodge")
g3 <- g3 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g3 <- g3 + coord_cartesian(ylim = c(-100, 250))
g3 <- g3 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g3 <- g3 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(g1)

g <- grid.arrange(arrangeGrob(g1 + labs(title = "A) ZDV vs ZDV + DID"), 
                              g2 + labs(title = "B) ZDV vs ZDV + ZAL"), 
                              g3 + labs(title = "C) ZDV vs DID"), 
                              nrow=1))

ggsave(paste0(dir,"/Table_Figure_support/sup_section3/s-Fig19-a.eps"),g , width = 16, height = 6)

########################################
### Figure 19-c  ctree + group lasso ###
########################################

### Fit the model ### 
fit <- mcre_agenet(X = X, Y = Y, W = W, W.hat = W.hat, gen.tree = "ctree", tune.init = "CV", tune.fin = "CV", seed = seed, nfolds = 5, learnrate = 0.001, ntrees = 333, depth = 3)
pred <- predict.mcre(fit, newdata = data.frame(Y = Y, W = as.factor(W), X))
tau.beta <- pred$tau.first

### ZDV vs ZDV + DID -----------------------------------------------------------
dat1 <- data.frame(Y = Y, W = W, tau = tau.beta[,1])
dat1 <- dat1[which(dat1$W == 0 | dat1$W == 1),]
dat1 <- dat1[order(dat1$tau),]
id1 <- 1:nrow(dat1)
split_vec1 <- split(id1, cut(seq_along(id1), 5, labels = FALSE))
ate1 <- c()
ate1.est <- c()
for(l in 1:length(split_vec1)){
  dat.sub1 <- dat1[split_vec1[[l]],]
  ate1.t <- mean(dat.sub1$Y[dat.sub1$W == 1]) - mean(dat.sub1$Y[dat.sub1$W == 0])
  ate1.est.t <- mean(dat.sub1$tau)
  ate1 <- c(ate1, ate1.t)
  ate1.est <- c(ate1.est, ate1.est.t)
}
ate1.dat <- data.frame(senario = "ZDV vs ZDV + DID", 
                       ate = c(ate1, ate1.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p1 <- clinfun::jonckheere.test(ate1.dat$ate[1:5], ate1.dat$ate[6:10])$p.value
cor1 <- cor(ate1.dat$ate[1:5], ate1.dat$ate[6:10], method = "spearman")

### ZDV vs ZDV + ZAL -----------------------------------------------------------
dat2 <- data.frame(Y = Y, W = W, tau = tau.beta[,2])
dat2 <- dat2[which(dat2$W == 0 | dat2$W == 2),]
dat2 <- dat2[order(dat2$tau),]
id2 <- 1:nrow(dat2)
split_vec2 <- split(id2, cut(seq_along(id2), 5, labels = FALSE))
ate2 <- c()
ate2.est <- c()
for(l in 1:length(split_vec2)){
  dat.sub2 <- dat2[split_vec2[[l]],]
  ate2.t <- mean(dat.sub2$Y[dat.sub2$W == 2]) - mean(dat.sub2$Y[dat.sub2$W == 0])
  ate2.est.t <- mean(dat.sub2$tau)
  ate2 <- c(ate2, ate2.t)
  ate2.est <- c(ate2.est, ate2.est.t)
}
ate2.dat <- data.frame(senario = "ZDV vs ZDV + ZAL", 
                       ate = c(ate2, ate2.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p2 <- clinfun::jonckheere.test(ate2.dat$ate[1:5], ate2.dat$ate[6:10])$p.value
cor2 <- cor(ate2.dat$ate[1:5], ate2.dat$ate[6:10], method = "spearman")

### ZDV vs DID -----------------------------------------------------------
dat3 <- data.frame(Y = Y, W = W, tau = tau.beta[,3])
dat3 <- dat3[which(dat3$W == 0 | dat3$W == 3),]
dat3 <- dat3[order(dat3$tau),]
id3 <- 1:nrow(dat3)
split_vec3 <- split(id3, cut(seq_along(id3), 5, labels = FALSE))
ate3 <- c()
ate3.est <- c()
for(l in 1:length(split_vec3)){
  dat.sub3 <- dat3[split_vec3[[l]],]
  ate3.t <- mean(dat.sub3$Y[dat.sub3$W == 3]) - mean(dat.sub3$Y[dat.sub3$W == 0])
  ate3.est.t <- mean(dat.sub3$tau)
  ate3 <- c(ate3, ate3.t)
  ate3.est <- c(ate3.est, ate3.est.t)
}
ate3.dat <- data.frame(senario = "ZDV vs DID", 
                       ate = c(ate3, ate3.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p3 <- clinfun::jonckheere.test(ate3.dat$ate[1:5], ate3.dat$ate[6:10])$p.value
cor3 <- cor(ate3.dat$ate[1:5], ate3.dat$ate[6:10], method = "spearman")

ate.dat <- rbind(ate1.dat, ate2.dat, ate3.dat)
ate.dat$senario <- factor(ate.dat$senario, levels = unique(ate.dat$senario))


library(ggplot2)
library(gridExtra)
g1 <- ggplot(ate.dat[1:10,], aes(x = sub, y = ate, fill = sub.c))
g1 <- g1 + geom_bar(stat = "identity", position = "dodge")
g1 <- g1 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g1 <- g1 + coord_cartesian(ylim = c(-100, 250))
g1 <- g1 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g1 <- g1 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

g2 <- ggplot(ate.dat[11:20,], aes(x = sub, y = ate, fill = sub.c))
g2 <- g2 + geom_bar(stat = "identity", position = "dodge")
g2 <- g2 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g2 <- g2 + coord_cartesian(ylim = c(-100, 250))
g2 <- g2 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g2 <- g2 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

g3 <- ggplot(ate.dat[21:30,], aes(x = sub, y = ate, fill = sub.c))
g3 <- g3 + geom_bar(stat = "identity", position = "dodge")
g3 <- g3 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g3 <- g3 + coord_cartesian(ylim = c(-100, 250))
g3 <- g3 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g3 <- g3 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(g1)

g <- grid.arrange(arrangeGrob(g1 + labs(title = "A) ZDV vs ZDV + DID"), 
                              g2 + labs(title = "B) ZDV vs ZDV + ZAL"), 
                              g3 + labs(title = "C) ZDV vs DID"), 
                              nrow=1))

ggsave(paste0(dir,"/Table_Figure_support/sup_section3/s-Fig19-c.eps"),g , width = 16, height = 6)

#################################################
### Figure 19-d  ctree + adaptive group lasso ###
#################################################

### Fit the model ### 
fit <- mcre_agenet(X = X, Y = Y, W = W, W.hat = W.hat, gen.tree = "ctree", tune.init = "CV", tune.fin = "CV", seed = seed, nfolds = 5, learnrate = 0.1, ntrees = 1000, depth = 3)
pred <- predict.mcre(fit, newdata = data.frame(Y = Y, W = as.factor(W), X))
tau.beta <- pred$tau

### ZDV vs ZDV + DID -----------------------------------------------------------
dat1 <- data.frame(Y = Y, W = W, tau = tau.beta[,1])
dat1 <- dat1[which(dat1$W == 0 | dat1$W == 1),]
dat1 <- dat1[order(dat1$tau),]
id1 <- 1:nrow(dat1)
split_vec1 <- split(id1, cut(seq_along(id1), 5, labels = FALSE))
ate1 <- c()
ate1.est <- c()
for(l in 1:length(split_vec1)){
  dat.sub1 <- dat1[split_vec1[[l]],]
  ate1.t <- mean(dat.sub1$Y[dat.sub1$W == 1]) - mean(dat.sub1$Y[dat.sub1$W == 0])
  ate1.est.t <- mean(dat.sub1$tau)
  ate1 <- c(ate1, ate1.t)
  ate1.est <- c(ate1.est, ate1.est.t)
}
ate1.dat <- data.frame(senario = "ZDV vs ZDV + DID", 
                       ate = c(ate1, ate1.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p1 <- clinfun::jonckheere.test(ate1.dat$ate[1:5], ate1.dat$ate[6:10])$p.value
cor1 <- cor(ate1.dat$ate[1:5], ate1.dat$ate[6:10], method = "spearman")

### ZDV vs ZDV + ZAL -----------------------------------------------------------
dat2 <- data.frame(Y = Y, W = W, tau = tau.beta[,2])
dat2 <- dat2[which(dat2$W == 0 | dat2$W == 2),]
dat2 <- dat2[order(dat2$tau),]
id2 <- 1:nrow(dat2)
split_vec2 <- split(id2, cut(seq_along(id2), 5, labels = FALSE))
ate2 <- c()
ate2.est <- c()
for(l in 1:length(split_vec2)){
  dat.sub2 <- dat2[split_vec2[[l]],]
  ate2.t <- mean(dat.sub2$Y[dat.sub2$W == 2]) - mean(dat.sub2$Y[dat.sub2$W == 0])
  ate2.est.t <- mean(dat.sub2$tau)
  ate2 <- c(ate2, ate2.t)
  ate2.est <- c(ate2.est, ate2.est.t)
}
ate2.dat <- data.frame(senario = "ZDV vs ZDV + ZAL", 
                       ate = c(ate2, ate2.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p2 <- clinfun::jonckheere.test(ate2.dat$ate[1:5], ate2.dat$ate[6:10])$p.value
cor2 <- cor(ate2.dat$ate[1:5], ate2.dat$ate[6:10], method = "spearman")

### ZDV vs DID -----------------------------------------------------------
dat3 <- data.frame(Y = Y, W = W, tau = tau.beta[,3])
dat3 <- dat3[which(dat3$W == 0 | dat3$W == 3),]
dat3 <- dat3[order(dat3$tau),]
id3 <- 1:nrow(dat3)
split_vec3 <- split(id3, cut(seq_along(id3), 5, labels = FALSE))
ate3 <- c()
ate3.est <- c()
for(l in 1:length(split_vec3)){
  dat.sub3 <- dat3[split_vec3[[l]],]
  ate3.t <- mean(dat.sub3$Y[dat.sub3$W == 3]) - mean(dat.sub3$Y[dat.sub3$W == 0])
  ate3.est.t <- mean(dat.sub3$tau)
  ate3 <- c(ate3, ate3.t)
  ate3.est <- c(ate3.est, ate3.est.t)
}
ate3.dat <- data.frame(senario = "ZDV vs DID", 
                       ate = c(ate3, ate3.est), 
                       sub = rep(c("S1", "S2", "S3", "S4", "S5"), 2),
                       sub.c = c("Actual", "Actual", "Actual", "Actual", "Actual", "Estimate", "Estimate", "Estimate", "Estimate", "Estimate")
)

p3 <- clinfun::jonckheere.test(ate3.dat$ate[1:5], ate3.dat$ate[6:10])$p.value
cor3 <- cor(ate3.dat$ate[1:5], ate3.dat$ate[6:10], method = "spearman")

ate.dat <- rbind(ate1.dat, ate2.dat, ate3.dat)
ate.dat$senario <- factor(ate.dat$senario, levels = unique(ate.dat$senario))


library(ggplot2)
library(gridExtra)
g1 <- ggplot(ate.dat[1:10,], aes(x = sub, y = ate, fill = sub.c))
g1 <- g1 + geom_bar(stat = "identity", position = "dodge")
g1 <- g1 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g1 <- g1 + coord_cartesian(ylim = c(-100, 250))
g1 <- g1 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g1 <- g1 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

g2 <- ggplot(ate.dat[11:20,], aes(x = sub, y = ate, fill = sub.c))
g2 <- g2 + geom_bar(stat = "identity", position = "dodge")
g2 <- g2 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g2 <- g2 + coord_cartesian(ylim = c(-100, 250))
g2 <- g2 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g2 <- g2 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

g3 <- ggplot(ate.dat[21:30,], aes(x = sub, y = ate, fill = sub.c))
g3 <- g3 + geom_bar(stat = "identity", position = "dodge")
g3 <- g3 + ylab("Heterogeneous Treatment Effect") + xlab("The subgroups ordered based on HTE (S1 - S5)")
g3 <- g3 + coord_cartesian(ylim = c(-100, 250))
g3 <- g3 + theme_bw() + theme(legend.justification = c(1.25,-0.5), legend.position = c(1,0), text = element_text(size = 15))
g3 <- g3 + geom_hline(yintercept = 0, color = "gray") + labs(fill = "")

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(g1)

g <- grid.arrange(arrangeGrob(g1 + labs(title = "A) ZDV vs ZDV + DID"), 
                              g2 + labs(title = "B) ZDV vs ZDV + ZAL"), 
                              g3 + labs(title = "C) ZDV vs DID"), 
                              nrow=1))

ggsave(paste0(dir,"/Table_Figure_support/sup_section3/s-Fig19-d.eps"),g , width = 16, height = 6)













