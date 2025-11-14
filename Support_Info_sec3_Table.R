library(speff2trial)
library(doParallel)
#########################

data(ACTG175)
dat <- ACTG175

###############
dir <- getwd()
source(paste0(dir, "/mcre_agenet.R"))

### Set seed ------------------------------------------------

sub.num <- 5 # divide into five subgroups

### Real data -----------------------------------------------

set.seed(1111)

Y <- as.numeric(dat$cd420 - dat$cd40)
W <- dat$arm
X <- dat[, c("age", "wtkg", "hemo", "homo", "drugs", "karnof", "race", "gender", "str2", "symptom",  "cd40", "cd80")]
dat <- data.frame(Y = Y, W = W, X)
dat <- dat[dat$cd40>= 200 & dat$cd40 <= 500,]

### Parameter candidate -------------------------------------

trees.can <- c(333, 666, 1000)
depth.can <- c(2, 3, 4)
learn.can <- c(0.1, 0.01, 0.001)

para.all <- NULL
for(t in trees.can){
  for(d in depth.can){
    for(l in learn.can){
      para.set <- c(t, d, l) 
      para.all <- rbind(para.all, para.set)
    }  
  }
}
rownames(para.all) <- 1:nrow(para.all)

### Parameter tuning ------------------------------------

mtr.all <- NULL

tune_fun <- function(seed){

  for(para in 1:nrow(para.all)){
  
    dat.t  <- dat.te <- dat
    
    ### Estimate the propensity score ----------------------------------------------
    PS.RF <-  nnet::multinom(W~., data = data.frame(W = as.factor(dat.t$W), dat.t[,-c(1:2)]))
    W.hat <- predict(PS.RF, data =  dat.t[,-c(1:2)], type = "prob")
    
    ### Fit the model CART -------------------------------------------------------
    fit0 <- mcre_agenet(X = dat.t[,-c(1:2)], Y = dat.t$Y, W = dat.t$W, W.hat = W.hat, gen.tree = "gbm", tune.init = "CV", tune.fin = "CV", nfolds = 5,
                        ntrees = para.all[para, 1], depth = para.all[para, 2], learnrate = para.all[para, 3], 
                        seed = seed)
    pred0 <- predict.mcre(fit0, newdata = data.frame(Y = dat.te$Y, W = as.factor(dat.te$W), dat.te[,-c(1:2)]))
    tau.gbm.gl <- pred0$tau.first
    tau.gbm.agl <- pred0$tau
    
    ### Fit the model ctree  -----------------------------------------------------
    fit1 <- mcre_agenet(X = dat.t[,-c(1:2)], Y = dat.t$Y, W = dat.t$W, W.hat = W.hat, gen.tree = "ctree", tune.init = "CV", tune.fin = "CV", nfolds = 5,
                        ntrees = para.all[para, 1], depth = para.all[para, 2], learnrate = para.all[para, 3], 
                        seed = seed)
    pred1 <- predict.mcre(fit1, newdata = data.frame(Y = dat.te$Y, W = as.factor(dat.te$W), dat.te[,-c(1:2)]))
    tau.ctree.gl <- pred1$tau.first
    tau.ctree.agl <- pred1$tau
    
    mtr1.can <- mtr2.can <- mtr3.can <- mtr4.can <- c()
    
    for(ind in unique(sort(dat.te$W))[-1]){
  
      tau.est <- cbind(gbm.gl = tau.gbm.gl[,ind], gbm.agl = tau.gbm.agl[,ind], ctree.gl = tau.ctree.gl[,ind], ctree.agl = tau.ctree.agl[,ind])
      dat.evl <- data.frame(Y = dat.te$Y, W = dat.te$W, tau = tau.est)
      dat.evl <- dat.evl[which(dat.te$W == 0 | dat.te$W == ind),]
      
      dat.evl1 <- dat.evl[order(dat.evl$tau.gbm.gl),]
      dat.evl2 <- dat.evl[order(dat.evl$tau.gbm.agl),]
      dat.evl3 <- dat.evl[order(dat.evl$tau.ctree.gl),]
      dat.evl4 <- dat.evl[order(dat.evl$tau.ctree.agl),]
      
      split_vec <- split(1:nrow(dat.evl), cut(seq_along(1:nrow(dat.evl)), sub.num, labels = FALSE))
      
      ate1.can <- ate2.can <- ate3.can <- ate4.can <- c()
      
      for(l in 1:length(split_vec)){
        
        dat.sub1 <- dat.evl1[split_vec[[l]],]
        dat.sub2 <- dat.evl2[split_vec[[l]],]
        dat.sub3 <- dat.evl3[split_vec[[l]],]
        dat.sub4 <- dat.evl4[split_vec[[l]],]
        
        ate1.t <- mean(dat.sub1$Y[dat.sub1$W == ind]) - mean(dat.sub1$Y[dat.sub1$W == 0])
        ate2.t <- mean(dat.sub2$Y[dat.sub2$W == ind]) - mean(dat.sub2$Y[dat.sub2$W == 0])
        ate3.t <- mean(dat.sub3$Y[dat.sub3$W == ind]) - mean(dat.sub3$Y[dat.sub3$W == 0])
        ate4.t <- mean(dat.sub4$Y[dat.sub4$W == ind]) - mean(dat.sub4$Y[dat.sub4$W == 0])
        
        ate1.est.t <- mean(dat.sub1$tau.gbm.gl)
        ate2.est.t <- mean(dat.sub2$tau.gbm.agl)
        ate3.est.t <- mean(dat.sub3$tau.ctree.gl)
        ate4.est.t <- mean(dat.sub4$tau.ctree.agl)
        
        ate1.can <- rbind(ate1.can, c(ate1.t, ate1.est.t))
        ate2.can <- rbind(ate2.can, c(ate2.t, ate2.est.t))
        ate3.can <- rbind(ate3.can, c(ate3.t, ate3.est.t))
        ate4.can <- rbind(ate4.can, c(ate4.t, ate4.est.t))
        
      }
      
      mtr1 <-  abs(mean(abs(ate1.can[,2] - ate1.can[,1]))/(cor(ate1.can[,2], ate1.can[,1])*mean(sign(ate1.can[,2]) ==  sign(ate1.can[,1])))) 
      mtr2 <-  abs(mean(abs(ate2.can[,2] - ate2.can[,1]))/(cor(ate2.can[,2], ate2.can[,1])*mean(sign(ate2.can[,2]) ==  sign(ate2.can[,1]))))  
      mtr3 <-  abs(mean(abs(ate3.can[,2] - ate3.can[,1]))/(cor(ate3.can[,2], ate3.can[,1])*mean(sign(ate3.can[,2]) ==  sign(ate3.can[,1])))) 
      mtr4 <-  abs(mean(abs(ate4.can[,2] - ate4.can[,1]))/(cor(ate4.can[,2], ate4.can[,1])*mean(sign(ate4.can[,2]) ==  sign(ate4.can[,1])))) 
      
      mtr1.can <- c(mtr1.can, mtr1)
      mtr2.can <- c(mtr2.can, mtr2)
      mtr3.can <- c(mtr3.can, mtr3)
      mtr4.can <- c(mtr4.can, mtr4)
      
    }    
    
    mtr.res <- c(para.all[para,], mean(mtr1.can), mean(mtr2.can), mean(mtr3.can), mean(mtr4.can))
    mtr.all <- rbind(mtr.all, mtr.res)
  }    
  
  return(mtr.all)
}

### Parallel calculation ------------------------------------------------------
cl <- makeCluster(detectCores()-7, outfile="")
registerDoParallel(cl)
tune.res <- foreach(seed = 1:10,.combine=rbind,.export=ls(envir = parent.frame())) %dopar% {tune_fun(seed)}
stopCluster(cl)

### Summary the results --------------------------------------------------------
res1 <- 0
id.can <- nrow(para.all)*c(0:9) + 1
for(r in 1:10){
  tune_ite <- tune.res[id.can[r]:(id.can[r] + nrow(para.all)- 1),]
  save.m <- paste0(dir,"/add_simulation/Section3_real_data_parameter_tune/tune_ite",r,".csv")
  write.csv(tune_ite[,4:7], save.m, row.names = FALSE)
  
  res0 <- read.csv(save.m)
  para.n <- tune_ite[,1:3]
  res1 <- res1 + res0
}
res1 <- round(res1/10,2)
res1 <- cbind(para.n, res1)

colnames(res1) <- c("Number of trees", "Depth of tree", "Learner rate", "gbm.gl", "gbm.agl", "ctree.gl", "ctree.agl")
write.csv(res1, paste0(dir,"/Table_Figure_support/sup_section3/Table9.csv"), row.names = FALSE)




  