library(doParallel)
#########################

dir <- getwd()
##########################

source(paste0(dir, "/data_gen.R"))
source(paste0(dir, "/mcre_agenet.R"))

### DGP parameters ###
n <- 2000
p <- 10
model.can <- expand.grid(1:3, 2:4)
model.can <- cbind(sim = 1:9, model.can)
Beta.can <- 2*rbind(c(1, 1, 1), c(-0.5, 1.0, 2.0), c(1.5, 1.5, -0.5), c(-1.5, 1.5, 0.5), c(-0.5, 2.0, 0.5))
group.can <- c(3,4,5)
ite.all <- 100

############################################################################################################################
model.all <- NULL
for(g in group.can){
  model.set <- cbind(model.can, g) 
  model.all <- rbind(model.all, model.set)
}

sim.all <- NULL

for(m in 1:nrow(model.all)){
  
  mod <- model.all[m,]
  
  sim_fun <- function(seed){
    
    set.seed(seed)
    
    ### Create the dataset ###
    sim.train <- data_gen(n = n, p = p, Beta = Beta.can[1:mod$g,], model = mod[2:3], seed = seed)
    sim.test <- data_gen(n = n, p = p, Beta = Beta.can[1:mod$g,], model = mod[2:3], seed = seed + 100)
    tau.test <- sim.test$tau
  
    PS.RF <- nnet::multinom(W~., data = data.frame(W = as.factor(sim.train$W), sim.train$X))
    W.hat <- predict(PS.RF, newdata = data.frame(sim.train$X), type = "prob")
    
    #####################
    ### Causal Forest ###
    #####################
    
    mcf <- grf::multi_arm_causal_forest(X = sim.train$X, W = as.factor(sim.train$W), Y = sim.train$Y, W.hat = W.hat)
    c.pred <- predict(mcf, sim.test$X, estimate.variance = TRUE)
    
    cate <- data.frame(c.pred$predictions)
    se   <- sqrt(data.frame(c.pred$variance.estimates))

    bound.cf <- list()
    cf.cover <- c()
    cf.width <- c()
    for(c in 1:(mod$g-1)){
      ci_lower <- cate[,c] - 1.96 * se[,c]
      ci_upper <- cate[,c] + 1.96 * se[,c]
      bound.cf[[c]] <- cbind(ci_lower, ci_upper)
      cf.cover[c] <- mean(tau.test[,c] >= ci_lower & tau.test[,c] <= ci_upper, na.rm = TRUE) 
      cf.width[c] <- mean(ci_upper - ci_lower)
    }
    
    cf.cover.res <- mean(cf.cover)
    cf.width.res <- mean(cf.width)
    
    ####################################################
    ### Proposed Method (GBM based Rule generation ) ###
    ####################################################
    
    trainprop <- 0.75
    alpha <- 0.05
    
    Y <- sim.train$Y
    W <- sim.train$W
    X <- sim.train$X
    
    ### Step 1 : Divide the data into train set and valid set (Balanced split)###
    
    trainid <- c()
    inds.train <- list()
    for(g in 0:(mod$g-1)){
      inds <- which(W == g); n <- length(inds)
      trainid_inds <- sample(n, floor(n * trainprop))
      trainid <- c(trainid, inds[trainid_inds])
      inds.train[[g+1]] <- inds
    }
    
    Xtrain <- X[trainid, ,drop=FALSE]; Wtrain <- W[trainid]; Ytrain <- Y[trainid]
    Xval <- X[-trainid, ,drop=FALSE]; Wval <- W[-trainid]; Yval <- Y[-trainid]
    
    ### Step 2 : Build model using train set ###
    
    ### Estimate the propensity score ###
    
    PS.RF <- nnet::multinom(W~., data = data.frame(W = as.factor(Wtrain), Xtrain))
    W.hat <- predict(PS.RF, newdata = data.frame(Xtrain), type = "prob")
    ### Fit the model ### 
    fit <- mcre_agenet(X = Xtrain, Y = Y[trainid], W = Wtrain, W.hat = W.hat, gen.tree = "gbm", tune.init = "CV", tune.fin = "CV", learnrate = 0.01, seed = seed)
    
    ### Step 3 : Compute the nonconformity measures ###
    pred <- predict.mcre(fit, newdata = data.frame(Y = Yval, W = as.factor(Wval), Xval))
    Wval_dum <- fastDummies::dummy_cols(Wval)[,-1]
    E.first <- abs(Yval - pred$y.first*Wval_dum); E <- abs(Yval - pred$y*Wval_dum)
    
    ### Step 4 : Compute the conformal interval for treatment estimate
    confint.first <- list()
    confint <- list()
    pred.all <- predict.mcre(fit, newdata = data.frame(sim.test$X))
    for(g in 1:length(unique(Wval))){
      confint.first[[g]] <- cbind(pred.all$y.first[,g] - quantile(E.first[Wval == (g-1),g], alpha/2), pred.all$y.first[,g] + quantile(E.first[Wval == (g-1),g], (1-alpha/2)))
      confint[[g]] <- cbind(pred.all$y[,g] - quantile(E[Wval == (g-1),g], alpha/2), pred.all$y[,g] + quantile(E[Wval == (g-1),g], (1-alpha/2)))
    }
    
    ### Step 5: Compute the conformal interval for estimated HTE
    confint_hte.first <-list()
    confint_hte <- list()
    pred.all <- predict.mcre(fit, newdata = data.frame(sim.test$X))
    for(g in 2:length(unique(Wval))){
      confint_hte.first[[g-1]] <- cbind(confint.first[[g]][,1] - confint.first[[1]][,2], confint.first[[g]][,2] - confint.first[[1]][,1])
      confint_hte[[g-1]] <- cbind(confint[[g]][,1] - confint[[1]][,2], confint[[g]][,2] - confint[[1]][,1])
    }  
    
    prop.cover <- c(); prop.cover.first <- c()
    prop.width <- c(); prop.width.first <- c()
    for(g in 1:(mod$g-1)){
      ci_lower.first <- confint_hte.first[[g]][,1];  ci_lower <- confint_hte[[g]][,1]
      ci_upper.first <- confint_hte.first[[g]][,2];  ci_upper <- confint_hte[[g]][,2]
      prop.cover.first[g] <- mean(tau.test[,g] >= ci_lower.first & tau.test[,g] <= ci_upper.first, na.rm = TRUE) 
      prop.cover[g] <- mean(tau.test[,g] >= ci_lower & tau.test[,g] <= ci_upper, na.rm = TRUE) 
      prop.width.first[g] <- mean(ci_upper.first - ci_lower.first)
      prop.width[g] <- mean(ci_upper - ci_lower)
    }
    
    prop.gbm.cover.res <- c(mean(prop.cover.first), mean(prop.cover)) 
    prop.gbm.width.res <- c(mean(prop.width.first), mean(prop.width)) 
    
    ######################################################
    ### Proposed Method (ctree based Rule generation ) ###
    ######################################################
    
    trainprop <- 0.75
    alpha <- 0.05
    
    Y <- sim.train$Y
    W <- sim.train$W
    X <- sim.train$X
    
    ### Step 1 : Divide the data into train set and valid set (Balanced split)###
    
    trainid <- c()
    inds.train <- list()
    for(g in 0:(mod$g-1)){
      inds <- which(W == g); n <- length(inds)
      trainid_inds <- sample(n, floor(n * trainprop))
      trainid <- c(trainid, inds[trainid_inds])
      inds.train[[g+1]] <- inds
    }
    
    Xtrain <- X[trainid, ,drop=FALSE]; Wtrain <- W[trainid]; Ytrain <- Y[trainid]
    Xval <- X[-trainid, ,drop=FALSE]; Wval <- W[-trainid]; Yval <- Y[-trainid]
    
    ### Step 2 : Build model using train set ###
    
    ### Estimate the propensity score ###
    
    PS.RF <- nnet::multinom(W~., data = data.frame(W = as.factor(Wtrain), Xtrain))
    W.hat <- predict(PS.RF, newdata = data.frame(Xtrain), type = "prob")
    ### Fit the model ### 
    fit <- mcre_agenet(X = Xtrain, Y = Y[trainid], W = Wtrain, W.hat = W.hat, gen.tree = "ctree", tune.init = "CV", tune.fin = "CV", learnrate = 0.01, seed = seed)
    
    ### Step 3 : Compute the nonconformity measures ###
    pred <- predict.mcre(fit, newdata = data.frame(Y = Yval, W = as.factor(Wval), Xval))
    Wval_dum <- fastDummies::dummy_cols(Wval)[,-1]
    E.first <- abs(Yval - pred$y.first*Wval_dum); E <- abs(Yval - pred$y*Wval_dum)
    
    ### Step 4 : Compute the conformal interval for treatment estimate
    confint.first <- list()
    confint <- list()
    pred.all <- predict.mcre(fit, newdata = data.frame(sim.test$X))
    for(g in 1:length(unique(Wval))){
      confint.first[[g]] <- cbind(pred.all$y.first[,g] - quantile(E.first[Wval == (g-1),g], alpha/2), pred.all$y.first[,g] + quantile(E.first[Wval == (g-1),g], (1-alpha/2)))
      confint[[g]] <- cbind(pred.all$y[,g] - quantile(E[Wval == (g-1),g], alpha/2), pred.all$y[,g] + quantile(E[Wval == (g-1),g], (1-alpha/2)))
    }
    
    ### Step 5: Compute the conformal interval for estimated HTE
    confint_hte.first <-list()
    confint_hte <- list()
    pred.all <- predict.mcre(fit, newdata = data.frame(sim.test$X))
    for(g in 2:length(unique(Wval))){
      confint_hte.first[[g-1]] <- cbind(confint.first[[g]][,1] - confint.first[[1]][,2], confint.first[[g]][,2] - confint.first[[1]][,1])
      confint_hte[[g-1]] <- cbind(confint[[g]][,1] - confint[[1]][,2], confint[[g]][,2] - confint[[1]][,1])
    }  
    
    prop.cover <- c(); prop.cover.first <- c()
    prop.width <- c(); prop.width.first <- c()
    for(g in 1:(mod$g-1)){
      ci_lower.first <- confint_hte.first[[g]][,1];  ci_lower <- confint_hte[[g]][,1]
      ci_upper.first <- confint_hte.first[[g]][,2];  ci_upper <- confint_hte[[g]][,2]
      prop.cover.first[g] <- mean(tau.test[,g] >= ci_lower.first & tau.test[,g] <= ci_upper.first, na.rm = TRUE) 
      prop.cover[g] <- mean(tau.test[,g] >= ci_lower & tau.test[,g] <= ci_upper, na.rm = TRUE) 
      prop.width.first[g] <- mean(ci_upper.first - ci_lower.first)
      prop.width[g] <- mean(ci_upper - ci_lower)
    }
    
    prop.ctree.cover.res <- c(mean(prop.cover.first), mean(prop.cover)) 
    prop.ctree.width.res <- c(mean(prop.width.first), mean(prop.width)) 
    
    cover.res <- data.frame(mcf = cf.cover.res, gbm.gl = prop.gbm.cover.res[1], gbm.agl = prop.gbm.cover.res[2],
                            ctree.gl = prop.ctree.cover.res[1], ctree.agl = prop.ctree.cover.res[2], evl = "cover_r")
    
    width.res <- data.frame(mcf = cf.width.res, gbm.gl = prop.gbm.width.res[1], gbm.agl = prop.gbm.width.res[2],
                            ctree.gl = prop.ctree.width.res[1], ctree.agl = prop.ctree.width.res[2], evl = "width")
    
    sim.evl <- rbind(cover.res, width.res)
    sim.evl$ite <- seed
    
    return(sim.evl)
  }
  
  
  ### Parallel ###
  cl <- makeCluster(detectCores()-7, outfile="")
  registerDoParallel(cl)
  sim.res <- foreach(seed = 1:ite.all,.combine=rbind,.export=ls(envir = parent.frame())) %dopar% {sim_fun(seed)}
  stopCluster(cl)
  
  sim.res$fun1 <- mod$Var1
  sim.res$fun2 <- mod$Var2
  sim.res$g <- mod$g
  
  save.res <- paste0(dir, "/add_simulation/Section4_conformal_prediction/cf_sim",mod$sim,"_n",n,"_p",p,"_",mod$prob,"_",mod$g,"g",".csv")
  write.csv(sim.res, save.res, row.names = FALSE)
}  

   