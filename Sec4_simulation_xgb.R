library(doParallel)
#########################

dir <- dir.save <- getwd()
##########################

source(paste0(dir, "/data_gen.R"))
source(paste0(dir, "/sim_method_xgb.R"))

### DGP parameters ###
n <- 2000
p <- 10
model.can <- expand.grid(1:3, 2:4)
model.can <- cbind(sim = 1:9, model.can)
Beta.can <- 2*rbind(c(1, 1, 1), c(-0.5, 1.0, 2.0), c(1.5, 1.5, -0.5), c(-1.5, 1.5, 0.5), c(-0.5, 2.0, 0.5))
prob.can <- c("rct", "obv")
group.can <- c(3,4,5)
ite.all <- 100

############################################################################################################################
model.all <- NULL
for(g in group.can){
  for(prob in prob.can){
    model.set <- cbind(model.can, prob, g) 
    model.all <- rbind(model.all, model.set)
  }
}

sim.all <- NULL

for(m in 1:nrow(model.all)){
  
  mod <- model.all[m,]
  
  sim_fun <- function(seed){
    
    ### Create the dataset ###
    sim.train <- data_gen(n = n, p = p, Beta = Beta.can[1:mod$g,], prob = mod$prob, model = mod[2:3], seed = seed)
    sim.test <- data_gen(n = n, p = p, Beta = Beta.can[1:mod$g,], prob = mod$prob, model = mod[2:3], seed = seed + 100)
    tau.test <- sim.test$tau
    
    ### Estimate the propensity score ###
    set.seed(seed)
    PS.fit <- nnet::multinom(W ~., data = data.frame(W = as.factor(sim.train$W), sim.train$X))
    W.hat <- predict(PS.fit, data = data.frame(sim.train$X), type = "probs")
    W.hat.test <- predict(PS.fit, data = data.frame(sim.test$X), type = "probs")
    
    #####################################a 
    t.tau.hat <- tboost(sim.train, sim.test, seed = seed)
    pred_t <- t.tau.hat

    s.tau.hat <- sboost(sim.train, sim.test, seed = seed)
    pred_s <- s.tau.hat
    
    x.tau.hat <- xboost(sim.train, sim.test, W.hat, seed = seed)
    pred_x <- x.tau.hat
    
    m.tau.hat <- mboost(sim.train, sim.test, W.hat, seed = seed)
    pred_m <-  m.tau.hat
    
    dr.tau.hat <- drboost(sim.train, sim.test, W.hat, seed = seed)
    pred_dr <- dr.tau.hat
    
    rfboost.tau.hat <- rfboost(sim.train, sim.test, W.hat, W.hat.test, seed = seed)
    pred_rf <- rfboost.tau.hat
    
    rboost.tau.hat <- rboost(sim.train, sim.test, W.hat, seed = seed)
    pred_r <- rboost.tau.hat
    
    ### TRUE & Estimated HTE ###
    
    res <- data.frame(tau = tau.test, 
                     sboost = pred_s$tau, 
                     tboost = pred_t$tau, 
                     xboost = pred_x$tau,
                     mboost = pred_m$tau,
                     drboost = pred_dr$tau,
                     rfboost = pred_rf$tau,
                     rboost = pred_r$tau)
    
    save.res <- paste0(dir.save, "/main_simulation/boost/prev_boost",mod$sim,"_n",n,"_p",p,"_",mod$prob,"_",mod$g,"g","_ite",seed,".csv")
    write.csv(res, save.res, row.names = FALSE)
    
    
    ## Evaluation Result -------------------------------------------------------
    mPEHE <- data.frame(sboost = mean(sqrt(apply((pred_s$tau - tau.test)^2, 2, mean))),
                        tboost = mean(sqrt(apply((pred_t$tau - tau.test)^2, 2, mean))),
                        xboost = mean(sqrt(apply((pred_x$tau - tau.test)^2, 2, mean))),
                        mboost = mean(sqrt(apply((pred_m$tau - tau.test)^2, 2, mean))),
                        drboost = mean(sqrt(apply((pred_dr$tau - tau.test)^2, 2, mean))),
                        rfboost = mean(sqrt(apply((pred_rf$tau - tau.test)^2, 2, mean))),
                        rboost = mean(sqrt(apply((pred_r$tau - tau.test)^2, 2, mean))))
    
    mbias <- data.frame(sboost = mean(abs(apply((pred_s$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        tboost = mean(abs(apply((pred_t$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        xboost = mean(abs(apply((pred_x$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        mboost = mean(abs(apply((pred_m$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        drboost = mean(abs(apply((pred_dr$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        rfboost = mean(abs(apply((pred_rf$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        rboost = mean(abs(apply((pred_r$tau - tau.test), 2, mean)/apply(tau.test,2,mean))))
    
    mcor <- data.frame(sboost = mean(diag(cor(pred_s$tau, tau.test, method = "spearman"))),
                       tboost = mean(diag(cor(pred_t$tau, tau.test, method = "spearman"))),
                       xboost = mean(diag(cor(pred_x$tau, tau.test, method = "spearman"))),
                       mboost = mean(diag(cor(pred_m$tau, tau.test, method = "spearman"))),
                       drboost = mean(diag(cor(pred_dr$tau, tau.test, method = "spearman"))),
                       rfboost = mean(diag(cor(pred_rf$tau, tau.test, method = "spearman"))),
                       rboost = mean(diag(cor(pred_r$tau, tau.test, method = "spearman"))))
    
    kappa <- data.frame(sboost = DescTools::CohenKappa(pred_s$rank, sim.test$rank),
                        tboost = DescTools::CohenKappa(pred_t$rank, sim.test$rank),
                        xboost = DescTools::CohenKappa(pred_x$rank, sim.test$rank),
                        mboost = DescTools::CohenKappa(pred_m$rank, sim.test$rank),
                        drboost = DescTools::CohenKappa(pred_dr$rank, sim.test$rank),
                        rfboost = DescTools::CohenKappa(pred_rf$rank, sim.test$rank),
                        rboost = DescTools::CohenKappa(pred_r$rank, sim.test$rank))
    
    mcor_ord <- data.frame(sboost = mean(diag(cor(t(pred_s$order), t(sim.test$order), method = "spearman"))),
                           tboost = mean(diag(cor(t(pred_t$order), t(sim.test$order), method = "spearman"))),
                           xboost = mean(diag(cor(t(pred_x$order), t(sim.test$order), method = "spearman"))),
                           mboost = mean(diag(cor(t(pred_m$order), t(sim.test$order), method = "spearman"))),
                           drboost = mean(diag(cor(t(pred_dr$order), t(sim.test$order), method = "spearman"))),
                           rfboost = mean(diag(cor(t(pred_rf$order), t(sim.test$order), method = "spearman"))),
                           rboost = mean(diag(cor(t(pred_r$order), t(sim.test$order), method = "spearman"))))
    
    time <- data.frame(sboost = pred_s$time, 
                      tboost = pred_t$time, 
                      xboost = pred_x$time,
                      mboost = pred_m$time,
                      drboost = pred_dr$time,
                      rfboost = pred_rf$time,
                      rboost = pred_r$time)
    
    
    sim.evl <- rbind(mPEHE = mPEHE, mbias = mbias, mcor = mcor, kappa = kappa, mcor_ord = mcor_ord, time = time)
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
  
  save.res <- paste0(dir.save, "/main_simulation/boost/boost_evl_sim",mod$sim,"_n",n,"_p",p,"_",mod$prob,"_",mod$g,"g",".csv")
  write.csv(sim.res, save.res, row.names = FALSE)
}
