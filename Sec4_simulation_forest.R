library(doParallel)
#########################

dir <- dir.save <- getwd()
##########################

source(paste0(dir, "/data_gen.R"))
source(paste0(dir, "/sim_method_forest.R"))

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
    
    #####################################
    t.tau.hat <- tforest(sim.train, sim.test, seed = seed)
    pred_t <- t.tau.hat

    s.tau.hat <- sforest(sim.train, sim.test, seed = seed)
    pred_s <- s.tau.hat
    
    x.tau.hat <- xforest(sim.train, sim.test, W.hat, seed = seed)
    pred_x <- x.tau.hat
    
    m.tau.hat <- mforest(sim.train, sim.test, W.hat, seed = seed)
    pred_m <-  m.tau.hat
    
    dr.tau.hat <- drforest(sim.train, sim.test, W.hat, seed = seed)
    pred_dr <- dr.tau.hat
    
    rfforest.tau.hat <- rfforest(sim.train, sim.test, W.hat, W.hat.test, seed = seed)
    pred_rf <- rfforest.tau.hat
    
    rforest.tau.hat <- rforest(sim.train, sim.test, W.hat, seed = seed)
    pred_r <- rforest.tau.hat
    
    ### TRUE & Estimated HTE ###
    res <- data.frame(tau = tau.test, 
                     sforest = pred_s$tau, 
                     tforest = pred_t$tau, 
                     xforest = pred_x$tau,
                     mforest = pred_m$tau,
                     drforest = pred_dr$tau,
                     rfforest = pred_rf$tau,
                     rforest = pred_r$tau)
    
    save.res <- paste0(dir.save, "/main_simulation/forest/prev_forest",mod$sim,"_n",n,"_p",p,"_",mod$prob,"_",mod$g,"g","_ite",seed,".csv")
    write.csv(res, save.res, row.names = FALSE)
    
    
    ## Evaluation Result -------------------------------------------------------
    mPEHE <- data.frame(sforest = mean(sqrt(apply((pred_s$tau - tau.test)^2, 2, mean))),
                        tforest = mean(sqrt(apply((pred_t$tau - tau.test)^2, 2, mean))),
                        xforest = mean(sqrt(apply((pred_x$tau - tau.test)^2, 2, mean))),
                        mforest = mean(sqrt(apply((pred_m$tau - tau.test)^2, 2, mean))),
                        drforest = mean(sqrt(apply((pred_dr$tau - tau.test)^2, 2, mean))),
                        rfforest = mean(sqrt(apply((pred_rf$tau - tau.test)^2, 2, mean))),
                        rforest = mean(sqrt(apply((pred_r$tau - tau.test)^2, 2, mean))))
    
    
    mbias <- data.frame(sforest = mean(abs(apply((pred_s$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        tforest = mean(abs(apply((pred_t$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        xforest = mean(abs(apply((pred_x$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        mforest = mean(abs(apply((pred_m$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        drforest = mean(abs(apply((pred_dr$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        rfforest = mean(abs(apply((pred_rf$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        rforest = mean(abs(apply((pred_r$tau - tau.test), 2, mean)/apply(tau.test,2,mean))))
    
    mcor <- data.frame(sforest = mean(diag(cor(pred_s$tau, tau.test, method = "spearman"))),
                       tforest = mean(diag(cor(pred_t$tau, tau.test, method = "spearman"))),
                       xforest = mean(diag(cor(pred_x$tau, tau.test, method = "spearman"))),
                       mforest = mean(diag(cor(pred_m$tau, tau.test, method = "spearman"))),
                       drforest = mean(diag(cor(pred_dr$tau, tau.test, method = "spearman"))),
                       rfforest = mean(diag(cor(pred_rf$tau, tau.test, method = "spearman"))),
                       rforest = mean(diag(cor(pred_r$tau, tau.test, method = "spearman"))))
    
    kappa <- data.frame(sforest = DescTools::CohenKappa(pred_s$rank, sim.test$rank),
                        tforest = DescTools::CohenKappa(pred_t$rank, sim.test$rank),
                        xforest = DescTools::CohenKappa(pred_x$rank, sim.test$rank),
                        mforest = DescTools::CohenKappa(pred_m$rank, sim.test$rank),
                        drforest = DescTools::CohenKappa(pred_dr$rank, sim.test$rank),
                        rfforest = DescTools::CohenKappa(pred_rf$rank, sim.test$rank),
                        rforest = DescTools::CohenKappa(pred_r$rank, sim.test$rank))
    
    mcor_ord <- data.frame(sforest = mean(diag(cor(t(pred_s$order), t(sim.test$order), method = "spearman"))),
                           tforest = mean(diag(cor(t(pred_t$order), t(sim.test$order), method = "spearman"))),
                           xforest = mean(diag(cor(t(pred_x$order), t(sim.test$order), method = "spearman"))),
                           mforest = mean(diag(cor(t(pred_m$order), t(sim.test$order), method = "spearman"))),
                           drforest = mean(diag(cor(t(pred_dr$order), t(sim.test$order), method = "spearman"))),
                           rfforest = mean(diag(cor(t(pred_rf$order), t(sim.test$order), method = "spearman"))),
                           rforest = mean(diag(cor(t(pred_r$order), t(sim.test$order), method = "spearman"))))
    
    time <- data.frame(sforest = pred_s$time, 
                      tforest = pred_t$time, 
                      xforest = pred_x$time,
                      mforest = pred_m$time,
                      drforest = pred_dr$time,
                      rfforest = pred_rf$time,
                      rforest = pred_r$time)
    
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
  
  save.res <- paste0(dir.save, "/main_simulation/forest/forest_evl_sim",mod$sim,"_n",n,"_p",p,"_",mod$prob,"_",mod$g,"g",".csv")
  write.csv(sim.res, save.res, row.names = FALSE)
}
