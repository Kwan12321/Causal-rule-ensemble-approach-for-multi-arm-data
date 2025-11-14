library(doParallel)
#########################

dir <- dir.save <- getwd()

##########################

source(paste0(dir, "/data_gen.R"))
source(paste0(dir, "/sim_method_bart.R"))

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

for(m in 42:nrow(model.all)){
  
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
    t.tau.hat <- tbart(sim.train, sim.test, seed = seed)
    pred_t <- t.tau.hat
                                                                                                                                                                                                                                                                                                                                                         
    s.tau.hat <- sbart(sim.train, sim.test, seed = seed)
    pred_s <- s.tau.hat
    
    x.tau.hat <- xbart(sim.train, sim.test, W.hat, seed = seed)
    pred_x <- x.tau.hat
    
    m.tau.hat <- mbart(sim.train, sim.test, W.hat, seed = seed)
    pred_m <-  m.tau.hat
    
    dr.tau.hat <- drbart(sim.train, sim.test, W.hat, seed = seed)
    pred_dr <- dr.tau.hat
    
    rfbart.tau.hat <- rfbart(sim.train, sim.test, W.hat, W.hat.test, seed = seed)
    pred_rf <- rfbart.tau.hat
    
    rbart.tau.hat <- rbart(sim.train, sim.test, W.hat, seed = seed)
    pred_r <- rbart.tau.hat
    
    ### TRUE & Estimated HTE ###
    res <- data.frame(tau = tau.test, 
                     sbart = pred_s$tau, 
                     tbart = pred_t$tau, 
                     xbart = pred_x$tau,
                     mbart = pred_m$tau,
                     drbart = pred_dr$tau,
                     rfbart = pred_rf$tau,
                     rbart = pred_r$tau)
    
    save.res <- paste0(dir.save, "/main_simulation/bart/prev_bart",mod$sim,"_n",n,"_p",p,"_",mod$prob,"_",mod$g,"g","_ite",seed,".csv")
    write.csv(res, save.res, row.names = FALSE)
    
    ## Evaluation Result -------------------------------------------------------
    
    mPEHE <- data.frame(sbart = mean(sqrt(apply((pred_s$tau - tau.test)^2, 2, mean))),
                        tbart = mean(sqrt(apply((pred_t$tau - tau.test)^2, 2, mean))),
                        xbart = mean(sqrt(apply((pred_x$tau - tau.test)^2, 2, mean))),
                        mbart = mean(sqrt(apply((pred_m$tau - tau.test)^2, 2, mean))),
                        drbart = mean(sqrt(apply((pred_dr$tau - tau.test)^2, 2, mean))),
                        rfbart = mean(sqrt(apply((pred_rf$tau - tau.test)^2, 2, mean))),
                        rbart = mean(sqrt(apply((pred_r$tau - tau.test)^2, 2, mean))))
    
    
    mbias <- data.frame(sbart = mean(abs(apply((pred_s$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        tbart = mean(abs(apply((pred_t$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        xbart = mean(abs(apply((pred_x$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        mbart = mean(abs(apply((pred_m$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        drbart = mean(abs(apply((pred_dr$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        rfbart = mean(abs(apply((pred_rf$tau - tau.test), 2, mean)/apply(tau.test,2,mean))),
                        rbart = mean(abs(apply((pred_r$tau - tau.test), 2, mean)/apply(tau.test,2,mean))))
    
    mcor <- data.frame(sbart = mean(diag(cor(pred_s$tau, tau.test, method = "spearman"))),
                       tbart = mean(diag(cor(pred_t$tau, tau.test, method = "spearman"))),
                       xbart = mean(diag(cor(pred_x$tau, tau.test, method = "spearman"))),
                       mbart = mean(diag(cor(pred_m$tau, tau.test, method = "spearman"))),
                       drbart = mean(diag(cor(pred_dr$tau, tau.test, method = "spearman"))),
                       rfbart = mean(diag(cor(pred_rf$tau, tau.test, method = "spearman"))),
                       rbart = mean(diag(cor(pred_r$tau, tau.test, method = "spearman"))))
    
    kappa <- data.frame(sbart = DescTools::CohenKappa(pred_s$rank, sim.test$rank),
                        tbart = DescTools::CohenKappa(pred_t$rank, sim.test$rank),
                        xbart = DescTools::CohenKappa(pred_x$rank, sim.test$rank),
                        mbart = DescTools::CohenKappa(pred_m$rank, sim.test$rank),
                        drbart = DescTools::CohenKappa(pred_dr$rank, sim.test$rank),
                        rfbart = DescTools::CohenKappa(pred_rf$rank, sim.test$rank),
                        rbart = DescTools::CohenKappa(pred_r$rank, sim.test$rank))
    
    mcor_ord <- data.frame(sbart = mean(diag(cor(t(pred_s$order), t(sim.test$order), method = "spearman"))),
                           tbart = mean(diag(cor(t(pred_t$order), t(sim.test$order), method = "spearman"))),
                           xbart = mean(diag(cor(t(pred_x$order), t(sim.test$order), method = "spearman"))),
                           mbart = mean(diag(cor(t(pred_m$order), t(sim.test$order), method = "spearman"))),
                           drbart = mean(diag(cor(t(pred_dr$order), t(sim.test$order), method = "spearman"))),
                           rfbart = mean(diag(cor(t(pred_rf$order), t(sim.test$order), method = "spearman"))),
                           rbart = mean(diag(cor(t(pred_r$order), t(sim.test$order), method = "spearman"))))
    
    time <- data.frame(sbart = pred_s$time, 
                      tbart = pred_t$time, 
                      xbart = pred_x$time,
                      mbart = pred_m$time,
                      drbart = pred_dr$time,
                      rfbart = pred_rf$time,
                      rbart = pred_r$time)
    
    
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
  
  save.res <- paste0(dir.save, "/main_simulation/bart/bart_evl_sim",mod$sim,"_n",n,"_p",p,"_",mod$prob,"_",mod$g,"g",".csv")
  write.csv(sim.res, save.res, row.names = FALSE)
}
