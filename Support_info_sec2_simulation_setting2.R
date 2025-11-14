library(doParallel)
#########################

dir <- dir.save <- getwd()
##########################

source(paste0(dir, "/data_gen_overlap2.R"))
source(paste0(dir, "/mcre_agenet.R"))
source(paste0(dir, "/sim_method_proposed.R"))

### DGP parameters ###
n <- 2000
p <- 10
model.can <- expand.grid(1:3, 2:4)
model.can <- cbind(sim = 1:9, model.can)
Beta.can <- 2*rbind(c(1, 1, 1), c(-0.5, 1.0, 2.0), c(1.5, 1.5, -0.5), c(-1.5, 1.5, 0.5), c(-0.5, 2.0, 0.5))
overlap.can <- c("strong", "moderate","weak")
group.can <- c(3,4,5)
ite.all <- 100

############################################################################################################################
model.all <- NULL
for(g in group.can){
  for(prob in overlap.can){
    model.set <- cbind(model.can, prob, g) 
    model.all <- rbind(model.all, model.set)
  }
}

sim.all <- NULL

for(m in 1:nrow(model.all)){
#for(m in 19:27){  
  
  mod <- model.all[m,]
  
  sim_fun <- function(seed){
    
    ### Create the dataset ###
    sim.train <- data_gen_overlap2(n = n, p = p, Beta = Beta.can[1:mod$g,], model = mod[2:3], overlap = mod$prob, seed = seed)
    sim.test <- data_gen_overlap2(n = n, p = p, Beta = Beta.can[1:mod$g,], model = mod[2:3], overlap = mod$prob, seed = seed + 100)
    tau.test <- sim.test$tau
    
    ### Estimate the propensity score ###
    set.seed(seed)
    PS.fit <- nnet::multinom(W ~., data = data.frame(W = as.factor(sim.train$W), sim.train$X))
    W.hat <- predict(PS.fit, data = data.frame(sim.train$X), type = "probs")
    W.hat.test <- predict(PS.fit, data = data.frame(sim.test$X), type = "probs")
    
    ### Proposed method (Origin) ###
    gbm.tau.hat <- gbm(sim.train, sim.test, W.hat, seed = seed)
    pred_gbm <- gbm.tau.hat
    
    ctree.tau.hat <- ctree(sim.train, sim.test, W.hat, seed = seed)
    pred_ctree <- ctree.tau.hat
    
    ### Proposed method (Covariate Balancing propensity score; Fong et al.,2018) ###
    gbm.b.tau.hat <- gbm.b(sim.train, sim.test, W.hat, seed = seed)
    pred_gbm.b <- gbm.b.tau.hat
    
    ctree.b.tau.hat <- ctree.b(sim.train, sim.test, W.hat, seed = seed)
    pred_ctree.b <- ctree.b.tau.hat
    
    ### DR - Learner + RuleFit ###
    gbm.dr.tau.hat <- try(gbm.dr(sim.train, sim.test, W.hat, seed = seed))
    if(class(gbm.dr.tau.hat) != "try-error"){
      pred_gbm.dr <- gbm.dr.tau.hat
    }else{
      pred_gbm.dr <- list(tau.first = matrix(NA,n,mod$g-1), tau = matrix(NA,n,mod$g-1),
                          rank.first = rep(NA,n), rank = rep(NA,n),
                          order.first = matrix(NA,n,mod$g), order = matrix(NA,n,mod$g), 
                          time.first = NA, time = NA,
                          term.first = NA, term = NA)
    }  
    
    ctree.dr.tau.hat <- try(ctree.dr(sim.train, sim.test, W.hat, seed = seed))
    if(class(ctree.dr.tau.hat) != "try-error"){
      pred_ctree.dr <- ctree.dr.tau.hat
    }else{
      pred_ctree.dr <- list(tau.first = matrix(NA,n,mod$g-1), tau = matrix(NA,n,mod$g-1),
                            rank.first = rep(NA,n), rank = rep(NA,n),
                            order.first = matrix(NA,n,mod$g), order = matrix(NA,n,mod$g), 
                            time.first = NA, time = NA,
                            term.first = NA, term = NA)
    }  
    
    ### TRUE & Estimated HTE ###
    res <- data.frame(tau = tau.test, 
                      gbm.first = pred_gbm$tau.first, 
                      gbm = pred_gbm$tau,
                      ctree.first = pred_ctree$tau.first,
                      ctree = pred_ctree$tau,
                      gbm.b.first = pred_gbm.b$tau.first, 
                      gbm.b = pred_gbm.b$tau,
                      ctree.b.first = pred_ctree.b$tau.first,
                      ctree.b = pred_ctree.b$tau,
                      gbm.dr.first = pred_gbm.dr$tau.first,
                      ctree.dr.first = pred_ctree.dr$tau)
    
    save.res <- paste0(dir.save, "/add_simulation/Section2_setting2/prop",mod$sim,"_n",n,"_p",p,"_",mod$prob,"_",mod$g,"g","_ite",seed,".csv")
    write.csv(res, save.res, row.names = FALSE)
    
    
    ## Evaluation Result -------------------------------------------------------
    mPEHE <- data.frame(gbm.first = mean(sqrt(apply((pred_gbm$tau.first - tau.test)^2, 2, mean, na.rm = TRUE))),
                        gbm = mean(sqrt(apply((pred_gbm$tau - tau.test)^2, 2, mean, na.rm = TRUE))),
                        ctree.first = mean(sqrt(apply((pred_ctree$tau.first - tau.test)^2, 2, mean, na.rm = TRUE))),
                        ctree = mean(sqrt(apply((pred_ctree$tau - tau.test)^2, 2, mean, na.rm = TRUE))),
                        gbm.b.first = mean(sqrt(apply((pred_gbm.b$tau.first - tau.test)^2, 2, mean, na.rm = TRUE))),
                        gbm.b = mean(sqrt(apply((pred_gbm.b$tau - tau.test)^2, 2, mean, na.rm = TRUE))),
                        ctree.b.first = mean(sqrt(apply((pred_ctree.b$tau.first - tau.test)^2, 2, mean, na.rm = TRUE))),
                        ctree.b = mean(sqrt(apply((pred_ctree.b$tau - tau.test)^2, 2, mean, na.rm = TRUE))),
                        gbm.dr.first = mean(sqrt(apply((pred_gbm.dr$tau.first - tau.test)^2, 2, mean, na.rm = TRUE))),
                        gbm.dr = mean(sqrt(apply((pred_gbm.dr$tau - tau.test)^2, 2, mean, na.rm = TRUE))),
                        ctree.dr.first = mean(sqrt(apply((pred_ctree.dr$tau.first - tau.test)^2, 2, mean, na.rm = TRUE))),
                        ctree.dr = mean(sqrt(apply((pred_ctree.dr$tau - tau.test)^2, 2, mean, na.rm = TRUE))))
                        
    
    mbias <- data.frame(gbm.first = mean(abs(apply((pred_gbm$tau.first - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))),
                        gbm = mean(abs(apply((pred_gbm$tau - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))),
                        ctree.first = mean(abs(apply((pred_ctree$tau.first - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))),
                        ctree = mean(abs(apply((pred_ctree$tau - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))),
                        gbm.b.first = mean(abs(apply((pred_gbm.b$tau.first - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))),
                        gbm.b = mean(abs(apply((pred_gbm.b$tau - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))),
                        ctree.b.first = mean(abs(apply((pred_ctree.b$tau.first - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))),
                        ctree.b = mean(abs(apply((pred_ctree.b$tau - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))),
                        gbm.dr.first = mean(abs(apply((pred_gbm.dr$tau.first - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))),
                        gbm.dr = mean(abs(apply((pred_gbm.dr$tau - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))),
                        ctree.dr.first = mean(abs(apply((pred_ctree.dr$tau.first - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))),
                        ctree.dr = mean(abs(apply((pred_ctree.dr$tau - tau.test), 2, mean, na.rm = TRUE)/apply(tau.test,2,mean))))
                        
                        
    mcor <- data.frame(gbm.first = mean(diag(cor(pred_gbm$tau.first, tau.test, method = "spearman"))),
                       gbm = mean(diag(cor(pred_gbm$tau, tau.test, method = "spearman"))),
                       ctree.first = mean(diag(cor(pred_ctree$tau.first, tau.test, method = "spearman"))),
                       ctree = mean(diag(cor(pred_ctree$tau, tau.test, method = "spearman"))),
                       gbm.b.first = mean(diag(cor(pred_gbm.b$tau.first, tau.test, method = "spearman"))),
                       gbm.b = mean(diag(cor(pred_gbm.b$tau, tau.test, method = "spearman"))),
                       ctree.b.first = mean(diag(cor(pred_ctree.b$tau.first, tau.test, method = "spearman"))),
                       ctree.b = mean(diag(cor(pred_ctree.b$tau, tau.test, method = "spearman"))),
                       gbm.dr.first = mean(diag(cor(pred_gbm.dr$tau.first, tau.test, method = "spearman"))),
                       gbm.dr = mean(diag(cor(pred_gbm.dr$tau, tau.test, method = "spearman"))),
                       ctree.dr.first = mean(diag(cor(pred_ctree.dr$tau.first, tau.test, method = "spearman"))),
                       ctree.dr = mean(diag(cor(pred_ctree.dr$tau, tau.test, method = "spearman"))))
                       
    kappa <- data.frame(gbm.first = DescTools::CohenKappa(pred_gbm$rank.first, sim.test$rank),
                        gbm = DescTools::CohenKappa(pred_gbm$rank, sim.test$rank),
                        ctree.first = DescTools::CohenKappa(pred_ctree$rank.first, sim.test$rank),
                        ctree = DescTools::CohenKappa(pred_ctree$rank, sim.test$rank),
                        gbm.b.first = DescTools::CohenKappa(pred_gbm.b$rank.first, sim.test$rank),
                        gbm.b = DescTools::CohenKappa(pred_gbm.b$rank, sim.test$rank),
                        ctree.b.first = DescTools::CohenKappa(pred_ctree.b$rank.first, sim.test$rank),
                        ctree.b = DescTools::CohenKappa(pred_ctree.b$rank, sim.test$rank),
                        gbm.dr.first = DescTools::CohenKappa(pred_gbm.dr$rank.first, sim.test$rank),
                        gbm.dr = DescTools::CohenKappa(pred_gbm.dr$rank, sim.test$rank),
                        ctree.dr.first = DescTools::CohenKappa(pred_ctree.dr$rank.first, sim.test$rank),
                        ctree.dr = DescTools::CohenKappa(pred_ctree.dr$rank, sim.test$rank))
                        
    mcor_ord <- data.frame(gbm.first = mean(diag(cor(t(pred_gbm$order.first), t(sim.test$order), method = "spearman"))),
                           gbm = mean(diag(cor(t(pred_gbm$order), t(sim.test$order), method = "spearman"))),
                           ctree.first = mean(diag(cor(t(pred_ctree$order.first), t(sim.test$order), method = "spearman"))),
                           ctree = mean(diag(cor(t(pred_ctree$order), t(sim.test$order), method = "spearman"))),
                           gbm.b.first = mean(diag(cor(t(pred_gbm.b$order.first), t(sim.test$order), method = "spearman"))),
                           gbm.b = mean(diag(cor(t(pred_gbm.b$order), t(sim.test$order), method = "spearman"))),
                           ctree.b.first = mean(diag(cor(t(pred_ctree.b$order.first), t(sim.test$order), method = "spearman"))),
                           ctree.b = mean(diag(cor(t(pred_ctree.b$order), t(sim.test$order), method = "spearman"))),
                           gbm.dr.first = mean(diag(cor(t(pred_gbm.dr$order.first), t(sim.test$order), method = "spearman"))),
                           gbm.dr = mean(diag(cor(t(pred_gbm.dr$order), t(sim.test$order), method = "spearman"))),
                           ctree.dr.first = mean(diag(cor(t(pred_ctree.dr$order.first), t(sim.test$order), method = "spearman"))),
                           ctree.dr = mean(diag(cor(t(pred_ctree.dr$order), t(sim.test$order), method = "spearman"))))
    
    term <- data.frame(gbm.first = pred_gbm$term.first, 
                       gbm = pred_gbm$term,
                       ctree.first = pred_ctree$term.first, 
                       ctree = pred_ctree$term, 
                       gbm.b.first = pred_gbm.b$term.first, 
                       gbm.b = pred_gbm.b$term,
                       ctree.b.first = pred_ctree.b$term.first, 
                       ctree.b = pred_ctree.b$term, 
                       gbm.dr.first = pred_gbm.dr$term.first,
                       gbm.dr = pred_gbm.dr$term,
                       ctree.dr.first = pred_ctree.dr$term.first,
                       ctree.dr = pred_ctree.dr$term)
                           
    time <- data.frame(gbm.first = pred_gbm$time.first, 
                       gbm = pred_gbm$time,
                       ctree.first = pred_ctree$time.first, 
                       ctree = pred_ctree$time, 
                       gbm.b.first = pred_gbm.b$time.first, 
                       gbm.b = pred_gbm.b$time,
                       ctree.b.first = pred_ctree.b$time.first, 
                       ctree.b = pred_ctree.b$time, 
                       gbm.dr.first = pred_gbm.dr$time.first,
                       gbm.dr = pred_gbm.dr$time,
                       ctree.dr.first = pred_ctree.dr$time.first,
                       ctree.dr = pred_ctree.dr$time)
    
    sim.evl <- rbind(mPEHE = mPEHE, mbias = mbias, mcor = mcor, kappa = kappa, mcor_ord = mcor_ord, term = term, time = time)
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
  
  save.res <- paste0(dir.save, "/add_simulation/Section2_setting2/prop_evl_sim",mod$sim,"_n",n,"_p",p,"_",mod$prob,"_",mod$g,"g",".csv")
  write.csv(sim.res, save.res, row.names = FALSE)
}
