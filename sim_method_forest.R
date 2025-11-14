
########################## T-learner (random forest)  ################################

tforest <- function(sim.train, sim.test, seed = 1){
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  X.test <- sim.test$X
  K.ind <- unique(sort(sim.train$W))
  
  start <- proc.time()
  
  res <- matrix(0, nrow(X.train), length(K.ind))
  
  for(w in K.ind){
    
    x.train <- X.train[which(W.train == w),]
    y.train <- Y.train[which(W.train == w)]

    fit = ranger::ranger(y = y.train, x = x.train)
    y_pred <- predict(fit, as.matrix(X.test))$prediction
    res[,which(K.ind == w)] <- y_pred
  }
  
  end <- (proc.time() - start)[3]
  
  tau.hat <- res[,-1] - res[,1]
  
  rank <- c()
  order.can <- NULL 
  for(t in 1:nrow(res)){
    rank[t] <- which.max(res[t,])
    order.can <- rbind(order.can, order(res[t,], decreasing = TRUE))
  }
  
  res.can <- list(tau = tau.hat, rank = rank, order = order.can,time = end)
  res.can
}

########################## S-learner (random forest)  ################################

sforest <- function(sim.train, sim.test, seed = 1){
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- cbind(W = sim.train$W, sim.train$X)
  X.test <- sim.test$X
  K.ind <- unique(sort(sim.train$W))
  
  start <- proc.time()
  
  res <- matrix(0, nrow(X.train), length(K.ind))
  
  fit = ranger::ranger(x = as.matrix(X.train), y = Y.train)
  
  
  for(w in K.ind){
    y_pred <- predict(fit, as.matrix(cbind(W = w, X.test)))$prediction
    res[,which(K.ind == w)] <- y_pred
  }
  
  end <- (proc.time() - start)[3]
  
  tau.hat <- res[,-1] - res[,1]
  
  rank <- c()
  order.can <- NULL 
  for(t in 1:nrow(res)){
    rank[t] <- which.max(res[t,])
    order.can <- rbind(order.can, order(res[t,], decreasing = TRUE))
  }
  
  res.can <- list(tau = tau.hat, rank = rank, order = order.can,time = end)
  res.can
}

########################## X-learner: naive extension of X-learner (random forest)  ##

xforest <- function(sim.train, sim.test, W.hat, seed = 1){
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  X.test <- sim.test$X
  K.ind <- unique(sort(sim.train$W))
  
  start <- proc.time()
  
  mu.can <- matrix(0, nrow(X.train), length(K.ind))
  
  for(w in K.ind){
    
    x.train <- X.train[which(W.train == w),]
    y.train <- Y.train[which(W.train == w)]

    fit = ranger::ranger(x = as.matrix(x.train), y = y.train)
    mu_pred <- predict(fit, as.matrix(X.train[which(W.train!=w),]))$prediction

    mu.can[which(W.train != w), which(K.ind == w)] <- mu_pred
    mu.can[which(W.train == w), which(K.ind == w)] <- y.train
  }
  
  tau.hat <- matrix(0, nrow(X.train), length(K.ind)-1)
  
  for(g in 2:(length(K.ind))){
    
    dk <- mu.can[which(W.train == K.ind[g]), g] - mu.can[which(W.train == K.ind[g]), 1]
    d0 <- mu.can[which(W.train == 0), g] - mu.can[which(W.train == 0), 1]
    
    t_fit_k = ranger::ranger(x = as.matrix(X.train[which(W.train == K.ind[g]),]), y = dk)
    t_fit_0 = ranger::ranger(x = as.matrix(X.train[which(W.train == 0),]), y = d0)
    
    tau_k_pred <- predict(t_fit_k, as.matrix(X.test))$prediction
    tau_0_pred <- predict(t_fit_0, as.matrix(X.test))$prediction
    
    tau.hat[,g-1] <- (W.hat[,g]/(W.hat[,g] + W.hat[,1]))*tau_k_pred + (W.hat[,1]/(W.hat[,g] + W.hat[,1]))*tau_0_pred
  }
  
  end <- (proc.time() - start)[3]
  
  rank <- c()
  order.can <- NULL 
  for(t in 1:nrow(tau.hat)){
    rank[t] <- which.max(cbind(0,tau.hat)[t,]) 
    order.can <- rbind(order.can, order(cbind(0,tau.hat)[t,], decreasing = TRUE))
  }
  
  res.can <- list(tau = tau.hat, rank = rank, order = order.can,time = end)
  res.can
}

########################## m-learner (random forest)  ################################

mforest <- function(sim.train, sim.test, W.hat, seed = 1){
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  X.test <- sim.test$X
  K.ind <- unique(sort(sim.train$W))
  
  start <- proc.time()
    
  tau.can <- matrix(0, nrow(X.train), length(K.ind)-1)
  
  for(g in 1:(length(K.ind)-1)){
    #zk <- Y.train*as.numeric(W.train == K.ind[g])/W.hat[,g]
    zk <- Y.train*as.numeric(W.train == K.ind[g + 1])/W.hat[,g + 1] - Y.train*as.numeric(W.train == 0)/W.hat[,1]
    fit = ranger::ranger(x = as.matrix(X.train), y = zk)
    tau_pred <- predict(fit, as.matrix(X.test))$prediction
    tau.can[,g] <- tau_pred
  }
  
  end <- (proc.time() - start)[3]
  
  tau.hat <- tau.can
  
  rank <- c()
  order.can <- NULL 
  for(t in 1:nrow(tau.hat)){
    rank[t] <- which.max(cbind(0,tau.hat)[t,]) 
    order.can <- rbind(order.can, order(cbind(0,tau.hat)[t,], decreasing = TRUE))
  }
  
  res.can <- list(tau = tau.hat, rank = rank, order = order.can,time = end)
  res.can
}

######################### DR-learner (random forest)  ################################

drforest <- function(sim.train, sim.test, W.hat, seed = 1){
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  X.test <- sim.test$X
  K.ind <- unique(sort(sim.train$W))
  
  start <- proc.time()
  
  res <- matrix(0, nrow(X.train), length(K.ind))
  
  for(w in K.ind){
    
    x.train <- X.train[which(W.train == w),]
    y.train <- Y.train[which(W.train == w)]
    
    fit = ranger::ranger(x = as.matrix(x.train), y = y.train)
    y_pred <- predict(fit, as.matrix(X.train))$prediction
    
    zk <- (Y.train - y_pred)*as.numeric(W.train == w)/W.hat[,w + 1] + y_pred
    z_fit = ranger::ranger(x = as.matrix(X.train), y = zk)
    tau_pred <- predict(z_fit, as.matrix(X.test))$prediction
    
    res[,which(K.ind == w)] <- tau_pred
  }
  
  end <- (proc.time() - start)[3]
  
  tau.hat <- res[,-1] - res[,1]
  
  rank <- c()
  order.can <- NULL 
  for(t in 1:nrow(res)){
    rank[t] <- which.max(res[t,])
    order.can <- rbind(order.can, order(res[t,], decreasing = TRUE))
  }
  
  res.can <- list(tau = tau.hat, rank = rank, order = order.can,time = end)
  res.can
}

########### Reference-free R-learner #########################

rfforest <- function(sim.train, sim.test, W.hat, W.hat.test, seed = 1){
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  X.test <- sim.test$X
  K.ind <- unique(sort(sim.train$W))
  pi.train <- W.hat
  pi.test <- W.hat.test
  
  start <- proc.time()
  
  fit = ranger::ranger(x = as.matrix(X.train), y = Y.train)
  m.train <- predict(fit, as.matrix(X.train))$prediction
  
  ######### Code from Zou et al.(2023) ########
  
  nobs = dim(X.train)[1]
  pobs = dim(X.train)[2]
  
  ################## Simplex Method #####################
  k = length(K.ind)
  z = rep(1, k-1)
  e = diag(x = 1, k-1)
  W1=(k-1)^(-0.5)*z
  W = cbind(W1, (k/(k-1))^(0.5)*e - z*(1+sqrt(k))/(k-1)^1.5)
  
  x.whole = as.matrix(cbind(1, X.train))
  x.new = sapply(seq(nobs), function(i){
    as.vector(outer(W[,(W.train[i]+1)] ,x.whole[i,]))/pi.train[i,(W.train[i]+1)]
  })
  x.new = t(x.new)
  penalty_f = c(rep(0,k-1), rep(1, pobs*(k-1)))
  fit.tau = glmnet::cv.glmnet(x.new, Y.train-m.train, family = "gaussian",
                              penalty.factor = penalty_f, intercept=FALSE)
  ## estimate tau for testing data ##
  x.test.whole = as.matrix(cbind(1, X.test))
  best.beta = coef(fit.tau, s="lambda.min")
  best.beta = matrix(best.beta[-1], nrow = pobs+1, byrow = T)
  est.tau = x.test.whole %*% best.beta %*% W / pi.test
  
  end <- (proc.time() - start)[3]
  
  tau.hat <- est.tau[,-1] - est.tau[,1]
  
  rank <- c()
  order.can <- NULL 
  for(t in 1:nrow(est.tau)){
    rank[t] <- which.max(est.tau[t,])
    order.can <- rbind(order.can, order(est.tau[t,], decreasing = TRUE))
  }
  
  res.can <- list(tau = tau.hat, rank = rank, order = order.can,time = end)
  res.can
}  

#################### R-Learner #################################################

rforest <- function(sim.train, sim.test, W.hat, seed = 1){
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  X.test <- sim.test$X
  K.ind <- unique(sort(sim.train$W))
  pi.train <- W.hat
  
  start <- proc.time()
  
  fit = ranger::ranger(x = as.matrix(X.train), y = Y.train)
  m.train <- predict(fit, as.matrix(X.train))$prediction
  
  ######### Code from Zou et al.(2023) ########
  
  nobs = dim(X.train)[1]
  pobs = dim(X.train)[2]
  
  Trt = as.factor(W.train)
  Trt.name = levels(W.train)
  Trt = as.numeric(W.train)
  ps.train = pi.train
  K.grp = K.ind
  y.adj = Y.train - m.train
  
  # R will have multiple results with different reference group
  r <- 1
  pi.hat = ps.train[, -r]
  W.ind = sapply(K.grp[-r], function(z){as.numeric(Trt==z)})
  adj.pi = W.ind - pi.hat
  x.tilde = matrix(as.vector(apply(adj.pi, 2, function(j) j * cbind(1,as.matrix(X.train)))), nrow = nrow(X.train), byrow = F)
  penalty_factor = rep(c(0, array(1, pobs)), length(K.grp)-1)
  fit.R = glmnet::cv.glmnet(x.tilde, y.adj, penalty.factor = penalty_factor, family = "gaussian", maxit = 100000, intercept = F)
  gc(verbose = FALSE)
  
  ## test data
  Est.R = NULL
  for ( k in 1:(length(K.grp)-1) ) {
    ind = rep(0, length(K.grp)-1); ind[k] = 1
    newX = t(ind) %x% cbind(1,as.matrix(X.test))
    est.R = predict(fit.R, newx = newX, type = "response", s = "lambda.min")
    Est.R = cbind(Est.R, est.R)
  }
  colnames(Est.R) <- paste(Trt.name[-which(K.grp == r)], Trt.name[r], sep = "-")
  
  end <- (proc.time() - start)[3]
  
  ######################################################################################
  
  tau.hat <- Est.R
  
  rank <- c()
  order.can <- NULL 
  for(t in 1:nrow(tau.hat)){
    rank[t] <- which.max(cbind(0,tau.hat)[t,]) 
    order.can <- rbind(order.can, order(cbind(0,tau.hat)[t,], decreasing = TRUE))
  }
  
  res.can <- list(tau = tau.hat, rank = rank, order = order.can,time = end)
  res.can
}
