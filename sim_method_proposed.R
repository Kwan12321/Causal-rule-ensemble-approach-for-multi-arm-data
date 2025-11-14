
########################## GBM  ################################

gbm <- function(sim.train, sim.test, W.hat, seed = 1){
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  
  start <- proc.time()
  
  fit <- mcre_agenet(X = X.train, Y = Y.train, W = W.train, W.hat = W.hat, learnrate = 0.01,
                     gen.tree = "gbm", tune.init = "CV", tune.fin = "CV", nfolds = 5, seed = seed)
  tau.hat <- predict.mcre(fit, newdata = sim.test$data)
  
  time.first <- fit$time.first
  time <- fit$time 

  rank <- c()
  rank.first <- c()
  order.can <- NULL 
  order.can.first <- NULL 
  for(t in 1:nrow(tau.hat$tau)){
    rank[t] <- which.max(cbind(0, tau.hat$tau)[t,])
    rank.first[t] <- which.max(cbind(0, tau.hat$tau.first)[t,])
    order.can <- rbind(order.can, order(cbind(0, tau.hat$tau)[t,], decreasing = TRUE))
    order.can.first <- rbind(order.can.first, order(cbind(0, tau.hat$tau.first)[t,], decreasing = TRUE))
  }
  
  ## Correct Variable Selection Rate ---------------
  r1 <- fit$beta[which(fit$beta != 0)]
  r0 <- fit$beta.first[which(fit$beta.first != 0)]
  err.var <- paste0("X",6:ncol(X.train), collapse = "|")
  err.var <- paste0("(",err.var,")")
  cor.rate.all <- 1 - length(grep(err.var, fit$rules))/length(fit$rules)
  cor.rate <- 1 - length(grep(err.var, names(r1)))/length(r1)
  cor.rate.first <- 1 - length(grep(err.var, names(r0)))/length(r0)
  
  ## Number of terms --------------------------------
  term.all <- length(fit$rules)
  term <- sum(fit$beta != 0)/length(unique(W.train))
  term.first <- sum(fit$beta.first != 0)/length(unique(W.train))
  
  ## Median of interaction -----------------------------
  int <- median(unlist(lapply(strsplit(names(fit$beta)[fit$beta != 0],"&"), length)) - 1)
  int.first <- median(unlist(lapply(strsplit(names(fit$beta.first)[fit$beta.first != 0],"&"), length)) - 1)
  
  res.can <- list(tau = tau.hat$tau, tau.first = tau.hat$tau.first,
                  rank = rank, rank.first = rank.first,
                  order = order.can, order.first = order.can.first,
                  time = time, time.first = time.first,
                  term = term, term.first = term.first, term.all = term.all,
                  cor.rate = cor.rate, cor.rate.first = cor.rate.first, cor.rate.all = cor.rate.all,
                  int = int, int.first = int.first)
  res.can
}

############## Conditional inference tree (Ctree) ##############################

ctree <- function(sim.train, sim.test, W.hat, seed = 1){
 
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  
  start <- proc.time()
  
  fit <- mcre_agenet(X = X.train, Y = Y.train, W = W.train, W.hat = W.hat, learnrate = 0.01, 
                     gen.tree = "ctree", tune.init = "CV", tune.fin = "CV", nfolds = 5, seed = seed)
  tau.hat <- predict.mcre(fit, newdata = sim.test$data)
  
  time.first <- fit$time.first
  time <- fit$time 
  
  rank <- c()
  rank.first <- c()
  order.can <- NULL 
  order.can.first <- NULL 
  for(t in 1:nrow(tau.hat$tau)){
    rank[t] <- which.max(cbind(0, tau.hat$tau)[t,])
    rank.first[t] <- which.max(cbind(0, tau.hat$tau.first)[t,])
    order.can <- rbind(order.can, order(cbind(0, tau.hat$tau)[t,], decreasing = TRUE))
    order.can.first <- rbind(order.can.first, order(cbind(0, tau.hat$tau.first)[t,], decreasing = TRUE))
  }
  
  ## Correct Variable Selection Rate ---------------
  r1 <- fit$beta[which(fit$beta != 0)]
  r0 <- fit$beta.first[which(fit$beta.first != 0)]
  err.var <- paste0("X",6:ncol(X.train), collapse = "|")
  err.var <- paste0("(",err.var,")")
  cor.rate.all <- 1 - length(grep(err.var, fit$rules))/length(fit$rules)
  cor.rate <- 1 - length(grep(err.var, names(r1)))/length(r1)
  cor.rate.first <- 1 - length(grep(err.var, names(r0)))/length(r0)
  
  ## Number of terms --------------------------------
  term.all <- length(fit$rules)
  term <- sum(fit$beta != 0)/length(unique(W.train))
  term.first <- sum(fit$beta.first != 0)/length(unique(W.train))
  
  ## Median of interaction -----------------------------
  int <- median(unlist(lapply(strsplit(names(fit$beta)[fit$beta != 0],"&"), length)) - 1)
  int.first <- median(unlist(lapply(strsplit(names(fit$beta.first)[fit$beta.first != 0],"&"), length)) - 1)
  
  res.can <- list(tau = tau.hat$tau, tau.first = tau.hat$tau.first,
                  rank = rank, rank.first = rank.first,
                  order = order.can, order.first = order.can.first,
                  time = time, time.first = time.first,
                  term = term, term.first = term.first, term.all = term.all,
                  cor.rate = cor.rate, cor.rate.first = cor.rate.first, cor.rate.all = cor.rate.all,
                  int = int, int.first = int.first)
  res.can
}


########################## GBM - Balance ################################

gbm.b <- function(sim.train, sim.test, W.hat, seed = 1){
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  
  start <- proc.time()
  
  fit <- mcre_agenet(X = X.train, Y = Y.train, W = W.train, W.hat = W.hat, learnrate = 0.01, 
                     gen.tree = "gbm.b", tune.init = "CV", tune.fin = "CV", nfolds = 5, seed = seed)
  tau.hat <- predict.mcre(fit, newdata = sim.test$data)
  
  time.first <- fit$time.first
  time <- fit$time 
  
  rank <- c()
  rank.first <- c()
  order.can <- NULL 
  order.can.first <- NULL 
  for(t in 1:nrow(tau.hat$tau)){
    rank[t] <- which.max(cbind(0, tau.hat$tau)[t,])
    rank.first[t] <- which.max(cbind(0, tau.hat$tau.first)[t,])
    order.can <- rbind(order.can, order(cbind(0, tau.hat$tau)[t,], decreasing = TRUE))
    order.can.first <- rbind(order.can.first, order(cbind(0, tau.hat$tau.first)[t,], decreasing = TRUE))
  }
  
  ## Correct Variable Selection Rate ---------------
  r1 <- fit$beta[which(fit$beta != 0)]
  r0 <- fit$beta.first[which(fit$beta.first != 0)]
  err.var <- paste0("X",6:ncol(X.train), collapse = "|")
  err.var <- paste0("(",err.var,")")
  cor.rate.all <- 1 - length(grep(err.var, fit$rules))/length(fit$rules)
  cor.rate <- 1 - length(grep(err.var, names(r1)))/length(r1)
  cor.rate.first <- 1 - length(grep(err.var, names(r0)))/length(r0)
  
  ## Number of terms --------------------------------
  term.all <- length(fit$rules)
  term <- sum(fit$beta != 0)/length(unique(W.train))
  term.first <- sum(fit$beta.first != 0)/length(unique(W.train))
  
  ## Median of interaction -----------------------------
  int <- median(unlist(lapply(strsplit(names(fit$beta)[fit$beta != 0],"&"), length)) - 1)
  int.first <- median(unlist(lapply(strsplit(names(fit$beta.first)[fit$beta.first != 0],"&"), length)) - 1)
  
  res.can <- list(tau = tau.hat$tau, tau.first = tau.hat$tau.first,
                  rank = rank, rank.first = rank.first,
                  order = order.can, order.first = order.can.first,
                  time = time, time.first = time.first,
                  term = term, term.first = term.first, term.all = term.all,
                  cor.rate = cor.rate, cor.rate.first = cor.rate.first, cor.rate.all = cor.rate.all,
                  int = int, int.first = int.first)
  res.can
}

############## Conditional inference tree (Ctree - Balance) ##############################

ctree.b <- function(sim.train, sim.test, W.hat, seed = 1){
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  
  start <- proc.time()
  
  fit <- mcre_agenet(X = X.train, Y = Y.train, W = W.train, W.hat = W.hat, learnrate = 0.01,
                     gen.tree = "ctree.b", tune.init = "CV", tune.fin = "CV", nfolds = 5, seed = seed)
  tau.hat <- predict.mcre(fit, newdata = sim.test$data)
  
  time.first <- fit$time.first
  time <- fit$time 
  
  rank <- c()
  rank.first <- c()
  order.can <- NULL 
  order.can.first <- NULL 
  for(t in 1:nrow(tau.hat$tau)){
    rank[t] <- which.max(cbind(0, tau.hat$tau)[t,])
    rank.first[t] <- which.max(cbind(0, tau.hat$tau.first)[t,])
    order.can <- rbind(order.can, order(cbind(0, tau.hat$tau)[t,], decreasing = TRUE))
    order.can.first <- rbind(order.can.first, order(cbind(0, tau.hat$tau.first)[t,], decreasing = TRUE))
  }
  
  ## Correct Variable Selection Rate ---------------
  r1 <- fit$beta[which(fit$beta != 0)]
  r0 <- fit$beta.first[which(fit$beta.first != 0)]
  err.var <- paste0("X",6:ncol(X.train), collapse = "|")
  err.var <- paste0("(",err.var,")")
  cor.rate.all <- 1 - length(grep(err.var, fit$rules))/length(fit$rules)
  cor.rate <- 1 - length(grep(err.var, names(r1)))/length(r1)
  cor.rate.first <- 1 - length(grep(err.var, names(r0)))/length(r0)
  
  ## Number of terms --------------------------------
  term.all <- length(fit$rules)
  term <- sum(fit$beta != 0)/length(unique(W.train))
  term.first <- sum(fit$beta.first != 0)/length(unique(W.train))
  
  ## Median of interaction -----------------------------
  int <- median(unlist(lapply(strsplit(names(fit$beta)[fit$beta != 0],"&"), length)) - 1)
  int.first <- median(unlist(lapply(strsplit(names(fit$beta.first)[fit$beta.first != 0],"&"), length)) - 1)
  
  res.can <- list(tau = tau.hat$tau, tau.first = tau.hat$tau.first,
                  rank = rank, rank.first = rank.first,
                  order = order.can, order.first = order.can.first,
                  time = time, time.first = time.first,
                  term = term, term.first = term.first, term.all = term.all,
                  cor.rate = cor.rate, cor.rate.first = cor.rate.first, cor.rate.all = cor.rate.all,
                  int = int, int.first = int.first)
  res.can
}

#################################  GBM - DR #################################### 

gbm.dr <- function(sim.train, sim.test, W.hat, seed = 1){
  
  ### Default parameters
  ntrees <- 333; depth <- 2; sampfrac <- NULL; learnrate <- 0.01
  
  ##############################################################
  ### Step 1 : Calculate the pesudo-outcome using DR-Learner ###
  ##############################################################
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  X.test <- sim.test$X
  K.ind <- unique(sort(sim.train$W))
  
  start <- proc.time()
  
  Y.tilde <- matrix(0, nrow(X.train), length(K.ind))
  
  for(w in K.ind){
    
    x.train <- X.train[which(W.train == w),]
    y.train <- Y.train[which(W.train == w)]
    
    fit = gbm::gbm(y~., data = data.frame(as.matrix(x.train), y = y.train), n.trees = ntrees, n.minobsinnode = 1)
    y_pred <- predict(fit, data.frame(X.train))
    
    zk <- (Y.train - y_pred)*as.numeric(W.train == w)/W.hat[,w + 1] + y_pred
    z_fit = gbm::gbm(y~., data = data.frame(as.matrix(X.train), y = zk), n.trees = ntrees)
    tau_pred <- predict(z_fit, data.frame(X.train))
    
    Y.tilde[,which(K.ind == w)] <- tau_pred
  }
  
  
  ################################
  ### Step 2 : Generate Rules  ###
  ################################
  
  
  Y.names <- paste0("Y",1:ncol(Y.tilde))
  colnames(Y.tilde) <- Y.names
  X.names <- colnames(X.train)
  
  # Train outcome & Formula
  data_with_y_learn <- data.learn <- cbind(Y.tilde, X.train)
  formula <- as.formula(paste0(paste(Y.names, collapse = "+"),"~",paste(X.names, collapse = "+")))
  
  # Determine the sample fraction size for training the base learner (Default using the setting in Friedman and Popescu (2008))
  n <- nrow(X.train)
  if (is.null(sampfrac)) {
    size <- min(n / 2, 100 + 6 * sqrt(n))
  } else {
    size <- sampfrac * n
  }
  
  # Create the index for training data for each base learner
  subsample <- list()
  subsample <- mapply(function(i) {
    if (size == n) {
      subsample[[i]] <- sample(1:n, size = size, replace = FALSE)
    } 
    else if (size < n) {
      subsample[[i]] <- sample(1:n, size = round(size), replace = FALSE)
    }
  }, i = 1:ntrees, SIMPLIFY = FALSE)
  
  # Calculate the number of terminal nodes for each tree using an exponential distribution
  maxdepth <- ceiling(log(2 + floor(rexp(ntrees, rate = 1 / (2^depth - 2))), base = 2))
  
  # Initialize estimates
  y <- data.learn[ , Y.names]
  eta_0 <- apply(y, 2, weighted.mean, weights = rep(1, nrow(y)))
  eta <- t(replicate(n = nrow(y), expr = eta_0))
  
  # Calculate the pseudo-response variable using the negative gradient function
  data_with_y_learn[, Y.names] <- y - eta
  
  # Initialize an empty vector to store the rules
  rules <- c()
  #var.imp <- numeric(ncol(X))
  #names(var.imp) <- X.names
  
  for(i in 1:ntrees) {
    
    # Set up controls for rpart tree
    tree.control <- rpart::rpart.control(maxdepth = maxdepth[i])
    
    # Fit the individual tree
    tree <- rpart::rpart(formula, control = tree.control, data = data_with_y_learn[subsample[[i]], ])
    
    # Decompose the tree into rules
    paths <- rpart::path.rpart(tree, nodes = rownames(tree$frame), print.it = FALSE, pretty = 0)
    paths <- unname(sapply(sapply(paths, `[`, index = -1), paste, collapse = " & ")[-1])
    
    # Append new rules to existing rules
    rules <- c(rules, paths)
    
    # Update eta using predictions from the new tree
    eta <- eta + learnrate * predict(tree, newdata = data_with_y_learn)
    
    # Update pseudo-response variable using the new eta
    data_with_y_learn[ ,Y.names] <- y - eta
    
  }
  
  rules <- unique(rules)
  print(paste("Number of unique rules:", length(rules)))
  
  # Remove duplicate and complementary rules
  rules <- remove_complement(X.train, rules)[[1]]
  if (length(rules) == 0) {
    stop("Error after removing complements: no rules remain.")
  }
  
  ################################
  ### Step 3 : Rule Ensemble   ###
  ################################
  
  data <- data.frame(Y = Y.train, W = W.train, X.train)

  # Modified version of linear terms ("Winsorized" linear terms)
  var.trunc <- winsorize(data)
  
  # Transform rules to binary matrix
  rulevars <- tran_rules(dat = data, rules = rules)
  
  # Combine the rule terms and linear terms
  base.fun <- cbind(rulevars, var.trunc)
  base.names <- c(rules, colnames(var.trunc))
  colnames(base.fun) <- base.names
  
  ## Scale of the base function (Freidman & Popescu, 2008)
  var.sd <- apply(data[,-c(1:2)], 2, sd)
  x.scale <- c(rep(1,ncol(rulevars)), var.sd/0.4)
  
  # Normalized the base functions 
  base.fun <- scale(base.fun, scale = x.scale, center = FALSE) 
  
  lasso.fit <- glmnet::cv.glmnet(y = Y.tilde, x = base.fun, family="mgaussian", nfolds = 5)
  beta.first0 <- predict(lasso.fit, base.fun, type = "coefficients")
  beta.first <- NULL
  for(b in 1:length(beta.first0)){
    beta.first0[[b]][-1] <- beta.first0[[b]][-1]/x.scale
    beta.first <- cbind(beta.first, beta.first0[[b]])
  }
  intercept.first <- beta.first[1,] 
  beta.first <- beta.first[-1,]

  end1 <- (proc.time() - start)[3]
  
  ### Prediction Using Test data ###
  X.test <- tran_rules(sim.test$data, rules)
  X.test <- cbind(X.test, sim.test$data[,-c(1,2)])
  
  tau.hat.first <- matrix(rep(intercept.first, nrow(sim.test$data)), nrow(sim.test$data), length(unique(sim.test$W)), byrow = TRUE) + as.matrix(X.test)%*%as.matrix(beta.first) 
  tau.hat.first <- tau.hat.first[,-1] - tau.hat.first[,1]
  
  ########################################
  ### Step4 : Rule Ensemble (Adaptive) ###
  ########################################
  
  penalty.factor <- rowMeans(as.matrix(beta.first))
  penalty.factor <- 1/abs(penalty.factor + .Machine$double.eps*2)
  
  adlasso.fit <- glmnet::cv.glmnet(y = Y.tilde, x = base.fun, family="mgaussian", nfolds = 5, penalty.factor = penalty.factor)
  
  beta0 <- predict(adlasso.fit, base.fun, type = "coefficients")
  beta <- NULL
  for(b in 1:length(beta0)){
    beta0[[b]][-1] <- beta0[[b]][-1]/x.scale
    beta <- cbind(beta, beta0[[b]])
  }
  intercept<- beta[1,] 
  beta <- beta[-1,]
  
  end2 <- (proc.time() - start)[3]
  
  tau.hat <- matrix(rep(intercept, nrow(sim.test$data)), nrow(sim.test$data), length(unique(sim.test$W)), byrow = TRUE) + as.matrix(X.test)%*%as.matrix(beta) 
  tau.hat <- tau.hat[,-1] - tau.hat[,1]
  
  ###################
  ### Evaluation  ###   
  ###################
  
  time.first <- end1
  time <- end2
  
  rank.first <- c()
  rank <- c()
  order.can.first <- NULL
  order.can <- NULL
  for(t in 1:nrow(tau.hat)){
    rank.first[t] <- which.max(cbind(0, tau.hat.first)[t,])
    rank[t] <- which.max(cbind(0, tau.hat)[t,])
    order.can.first <- rbind(order.can.first, order(cbind(0, tau.hat.first)[t,], decreasing = TRUE))
    order.can <- rbind(order.can, order(cbind(0, tau.hat)[t,], decreasing = TRUE))
  }
  
  ## Correct Variable Selection Rate ---------------
  r0 <- beta.first[apply(beta.first != 0, 1, sum) == 3,]
  r1 <- beta[apply(beta != 0, 1, sum) == 3,]
  err.var <- paste0("X",6:ncol(X.train), collapse = "|")
  err.var <- paste0("(",err.var,")")
  cor.rate.first <- 1 - length(grep(err.var, rownames(r0)))/length(r0)
  cor.rate <- 1 - length(grep(err.var, rownames(r1)))/length(r1)
  
  ## Number of terms --------------------------------
  term.first <- sum(beta.first != 0)/length(unique(W.train))
  term <- sum(beta != 0)/length(unique(W.train))
  
  res.can <- list(tau.first = tau.hat.first, tau = tau.hat,  
                  rank.first = rank.first, rank = rank,
                  order.first = order.can.first, order = order.can,
                  time.first = time.first, time = time,
                  term.first = term.first, term = term)
  res.can
  
}

#################################  CTREE - DR #################################### 

ctree.dr <- function(sim.train, sim.test, W.hat, seed = 1){
  
  ### Default parameters
  ntrees <- 333; depth <- 2; sampfrac <- NULL; learnrate <- 0.01
  
  ##############################################################
  ### Step 1 : Calculate the pesudo-outcome using DR-Learner ###
  ##############################################################
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  X.test <- sim.test$X
  K.ind <- unique(sort(sim.train$W))
  
  start <- proc.time()
  
  Y.tilde <- matrix(0, nrow(X.train), length(K.ind))
  
  for(w in K.ind){
    
    x.train <- X.train[which(W.train == w),]
    y.train <- Y.train[which(W.train == w)]
    
    fit = gbm::gbm(y~., data = data.frame(as.matrix(x.train), y = y.train), n.trees = ntrees, n.minobsinnode = 1)
    y_pred <- predict(fit, data.frame(X.train))
    
    zk <- (Y.train - y_pred)*as.numeric(W.train == w)/W.hat[,w + 1] + y_pred
    z_fit = gbm::gbm(y~., data = data.frame(as.matrix(X.train), y = zk), n.trees = ntrees)
    tau_pred <- predict(z_fit, data.frame(X.train))
    
    Y.tilde[,which(K.ind == w)] <- tau_pred
  }
  
  
  ################################
  ### Step 2 : Generate Rules  ###
  ################################
  
  
  Y.names <- paste0("Y",1:ncol(Y.tilde))
  colnames(Y.tilde) <- Y.names
  X.names <- colnames(X.train)
  
  # Train outcome & Formula
  data_with_y_learn <- data.learn <- cbind(Y.tilde, X.train)
  formula <- as.formula(paste0(paste(Y.names, collapse = "+"),"~",paste(X.names, collapse = "+")))
  
  # Determine the sample fraction size for training the base learner (Default using the setting in Friedman and Popescu (2008))
  n <- nrow(X.train)
  if (is.null(sampfrac)) {
    size <- min(n / 2, 100 + 6 * sqrt(n))
  } else {
    size <- sampfrac * n
  }
  
  # Create the index for training data for each base learner
  subsample <- list()
  subsample <- mapply(function(i) {
    if (size == n) {
      subsample[[i]] <- sample(1:n, size = size, replace = FALSE)
    } 
    else if (size < n) {
      subsample[[i]] <- sample(1:n, size = round(size), replace = FALSE)
    }
  }, i = 1:ntrees, SIMPLIFY = FALSE)
  
  # Calculate the number of terminal nodes for each tree using an exponential distribution
  maxdepth <- ceiling(log(2 + floor(rexp(ntrees, rate = 1 / (2^depth - 2))), base = 2))
  
  # Initialize estimates
  y <- data.learn[ , Y.names]
  eta_0 <- apply(y, 2, weighted.mean, weights = rep(1, nrow(y)))
  eta <- t(replicate(n = nrow(y), expr = eta_0))
  
  # Calculate the pseudo-response variable using the negative gradient function
  data_with_y_learn[, Y.names] <- y - eta
  
  # Initialize an empty vector to store the rules
  rules <- c()
  
  for(i in 1:ntrees) {
    
    # Set up controls for ctree
    ctree.control <- partykit::ctree_control(maxdepth = maxdepth[i])
    
    # Fit the individual tree
    tree <- partykit::ctree(formula = formula, control = ctree.control, data = data_with_y_learn[subsample[[i]], ])
    
    # Decompose tree into rules
    paths <- list.rules(tree)
    
    # Append new rules to existing rules
    rules <- c(rules, paths)
    
    # Update eta using predictions from the new tree
    eta <- eta + learnrate * predict(tree, newdata = data_with_y_learn)
    
    # Update pseudo-response variable using the new eta
    data_with_y_learn[ ,Y.names] <- y - eta
    
  }
  
  rules <- unique(rules)
  print(paste("Number of unique rules:", length(rules)))
  
  # Remove duplicate and complementary rules
  rules <- remove_complement(X.train, rules)[[1]]
  if (length(rules) == 0) {
    stop("Error after removing complements: no rules remain.")
  }
  
  ################################
  ### Step 3 : Rule Ensemble   ###
  ################################
  
  data <- data.frame(Y = Y.train, W = W.train, X.train)
  
  # Modified version of linear terms ("Winsorized" linear terms)
  var.trunc <- winsorize(data)
  
  # Transform rules to binary matrix
  rulevars <- tran_rules(dat = data, rules = rules)
  
  # Combine the rule terms and linear terms
  base.fun <- cbind(rulevars, var.trunc)
  base.names <- c(rules, colnames(var.trunc))
  colnames(base.fun) <- base.names
  
  ## Scale of the base function (Freidman & Popescu, 2008)
  var.sd <- apply(data[,-c(1:2)], 2, sd)
  x.scale <- c(rep(1,ncol(rulevars)), var.sd/0.4)
  
  # Normalized the base functions 
  base.fun <- scale(base.fun, scale = x.scale, center = FALSE) 
  
  lasso.fit <- glmnet::cv.glmnet(y = Y.tilde, x = base.fun, family="mgaussian", nfolds = 5)
  beta.first0 <- predict(lasso.fit, base.fun, type = "coefficients")
  beta.first <- NULL
  for(b in 1:length(beta.first0)){
    beta.first0[[b]][-1] <- beta.first0[[b]][-1]/x.scale
    beta.first <- cbind(beta.first, beta.first0[[b]])
  }
  
  end1 <- (proc.time() - start)[3]
  
  intercept.first <- beta.first[1,] 
  beta.first <- beta.first[-1,]
  
  ### Prediction Using Test data ###
  X.test <- tran_rules(sim.test$data, rules)
  X.test <- cbind(X.test, sim.test$data[,-c(1,2)])
  
  tau.hat.first <- matrix(rep(intercept.first, nrow(sim.test$data)), nrow(sim.test$data), length(unique(sim.test$W)), byrow = TRUE) + as.matrix(X.test)%*%as.matrix(beta.first) 
  tau.hat.first <- tau.hat.first[,-1] - tau.hat.first[,1]
  
  ########################################
  ### Step4 : Rule Ensemble (Adaptive) ###
  ########################################
  
  penalty.factor <- rowMeans(as.matrix(beta.first))
  penalty.factor <- 1/abs(penalty.factor + .Machine$double.eps*2)
  
  adlasso.fit <- glmnet::cv.glmnet(y = Y.tilde, x = base.fun, family="mgaussian", nfolds = 5, penalty.factor = penalty.factor)
  
  beta0 <- predict(adlasso.fit, base.fun, type = "coefficients")
  beta <- NULL
  for(b in 1:length(beta0)){
    beta0[[b]][-1] <- beta0[[b]][-1]/x.scale
    beta <- cbind(beta, beta0[[b]])
  }
  intercept<- beta[1,] 
  beta <- beta[-1,]
  
  end2 <- (proc.time() - start)[3]
  
  tau.hat <- matrix(rep(intercept, nrow(sim.test$data)), nrow(sim.test$data), length(unique(sim.test$W)), byrow = TRUE) + as.matrix(X.test)%*%as.matrix(beta) 
  tau.hat <- tau.hat[,-1] - tau.hat[,1]
  
  ###################
  ### Evaluation  ###   
  ###################
  
  time.first <- end1
  time <- end2
  
  rank.first <- c()
  rank <- c()
  order.can.first <- NULL
  order.can <- NULL
  for(t in 1:nrow(tau.hat)){
    rank.first[t] <- which.max(cbind(0, tau.hat.first)[t,])
    rank[t] <- which.max(cbind(0, tau.hat)[t,])
    order.can.first <- rbind(order.can.first, order(cbind(0, tau.hat.first)[t,], decreasing = TRUE))
    order.can <- rbind(order.can, order(cbind(0, tau.hat)[t,], decreasing = TRUE))
  }
  
  ## Correct Variable Selection Rate ---------------
  r0 <- beta.first[apply(beta.first != 0, 1, sum) == 3,]
  r1 <- beta[apply(beta != 0, 1, sum) == 3,]
  err.var <- paste0("X",6:ncol(X.train), collapse = "|")
  err.var <- paste0("(",err.var,")")
  cor.rate.first <- 1 - length(grep(err.var, rownames(r0)))/length(r0)
  cor.rate <- 1 - length(grep(err.var, rownames(r1)))/length(r1)
  
  ## Number of terms --------------------------------
  term.first <- sum(beta.first != 0)/length(unique(W.train))
  term <- sum(beta != 0)/length(unique(W.train))
  
  res.can <- list(tau.first = tau.hat.first, tau = tau.hat,  
                  rank.first = rank.first, rank = rank,
                  order.first = order.can.first, order = order.can,
                  time.first = time.first, time = time,
                  term.first = term.first, term = term)
  res.can
}

########################## Causal Forest (CF) ##################################

cf <- function(sim.train, sim.test, W.hat, seed = 1){
  
  set.seed(seed)
  
  Y.train <- as.numeric(sim.train$Y)
  W.train <- sim.train$W
  X.train <- sim.train$X
  
  start <- proc.time()
  
  fit <- mcre_agenet(X = X.train, Y = Y.train, W = W.train, W.hat = W.hat, 
                     gen.tree = "cf", tune.init = "CV", tune.fin = "CV", nfolds = 5, seed = seed)
  tau.hat <- predict.mcre(fit, newdata = sim.test$data)
  
  time.first <- fit$time.first
  time <- fit$time 
  
  rank <- c()
  rank.first <- c()
  order.can <- NULL 
  order.can.first <- NULL 
  for(t in 1:nrow(tau.hat$tau)){
    rank[t] <- which.max(cbind(0, tau.hat$tau)[t,])
    rank.first[t] <- which.max(cbind(0, tau.hat$tau.first)[t,])
    order.can <- rbind(order.can, order(cbind(0, tau.hat$tau)[t,], decreasing = TRUE))
    order.can.first <- rbind(order.can.first, order(cbind(0, tau.hat$tau.first)[t,], decreasing = TRUE))
  }
  
  ## Correct Variable Selection Rate ---------------
  r1 <- fit$beta[which(fit$beta != 0)]
  r0 <- fit$beta.first[which(fit$beta.first != 0)]
  err.var <- paste0("X",6:ncol(X.train), collapse = "|")
  err.var <- paste0("(",err.var,")")
  cor.rate.all <- 1 - length(grep(err.var, fit$rules))/length(fit$rules)
  cor.rate <- 1 - length(grep(err.var, names(r1)))/length(r1)
  cor.rate.first <- 1 - length(grep(err.var, names(r0)))/length(r0)
  
  ## Number of terms --------------------------------
  term.all <- length(fit$rules)
  term <- sum(fit$beta != 0)/length(unique(W.train))
  term.first <- sum(fit$beta.first != 0)/length(unique(W.train))
  
  ## Median of interaction -----------------------------
  int <- median(unlist(lapply(strsplit(names(fit$beta)[fit$beta != 0],"&"), length)) - 1)
  int.first <- median(unlist(lapply(strsplit(names(fit$beta.first)[fit$beta.first != 0],"&"), length)) - 1)
  
  res.can <- list(tau = tau.hat$tau, tau.first = tau.hat$tau.first,
                  rank = rank, rank.first = rank.first,
                  order = order.can, order.first = order.can.first,
                  time = time, time.first = time.first,
                  term = term, term.first = term.first, term.all = term.all,
                  cor.rate = cor.rate, cor.rate.first = cor.rate.first, cor.rate.all = cor.rate.all,
                  int = int, int.first = int.first)
  res.can
}

