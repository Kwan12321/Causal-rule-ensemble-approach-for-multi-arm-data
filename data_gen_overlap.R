##### DGP function ###
data_gen_overlap <- function(n, # Sample size
                     p, # Number of variables
                     Beta,  # Coefficient matrix
                     model, # Data generation model
                     rho = 0, # Correlation between continuous variables
                     bin = "b", # balance (b) and unbalance (ub) of binary valuables
                     overlap = "moderate",
                     seed = 1
){
  
  # n <- 1e4
  # p <- 10
  # Beta <- 2*rbind(c(1, 1, 1), c(-0.5, 1.0, 2.0), c(1.5, 1.5, -0.5), c(-1.5, 1.5, 0.5), c(-0.5, 2.0, 0.5))
  # rho <- 0
  # bin <- "b"
  # model <- c(1,2)
  # seed <- 1
  
  # Average for normal & prob for binomial
  if(overlap == "weak"){
    avg.can <- c(-1, 0, 1, 2, 3)
    prob.can <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  }else if(overlap == "moderate"){
    avg.can <- c(-0.5, 0, 0.5, 1, 1.5)
    prob.can <- c(0.4, 0.45, 0.5, 0.55, 0.6)
  }else if(overlap == "strong"){
    avg.can <- rep(0,5)
    prob.can <- rep(0.5,5)
  }  
    #avg.can <- rep(0,5)
    #prob.can <- c(0.5, 0.3, 0.4, 0.6, 0.7)
    #prob.can <- rep(0.5,5)
  
  # > apply(treat%*%t(Beta), 2, mean)
  # [1] 0.2874069 1.5086238 1.7418367 0.4871168 0.9876845
  
  treat_can <- NULL
  main_can <- NULL
  X_can <- list()
  
  #X_all <- NULL
  
  for(g in seq_len(nrow(Beta))){
  
    ## Covariates (odds: continuous; even: binary)
    set.seed(seed)
    X <- pnorm(MASS::mvrnorm(n, rep(0, p), stats::toeplitz(rho^seq(0, p - 1))))
    
    ## Different distribution for each group 
    X.cont <- matrix(rnorm(n*p, avg.can[g], 1), n, p)
    X.bin <- matrix(rbinom(n*p, 1, prob.can[g]), n, p)
    
    # ## weak, moderate & strong overlap settings
    # if(overlap == "weak"){
    #   X.cont0 <- X.cont
    #   X.bin0 <- X.bin
    # }else if(overlap == "moderate"){
    #   X.cont0[,c(1,3)] <- X.cont[,c(1,3)]
    #   X.bin0[,c(2,4)] <- X.bin[,c(2,4)]
    # }else if(overlap == "strong"){
    #   X.cont0 <- X.cont0
    #   X.bin0 <- X.bin0
    # }
      
    X[,(1:p)%%2==1] <- X.cont[,(1:p)%%2==1]
    X[,(1:p)%%2==0] <- X.bin[,(1:p)%%2==0]
    X <- data.frame(X)
    
    #### The main effect function ####
    if(model[1] == 1){
      main <- X[,1] + 1.5*X[,2] + X[,3] - 1.5*X[,4] + X[,5]
    }else if(model[1] == 2){  
      main <- 2*I(X[,1] > -1)*I(X[,3] < 1) - 2*I(X[,2] < 0.5) - 2*I(X[,3] > -1)*I(X[,5] < 1) + 2*I(X[,4] > 0.5)
    }else if(model[1] == 3){  
      main <- X[,1]^2 + 0.3*X[,2] - 2.3*X[,4] + X[,3]*X[,5]
    }
    
    #### The treatment effect function ####
    if(model[2] == 1){
      treat <- matrix(1,n,3)
    }else if(model[2] == 2){
      treat <- cbind((0.5*X[,1] + X[,2]),(0.5*X[,3] + X[,4]),(0.5*X[,5] + X[,2]))
    }else if(model[2] == 3){  
      treat <- cbind(0.7*(I(X[,1]>0) + I(X[,2]>0.5)), 0.7*(I(X[,3]>0) + I(X[,4]>0.5)), 0.7*(I(X[,5]>0) + I(X[,2]>0.5)))
    }else if(model[2] == 4){  
      treat <- cbind(0.75*sin(X[,1]) + X[,2], 0.75*sin(X[,3]) + X[,4], 0.75*sin(X[,5]) + X[,2]) 
    }
    
    treat_can <- cbind(treat_can, treat%*%Beta[g,])
    main_can <- cbind(main_can, main)
    X_can[[g]] <- X 
    
    #X_all <- cbind(X_all, X[,1])
  }
  
  Y_can <- main_can + treat_can
  
  
  # df <- data.frame(X_all)
  # df <- reshape2::melt(df)
  # 
  # g <- ggplot(df, aes(x=value, fill=variable)) + geom_density(alpha=.3)
  # 
  
  #### Create group indicators ####
  if(nrow(Beta) == 3){
    probs <- c(0.48, 0.29, 0.23)
  }else if(nrow(Beta) == 4){
    probs <- c(0.41, 0.25, 0.19, 0.15)
  }else if(nrow(Beta) == 5){
    probs <- c(0.38, 0.23, 0.18, 0.13, 0.08)
  }
  
  W <- numeric(n)
  W_dummy <- matrix(0,nrow(Beta), n)
  Y0 <- c(); W0 <- c(); X0 <- matrix(0, n, p)
  
  set.seed(seed)
  for(i in 1:n){
    W_dummy <- rmultinom(1,1,probs)
    W0[i] <- as.numeric(c(0:(nrow(Beta) - 1))%*%W_dummy)
    g.id <- W0[i] + 1
    X0[i,] <- as.matrix(X_can[[g.id]][i,,drop = FALSE])
    Y0[i] <- rnorm(1, Y_can[i,g.id], 1)
  }
  data <- data.frame(Y = Y0, W = W0, X0)

  #### True treatment effect #### 
  tau <- treat_can[,-1] - treat_can[,1]
  
  rank <- c()
  order.can <- NULL
  for(t in 1:nrow(treat_can)){
    rank[t] <- which.max(treat_can[t,])
    order.can <- rbind(order.can, order(treat_can[t,], decreasing = TRUE))
  }
  
  return(list(data = data,Y = Y0, W = W0, X = data.frame(X0), tau = tau, rank = rank, order = order.can))
}

################################################################################

# n <- 2000
# p <- 10
# Beta <- rbind(c(2.5,-2.5,0.5), c(-0.5, 1.0, 2.5), c(1.5, 2.5, -0.5), c(-1.0, 1.5, 0.5), c(-0.5, 2.0, 0.5))[1:3,]
# prob <- "obv"
# rho <- 0
# bin <- "b"
# model <- c(1,2)
# seed <- 1
# 
# train <- data_gen_overlap(n, p, Beta[1:3,],  model, rho = 0, bin = "b", overlap = "weak", seed = 1)
# test <- data_gen_overlap(n, p, Beta[1:3,],  model, rho = 0, bin = "b", overlap = "weak",seed = 11)
# 
# ## Inverse weights ###
# PS.fit <- nnet::multinom(W ~., data = data.frame(W = as.factor(train$W), train$X))
# W.hat <- predict(PS.fit, data = data.frame(sim.test$X), type = "probs")
# ## trimming the propensity score ()
# r_upper <- apply(W.hat, 2, quantile, 0.9)
# r_lower <- apply(W.hat, 2, quantile, 0.1)
# 
# for(c in 1:length(unique(train$W))){
#   W.hat[,c][W.hat[,c] >= r_upper[c]] <- r_upper[c]
#   W.hat[,c][W.hat[,c] <= r_lower[c]] <- r_lower[c]
# }
# 
# fit0 <- mbart(train, test, 1/W.hat)
# cc0 <- c()
# for(i in 1:2000){
#   cc0[i] <- cor(fit0$order[i,], test$order[i,])
# }
# ## Entropy Balancing weights ###
# d1 <- train$data
# d1[,2] <- as.factor(d1[,2])
# W1 <- WeightIt::weightit(W ~., data = d1[,-1], method = "cbps", estimand = "ATE")
# dum1 <- fastDummies::dummy_cols(d1[,2])[,-1]
# IPW1 <- W1$weights*dum1
# 
# ## trimming the propensity score ()
# r_upper <- apply(IPW1, 2, quantile, 0.9)
# r_lower <- apply(IPW1, 2, quantile, 0.1)
# 
# for(c in 1:length(unique(train$W))){
#   IPW1[,c][IPW1[,c] >= r_upper[c]] <- r_upper[c]
#   IPW1[,c][IPW1[,c] <= r_lower[c]] <- r_lower[c]
# }
# 
# 
# fit1 <- mbart(train, test, IPW1)
# # cor(fit1$tau, test$tau)
# # apply((fit1$tau - test$tau), 2, mean)
# # apply((fit1$tau - test$tau)^2, 2, mean)
# 
# cc1 <- c()
# for(i in 1:2000){
#   cc1[i] <- cor(fit1$order[i,], test$order[i,])
# }
# 
# 
# ## Stable Balancing weights ###
# d2 <- train$data
# d2[,2] <- as.factor(d2[,2])
# W2 <- WeightIt::weightit(W ~., data = d2[,-1], method = "optweight", estimand = "ATE")
# dum2 <- fastDummies::dummy_cols(d2[,2])[,-1]
# IPW2 <- W2$weights*dum2
# 
# ## trimming the propensity score ()
# r_upper <- apply(IPW2, 2, quantile, 0.9)
# r_lower <- apply(IPW2, 2, quantile, 0.1)
# 
# for(c in 1:length(unique(train$W))){
#   IPW2[,c][IPW2[,c] >= r_upper[c]] <- r_upper[c]
#   IPW2[,c][IPW2[,c] <= r_lower[c]] <- r_lower[c]
# }
# 
# fit2 <- mbart(train, test, IPW2)
# 
# cor(fit2$tau, test$tau)
# apply((fit2$tau - test$tau), 2, mean)
# apply((fit2$tau - test$tau)^2, 2, mean)
# 
# 
# cc2 <- c()
# for(i in 1:2000){
#   cc2[i] <- cor(fit2$order[i,], test$order[i,])
# }
# 
# ## Optweight Balancing weights ###
# d3 <- train$data
# d3[,2] <- as.factor(d3[,2])
# W3 <- WeightIt::weightit(W ~., data = d3[,-1], method = "glm", estimand = "ATE")
# dum3 <- fastDummies::dummy_cols(d3[,2])[,-1]
# IPW3 <- W3$weights*dum3
# 
# ## trimming the propensity score ()
# r_upper <- apply(IPW3, 2, quantile, 0.9)
# r_lower <- apply(IPW3, 2, quantile, 0.1)
# 
# for(c in 1:length(unique(train$W))){
#   IPW3[,c][IPW3[,c] >= r_upper[c]] <- r_upper[c]
#   IPW3[,c][IPW3[,c] <= r_lower[c]] <- r_lower[c]
# }
# 
# fit3 <- mbart(train, test, IPW3)
# 
# cc3 <- c()
# for(i in 1:2000){
#   cc3[i] <- cor(fit3$order[i,], test$order[i,])
# }
# 
# cor(fit3$tau, test$tau)
# apply((fit3$tau - test$tau), 2, mean)
# apply((fit3$tau - test$tau)^2, 2, mean)
# 
# 
# #########################################################
# 
# ## less or moderate overlap ##




