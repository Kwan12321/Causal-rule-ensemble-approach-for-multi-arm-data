##### DGP function ###
data_gen <- function(n, # Sample size
                     p, # Number of variables
                     Beta,  # Coefficient matrix
                     model, # Data generation model
                     prob = "rct", # Probability for each group
                     rho = 0, # Correlation between continuous variables
                     bin = "b", # balance (b) and unbalance (ub) of binary valuables
                     seed = 1
){
  
  #n <- 2000
  #p <- 10
  #Beta <- rbind(c(2.5,-2.5,0.5), c(-0.5, 1.0, 2.5), c(1.5, 2.5, -0.5), c(-1.0, 1.5, 0.5), c(-0.5, 2.0, 0.5))
  #Beta <- rbind(c(1, 1, 1), c(-0.5, 1.0, 2.0), c(1.5, 1.5, -0.5), c(-1.5, 1.5, 0.5), c(-0.5, 2.0, 0.5))
  #prob <- "obv"
  #rho <- 0
  #bin <- "b"
  #model <- c(3,3)
  #seed <- 1
  
  set.seed(seed)
  
  # Covariates (odds: continuous; even: binary)
  X <- pnorm(MASS::mvrnorm(n, rep(0, p), stats::toeplitz(rho^seq(0, p - 1))))
  X.cont <- matrix(rnorm(n * p), n, p)
  if(bin == "b"){
    X.bin <- matrix(rbinom(n*p, 1, 0.5), n, p)
  }else if(bin == "ub"){  
    X.bin <- t(matrix(rbinom(p*n, 1, sample(seq(0.1,0.9,0.1), p, replace = TRUE)), p, n))
  }
  
  X[,(1:p)%%2==1] <- X.cont[,(1:p)%%2==1]
  X[,(1:p)%%2==0] <- X.bin[,(1:p)%%2==0]
  X <- data.frame(X)
  
  #### The main effect function ####
  if(model[1] == 1){
    main <- 0.6*X[,1] + 0.9*X[,2] + 0.6*X[,3] - 0.9*X[,4] + 0.6*X[,5]
  }else if(model[1] == 2){  
    main <- 1.2*I(X[,1] > -1)*I(X[,3] < 1) - 1.2*I(X[,2] < 0.5) - 1.2*I(X[,3] > -1)*I(X[,5] < 1) + 1.2*I(X[,4] > 0.5)
  }else if(model[1] == 3){  
    main <- 0.6*X[,1]^2 + 0.5*X[,2]*X[,3] - 1.2*cos(pi*X[,4]*X[,5])
  }
  
  #### The treatment effect function ####
  if(model[2] == 1){
    treat <- matrix(0.5,n,3)
  }else if(model[2] == 2){
    treat <- cbind((0.5*X[,1] + X[,2]),(0.5*X[,3] + X[,4]),(0.5*X[,5] + X[,2]))
  }else if(model[2] == 3){  
    treat <- cbind(1.4*I(X[,1]> 0) - 0.3*I(X[,2]>0.5), 1.4*I(X[,3]> 0) - 0.3*I(X[,4]>0.5), 1.4*I(X[,5]>0) - 0.3*I(X[,2]>0.5))
  }else if(model[2] == 4){  
    treat <- cbind(0.75*sin(X[,1]) + X[,2], 0.75*sin(X[,3]) + X[,4], 0.75*sin(X[,5]) + X[,2]) 
  }
  
  #### True treatment effect #### 
  if(model[2] != 5){
    treat_can <- treat%*%t(Beta)
  }else if(model[2] == 5){
    Beta.can <- matrix(0, 5, 3)
    Beta.can[1:nrow(Beta),] <- Beta
    treat_can1 <- cbind((0.5*X[,1] + X[,2]),(0.5*X[,3] + X[,4]),(0.5*X[,5] + X[,2]))%*%t(Beta.can)[,1]
    treat_can2 <- cbind(1.4*I(X[,1]> 0) - 0.3*I(X[,2]>0.5), 1.4*I(X[,3]> 0) - 0.3*I(X[,4]>0.5), 1.4*I(X[,5]>0) - 0.3*I(X[,2]>0.5))%*%t(Beta.can)[,2]
    treat_can3 <- cbind(0.75*sin(X[,1]) + X[,2], 0.75*sin(X[,3]) + X[,4], 0.75*sin(X[,5]) + X[,2])%*%t(Beta.can)[,3]
    treat_can4 <- cbind(1.4*I(X[,1]> 0) - 0.3*I(X[,2]>0.5), 1.4*I(X[,3]> 0) - 0.3*I(X[,4]>0.5), 1.4*I(X[,5]>0) - 0.3*I(X[,2]>0.5))%*%t(Beta.can)[,4]
    treat_can5 <- cbind(0.75*sin(X[,1]) + X[,2], 0.75*sin(X[,3]) + X[,4], 0.75*sin(X[,5]) + X[,2])%*%t(Beta.can)[,5]
    treat_can <- cbind(treat_can1, treat_can2,treat_can3,treat_can4,treat_can5)[,1:nrow(Beta)]
  }
  tau <- treat_can[,-1] - treat_can[,1]
  
  rank <- c()
  order.can <- NULL
  for(t in 1:nrow(treat_can)){
    rank[t] <- which.max(treat_can[t,])
    order.can <- rbind(order.can, order(treat_can[t,], decreasing = TRUE))
  }
  
  #### Create group indicators #### 
  probs <- prob_gen(X, group_num = nrow(Beta), prob = prob)
  W <- numeric(n)
  W_dummy <- matrix(0,nrow(Beta), n)
  for(i in 1:n){
    W_dummy[,i] <- rmultinom(1,1,probs[i,])
    W[i] <- as.numeric(c(0:(nrow(Beta) - 1))%*%W_dummy[,i])
  }
  
  #### Outcome model ####
  Y <- main + rowSums(treat_can*t(W_dummy))
  Y <- rnorm(n, Y, 1)
  print(paste0("SNR = ", sd(main + rowSums(treat_can*t(W_dummy))), ":1"))
  
  #### simulation dataset ####
  data <- data.frame(Y = Y, W = W, X) 
  
  return(list(data = data,Y = Y, W = W, X = X, tau = tau, rank = rank, order = order.can, gps = probs))
}

### Sub: Probabilities generation using multinomial logistic regression

prob_gen <- function(X, group_num, prob = "rct"){
  
  if(prob == "rct"){
    probs <- matrix(1/group_num, n, group_num)
    
  }else if(prob == "obv"){
    
    # Relative Risk from control
    ex1 <- exp(-(0.50 + 0.1*X[,1] + 0.2*X[,2] + 0.3*X[,3] - 0.2*X[,4] + 0.7*X[,5]))
    ex2 <- exp(-(0.75 + 0.2*X[,1] + 0.4*X[,2] + 0.6*X[,3] - 0.4*X[,4] + 0.3*X[,5]))
    ex3 <- exp(-(1.00 + 0.2*X[,1] + 0.5*X[,2] + 0.5*X[,3] - 0.5*X[,4] + 0.3*X[,5]))
    ex4 <- exp(-(1.50 + 0.3*X[,1] + 0.4*X[,2] + 0.2*X[,3] - 0.4*X[,4] + 0.1*X[,5]))
    ex <- cbind(ex1, ex2, ex3, ex4)
    
    # TRUE GPS
    prob_t <- ex[,1:(group_num - 1)]/(1 + apply(ex[,1:(group_num - 1)], 1, sum))
    prob_c <- 1 - apply(prob_t, 1, sum)
    probs <- data.frame(cbind(prob_c,prob_t))
    colnames(probs) <- paste0("Group",1:group_num)
  }
  probs
}


