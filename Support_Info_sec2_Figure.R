
dir <- getwd()
# 
n <- 1e5
p <- 10
Beta <- 2*rbind(c(1, 1, 1), c(-0.5, 1.0, 2.0), c(1.5, 1.5, -0.5), c(-1.5, 1.5, 0.5), c(-0.5, 2.0, 0.5))[1:4,]
prob <- "obv"


seed <- 1
set.seed(seed)

### Sub: Probabilities generation using multinomial logistic regression

prob_gen <- function(X, group_num, prob = "rct", overlap = "strong"){
  
  if(prob == "rct"){
    
    probs <- matrix(1/group_num, n, group_num)
    
  }else if(prob == "obv"){
    
    if(overlap == "strong"){
      
      eta <- 1
      
    }else if(overlap == "weak"){
      
      if(group_num == 3){
        eta <- 7
      }else if(group_num == 4){
        eta <- 5
      }else if(group_num == 5){
        eta <- 3
      }
      
    }else if(overlap == "moderate"){
      
      if(group_num == 3){
        eta <- 3.5
      }else if(group_num == 4){
        eta <- 2.5
      }else if(group_num == 5){
        eta <- 1.5
      }
    }
    
    # Relative Risk from control
    ex1 <- exp(-(0.50*eta + 0.1*X[,1] + 0.2*X[,2] + 0.3*X[,3] - 0.2*X[,4] + 0.7*X[,5]))
    ex2 <- exp(-(0.75*eta + 0.2*X[,1] + 0.4*X[,2] + 0.6*X[,3] - 0.4*X[,4] + 0.3*X[,5]))
    ex3 <- exp(-(1.00*eta + 0.2*X[,1] + 0.5*X[,2] + 0.5*X[,3] - 0.5*X[,4] + 0.3*X[,5]))
    ex4 <- exp(-(1.50*eta + 0.3*X[,1] + 0.4*X[,2] + 0.2*X[,3] - 0.4*X[,4] + 0.1*X[,5]))
    ex <- cbind(ex1, ex2, ex3, ex4)
    
    # TRUE GPS
    prob_t <- ex[,1:(group_num - 1)]/(1 + apply(ex[,1:(group_num - 1)], 1, sum))
    prob_c <- 1 - apply(prob_t, 1, sum)
    probs <- data.frame(cbind(prob_c,prob_t))
    colnames(probs) <- paste0("Group",1:group_num)
  }
  probs
}

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

# ### Plot Distribute (N = 1e5) ### 

probs.melt.all <- NULL
overlap.can <- c("strong","moderate", "weak")
for(op in overlap.can){
for(gg in 3:5){

  probs <- prob_gen(X, group_num = nrow(Beta.can[1:gg,]), prob = prob, overlap = op)
  probs.melt <- reshape2::melt(probs)
  probs.melt$g <- paste0("Number of groups = ", gg)
  probs.melt$op <- op
  print(op)
  probs.melt.all <- rbind(probs.melt.all, probs.melt)
}
}
colnames(probs.melt.all)[1:2] <- c("Groups", "GPS")
probs.melt.all$op <- factor(probs.melt.all$op, levels = overlap.can)
g <- ggplot(probs.melt.all, aes(x=GPS, fill=Groups)) + geom_density(alpha=.3) + facet_grid(op~g, scales = "free")
g <- g + theme_bw() + theme(legend.position = "bottom", legend.direction = "horizontal", text=element_text(size= 20))
g <- g + xlab("Generalized propensity score")

save <- paste0(dir, "/Table_Figure_support/sup_section2/s-Fig17.pdf")
ggsave(save, g, width = 16, height = 7)
