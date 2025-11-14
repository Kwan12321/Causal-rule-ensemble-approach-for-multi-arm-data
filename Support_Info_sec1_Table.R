library(gridExtra)
library(ggplot2)
library(greekLetters)
library(tidyr)
library(ggpubr)
library(dplyr)
#########################

dir <- getwd()
### DGP parameters ###
n <- 2000
p <- 10
model.can <- expand.grid(1:3, 2:4)
model.can <- cbind(sim = 1:nrow(model.can), model.can)
prob.can <- c("rct", "obv")
group.can <- c(3,4,5)
ite.all <- 100

model.all <- NULL
for(g in group.can){
  for(prob in prob.can){
    model.set <- cbind(model.can, prob, g) 
    model.all <- rbind(model.all, model.set)
  }
}

Fun1.can <- c("Linear", "Stepwise", "Non-linear")
Fun2.can <- c("", "Linear", "Stepwise", "Non-linear")
Fun1.can <- paste(Fun1.can, greeks("mu"))
Fun2.can <- paste(Fun2.can, greeks("delta"))
Fun1.can <- factor(Fun1.can, levels = Fun1.can)
Fun2.can <- factor(Fun2.can, levels = Fun2.can)

mPEHE.can3 <- NULL
mbias.can3 <- NULL
kappa.can3 <- NULL
mcor_ord.can3 <- NULL
term.can3 <- NULL

for(method in c("bart","forest","boost", "prop")){
  
  mPEHE.list <- c()
  mbias.list <- c()
  kappa.list <- c()
  mcord.list <- c()
  
  mPEHE.can1 <- NULL
  mbias.can1 <- NULL
  kappa.can1 <- NULL
  mcor_ord.can1 <- NULL
  term.can1 <- NULL
  
  mPEHE.can2 <- NULL
  mbias.can2 <- NULL
  kappa.can2 <- NULL
  mcor_ord.can2 <- NULL
  term.can2 <- NULL
  
  
  for(prob in c("rct", "obv")){
    for(g in 3:5){
      for(s in 1:nrow(model.can)){
        
        mPEHE.can1 <- NULL
        mbias.can1 <- NULL
        kappa.can1 <- NULL
        mcor_ord.can1 <- NULL
        term.can1 <- NULL
        
        # if(method != "prop"){
        #   dir <- "Z:/MCRE/"
        #   save.res <- paste0(dir, method,"/",method,"_evl_sim",s,"_n",n,"_p",p,"_",prob,"_",g,"g",".csv")
        # }else if(method == "prop"){  
        #   dir <- "Z:/m-cre/"
        #   save.res <- paste0(dir, "sim.res/",method,"_evl_sim",s,"_n",n,"_p",p,"_",prob,"_",g,"g",".csv")
        # }  
        
        save.res <- paste0(dir, "/main_simulation/",method,"/",method,"_evl_sim",s,"_n",n,"_p",p,"_",prob,"_",g,"g",".csv")
        dat <- read.csv(save.res)
        
        mPEHE.can0 <- NULL
        mbias.can0 <- NULL
        kappa.can0 <- NULL
        mcor_ord.can0 <- NULL
        term.can0 <- NULL
        
        for(ite in 1:ite.all){
          if(method != "prop"){
            mPEHE <- dat[dat$ite == ite,][1, 1:7]
            mbias <- dat[dat$ite == ite,][2, 1:7]
            kappa <- dat[dat$ite == ite,][4, 1:7]
            mcor_ord <- dat[dat$ite == ite,][5, 1:7]
            
          }else{
            
            colnames(dat)[1:12] <- c("gbm.gl","gbm.agl","ctree.gl","ctree.agl", "gbm.gl.b","gbm.agl.b","ctree.gl.b","ctree.agl.b", "gbm.l.dr","gbm.al.dr","ctree.l.dr","ctree.al.dr")
            mPEHE <- dat[dat$ite == ite,][1, c(1:12)]
            mbias <- dat[dat$ite == ite,][2, c(1:12)]
            kappa <- dat[dat$ite == ite,][4, c(1:12)]
            mcor_ord <- dat[dat$ite == ite,][5, c(1:12)]
            mterms <- dat[dat$ite == ite,][6, c(1:12)]
            
            term.can0 <- rbind(term.can0, mterms)
          }  
          mPEHE.can0 <- rbind(mPEHE.can0, mPEHE)
          mbias.can0 <- rbind(mbias.can0, mbias)
          kappa.can0 <- rbind(kappa.can0, kappa)
          mcor_ord.can0 <- rbind(mcor_ord.can0, mcor_ord)
        }
        mPEHE.can1 <- cbind(mPEHE.can1, as.matrix(mPEHE.can0))
        mbias.can1 <- cbind(mbias.can1, as.matrix(mbias.can0))
        kappa.can1 <- cbind(kappa.can1, as.matrix(kappa.can0))
        mcor_ord.can1 <- cbind(mcor_ord.can1, as.matrix(mcor_ord.can0))
        # if(method == "prop"){
        #   term.can1 <- cbind(term.can1, as.matrix(term.can0))
        # }  
        
        mPEHE.can1 <- data.frame(mPEHE.can1)
        mbias.can1 <- data.frame(mbias.can1)
        kappa.can1 <- data.frame(kappa.can1)
        mcor_ord.can1 <- data.frame(mcor_ord.can1)
        term.can1 <- data.frame(term.can1)
        
        mPEHE.can1$fun1 <- Fun1.can[model.can[s,2]]; mPEHE.can1$fun2 <- Fun2.can[model.can[s,3]]; mPEHE.can1$g <- g
        mbias.can1$fun1 <- Fun1.can[model.can[s,2]]; mbias.can1$fun2 <- Fun2.can[model.can[s,3]]; mbias.can1$g <- g
        kappa.can1$fun1 <- Fun1.can[model.can[s,2]]; kappa.can1$fun2 <- Fun2.can[model.can[s,3]]; kappa.can1$g <- g
        mcor_ord.can1$fun1 <- Fun1.can[model.can[s,2]]; mcor_ord.can1$fun2 <- Fun2.can[model.can[s,3]]; mcor_ord.can1$g <- g
        #term.can1$fun1 <- Fun1.can[model.can[s,2]]; term.can1$fun2 <- Fun2.can[model.can[s,3]]; term.can1$g <- g
        
        mPEHE.can1$prob <- prob
        mbias.can1$prob <- prob
        kappa.can1$prob <- prob
        mcor_ord.can1$prob <- prob
        #term.can1$prob <- prob
        
        mPEHE.can2 <- rbind(mPEHE.can2, mPEHE.can1)
        mbias.can2 <- rbind(mbias.can2, mbias.can1)
        kappa.can2 <- rbind(kappa.can2, kappa.can1)
        mcor_ord.can2 <- rbind(mcor_ord.can2, mcor_ord.can1)
        #term.can2 <- rbind(term.can2, term.can1)
      }
    }
  }
  
  mPEHE.can_2 <- reshape2::melt(mPEHE.can2, id = c("fun1","fun2","g","prob"))
  mPEHE.can_2$value <- as.numeric(mPEHE.can_2$value)
  mPEHE.can_2 <- mPEHE.can_2 %>% group_by(prob, g, fun1, fun2, variable) %>% summarise_all(list(mean = mean, sd = sd))
  mPEHE_min <- mPEHE.can_2 %>% group_by(prob, g, fun1, fun2) %>% summarise(min = min(mean), .groups = "drop")
  mPEHE_winners <- mPEHE.can_2 %>% group_by(prob, g, fun1, fun2) %>% slice_min(mean, n = 1, with_ties = TRUE) %>%
    mutate(is_best = TRUE) %>% select(prob, g, fun1, fun2, variable, is_best)
  mPEHE.can_2 <- mPEHE.can_2 %>% left_join(mPEHE_winners, by = c("prob","g","fun1","fun2","variable")) %>% 
    mutate(is_best = replace_na(is_best, FALSE))%>%filter(is_best)
  
  
  mbias.can_2 <- reshape2::melt(mbias.can2, id = c("fun1","fun2","g","prob"))
  mbias.can_2$value <- as.numeric(mbias.can_2$value)
  mbias.can_2 <- mbias.can_2 %>% group_by(prob, g, fun1, fun2, variable) %>% summarise_all(list(mean = mean, sd = sd))
  mbias_min <- mbias.can_2 %>% group_by(prob, g, fun1, fun2) %>% summarise(min = min(mean), .groups = "drop")
  mbias_winners <- mbias.can_2 %>% group_by(prob, g, fun1, fun2) %>% slice_min(mean, n = 1, with_ties = TRUE) %>%
    mutate(is_best = TRUE) %>% select(prob, g, fun1, fun2, variable, is_best)
  mbias.can_2 <- mbias.can_2 %>% left_join(mbias_winners, by = c("prob","g","fun1","fun2","variable")) %>% 
    mutate(is_best = replace_na(is_best, FALSE))%>%filter(is_best)
  
  
  kappa.can_2 <- reshape2::melt(kappa.can2, id = c("fun1","fun2","g","prob"))
  kappa.can_2$value <- as.numeric(kappa.can_2$value)
  kappa.can_2 <- kappa.can_2 %>% group_by(prob, g, fun1, fun2, variable) %>% summarise_all(list(mean = mean, sd = sd))
  kappa_min <- kappa.can_2 %>% group_by(prob, g, fun1, fun2) %>% summarise(min = max(mean), .groups = "drop")
  kappa_winners <- kappa.can_2 %>% group_by(prob, g, fun1, fun2) %>% slice_max(mean, n = 1, with_ties = TRUE) %>%
    mutate(is_best = TRUE) %>% select(prob, g, fun1, fun2, variable, is_best)
  kappa.can_2 <- kappa.can_2 %>% left_join(kappa_winners, by = c("prob","g","fun1","fun2","variable")) %>% 
    mutate(is_best = replace_na(is_best, FALSE))%>%filter(is_best)
  
  
  mcor_ord.can_2 <- reshape2::melt(mcor_ord.can2, id = c("fun1","fun2","g","prob"))
  mcor_ord.can_2$value <- as.numeric(mcor_ord.can_2$value)
  mcor_ord.can_2 <- mcor_ord.can_2 %>% group_by(prob, g, fun1, fun2, variable) %>% summarise_all(list(mean = mean, sd = sd))
  mcor_ord_min <- mcor_ord.can_2 %>% group_by(prob, g, fun1, fun2) %>% summarise(min = max(mean), .groups = "drop")
  mcor_ord_winners <- mcor_ord.can_2 %>% group_by(prob, g, fun1, fun2) %>% slice_max(mean, n = 1, with_ties = TRUE) %>%
    mutate(is_best = TRUE) %>% select(prob, g, fun1, fun2, variable, is_best)
  mcor_ord.can_2 <- mcor_ord.can_2 %>% left_join(mcor_ord_winners, by = c("prob","g","fun1","fun2","variable")) %>% 
    mutate(is_best = replace_na(is_best, FALSE))%>%filter(is_best)
  
  
  mPEHE.can_2$mean <- format(round(mPEHE.can_2$mean, 2), nsmall = 3); mbias.can_2$mean <- format(round(mbias.can_2$mean, 2), nsmall = 3)
  kappa.can_2$mean <- format(round(kappa.can_2$mean, 2), nsmall = 3); mcor_ord.can_2$mean <- format(round(mcor_ord.can_2$mean, 2), nsmall = 3)
  
  mPEHE.can_2$sd <- format(round(mPEHE.can_2$sd, 2), nsmall = 3); mbias.can_2$sd <- format(round(mbias.can_2$sd, 2), nsmall = 3)
  kappa.can_2$sd <- format(round(kappa.can_2$sd, 2), nsmall = 3); mcor_ord.can_2$sd <- format(round(mcor_ord.can_2$sd, 2), nsmall = 3)
  
  scenario <- mPEHE.can_2[,1:4]
  
  mPEHE.can3 <- cbind(mPEHE.can3, as.matrix(mPEHE.can_2[,5:7]))
  mbias.can3 <- cbind(mbias.can3, as.matrix(mbias.can_2[,5:7]))
  kappa.can3 <- cbind(kappa.can3, as.matrix(kappa.can_2[,5:7]))
  mcor_ord.can3 <- cbind(mcor_ord.can3, as.matrix(mcor_ord.can_2[,5:7]))

}    
  
mPEHE.all <- cbind(scenario, mPEHE.can3); mbias.all <- cbind(scenario, mbias.can3)
kappa.all <- cbind(scenario, kappa.can3); mcor_ord.all <- cbind(scenario, mcor_ord.can3)
colnames(mPEHE.all) <- colnames(mbias.all) <- colnames(kappa.all) <- colnames(mcor_ord.all) <- c("prob","Number of /n groups", "Main effect /n function", "Treatment effect /n function","Method", "Mean", "SD", "Method", "Mean", "SD", "Method", "Mean", "SD", "Method", "Mean", "SD")

for(prob in c("rct", "obv")){
  write.csv(mPEHE.all[mPEHE.all$prob == prob,][,-1], paste0(dir,"/Table_Figure_support/sup_section1/Table1_",prob,".csv"), fileEncoding = "CP932", row.names = FALSE)
  write.csv(mbias.all[mbias.all$prob == prob,][,-1], paste0(dir,"/Table_Figure_support/sup_section1/Table2_",prob,".csv"), fileEncoding = "CP932", row.names = FALSE)
  write.csv(kappa.all[kappa.all$prob == prob,][,-1], paste0(dir,"/Table_Figure_support/sup_section1/Table3_",prob,".csv"), fileEncoding = "CP932", row.names = FALSE)
  write.csv(mcor_ord.all[mcor_ord.all$prob == prob,][,-1], paste0(dir,"/Table_Figure_support/sup_section1/Table4_",prob,".csv"), fileEncoding = "CP932", row.names = FALSE)
}
        