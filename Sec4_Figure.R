library(gridExtra)
library(ggplot2)
library(tidyverse)
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

Fun1.can <- c("Linear~mu", "Stepwise~mu", "Non-linear~mu")
Fun2.can <- c("", "Linear~delta", "Stepwise~delta", "Non-linear~delta")
Fun1.can <- factor(Fun1.can, levels = Fun1.can)
Fun2.can <- factor(Fun2.can, levels = Fun2.can)

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
      
      for(method in c("bart","prop")){
        save.res <- paste0(dir, "/main_simulation/", method,"/",method,"_evl_sim",s,"_n",n,"_p",p,"_",prob,"_",g,"g",".csv")
        dat <- read.csv(save.res)
        
        mPEHE.can0 <- NULL
        mbias.can0 <- NULL
        kappa.can0 <- NULL
        mcor_ord.can0 <- NULL
        term.can0 <- NULL
        
        for(ite in 1:ite.all){
          if(method != "prop"){
            mPEHE <- dat[dat$ite == ite,][1, c(1,3,7)]
            mbias <- dat[dat$ite == ite,][2, c(1,3,7)]
            kappa <- dat[dat$ite == ite,][4, c(1,3,7)]
            mcor_ord <- dat[dat$ite == ite,][5, c(1,3,7)]
            
            # mPEHE <- dat[dat$ite == ite,][1, c(1,2,3,5)]
            # mbias <- dat[dat$ite == ite,][2, c(1,2,3,5)]
            # kappa <- dat[dat$ite == ite,][4, c(1,2,3,5)]
            # mcor_ord <- dat[dat$ite == ite,][5, c(1,2,3,5)]
            
          }else{
            colnames(dat)[1:4] <- c("gbm.gl","gbm.agl","ctree.gl","ctree.agl")
            mPEHE <- dat[dat$ite == ite,][1, c(1:4)]
            mbias <- dat[dat$ite == ite,][2, c(1:4)]
            kappa <- dat[dat$ite == ite,][4, c(1:4)]
            mcor_ord <- dat[dat$ite == ite,][5, c(1:4)]
            mterms <- dat[dat$ite == ite,][6, c(1:4)]
            
            # mPEHE <- dat[dat$ite == ite,][1, c(1,7)]
            # mbias <- dat[dat$ite == ite,][2, c(1,7)]
            # kappa <- dat[dat$ite == ite,][4, c(1,7)]
            # mcor_ord <- dat[dat$ite == ite,][5, c(1,7)]
            # mterms <- dat[dat$ite == ite,][6, c(1,7)]
            
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
        if(method == "prop"){
          term.can1 <- cbind(term.can1, as.matrix(term.can0))
        }  
      }  
       
      mPEHE.can1 <- data.frame(mPEHE.can1)
      mbias.can1 <- data.frame(mbias.can1)
      kappa.can1 <- data.frame(kappa.can1)
      mcor_ord.can1 <- data.frame(mcor_ord.can1)
      term.can1 <- data.frame(term.can1)
      
      mPEHE.can1$fun1 <- Fun1.can[model.can[s,2]]; mPEHE.can1$fun2 <- Fun2.can[model.can[s,3]]; mPEHE.can1$g <- g
      mbias.can1$fun1 <- Fun1.can[model.can[s,2]]; mbias.can1$fun2 <- Fun2.can[model.can[s,3]]; mbias.can1$g <- g
      kappa.can1$fun1 <- Fun1.can[model.can[s,2]]; kappa.can1$fun2 <- Fun2.can[model.can[s,3]]; kappa.can1$g <- g
      mcor_ord.can1$fun1 <- Fun1.can[model.can[s,2]]; mcor_ord.can1$fun2 <- Fun2.can[model.can[s,3]]; mcor_ord.can1$g <- g
      term.can1$fun1 <- Fun1.can[model.can[s,2]]; term.can1$fun2 <- Fun2.can[model.can[s,3]]; term.can1$g <- g
      
      mPEHE.can1$prob <- prob
      mbias.can1$prob <- prob
      kappa.can1$prob <- prob
      mcor_ord.can1$prob <- prob
      term.can1$prob <- prob
      
      mPEHE.can2 <- rbind(mPEHE.can2, mPEHE.can1)
      mbias.can2 <- rbind(mbias.can2, mbias.can1)
      kappa.can2 <- rbind(kappa.can2, kappa.can1)
      mcor_ord.can2 <- rbind(mcor_ord.can2, mcor_ord.can1)
      term.can2 <- rbind(term.can2, term.can1)
    }
  }
}  
 


### Create user-defined function, which extracts legends from ggplots ###

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

#########################################################################

df <- mbias.can2
#colplate <- c("#FDD0A2FF", "#C6DBEFFF", "#C7E9C0FF", "#DADAEBFF","#ffb6db", "#7e6148ff", "#5c88daff","#42b540ff","#e377c2ff")
colplate <- c("#FDD0A2FF", "#C7E9C0FF","#ffb6db", "#7e6148ff", "#5c88daff","#42b540ff","#e377c2ff") 
y_range <- c(0,0.2)
y_lab <- "mPEHE"


pic_fig <- function(df, colplate, y_range, type){
  
  pic.can <- list()
  pic.can0 <- list()
  count <- 1
  count0 <- 1
  
  for(prob in c("rct","obv")){
    
    df0 <- df[,!is.element(colnames(df),c("prob"))]
    df0 <- df0[df$prob == prob,]
    
    df0 <- reshape2::melt(df0, id = c("fun1","fun2","g"))
    df0$value <- as.numeric(df0$value)
    
    df_0 <- group_by(df0, fun1, fun2, g, variable) %>% summarise_all(list(mean = mean, sd = sd))
    
    for(g in unique(df$g)){ 
      
      df1 <- df0[df0$g == g,]

      sum <- df1 %>% group_by(g, fun1, fun2, variable) %>% summarise_all(list(mean = mean, sd = sd))

      if(type == "mPEHE"| type == "Mean absolute relative bias"|type == "Number of base functions"){
        facet_min <- sum %>% group_by(fun1, fun2) %>% summarise(min_mbias = min(mean), .groups = "drop")
        winners <- sum %>% group_by(g, fun1, fun2) %>% slice_min(mean, n = 1, with_ties = TRUE) %>%
          mutate(is_best = TRUE) %>% select(g, fun1, fun2, variable, is_best)
      }else if(type == "Mean of Cohen's kappa"| type == "Mean of Spearman's rank correlation"){
        facet_min <- sum %>% group_by(fun1, fun2) %>% summarise(min_mbias = max(mean), .groups = "drop")
        winners <- sum %>% group_by(g, fun1, fun2) %>% slice_max(mean, n = 1, with_ties = TRUE) %>%
          mutate(is_best = TRUE) %>% select(g, fun1, fun2, variable, is_best)
      }
      
      if(type == "Number of base functions"){
        #shape.val <- c(4,5,6,7)
        shape.val <- rep(4,4)
        col <- colplate[4:7]
      }else{
        #shape.val <- c(0,1,2,4,5,6,7)
        shape.val <- c(rep(1,3), rep(4,4))
        col <- colplate
      }
      
      df_1 <- sum %>% left_join(winners, by = c("g","fun1","fun2","variable")) %>% mutate(is_best = replace_na(is_best, FALSE))

      gb <- ggplot(sum, aes(x = variable, y = mean, color = variable, shape = variable)) + geom_point(size = 1.5)
      gb <- gb + scale_color_manual(values = col) + scale_shape_manual(values = shape.val)
      gb <- gb + geom_hline(data = facet_min, aes(yintercept = min_mbias), color = "red", linetype = "dashed", linewidth = 0.5)
      gb <- gb + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(.2), linewidth = 0.5)
      gb <- gb + guides(color = "none") + facet_grid(fun1 ~ fun2, labeller = label_parsed) + theme_bw() + coord_cartesian(ylim= y_range)
      gb <- gb + theme(legend.position = "none", legend.direction = "horizontal", text=element_text(size= 10), panel.grid.major = element_line(size = 0.3),  # Major grid line width
                       panel.grid.minor = element_line(size = 0.25),panel.border = element_rect(size = 0.5))
      gb <- gb + xlab("") + ylab("")

      pic.can[[count]] <- gb
      count <- count + 1
    }
  }
  
  library(ggpubr)
  gg <-  ggarrange(pic.can[[1]] + labs(title = paste0("RCT setting: 3-Groups \n (",type,")"), color = "Method", shape = "Method")+ guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=15), legend.title = element_text(size=15), axis.text.x=element_blank(), text = element_text(size = 15)),

                   pic.can[[2]] + labs(title = paste0("RCT setting: 4-Groups \n (",type,")"), color = "Method", shape = "Method") + guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=15), legend.title = element_text(size=15), axis.text.x=element_blank(), text = element_text(size = 15)),

                   pic.can[[3]] + labs(title = paste0("RCT setting: 5-Groups \n (",type,")"), color = "Method", shape = "Method")+ guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=15), legend.title = element_text(size=15), axis.text.x=element_blank(), text = element_text(size = 15)),

                   pic.can[[4]] + labs(title = paste0("Observational study setting: 3-Groups \n (",type,")"), color = "Method", shape = "Method") + guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=15), legend.title = element_text(size=15), axis.text.x=element_blank(), text = element_text(size = 15)),

                   pic.can[[5]] + labs(title = paste0("Observational study setting: 4-Groups \n (",type,")"), color = "Method", shape = "Method")+ guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=15), legend.title = element_text(size=15), axis.text.x=element_blank(), text = element_text(size = 15)),

                   pic.can[[6]] + labs(title = paste0("Observational study setting: 5-Groups \n (",type,")"), color = "Method", shape = "Method")+ guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=15), legend.title = element_text(size=15), axis.text.x=element_blank(), text = element_text(size = 15)),
                   ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom")
}

mPEHE.pic <- pic_fig(mPEHE.can2, colplate, y_range = c(0.15, 0.8), type = "mPEHE")
ggsave(paste0(dir,"/Table_Figure/Section4/Fig2.eps"),mPEHE.pic , width = 18, height = 10)
mBias.pic <- pic_fig(mbias.can2, colplate, y_range = c(0.0, 0.35), type = "Mean absolute relative bias")
ggsave(paste0(dir,"/Table_Figure/Section4/Fig3.eps"),mBias.pic , width = 18, height = 10)
mkappa.pic <- pic_fig(kappa.can2, colplate, y_range = c(0.6, 1), type = "Mean of Cohen's kappa")
ggsave(paste0(dir,"/Table_Figure/Section4/Fig4.eps"),mkappa.pic , width = 18, height = 10)
mcor_ord.pic <- pic_fig(mcor_ord.can2, colplate, y_range = c(0.65, 1), type = "Mean of Spearman's rank correlation")
ggsave(paste0(dir,"/Table_Figure/Section4/Fig5.eps"),mcor_ord.pic , width = 18, height = 10)
rule_terms.pic <- pic_fig(term.can2, colplate, y_range = c(0, 120), type = "Number of base functions")
ggsave(paste0(dir,"/Table_Figure/Section4/Fig6.eps"),rule_terms.pic , width = 18, height = 10)


