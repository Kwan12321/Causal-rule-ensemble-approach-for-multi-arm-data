library(gridExtra)
library(ggplot2)
library(greekLetters)
library(tidyr)
library(ggpubr)
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


fig.num <- 1

for(method in c("bart","forest","boost","prop")){

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
      
      #save.res <- paste0(dir, method,"/",method,"_evl_sim",s,"_n",n,"_p",p,"_",prob,"_",g,"g",".csv")
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
          # mPEHE <- dat[dat$ite == ite,][1, c(1:4)]
          # mbias <- dat[dat$ite == ite,][2, c(1:4)]
          # kappa <- dat[dat$ite == ite,][4, c(1:4)]
          # mcor_ord <- dat[dat$ite == ite,][5, c(1:4)]
          # mterms <- dat[dat$ite == ite,][6, c(1:4)]
          
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
      if(method == "prop"){
        term.can1 <- cbind(term.can1, as.matrix(term.can0))
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
 


### Create user-defined function, which extracts legends from ggplots ###

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

#########################################################################

if(method != "prop"){
  colplate <- c("#8DD3C7", "#BC80BD", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69")
  shapes <- rep(1,length(colplate))
}else if(method == "prop"){
  colplate <- c( "#7e6148ff", "#5c88daff","#42b540ff","#e377c2ff","#FDD0A2FF", "#C6DBEFFF", "#C7E9C0FF", "#DADAEBFF","#ffb6db", "#FFED6F","#920000","#004949")
  shapes <- c(rep(4,length(colplate)-4), rep(1,4)) 
}

#########################################################################

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
    
    for(g in unique(df$g)){
      
     df1 <- df0[df0$g == g,]

     sum <- df1 %>% group_by(g, fun1, fun2, variable) %>% summarise_all(list(mean = mean, sd = sd))

     if(type == "mPEHE"| type == "Mean absolute relative bias"){
       facet_min <- sum %>% group_by(fun1, fun2) %>% summarise(min_mbias = min(mean), .groups = "drop")
       winners <- sum %>% group_by(g, fun1, fun2) %>% slice_min(mean, n = 1, with_ties = TRUE) %>%
         mutate(is_best = TRUE) %>% select(g, fun1, fun2, variable, is_best)
     }else if(type == "Mean Cohen's kappa"| type == "Mean Spearman's rank correlation"){
       facet_min <- sum %>% group_by(fun1, fun2) %>% summarise(min_mbias = max(mean), .groups = "drop")
       winners <- sum %>% group_by(g, fun1, fun2) %>% slice_max(mean, n = 1, with_ties = TRUE) %>%
         mutate(is_best = TRUE) %>% select(g, fun1, fun2, variable, is_best)
     }

     df_1 <- sum %>% left_join(winners, by = c("g","fun1","fun2","variable")) %>% mutate(is_best = replace_na(is_best, FALSE))

     gb <- ggplot(sum, aes(x = variable, y = mean, color = variable, shape = variable)) + geom_point(size = 1.5)
     gb <- gb + scale_color_manual(values = colplate) + scale_shape_manual(values = shapes)
     gb <- gb + geom_hline(data = facet_min, aes(yintercept = min_mbias), color = "red", linetype = "dashed", linewidth = 0.5)
     gb <- gb + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(.2), linewidth = 0.5)
     gb <- gb + guides(color = "none") + facet_grid(fun1 ~ fun2) + theme_bw() + coord_cartesian(ylim= y_range)
     gb <- gb + theme(legend.position = "none", legend.direction = "horizontal", text=element_text(size= 30))
     gb <- gb + xlab("") + ylab("")
      
     pic.can[[count]] <- gb
     count <- count + 1
    }
  }
  
  library(ggpubr)
  gg <-  ggarrange(pic.can[[1]] + labs(title = paste0("RCT setting: 3-Groups \n (",type,")"), color = "Method", shape = "Method")+ guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),

                   pic.can[[2]] + labs(title = paste0("RCT setting: 4-Groups \n (",type,")"), color = "Method", shape = "Method") + guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),

                   pic.can[[3]] + labs(title = paste0("RCT setting: 5-Groups \n (",type,")"), color = "Method", shape = "Method")+ guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),

                   pic.can[[4]] + labs(title = paste0("Observational study setting: 3-Groups \n (",type,")"), color = "Method", shape = "Method") + guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),

                   pic.can[[5]] + labs(title = paste0("Observational study setting: 4-Groups \n (",type,")"), color = "Method", shape = "Method")+ guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),

                   pic.can[[6]] + labs(title = paste0("Observational study setting: 5-Groups \n (",type,")"), color = "Method", shape = "Method")+ guides(color = guide_legend(nrow = 1))
                   + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),
                   ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom")
  
}
  
mPEHE.pic <- pic_fig(mPEHE.can2, colplate, y_range = c(0, 2.5), type = "mPEHE")
ggsave(paste0(dir,"/Table_Figure_support/sup_section1/s-Fig",fig.num,".eps"),mPEHE.pic , width = 20, height = 18)
fig.num <- fig.num + 1
mBias.pic <- pic_fig(mbias.can2, colplate, y_range = c(0.0, 0.5), type = "Mean absolute relative bias")
ggsave(paste0(dir,"/Table_Figure_support/sup_section1/s-Fig",fig.num,".eps"),mBias.pic , width = 20, height = 18)
fig.num <- fig.num + 1
mkappa.pic <- pic_fig(kappa.can2, colplate, y_range = c(0, 1), type = "Mean Cohen's kappa")
ggsave(paste0(dir,"/Table_Figure_support/sup_section1/s-Fig",fig.num,".eps"),mkappa.pic , width = 20, height = 18)
fig.num <- fig.num + 1
mcor_ord.pic <- pic_fig(mcor_ord.can2, colplate, y_range = c(0, 1), type = "Mean Spearman's rank correlation")
ggsave(paste0(dir,"/Table_Figure_support/sup_section1/s-Fig",fig.num,".eps"),mcor_ord.pic , width = 20, height = 18)
fig.num <- fig.num + 1
}

