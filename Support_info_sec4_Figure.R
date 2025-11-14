library(gridExtra)
library(ggplot2)
library(greekLetters)
library(tidyr)
library(ggpubr)
#########################

dir <- getwd()

######################################################
### Supporting information : Section 4 (Figure 20) ###
######################################################

### DGP parameters ###
n <- 2000
p <- 10
model.can <- expand.grid(1:3, 2:4)
model.can <- cbind(sim = 1:9, model.can)
Beta.can <- 2*rbind(c(1, 1, 1), c(-0.5, 1.0, 2.0), c(1.5, 1.5, -0.5), c(-1.5, 1.5, 0.5), c(-0.5, 2.0, 0.5))
group.can <- c(3,4,5)
ite.all <- 100

############################################################################################################################
model.all <- NULL
for(g in group.can){
  model.set <- cbind(model.can, g) 
  model.all <- rbind(model.all, model.set)
}

cover.all <- NULL
width.all <- NULL

for(m in 1:nrow(model.all)){
  
  mod <- model.all[m,]
  save.res <- paste0(dir, "/add_simulation/Section4_conformal_prediction/cf_sim",mod$sim,"_n",n,"_p",p,"_",mod$prob,"_",mod$g,"g",".csv")
  dat <- read.csv(save.res)
  
  cover <- dat[dat$evl == "cover_r",c(1:5,8:10)]
  width <- dat[dat$evl == "width",c(1:5,8:10)]
  
  cover.all <- rbind(cover.all, cover)
  width.all <- rbind(width.all, width)
  
}                     

cover.all <- reshape2::melt(cover.all, id = c("fun1","fun2","g"))
width.all <- reshape2::melt(width.all, id = c("fun1","fun2","g"))

Fun1.can <- c("Linear", "Stepwise", "Non-linear")
Fun2.can <- c("", "Linear", "Stepwise", "Non-linear")
Fun1.can <- paste(Fun1.can, greeks("mu"))
Fun2.can <- paste(Fun2.can, greeks("delta"))
Fun1.can <- factor(Fun1.can, levels = Fun1.can)
Fun2.can <- factor(Fun2.can, levels = Fun2.can)

cover.all[,"fun1"] <- Fun1.can[cover.all[,"fun1"]]; width.all[,"fun1"] <- Fun1.can[width.all[,"fun1"]]
cover.all[,"fun2"] <- Fun2.can[cover.all[,"fun2"]]; width.all[,"fun2"] <- Fun2.can[width.all[,"fun2"]]

colplate <- c("skyblue","#7e6148ff", "#5c88daff","#42b540ff","#e377c2ff")

gc1 <- ggplot(cover.all[cover.all$g == 3,], aes(x = variable, y = value, fill = variable)) + geom_boxplot()
gc1 <- gc1 + scale_color_manual(values = colplate) 
gc1 <- gc1 + guides(color = "none") + facet_grid(fun1 ~ fun2) + theme_bw() + coord_cartesian(ylim= c(0.5,1))
gc1 <- gc1 + theme(legend.position = "none", legend.direction = "horizontal", text=element_text(size= 30))
gc1 <- gc1 + xlab("") + ylab("")

gc2 <- ggplot(cover.all[cover.all$g == 4,], aes(x = variable, y = value, fill = variable)) + geom_boxplot()
gc2 <- gc2 + scale_color_manual(values = colplate) 
gc2 <- gc2 + guides(color = "none") + facet_grid(fun1 ~ fun2) + theme_bw() + coord_cartesian(ylim= c(0.5,1))
gc2 <- gc2 + theme(legend.position = "none", legend.direction = "horizontal", text=element_text(size= 30))
gc2 <- gc2 + xlab("") + ylab("")

gc3 <- ggplot(cover.all[cover.all$g == 5,], aes(x = variable, y = value, fill = variable)) + geom_boxplot()
gc3 <- gc3 + scale_color_manual(values = colplate) 
gc3 <- gc3 + guides(color = "none") + facet_grid(fun1 ~ fun2) + theme_bw() + coord_cartesian(ylim= c(0.5,1))
gc3 <- gc3 + theme(legend.position = "none", legend.direction = "horizontal", text=element_text(size= 30))
gc3 <- gc3 + xlab("") + ylab("")

gw1 <- ggplot(width.all[width.all$g == 3,], aes(x = variable, y = value, fill = variable)) + geom_boxplot()
gw1 <- gw1 + scale_color_manual(values = colplate) 
gw1 <- gw1 + guides(color = "none") + facet_grid(fun1 ~ fun2) + theme_bw() #+ coord_cartesian(ylim= c(0,6))
gw1 <- gw1 + theme(legend.position = "none", legend.direction = "horizontal", text=element_text(size= 30))
gw1 <- gw1 + xlab("") + ylab("")

gw2 <- ggplot(width.all[width.all$g == 4,], aes(x = variable, y = value, fill = variable)) + geom_boxplot()
gw2 <- gw2 + scale_color_manual(values = colplate) 
gw2 <- gw2 + guides(color = "none") + facet_grid(fun1 ~ fun2) + theme_bw() + coord_cartesian(ylim= c(0,6))
gw2 <- gw2 + theme(legend.position = "none", legend.direction = "horizontal", text=element_text(size= 30))
gw2 <- gw2 + xlab("") + ylab("")

gw3 <- ggplot(width.all[width.all$g == 5,], aes(x = variable, y = value, fill = variable)) + geom_boxplot()
gw3 <- gw3 + scale_color_manual(values = colplate) 
gw3 <- gw3 + guides(color = "none") + facet_grid(fun1 ~ fun2) + theme_bw() + coord_cartesian(ylim= c(0,6))
gw3 <- gw3 + theme(legend.position = "none", legend.direction = "horizontal", text=element_text(size= 30))
gw3 <- gw3 + xlab("") + ylab("")


gg <-  ggarrange(gc1 + labs(title = paste0("Coverage rate: 3-Groups"), fill = "Method")+ guides(color = guide_legend(nrow = 1))
                 + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),
                 
                 gc2 + labs(title = paste0("Coverage rate: 4-Groups"), fill = "Method") + guides(color = guide_legend(nrow = 1))
                 + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),
                 
                 gc3 + labs(title = paste0("Coverage rate: 5-Groups"), fill = "Method")+ guides(color = guide_legend(nrow = 1))
                 + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),
                 
                 gw1 + labs(title = paste0("Conformal interval width: 3-Groups"), fill = "Method") + guides(color = guide_legend(nrow = 1))
                 + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),
                 
                 gw2 + labs(title = paste0("Conformal interval width: 4-Groups"), fill = "Method")+ guides(color = guide_legend(nrow = 1))
                 + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),
                 
                 gw3 + labs(title = paste0("Conformal interval width: 5-Groups"), fill = "Method")+ guides(color = guide_legend(nrow = 1))
                 + theme(legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x=element_blank(), text = element_text(size = 15)),
                 ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom")

ggsave(paste0(dir,"/Table_Figure_support/sup_section4/s-Fig20.eps"),gg , width = 20, height = 18)



  