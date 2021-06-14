
rm(list=ls())
library(here)
library(SuperLearner)
library(tidyverse)

tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")
tableau11 <- c("Black","#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")


#simulate data

set.seed(12345)
U=rnorm(100)
X=1:100
X2=ifelse(X>50,1,0)
Y= X^2 - 18*X +1000*X2 + 100*U + 100*U^3 + 100000
summary(Y)
Y=(scale(Y, center=F)-.9)/10
plotdf=data.frame(Y,X,U)



# lib=c("SL.glm", "SL.gam", "SL.randomForest")
# res.SL= SuperLearner(Y=plotdf$Y,X=plotdf %>% select(X), SL.library = lib)
# res.SL

CV.glm= CV.SuperLearner(Y=plotdf$Y,X=plotdf %>% select(X), SL.library = "SL.glm",
                     cvControl = list(V = 5, shuffle = FALSE),
                     innerCvControl = list(list(V = 5)))
perfglm=summary(CV.glm)
perfglm$Table[1,2]

res.glm= SuperLearner(Y=plotdf$Y,X=plotdf %>% select(X), SL.library = "SL.glm")
glm.pred=data.frame(X=1:100, predY=predict(res.glm)$pred)


p <- ggplot(plotdf, aes(x=X, y=Y)) + geom_point() +
  xlab("Diabetes duration at baseline") + ylab("Dementia risk") +
  geom_path(aes(x=X, y=predY), data=glm.pred, color=tableau10[1], size=1.25, alpha=0) + 
  theme_minimal() + theme( panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank())
p
p1 <- ggplot(plotdf, aes(x=X, y=Y)) + geom_point() +
  geom_path(aes(x=X, y=predY), data=glm.pred, color=tableau10[1], size=1.25) + 
  xlab("Diabetes duration at baseline") + ylab("Dementia risk") +
  theme_minimal() + theme( panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank())



lib=c("SL.gam")
CV.gam= CV.SuperLearner(Y=plotdf$Y,X=plotdf %>% select(X), SL.library = lib,
                        cvControl = list(V = 5, shuffle = FALSE),
                        innerCvControl = list(list(V = 5)))
perfgam=summary(CV.gam)
perfgam$Table[1,2]

res.gam= SuperLearner(Y=plotdf$Y,X=plotdf %>% select(X), SL.library = lib)
gam.pred=data.frame(X=1:100, predY=predict(res.gam)$pred)
p2 <- ggplot(plotdf, aes(x=X, y=Y)) + geom_point() +
  geom_path(aes(x=X, y=predY), data=glm.pred, color=tableau10[1]) + 
  geom_path(aes(x=X, y=predY), data=gam.pred, color=tableau10[2], size=1.25) + 
  xlab("Diabetes duration at baseline") + ylab("Dementia risk") +
  theme_minimal() + theme( panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank())


lib=c("SL.randomForest")
CV.RF= CV.SuperLearner(Y=plotdf$Y,X=plotdf %>% select(X), SL.library = lib,
                        cvControl = list(V = 5, shuffle = FALSE),
                        innerCvControl = list(list(V = 5)))
perfRF=summary(CV.RF)
perfRF$Table[1,2]

res.RF= SuperLearner(Y=plotdf$Y,X=plotdf %>% select(X), SL.library = lib)
RF.pred=data.frame(X=1:100, predY=predict(res.RF)$pred)
p3 <- ggplot(plotdf, aes(x=X, y=Y)) + geom_point() +
  geom_path(aes(x=X, y=predY), data=glm.pred, color=tableau10[1]) + 
  geom_path(aes(x=X, y=predY), data=gam.pred, color=tableau10[2]) + 
  geom_path(aes(x=X, y=predY), data=RF.pred, color=tableau10[3], size=1.25) + 
  xlab("Diabetes duration at baseline") + ylab("Dementia risk") +
  theme_minimal() + theme( panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank())

lib=c("SL.glm","SL.gam", "SL.randomForest")
CV.SL= CV.SuperLearner(Y=plotdf$Y,X=plotdf %>% select(X), SL.library = lib,
                       cvControl = list(V = 5, shuffle = FALSE),
                       innerCvControl = list(list(V = 5)))
coef(CV.SL)
perfSL=summary(CV.SL)
perfSL$Table[1,2]

res.SL= SuperLearner(Y=plotdf$Y,X=plotdf %>% select(X), SL.library = lib)
SL.pred=data.frame(X=1:100, predY=predict(res.SL)$pred)
p4 <- ggplot(plotdf, aes(x=X, y=Y)) + geom_point() +
  geom_path(aes(x=X, y=predY), data=glm.pred, color=tableau10[1]) + 
  geom_path(aes(x=X, y=predY), data=gam.pred, color=tableau10[2]) + 
  geom_path(aes(x=X, y=predY), data=RF.pred, color=tableau10[3]) + 
  geom_path(aes(x=X, y=predY), data=SL.pred, color=tableau10[4], size=1.25) + 
  xlab("Diabetes duration at baseline") + ylab("Dementia risk") +
  theme_minimal() + theme( panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank())
p4

perfglm$Table[1,2]
perfgam$Table[1,2]
perfRF$Table[1,2]
perfSL$Table[1,2]

w=3
h=2
ggsave(p, file=here("figures/SLexample1.png"),  width = w, height = h)
ggsave(p1, file=here("figures/SLexample2.png"),  width = w, height = h)
ggsave(p2, file=here("figures/SLexample3.png"),  width = w, height = h)
ggsave(p3, file=here("figures/SLexample4.png"),  width = w, height = h)
ggsave(p4, file=here("figures/SLexample5.png"),  width = w, height = h)
