## ############################################
## Summarizing the Simulation Results
## ############################################
library("tidyverse")
library("gridExtra")


## location to save plots:
plot_path <- "/home/dom/ownCloud/PoIGLMShared/Manuscript/Figures"
##
B         <- 1000
DGP       <- 3
N.seq     <- c(100,200,500,1000,3000)
p.seq     <- c(100,500,1000)

## ##########################################
Sim_DGP3_df <- expand.grid(
  DGP             = DGP, 
  p               = p.seq,
  N               = N.seq,
  repet           = seq(B),
  ##
  "beta0.hat.PoI" = NA,  
  "beta1.hat.PoI" = NA,  
  "beta2.hat.PoI" = NA,  
  "beta3.hat.PoI" = NA,
  "beta4.hat.PoI" = NA,  
  "tau1.hat.PoI"  = NA, 
  "tau2.hat.PoI"  = NA, 
  "tau3.hat.PoI"  = NA, 
  "tau4.hat.PoI"  = NA, 
  "S.hat.PoI"     = NA, 
  "P.hat.PoI"	= NA, #  
  "beta0.hat.TRH" = NA, 
  "beta1.hat.TRH" = NA, 
  "beta2.hat.TRH" = NA, 
  "beta3.hat.TRH" = NA, 
  "beta4.hat.TRH" = NA, 
  "tau1.hat.TRH"  = NA,
  "tau2.hat.TRH"  = NA,
  "tau3.hat.TRH"  = NA,
  "tau4.hat.TRH"  = NA,
  "S.hat.TRH"     = NA,
  "P.hat.TRH"	= NA, #
  "beta0"         = NA,
  "beta1"         = NA,
  "beta2"         = NA,
  "beta3"         = NA,
  "beta4"         = NA,
  "tau1"          = NA,
  "tau2"          = NA,
  "tau3"          = NA,
  "tau4"          = NA,
  "S"             = NA
)


for(N in N.seq){
  for(p in p.seq){
    if(DGP==3){
      beta0         <- 1
      beta          <- c(-6, 6, -5, 5)
      S             <- length(beta)
      tau.true      <- c(1/6, 2/6, 4/6, 5/6) 
      t.grid        <- (1:p-1)/(p-1)
      tau.ind.true  <- rep(NA,S)
      for(s in 1:S){
        tau.ind.true[s]  <- which.min(abs(t.grid - tau.true[s]))
      }
      ##
      load(paste0("Simulation_Results/Sim_Results_DGP=", DGP, "_N=", N, "_p=", p, ".RData"))
      ##
      slct <- Sim_DGP3_df$DGP==DGP & Sim_DGP3_df$N==N & Sim_DGP3_df$p==p 
      ##
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta0.hat.PoI"] <- sim.results[ , colnames(sim.results)=="beta0.hat.PoI"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta1.hat.PoI"] <- sim.results[ , colnames(sim.results)=="beta1.hat.PoI"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta2.hat.PoI"] <- sim.results[ , colnames(sim.results)=="beta2.hat.PoI"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta3.hat.PoI"] <- sim.results[ , colnames(sim.results)=="beta3.hat.PoI"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta4.hat.PoI"] <- sim.results[ , colnames(sim.results)=="beta4.hat.PoI"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau1.hat.PoI"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau1.ind.hat.PoI"]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau2.hat.PoI"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau2.ind.hat.PoI"]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau3.hat.PoI"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau3.ind.hat.PoI"]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau4.hat.PoI"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau4.ind.hat.PoI"]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="S.hat.PoI"]     <- sim.results[ , colnames(sim.results)=="S.hat.PoI"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="P.hat.PoI"]     <- mean(sim.results[ , colnames(sim.results)=="S.hat.PoI"]==4,na.rm=T)      
      ##
      
      ##
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta0.hat.TRH"] <- sim.results[ , colnames(sim.results)=="beta0.hat.TRH"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta1.hat.TRH"] <- sim.results[ , colnames(sim.results)=="beta1.hat.TRH"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta2.hat.TRH"] <- sim.results[ , colnames(sim.results)=="beta2.hat.TRH"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta3.hat.TRH"] <- sim.results[ , colnames(sim.results)=="beta3.hat.TRH"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta4.hat.TRH"] <- sim.results[ , colnames(sim.results)=="beta4.hat.TRH"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau1.hat.TRH"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau1.ind.hat.TRH"]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau2.hat.TRH"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau2.ind.hat.TRH"]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau3.hat.TRH"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau3.ind.hat.TRH"]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau4.hat.TRH"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau4.ind.hat.TRH"]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="S.hat.TRH"]     <- sim.results[ , colnames(sim.results)=="S.hat.TRH"]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="P.hat.TRH"]     <- mean(sim.results[ , colnames(sim.results)=="S.hat.TRH"]==4,na.rm=T)
      ##
      
      ##
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta0"]          <-  beta0
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta1"]          <-  beta[1]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta2"]          <-  beta[2]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta3"]          <-  beta[3]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="beta4"]          <-  beta[4]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau1"]           <-  t.grid[tau.ind.true[1]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau2"]           <-  t.grid[tau.ind.true[2]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau3"]           <-  t.grid[tau.ind.true[3]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="tau4"]           <-  t.grid[tau.ind.true[4]]
      Sim_DGP3_df[slct,  colnames(Sim_DGP3_df)=="S"]              <-  S
      rm(sim.results)
    }
  }
}





#####################################################################
#####################################################################
xxx<-Sim_DGP3_df[Sim_DGP3_df[,1]==3,]
head(xxx)
xxx %>% group_by(N,p) %>% summarize(mean=mean(S.hat.TRH), na.rm=T)
xxx %>% group_by(N,p) %>% summarize(mean=mean(S.hat.PoI), na.rm=T)
xxx %>% group_by(N,p) %>% summarize(mean=mean(S.hat.PoI==S), na.rm=T)
xxx %>% group_by(N,p) %>% summarize(mean=mean(S.hat.TRH==S), na.rm=T)
#####################################################################
#####################################################################
## #################################################
SimRes_DGP3 <- tibble::as_tibble(Sim_DGP3_df) %>% 
  dplyr::mutate(p   = as.factor(p)) %>% 
  dplyr::mutate(N   = as.factor(N)) %>% 
  tidyr::gather(key = Estim_Types, value = Estim_Values, beta0.hat.PoI:P.hat.TRH, factor_key=TRUE) %>% 
  dplyr::mutate(Targets = fct_recode(Estim_Types, 
                                     "beta[0]==1" = "beta0.hat.PoI", "beta[0]==1" = "beta0.hat.TRH", 
                                     "beta[1]==-6"= "beta1.hat.PoI", "beta[1]==-6"= "beta1.hat.TRH", 
                                     "beta[2]==6" = "beta2.hat.PoI", "beta[2]==6" = "beta2.hat.TRH", 
                                     "beta[3]==-5"= "beta3.hat.PoI", "beta[3]==-5"= "beta3.hat.TRH", 
                                     "beta[4]==5" = "beta4.hat.PoI", "beta[4]==5" = "beta4.hat.TRH", 
                                     "tau[1]==1/6"= "tau1.hat.PoI",  "tau[1]==1/6"= "tau1.hat.TRH", 
                                     "tau[2]==2/6"= "tau2.hat.PoI",  "tau[2]==2/6"= "tau2.hat.TRH", 
                                     "tau[3]==4/6"= "tau3.hat.PoI",  "tau[3]==4/6"= "tau3.hat.TRH", 
                                     "tau[4]==5/6"= "tau4.hat.PoI",  "tau[4]==5/6"= "tau4.hat.TRH", 
                                     "S==4"       = "S.hat.PoI",     "S==4"       = "S.hat.TRH",
                                     "Phat==1"	  = "P.hat.PoI",	   "Phat==1"    = "P.hat.TRH")) %>% 
  dplyr::mutate(Estimators = fct_recode(Estim_Types, 
                                        "PoI" = "beta0.hat.PoI", "PoI" = "S.hat.PoI",
                                        "PoI" = "P.hat.PoI",
                                        "PoI" = "beta1.hat.PoI", "PoI" = "beta2.hat.PoI", 
                                        "PoI" = "beta3.hat.PoI", "PoI" = "beta4.hat.PoI", 
                                        "PoI" = "tau1.hat.PoI",  "PoI" = "tau2.hat.PoI", 
                                        "PoI" = "tau3.hat.PoI",  "PoI" = "tau4.hat.PoI", 
                                        "TRH" = "beta0.hat.TRH", "TRH" = "S.hat.TRH",
                                        "TRH" = "P.hat.TRH",
                                        "TRH" = "beta1.hat.TRH", "TRH" = "beta2.hat.TRH", 
                                        "TRH" = "beta3.hat.TRH", "TRH" = "beta4.hat.TRH", 
                                        "TRH" = "tau1.hat.TRH",  "TRH" = "tau2.hat.TRH", 
                                        "TRH" = "tau3.hat.TRH",  "TRH" = "tau4.hat.TRH")) %>% 
  group_by(p) %>% 
  dplyr::mutate(True_Values = case_when(
    Targets == "beta[0]==1"  ~ beta0,
    Targets == "beta[1]==-6" ~ beta1,
    Targets == "beta[2]==6"  ~ beta2,
    Targets == "beta[3]==-5" ~ beta3,
    Targets == "beta[4]==5"  ~ beta4,
    Targets == "tau[1]==1/6" ~ tau1,
    Targets == "tau[2]==2/6" ~ tau2,
    Targets == "tau[3]==4/6" ~ tau3,
    Targets == "tau[4]==5/6" ~ tau4,
    Targets == "S==4"        ~ as.double(S),
    Targets == "Phat==1"        ~ 1)) %>% 
  ungroup() %>% 
  dplyr::select(N, p, Targets, Estimators, Estim_Values, True_Values)



SimRes_DGP3%>%dplyr::select(Phat==1)



###############################################################

## #################
## DGP 3
## #################

#Boxplots: whiskers|-|==box==|----| whiser
fff <- function(x) {
  r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9),na.rm=T)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


#legend below graphs:
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}




########################
#DATA:
SimResSum_DGP3 <- SimRes_DGP3 %>% dplyr::filter(DGP==3)


########################
#beta0
########################
Plot_DGP3_1 <- ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "beta[0]==1"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5)+
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")
#get ylim for coord caretesion    
ylim.aux<-max(abs(c(ggplot_build(Plot_DGP3_1)$data[[2]]$ymin,ggplot_build(Plot_DGP3_1)$data[[2]]$ymax)))
ylim.val<-c(-ylim.aux,ylim.aux)
#
Plot_DGP3_1<-Plot_DGP3_1+
  coord_cartesian(ylim = ylim.val)  +  
  #    coord_cartesian(ylim = c(-2, 2))  +  
  labs(x = "", y="", title=expression(paste(paste(hat(alpha)-alpha*phantom(hat(beta)[1]))))) +
  scale_color_manual(name="",
                     labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                     values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  theme_bw() +
  theme(plot.margin     = margin(t =0, r = 0.5, b = 0, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))



Plot_DGP3_1
#########################################################

#########################################################
#beta1
#########################################################

Plot_DGP3_2<-ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "beta[1]==-6"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5)+
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")

#get ylim for coord caretesion    
ylim.aux<-max(abs(c(ggplot_build(Plot_DGP3_2)$data[[2]]$ymin,ggplot_build(Plot_DGP3_2)$data[[2]]$ymax)))
ylim.val<-c(-ylim.aux,ylim.aux)
# recontinue plotting
Plot_DGP3_2<-Plot_DGP3_2+
  coord_cartesian(ylim = ylim.val)  +  
  labs(x = "", y="", title=expression(paste(hat(beta)[1]-beta[1]))) +
  scale_color_manual(name="",
                     labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                     values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  theme_bw() +
  theme(plot.margin     = margin(t = 0, r = 0.5, b = 0, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))

Plot_DGP3_2

#########################################################
#beta2
#########################################################


Plot_DGP3_3<-ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "beta[2]==6"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5)+
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")

#get ylim for coord caretesion    
ylim.aux<-max(abs(c(ggplot_build(Plot_DGP3_3)$data[[2]]$ymin,ggplot_build(Plot_DGP3_3)$data[[2]]$ymax)))
ylim.val<-c(-ylim.aux,ylim.aux)
# recontinue plotting
Plot_DGP3_3<-Plot_DGP3_3+
  coord_cartesian(ylim = ylim.val)  +  
  labs(x = "", y="", title=expression(paste(hat(beta)[2]-beta[2]))) +
  scale_color_manual(name="",
                     labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                     values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  theme_bw() +
  theme(plot.margin     = margin(t = 0, r = 0.5, b = 0, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))

Plot_DGP3_2



#########################################################
#beta3
#########################################################


Plot_DGP3_4<-ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "beta[3]==-5"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5)+
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")

#get ylim for coord caretesion    
ylim.aux<-max(abs(c(ggplot_build(Plot_DGP3_4)$data[[2]]$ymin,ggplot_build(Plot_DGP3_4)$data[[2]]$ymax)))
ylim.val<-c(-ylim.aux,ylim.aux)
# recontinue plotting
Plot_DGP3_4<-Plot_DGP3_4+
  coord_cartesian(ylim = ylim.val)  +  
  labs(x = "", y="", title=expression(paste(hat(beta)[3]-beta[3]))) +
  scale_color_manual(name="",
                     labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                     values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  theme_bw() +
  theme(plot.margin     = margin(t = 0, r = 0.5, b = 0, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))

Plot_DGP3_4




#########################################################
#beta4
#########################################################


Plot_DGP3_5<-ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "beta[4]==5"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5)+
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")

#get ylim for coord caretesion    
ylim.aux<-max(abs(c(ggplot_build(Plot_DGP3_5)$data[[2]]$ymin,ggplot_build(Plot_DGP3_5)$data[[2]]$ymax)))
ylim.val<-c(-ylim.aux,ylim.aux)
# recontinue plotting
Plot_DGP3_5<-Plot_DGP3_5+
  coord_cartesian(ylim = ylim.val)  +  
  labs(x = "", y="", title=expression(paste(hat(beta)[4]-beta[4]))) +
  scale_color_manual(name="",
                     labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                     values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  theme_bw() +
  theme(plot.margin     = margin(t = 0, r = 0.5, b = 0, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))

Plot_DGP3_5


#########################################################
#tau1
#########################################################
Plot_DGP3_6 <- ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "tau[1]==1/6"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5) +
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")

#get ylim for coord caretesion    
ylim.aux<-max(abs(c(ggplot_build(Plot_DGP3_6)$data[[2]]$ymin,ggplot_build(Plot_DGP3_6)$data[[2]]$ymax)))
ylim.val<-c(-ylim.aux,ylim.aux)
# recontinue plotting
Plot_DGP3_6<-Plot_DGP3_6+
  coord_cartesian(ylim = ylim.val)  +  
  labs(x = "", y="", title=expression(paste(paste(hat(tau)[1]-tau[1]*phantom(hat(beta)[1]))))) +
  scale_color_manual(name="",
                     labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                     values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  theme_bw() +
  theme(plot.margin     = margin(t = 2, r = 0.5, b = -2, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))


Plot_DGP3_6 


#########################################################
#tau2
#########################################################
Plot_DGP3_7 <- ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "tau[2]==2/6"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5) +
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")

#get ylim for coord caretesion    
ylim.aux<-max(abs(c(ggplot_build(Plot_DGP3_7)$data[[2]]$ymin,ggplot_build(Plot_DGP3_7)$data[[2]]$ymax)))
ylim.val<-c(-ylim.aux,ylim.aux)
# recontinue plotting
Plot_DGP3_7<-Plot_DGP3_7+
  coord_cartesian(ylim = ylim.val)  +  
  labs(x = "", y="", title=expression(paste(paste(hat(tau)[2]-tau[2]*phantom(hat(beta)[1]))))) +
  scale_color_manual(name="",
                     labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                     values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  theme_bw() +
  theme(plot.margin     = margin(t = 2, r = 0.5, b = -2, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))



Plot_DGP3_7 




#########################################################
#tau3
#########################################################
Plot_DGP3_8 <- ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "tau[3]==4/6"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5) +
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")

#get ylim for coord caretesion    
ylim.aux<-max(abs(c(ggplot_build(Plot_DGP3_8)$data[[2]]$ymin,ggplot_build(Plot_DGP3_8)$data[[2]]$ymax)))
ylim.val<-c(-ylim.aux,ylim.aux)
# recontinue plotting
Plot_DGP3_8<-Plot_DGP3_8+
  coord_cartesian(ylim = ylim.val)  +  
  labs(x = "", y="", title=expression(paste(paste(hat(tau)[3]-tau[3]*phantom(hat(beta)[1]))))) +
  scale_color_manual(name="",
                     labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                     values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  theme_bw() +
  theme(plot.margin     = margin(t = 2, r = 0.5, b = -2, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))



Plot_DGP3_8



#########################################################
#tau4
#########################################################
Plot_DGP3_9 <- ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "tau[4]==5/6"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5) +
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")

#get ylim for coord caretesion    
ylim.aux<-max(abs(c(ggplot_build(Plot_DGP3_9)$data[[2]]$ymin,ggplot_build(Plot_DGP3_9)$data[[2]]$ymax)))
ylim.val<-c(-ylim.aux,ylim.aux)
# recontinue plotting
Plot_DGP3_9<-Plot_DGP3_9+
  coord_cartesian(ylim = ylim.val)  +  
  labs(x = "", y="", title=expression(paste(paste(hat(tau)[4]-tau[4]*phantom(hat(beta)[1]))))) +
  scale_color_manual(name="",
                     labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                     values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  theme_bw() +
  theme(plot.margin     = margin(t = 2, r = 0.5, b = -2, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))



Plot_DGP3_9






#########################################################
#Shat
#########################################################
Plot_DGP3_10 <- ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "S==4"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5) +
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")

#get ylim for coord caretesion    
ylim.aux<-max(abs(c(ggplot_build(Plot_DGP3_10)$data[[2]]$ymin,ggplot_build(Plot_DGP3_10)$data[[2]]$ymax)))
ylim.val<-c(-ylim.aux,ylim.aux)
# recontinue plotting
Plot_DGP3_10<-Plot_DGP3_10+
  coord_cartesian(ylim = ylim.val)  + 
  labs(x = "", y="", title=expression(paste(paste(hat(S)-S*phantom(hat(beta)[1]))))) +
  scale_color_manual(name="",
                     labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                     values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  theme_bw() +
  theme(plot.margin     = margin(t = 2, r = 0.5, b = -2, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))


Plot_DGP3_10



#########################################################
#Phat(Shat=S)
#########################################################

Plot_DGP3_11 <- ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "Phat==1"),
  mapping = aes(x = N, y = Estim_Values, group=Estimators:p, color = p)) +
  geom_point()+
  geom_line(aes(lty=Estimators))+scale_color_manual(name="",
                                                    labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                                                    values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  geom_hline(yintercept = 1,linetype=5)+
  coord_cartesian(ylim = c(0, 1))+
  labs(x = "", y="", title=expression(paste(paste(hat(P)(hat(S)==S)*phantom(hat(beta)[1]))))) +
  theme_bw() +
  theme(plot.margin     = margin(t = 2, r = 0.5, b = -2, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +guides(fill=guide_legend(ncol=2))+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))

#    theme(legend.justification = c(1, 0), legend.position = c(1, 1), )+theme(legend.direction="horizontal",legend.box = "horizontal")
Plot_DGP3_11

#####################################################################
##################get legend:

Plot_DGP3_2.aux<-ggplot(
  data = SimResSum_DGP3 %>% dplyr::filter(Targets == "beta[1]==-6"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5)+
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")

#get ylim for coord caretesion    
ylim.aux<-max(abs(c(ggplot_build(Plot_DGP3_2)$data[[2]]$ymin,ggplot_build(Plot_DGP3_2)$data[[2]]$ymax)))
ylim.val<-c(-ylim.aux,ylim.aux)
# recontinue plotting
Plot_DGP3_2.aux<-Plot_DGP3_2.aux+
  coord_cartesian(ylim = ylim.val)  +  
  labs(x = "", y="", title=expression(paste(hat(beta)[1]-beta[1]))) +
  scale_color_manual(name="",
                     labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                     values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
  theme_bw() +
  theme(plot.margin     = margin(t = 2, r = 0.5, b = 0, l = -3, unit = "mm")) +
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "bottom")+
  theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))+
  theme(legend.text=element_text(size=13), legend.title=element_text(size=13))+
  theme(legend.key.size = unit(1.5,"line"))



Plot_DGP3_2.aux
mylegend <- g_legend(Plot_DGP3_2.aux)





Plot_DGP3<-grid.arrange(arrangeGrob(Plot_DGP3_6, Plot_DGP3_7, Plot_DGP3_8, Plot_DGP3_9,Plot_DGP3_10,
                                    Plot_DGP3_1, Plot_DGP3_2, Plot_DGP3_3,Plot_DGP3_4, Plot_DGP3_5,
                                    nrow=2#,
                                    # top = ""
),mylegend, nrow=2,heights=c(6,0.2))

ggsave(filename = "Plot_DGP3BOX.pdf", plot = Plot_DGP3, device = "pdf", path = plot_path, width = 11.5, height = 8)





Plot_DGP3a<-grid.arrange(arrangeGrob(Plot_DGP3_6, Plot_DGP3_7, Plot_DGP3_8, Plot_DGP3_9,Plot_DGP3_11,
                                     Plot_DGP3_1, Plot_DGP3_2, Plot_DGP3_3,Plot_DGP3_4, Plot_DGP3_5,
                                     nrow=2),#,
                         # top = "Estimation errors for different sample sizes n")
                         mylegend, nrow=2,heights=c(6,0.2))


ggsave(filename = "Plot_DGP_3.pdf", plot = Plot_DGP3a, device = "pdf", path = plot_path, width=11.5, height=8)   #11.5,8


