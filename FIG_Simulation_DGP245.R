## ############################################
## Summarizing the Simulation Results
## ############################################
library("tidyverse")
library("gridExtra")

# location to save plots:
plot_path <- "/home/dom/ownCloud/PoIGLMShared/Manuscript/Figures"


B       <- 1000
DGP.aux <- c(2,4,5)
N.seq   <- c(100,200,500,1000,3000)
p.seq   <- c(100,500,1000)

## ##########################################
Sim_DGP245_df <- expand.grid(
  DGP             = DGP.aux,
  p               = p.seq,
  N               = N.seq,
  repet           = seq(B),
  ##
  "beta0.hat.PoI" = NA,  
  "beta1.hat.PoI" = NA,  
  "beta2.hat.PoI" = NA,  
  "tau1.hat.PoI"  = NA, 
  "tau2.hat.PoI"  = NA, 
  "S.hat.PoI"     = NA, 
  "P.hat.PoI"	= NA, #  
  "beta0.hat.TRH" = NA, 
  "beta1.hat.TRH" = NA, 
  "beta2.hat.TRH" = NA, 
  "tau1.hat.TRH"  = NA,
  "tau2.hat.TRH"  = NA,
  "S.hat.TRH"     = NA,
  "P.hat.TRH"	= NA, #
  "beta0"         = NA,
  "beta1"         = NA,
  "beta2"         = NA,
  "tau1"          = NA,
  "tau2"          = NA,
  "S"             = NA
)


for(DGP in DGP.aux ){
  #for(DGP in c(2)){
  for(N in N.seq){ #N<-100;p<-100
    for(p in p.seq){
      if(any( DGP==c(2,4,5) )){
        beta0         <- 1
        beta          <- c(-6, 5)
        #if(AUX==FALSE){if(DGP==4){beta <- c(-5,4)} }  #
        S             <- length(beta)
        tau.true      <- c(1/3,2/3) 
        t.grid        <- (1:p-1)/(p-1)
        tau.ind.true  <- rep(NA,S)
        for(s in 1:S){
          tau.ind.true[s]  <- which.min(abs(t.grid - tau.true[s]))
        }
      }
      ##
      load(paste0("Simulation_Results/Sim_Results_DGP=", DGP, "_N=", N, "_p=", p, ".RData"))
      ##
      slct <- Sim_DGP245_df$DGP==DGP & Sim_DGP245_df$N==N & Sim_DGP245_df$p==p 
      ##
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="beta0.hat.PoI"] <- sim.results[ , colnames(sim.results)=="beta0.hat.PoI"]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="beta1.hat.PoI"] <- sim.results[ , colnames(sim.results)=="beta1.hat.PoI"]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="beta2.hat.PoI"] <- sim.results[ , colnames(sim.results)=="beta2.hat.PoI"]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="tau1.hat.PoI"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau1.ind.hat.PoI"]]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="tau2.hat.PoI"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau2.ind.hat.PoI"]]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="S.hat.PoI"]     <- sim.results[ , colnames(sim.results)=="S.hat.PoI"]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="P.hat.PoI"]     <- mean(sim.results[ , colnames(sim.results)=="S.hat.PoI"]==2,na.rm=T)      
      ##
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="beta0.hat.TRH"] <- sim.results[ , colnames(sim.results)=="beta0.hat.TRH"]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="beta1.hat.TRH"] <- sim.results[ , colnames(sim.results)=="beta1.hat.TRH"]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="beta2.hat.TRH"] <- sim.results[ , colnames(sim.results)=="beta2.hat.TRH"]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="tau1.hat.TRH"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau1.ind.hat.TRH"]]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="tau2.hat.TRH"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau2.ind.hat.TRH"]]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="S.hat.TRH"]     <- sim.results[ , colnames(sim.results)=="S.hat.TRH"]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="P.hat.TRH"]     <- mean(sim.results[ , colnames(sim.results)=="S.hat.TRH"]==2,na.rm=T)
      ##
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="beta0"]          <-  beta0
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="beta1"]          <-  beta[1]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="beta2"]          <-  beta[2]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="tau1"]           <-  t.grid[tau.ind.true[1]]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="tau2"]           <-  t.grid[tau.ind.true[2]]
      Sim_DGP245_df[slct,  colnames(Sim_DGP245_df)=="S"]              <-  S
      rm(sim.results)
    }
  }
}


#####################################################################
#####################################################################
xxx <- Sim_DGP245_df[Sim_DGP245_df[,1]==4,]
head(xxx)
xxx %>% group_by(N,p) %>% summarize(mean=mean(S.hat.TRH), na.rm=T)
xxx %>% group_by(N,p) %>% summarize(mean=mean(S.hat.PoI), na.rm=T)
xxx %>% group_by(N,p) %>% summarize(mean=mean(S.hat.PoI==S), na.rm=T)
xxx %>% group_by(N,p) %>% summarize(mean=mean(S.hat.TRH==S), na.rm=T)
#####################################################################
#####################################################################


SimRes_DGP245 <- tibble::as_tibble(Sim_DGP245_df) %>% 
  dplyr::mutate(p   = as.factor(p)) %>% 
  dplyr::mutate(N   = as.factor(N)) %>% 
  dplyr::mutate(DGP = as.factor(DGP)) %>% 
  tidyr::gather(key = Estim_Types, value = Estim_Values, beta0.hat.PoI:P.hat.TRH, factor_key=TRUE) %>% 
  dplyr::mutate(Targets = fct_recode(Estim_Types, 
                                     "beta[0]==1" = "beta0.hat.PoI", "beta[0]==1" = "beta0.hat.TRH", 
                                     "beta[1]==-6"= "beta1.hat.PoI", "beta[1]==-6"= "beta1.hat.TRH", 
                                     "beta[2]==5" = "beta2.hat.PoI", "beta[2]==5" = "beta2.hat.TRH", 
                                     "tau[1]==1/3"= "tau1.hat.PoI",  "tau[1]==1/3"= "tau1.hat.TRH", 
                                     "tau[2]==2/3"= "tau2.hat.PoI",  "tau[2]==2/3"= "tau2.hat.TRH", 
                                     "S==2"       = "S.hat.PoI",     "S==2"       = "S.hat.TRH",
                                     "Phat==1"	  = "P.hat.PoI",	   "Phat==1"    = "P.hat.TRH")) %>% 
  dplyr::mutate(Estimators = fct_recode(Estim_Types, 
                                        "PoI" = "beta0.hat.PoI", "PoI" = "S.hat.PoI",
                                        "PoI" = "P.hat.PoI",
                                        "PoI" = "beta1.hat.PoI", "PoI" = "beta2.hat.PoI", 
                                        "PoI" = "tau1.hat.PoI",  "PoI" = "tau2.hat.PoI", 
                                        "TRH" = "beta0.hat.TRH", "TRH" = "S.hat.TRH",
                                        "TRH" = "P.hat.TRH",
                                        "TRH" = "beta1.hat.TRH", "TRH" = "beta2.hat.TRH", 
                                        "TRH" = "tau1.hat.TRH",  "TRH" = "tau2.hat.TRH")) %>% 
  group_by(p) %>% 
  dplyr::mutate(True_Values = case_when(
    Targets == "beta[0]==1"  ~ beta0,
    Targets == "beta[1]==-6" ~ beta1,
    Targets == "beta[2]==5"  ~ beta2,
    Targets == "tau[1]==1/3" ~ tau1,
    Targets == "tau[2]==2/3" ~ tau2,
    Targets == "S==2"        ~ as.double(S),
    Targets == "Phat==1"        ~ 1 )) %>% 
  ungroup() %>% 
  dplyr::select(N, p, DGP, Targets, Estimators, Estim_Values, True_Values)


# Boxplots: whiskers|-|==box==|----| whiser
fff <- function(x) {
  r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9),na.rm=T)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


# legend below graphs:
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}




for(DGP.sel in DGP.aux){
  #select data to plot, all DGPS used here have 2 poi - use same generator!
  #SimResSum_DGP2 <- SimRes_DGP245 %>% dplyr::filter(DGP==2)
  #SimResSum_DGP2 <- SimRes_DGP245 %>% dplyr::filter(DGP==5)
  SimResSum_DGP2 <- SimRes_DGP245 %>% dplyr::filter(DGP==DGP.sel)
  
  #beta0
  Plot_DGP2_1 <- ggplot(
    data = SimResSum_DGP2 %>% dplyr::filter(Targets == "beta[0]==1"),
    mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
    geom_hline(yintercept = 0,linetype=5)+
    stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
    stat_summary(fun.data = fff, geom="boxplot",position="dodge")
  #get ylim for coord caretesion    
  ylim.aux<-max(abs(c(ggplot_build(Plot_DGP2_1)$data[[2]]$ymin,ggplot_build(Plot_DGP2_1)$data[[2]]$ymax)))
  ylim.val<-c(-ylim.aux,ylim.aux)
  #
  Plot_DGP2_1<-Plot_DGP2_1+
    coord_cartesian(ylim = ylim.val)  +  
    #    coord_cartesian(ylim = c(-2, 2))  +  
    labs(x = "", y="", title=expression(paste(paste(hat(alpha)-alpha*phantom(hat(beta)[1]))))) +
    scale_color_manual(name="",
                       labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                       values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
    theme_bw() +
    theme(plot.margin     = margin(t = -2, r = 0.5, b = 0, l = -3, unit = "mm")) +
    theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
    theme(legend.position = "none")+
    theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))
  
  
  
  Plot_DGP2_1
  
  
  #beta1
  Plot_DGP2_2<-ggplot(
    data = SimResSum_DGP2 %>% dplyr::filter(Targets == "beta[1]==-6"),
    mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
    geom_hline(yintercept = 0,linetype=5)+
    stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
    stat_summary(fun.data = fff, geom="boxplot",position="dodge")
  
  #get ylim for coord caretesion    
  ylim.aux<-max(abs(c(ggplot_build(Plot_DGP2_2)$data[[2]]$ymin,ggplot_build(Plot_DGP2_2)$data[[2]]$ymax)))
  ylim.val<-c(-ylim.aux,ylim.aux)
  # recontinue plotting
  Plot_DGP2_2<-Plot_DGP2_2+
    coord_cartesian(ylim = ylim.val)  +  
    labs(x = "", y="", title=expression(paste(hat(beta)[1]-beta[1]))) +
    scale_color_manual(name="",
                       labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                       values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
    theme_bw() +
    theme(plot.margin     = margin(t = -2, r = 0.5, b = 0, l = -3, unit = "mm")) +
    theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
    theme(legend.position = "none")+
    theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))
  
  
  
  Plot_DGP2_2
  
  
  #beta2
  Plot_DGP2_3<-ggplot(
    data = SimResSum_DGP2 %>% dplyr::filter(Targets == "beta[2]==5"),
    mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
    geom_hline(yintercept = 0,linetype=5)+
    stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
    stat_summary(fun.data = fff, geom="boxplot",position="dodge")
  
  #get ylim for coord caretesion    
  ylim.aux<-max(abs(c(ggplot_build(Plot_DGP2_3)$data[[2]]$ymin,ggplot_build(Plot_DGP2_3)$data[[2]]$ymax)))
  ylim.val<-c(-ylim.aux,ylim.aux)
  # recontinue plotting
  Plot_DGP2_3<-Plot_DGP2_3+
    coord_cartesian(ylim = ylim.val)  +  
    labs(x = "", y="", title=expression(paste(hat(beta)[2]-beta[2]))) +
    scale_color_manual(name="",
                       labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                       values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
    theme_bw() +
    theme(plot.margin     = margin(t = -2, r = 0.5, b = 0, l = -3, unit = "mm")) +
    theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
    theme(legend.position = "none")+
    theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))
  
  
  
  Plot_DGP2_3

  
  #tau1 
  Plot_DGP2_4 <- ggplot(
    data = SimResSum_DGP2 %>% dplyr::filter(Targets == "tau[1]==1/3"),
    mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
    geom_hline(yintercept = 0,linetype=5) +
    stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
    stat_summary(fun.data = fff, geom="boxplot",position="dodge")
  
  #get ylim for coord caretesion    
  ylim.aux<-max(abs(c(ggplot_build(Plot_DGP2_4)$data[[2]]$ymin,ggplot_build(Plot_DGP2_4)$data[[2]]$ymax)))
  ylim.val<-c(-ylim.aux,ylim.aux)
  # recontinue plotting
  Plot_DGP2_4<-Plot_DGP2_4+
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
  
  
  
  Plot_DGP2_4
  
  #tau2 
  Plot_DGP2_5 <- ggplot(
    data = SimResSum_DGP2 %>% dplyr::filter(Targets == "tau[2]==2/3"),
    mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
    geom_hline(yintercept = 0,linetype=5) +
    stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
    stat_summary(fun.data = fff, geom="boxplot",position="dodge")
  
  #get ylim for coord caretesion    
  ylim.aux<-max(abs(c(ggplot_build(Plot_DGP2_5)$data[[2]]$ymin,ggplot_build(Plot_DGP2_5)$data[[2]]$ymax)))
  ylim.val<-c(-ylim.aux,ylim.aux)
  # recontinue plotting
  Plot_DGP2_5<-Plot_DGP2_5+
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
  
  Plot_DGP2_5
  
  #Shat
  Plot_DGP2_6 <- ggplot(
    data = SimResSum_DGP2 %>% dplyr::filter(Targets == "S==2"),
    mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
    geom_hline(yintercept = 0,linetype=5) +
    stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
    stat_summary(fun.data = fff, geom="boxplot",position="dodge")
  
  #get ylim for coord caretesion    
  ylim.aux<-max(abs(c(ggplot_build(Plot_DGP2_6)$data[[2]]$ymin,ggplot_build(Plot_DGP2_6)$data[[2]]$ymax)))
  ylim.val<-c(-ylim.aux,ylim.aux)
  # recontinue plotting
  Plot_DGP2_6<-Plot_DGP2_6+
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
  
  
  Plot_DGP2_6
  
  
  Plot_DGP2_7 <- 
    ggplot(
      data = SimResSum_DGP2 %>% dplyr::filter(Targets == "Phat==1"),
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
  Plot_DGP2_7
  
  ##################get legend:
  Plot_DGP1_2.aux<-ggplot(
    data = SimResSum_DGP2 %>% dplyr::filter(Targets == "beta[2]==5"),
    mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
    geom_hline(yintercept = 0,linetype=5)+
    stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
    stat_summary(fun.data = fff, geom="boxplot",position="dodge")+
    #   geom_errorbar(aes(ymin=min, ymax= max), position ="dodge", width=0.5)
    # geom_boxplot(outlier.alpha = 0)
    coord_cartesian(ylim = c(-3, 3))  +  
    #    scale_y_continuous(limits = c(-2, 2)) +
    labs(x = "", y="", title=expression(paste(hat(beta)[1]-beta[1]))) +
    scale_color_manual(name="",
                       labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                       values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
    #    scale_y_continuous(limits = c(-0.01, 0.01)) +
    theme_bw() +
    theme(plot.margin     = margin(t = 2, r = 0.5, b = 0, l = 0, unit = "mm")) +
    theme(plot.title = element_text(size=9,)) +
    theme(legend.position = "bottom")+
    theme(legend.text=element_text(size=13), legend.title=element_text(size=13))+
    theme(legend.key.size = unit(1.5,"line"))
  
  
  Plot_DGP1_2.aux
  #Plot_DGP1_2.aux + theme(legend.key.size = unit(1.5,"line"))
  mylegend<-g_legend(Plot_DGP1_2.aux)
  #################
  #Plot_DGP2<-grid.arrange(arrangeGrob(Plot_DGP2_4, Plot_DGP2_5, Plot_DGP2_6,
  #						Plot_DGP2_1, Plot_DGP2_2, Plot_DGP2_3,#
  #											 nrow=2,
  # top = "Estimation errors for different sample sizes n"),mylegend, nrow=2,heights=c(6,0.2))
  
  Plot_DGP2<-grid.arrange(arrangeGrob(Plot_DGP2_4, Plot_DGP2_5, Plot_DGP2_6,
                                      Plot_DGP2_1, Plot_DGP2_2, Plot_DGP2_3,
                                      nrow=2),mylegend, nrow=2,heights=c(6,0.2))
  
  
  
  ggsave(filename = paste0("Plot_DGP",DGP.sel,"BOX.pdf"), plot = Plot_DGP2, device = "pdf", path = plot_path, width = 11.5, height = 8)
  ###########
  
  #Plot_DGP2A<-grid.arrange(arrangeGrob(Plot_DGP2_4, Plot_DGP2_5, Plot_DGP2_7,
  #						Plot_DGP2_1, Plot_DGP2_2, Plot_DGP2_3,
  #nrow=2,
  #top = "Estimation errors for different sample sizes n"),mylegend, nrow=2,heights=c(6,0.2))
  
  Plot_DGP2A <- grid.arrange(arrangeGrob(Plot_DGP2_4, Plot_DGP2_5, Plot_DGP2_7,
                                         Plot_DGP2_1, Plot_DGP2_2, Plot_DGP2_3,
                                         nrow=2),mylegend, nrow=2,heights=c(6,0.2))
  
  
  
  
  #ggsave(filename = "Plot_DGP2NEULINES.pdf", plot = Plot_DGP2A, device = "pdf", path = plot_path, width = 8, height = 6)
  ggsave(filename = paste0("Plot_DGP_",DGP.sel,".pdf"), plot = Plot_DGP2A, device = "pdf", path = plot_path, width = 11.5, height = 8)  #8 und 6
}


# 
