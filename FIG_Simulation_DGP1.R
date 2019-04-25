## ############################################
## Summarizing the Simulation Results
## ############################################
library("tidyverse")
library("gridExtra")

# location to save plots:
plot_path <- "/home/dom/ownCloud/PoIGLMShared/Manuscript/Figures"

B     <- 1000
DGP   <- 1
N.seq <- c(100,200,500,1000,3000)
p.seq <- c(100,500,1000)

# N<-100
# p<-100


## ##########################################
Sim_DGP1_df <- expand.grid(
  DGP             = 1, 
  p               = p.seq,
  N               = N.seq,
  repet           = seq(B),
  ##
  "beta0.hat.PoI" = NA,  
  "beta1.hat.PoI" = NA,  
  "tau1.hat.PoI"  = NA, 
  "beta0.hat.LMcK"= NA, 
  "beta1.hat.LMcK"= NA, 
  "tau1.hat.LMcK" = NA,
  "beta0"         = NA,
  "beta1"         = NA,
  "tau1"          = NA
)

for(N in N.seq){
  for(p in p.seq){
    if(DGP==1){
      beta0         <- 1
      beta          <- c(4)
      S             <- length(beta)
      tau.true      <- c(1/2) 
      t.grid        <- (1:p-1)/(p-1)
      tau.ind.true  <- rep(NA,S)
      for(s in 1:S){tau.ind.true[s] <- which.min(abs(t.grid - tau.true[s]))}
    }
    ##
    load(paste0("Simulation_Results/Sim_Results_DGP=", DGP, "_N=", N, "_p=", p, ".RData"))
    ##
    slct <- Sim_DGP1_df$DGP==DGP & Sim_DGP1_df$N==N & Sim_DGP1_df$p==p 
    ##
    Sim_DGP1_df[slct,  colnames(Sim_DGP1_df)=="beta0.hat.PoI"]  <- sim.results[ , colnames(sim.results)=="beta0.hat.PoI"]
    Sim_DGP1_df[slct,  colnames(Sim_DGP1_df)=="beta1.hat.PoI"]  <- sim.results[ , colnames(sim.results)=="beta1.hat.PoI"]
    Sim_DGP1_df[slct,  colnames(Sim_DGP1_df)=="tau1.hat.PoI"]   <- t.grid[sim.results[ , colnames(sim.results)=="tau1.ind.hat.PoI"]]
    ##
    Sim_DGP1_df[slct,  colnames(Sim_DGP1_df)=="beta0.hat.LMcK"] <- sim.results[ , colnames(sim.results)=="beta0.hat.LMcK"]
    Sim_DGP1_df[slct,  colnames(Sim_DGP1_df)=="beta1.hat.LMcK"] <- sim.results[ , colnames(sim.results)=="beta1.hat.LMcK"]
    Sim_DGP1_df[slct,  colnames(Sim_DGP1_df)=="tau1.hat.LMcK"]  <- t.grid[sim.results[ , colnames(sim.results)=="tau1.ind.hat.LMcK"]]
    ##
    Sim_DGP1_df[slct,  colnames(Sim_DGP1_df)=="beta0"]          <- beta0
    Sim_DGP1_df[slct,  colnames(Sim_DGP1_df)=="beta1"]          <- beta
    Sim_DGP1_df[slct,  colnames(Sim_DGP1_df)=="tau1"]           <- t.grid[tau.ind.true]
    rm(sim.results)
  }
}

## #################################################
SimRes_DGP1 <- tibble::as_tibble(Sim_DGP1_df) %>% 
  dplyr::mutate(p=as.factor(p)) %>% 
  dplyr::mutate(N=as.factor(N)) %>% 
  tidyr::gather(key = Estim_Types, value = Estim_Values, beta0.hat.PoI:tau1.hat.LMcK, factor_key=TRUE) %>% 
  dplyr::mutate(Targets = fct_recode(Estim_Types,
                                     "beta[0]==1"  = "beta0.hat.PoI", "beta[0]==1" = "beta0.hat.LMcK",
                                     "beta[1]==4"  = "beta1.hat.PoI", "beta[1]==4" = "beta1.hat.LMcK",
                                     "tau[1]==0.5" = "tau1.hat.PoI",  "tau[1]==0.5"= "tau1.hat.LMcK")) %>%
  dplyr::mutate(Estimators = fct_recode(Estim_Types, 
                                        "PoI"  = "beta0.hat.PoI",  "PoI"  = "beta1.hat.PoI",  "PoI"  = "tau1.hat.PoI", 
                                        "LMcK" = "beta0.hat.LMcK", "LMcK" = "beta1.hat.LMcK", "LMcK" = "tau1.hat.LMcK")) %>% 
  group_by(p) %>% 
  dplyr::mutate(True_Values = case_when(
    Targets == "beta[0]==1"  ~ beta0,
    Targets == "beta[1]==4"  ~ beta1,
    Targets == "tau[1]==0.5" ~ tau1)) %>% 
  ungroup() %>% 
  dplyr::select(N, p, Estimators, Targets, Estim_Values, True_Values)


## #################
## DGP 1
## #################

# Boxplots
fff <- function(x) {
  r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9),na.rm=T)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


# alpha
Plot_DGP1_1 <- ggplot(
    data = SimRes_DGP1 %>% dplyr::filter(Targets == "beta[0]==1"),
    mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
    geom_hline(yintercept = 0,linetype=5)+
    stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
    stat_summary(fun.data = fff, geom="boxplot",position="dodge")
# get ylim for coord caretesion    
ylim.aux <- max(abs(c(ggplot_build(Plot_DGP1_1)$data[[2]]$ymin,ggplot_build(Plot_DGP1_1)$data[[2]]$ymax)))
ylim.val <- c(-ylim.aux,ylim.aux)


Plot_DGP1_1 <- Plot_DGP1_1+
   coord_cartesian(ylim = ylim.val)  +  
    labs(x = "", y="", title=expression(paste(paste(hat(alpha)-alpha*phantom(hat(beta)[1]))))) +
    scale_color_manual(name="",
                       labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                       values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
    theme_bw() +
    theme(plot.margin     = margin(t = 2, r = 0.5, b = 0, l = -3, unit = "mm")) +
    theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
    theme(legend.position = "none")+ 
theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))



Plot_DGP1_1


# beta1
Plot_DGP1_2 <- ggplot(
    data = SimRes_DGP1 %>% dplyr::filter(Targets == "beta[1]==4"),
    mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
    geom_hline(yintercept = 0,linetype=5)+
    stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
    stat_summary(fun.data = fff, geom="boxplot",position="dodge")

# get ylim for coord caretesion    
ylim.aux <- max(abs(c(ggplot_build(Plot_DGP1_2)$data[[2]]$ymin,ggplot_build(Plot_DGP1_2)$data[[2]]$ymax)))
ylim.val <- c(-ylim.aux,ylim.aux)
#
Plot_DGP1_2 <- Plot_DGP1_2 +
   coord_cartesian(ylim = ylim.val)  +  
    labs(x = "", y="", title=expression(paste(hat(beta)[1]-beta[1]*phantom(hat(beta)[1])))) +
    scale_color_manual(name="",
                       labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                        values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
    theme_bw() +
    theme(plot.margin     = margin(t = 2, r = 0.5, b = 0, l = -3, unit = "mm")) +
    theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
    theme(legend.position = "none")+
theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))

    
Plot_DGP1_2


# tau
Plot_DGP1_3 <- ggplot(
    data = SimRes_DGP1 %>% dplyr::filter(Targets == "tau[1]==0.5"),
    mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
    geom_hline(yintercept = 0,linetype=5) +
    stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
    stat_summary(fun.data = fff, geom="boxplot",position="dodge")

# get ylim for coord caretesion    
ylim.aux <- max(abs(c(ggplot_build(Plot_DGP1_3)$data[[2]]$ymin,ggplot_build(Plot_DGP1_3)$data[[2]]$ymax)))
ylim.val <- c(-ylim.aux,ylim.aux)
#
Plot_DGP1_3 <- Plot_DGP1_3+
   coord_cartesian(ylim = ylim.val)  +  
coord_cartesian(ylim = c(-0.005, 0.005))+
    labs(x = "", y="", title=expression(paste(paste(hat(tau)[1]-tau[1]*phantom(hat(beta)[1]))))) +
    scale_color_manual(name="",
                       labels = c("100"="p=100", "500"="p=500", "1000"="p=1000"),
                       values = c("100"=gray(.6), "500"=gray(.4), "1000"=gray(0))) +
    theme_bw() +
    theme(plot.margin     = margin(t = 2, r = 0.5, b = 0, l = -3, unit = "mm")) +
    theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
    theme(legend.position = "none")+
theme(axis.text.y = element_text(angle = 90,hjust=0.5 ))

Plot_DGP1_3



################## 
# get legend:
Plot_DGP1_2.aux <- ggplot(
  data = SimRes_DGP1 %>% dplyr::filter(Targets == "beta[1]==4"),
  mapping = aes(x = N, y = Estim_Values-True_Values,  linetype = Estimators, color = p)) +
  geom_hline(yintercept = 0,linetype=5)+
  stat_summary(fun.data=fff, position="dodge",geom="errorbar")+
  stat_summary(fun.data = fff, geom="boxplot",position="dodge")+
  # geom_errorbar(aes(ymin=min, ymax= max), position ="dodge", width=0.5)
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
  theme(plot.title = element_text(size=12, margin = margin(b = -1))) +
  theme(legend.position = "bottom")+
  theme(legend.text=element_text(size=13), legend.title=element_text(size=13))+
  theme(legend.key.size = unit(1.5,"line"))


Plot_DGP1_2.aux
mylegend <- g_legend(Plot_DGP1_2.aux)
#################

Plot_DGP1 <- grid.arrange(arrangeGrob(Plot_DGP1_3, Plot_DGP1_1, Plot_DGP1_2, nrow=1), mylegend, nrow=2, heights=c(6,0.5))

ggsave(filename = "Plot_DGP_1.pdf", plot = Plot_DGP1, device = "pdf", path = plot_path, width = 10, height = 4) #8,6
