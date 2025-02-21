#rm(list=ls())
# Please confirm that the following packages have been installed.
require(readxl)
require(MCMCpack)
require(progress)
require(coda)
require(ggplot2)
require(tidyverse)
setwd("working example for illustration")

#### This is a working example for illustration using simulation data, which allows to run the core code of Bayesian model.####

#0. data import----

# simulation data
df <- read_excel("simulation data.xlsx", sheet = "Sheet1") 

# 1.1 Estimating RSV hospitalisation by chronological month of age----

set.seed(1114)
n.iteration <- 600 # Note this is deliberately set to a relatively small number.
n.burnin <- 100
n.thin <- 5
source("example_function.R")

#### Age distribution stratified by income level (income level in the median year of the study)
# This step may take a few minutes!
paramSummary_Income <- do.call(rbind, by(df, df$Income_level, function(data) {
  genRes(inputdata = data, Group = data$Income_level[1], n.iteration = n.iteration, n.burnin = n.burnin, n.thin = n.thin)
}))

paramSummary_Income$Group <- factor(paramSummary_Income$Group,levels =c("H","UM","LM","L"),labels=c("HICs","UMICs","LMICs/LICs","LICs"))

genAgePlots(inputdata = paramSummary_Income,type ="Income",analysis = "main") #plots for age distribution

#plot for age distribution of cumulative proportion
ggsave(ggplot(data=paramSummary_Income, aes(x=parameter, y=param_cumsum.est, ymin=param_cumsum.lci, ymax=param_cumsum.uci, fill=Group, linetype=Group)) + 
         geom_line()+
         geom_ribbon(alpha=0.2)+
         labs(x="Age",y="Cumulative Proportion",fill="Income",linetype="Income")+
         scale_x_continuous(breaks = seq(1,12,1),labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                                            "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
         theme(legend.position = "bottom",
               legend.justification = "center",
               text = element_text(size = 20),
               panel.background = element_blank(),
               axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
               axis.line.y = element_line(linetype=1,color="black",linewidth=0.25)),
       filename = paste("results/","AgeDistribution_Income_cum",".pdf",sep = ""), width = 14, height = 7)


