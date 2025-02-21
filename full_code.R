#rm(list=ls())
# Please confirm that the following  packages have been installed.
library(readxl)
library(MCMCpack)
library(progress)
library(coda)
library(ggplot2)
library(tidyverse)
library(rgdal)
library(ggpubr)

#### Considering the data licences, only the full code is shown here without the corresponding data. 
#### We have given simulation data in the 'working example for illustration' file, which allows to run the core code of Bayesian model.

#0. data import----

# Exclude studies reported pooled results from multiple countries across different income levels
# and studies that report data after 2019
mixedID <- c("L61501", "L68601") 
after2019 <-c("L58201","L64201")
metaData <- read_excel("SCEQ_v4.2.xlsx", sheet = "Cataloguing")%>% 
  filter(!is.na(Study_ID)& !(Study_ID %in% mixedID )& !(Study_ID %in% after2019 ))

#data for main analysis
df <- read_excel("SCEQ_v4.2.xlsx", sheet = "Extraction") 
df <- df %>% filter(!is.na(Study_ID) & !(Study_ID %in% mixedID )) %>% left_join(metaData,by=c("Study_ID","m_analysis"))
df.IP <- df %>% filter(Setting=="Hospital"& m_analysis==1) 

#Country and latitude/longitude information for each study
LocationInfo <- read_excel("LocationInfo.xlsx", sheet = "LocationInfo")
LocationInfo <- LocationInfo %>% inner_join(metaData[,c("Study_ID", "Income")]) #Country income levels for the median year of the study period

LocationInfo$lat_group <- NA
for (i in 1:nrow(LocationInfo)) {
  if (!is.na(LocationInfo[i, "lat"])) {
    if (LocationInfo[i, "lat"] > 35) {
      LocationInfo[i, "lat_group"] <- "Temperate_region"
    } else if (LocationInfo[i, "lat"] < 35 & LocationInfo[i, "lat"] > 23.5) {
      LocationInfo[i, "lat_group"] <- "Subtropical_region"
    } else if (LocationInfo[i, "lat"] < 23.5 & LocationInfo[i, "lat"] > -23.5) {
      LocationInfo[i, "lat_group"] <- "Tropical_region"
    } else if (LocationInfo[i, "lat"] < -23.5& LocationInfo[i, "lat"] > -35) {
      LocationInfo[i, "lat_group"] <- "Subtropical_region"
    }else if (LocationInfo[i, "lat"] < -35) {
      LocationInfo[i, "lat_group"] <- "Temperate_region"
    }
  }
}
df.IP <- df.IP %>% left_join(LocationInfo[,c("Study_ID","lat_group")],multiple = "first")
df.IP$Income_lat <- paste(df.IP$Income,df.IP$lat_group,sep = "_")

# 1.1 Estimating RSV hospitalisation by chronological month of age----

set.seed(1114)
n.iteration <- 6000
n.burnin <- 1000
n.thin <- 10
source("full_function.R")

#### Age distribution stratified by income level (income level in the median year of the study)

paramSummary_Income <- do.call(rbind, by(df.IP, df.IP$Income, function(data) {
  genRes(inputdata = data, Group = data$Income[1], n.iteration = n.iteration, n.burnin = n.burnin, n.thin = n.thin)
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

#### Age distribution stratified by latitude
set.seed(1115)
paramSummary_Lat <- do.call(rbind,by(df.IP,df.IP$lat_group,function(data) {
  genRes(inputdata = data, Group = data$lat_group[1], n.iteration = n.iteration, n.burnin = n.burnin, n.thin = n.thin)
}))

paramSummary_Lat$Group <- factor(paramSummary_Lat$Group,levels =
                                   c("Temperate_region","Subtropical_region","Tropical_region"),
                                 labels=c("Temperate region","Subtropical region","Tropical region"))

genAgePlots(inputdata = paramSummary_Lat,type ="regions",analysis = "main") #plots for age distribution

#### Simultaneously by country income and latitudinal group
set.seed(1116)
paramSummary_Income.Lat <- do.call(rbind,by(df.IP,df.IP$Income_lat,function(data) {
  genRes(inputdata = data, Group = data$Income_lat[1], n.iteration = n.iteration, n.burnin = n.burnin, n.thin = n.thin)
}))

paramSummary_Income.Lat$Income <- sapply(strsplit(as.character(paramSummary_Income.Lat$Group), "_"), function(x) x[1])
paramSummary_Income.Lat$Income <- factor(paramSummary_Income.Lat$Income,levels =c("H","UM","LM"),labels=c("HICs","UMICs","LMICs/LICs"))

paramSummary_Income.Lat$Latitude <- sapply(strsplit(as.character(paramSummary_Income.Lat$Group), "_"), function(x) paste(x[2], x[3], sep = " "))
paramSummary_Income.Lat$Latitude <- factor(paramSummary_Income.Lat$Latitude,levels =
                                   c("Temperate region","Subtropical region","Tropical region"),
                                 labels=c("Temperate region","Subtropical region","Tropical region"))

paramSummary_Income.Lat$Group <- factor(paramSummary_Income.Lat$Group,levels =
                                          c("H_Temperate_region","H_Tropical_region","H_Subtropical_region",
                                            "UM_Temperate_region","UM_Tropical_region","UM_Subtropical_region",
                                            "LM_Temperate_region","LM_Tropical_region","LM_Subtropical_region"),labels=
                                          c("HIC (Temperate region)","HIC (Tropical region)","HIC (Subtropical region)",
                                            "UMIC (Temperate region)","UMIC (Tropical region)","UMIC (Subtropical region)",
                                            "LMIC/LIC (Temperate region)","LMIC/LIC (Tropical region)","LMIC/LIC (Subtropical region)"))

ggsave(ggplot(data=paramSummary_Income.Lat, aes(x=parameter, y=param.est, ymin=param.lci, ymax=param.uci, fill=Latitude, linetype=Latitude)) + 
         geom_line()+
         geom_ribbon(alpha=0.2)+
         labs(x="Age",y="Proportion",fill="regions",linetype="regions")+
         scale_x_continuous(breaks = seq(1,12,1),labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                                            "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
         theme(legend.position = "bottom",
               legend.justification = "center",
               text = element_text(size = 20),
               panel.background = element_blank(),
               axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
               axis.line.y = element_line(linetype=1,color="black",linewidth=0.25))+
       facet_wrap(~paramSummary_Income.Lat$Income,nrow = 3),
       filename = paste("results/","AgeDistribution","_","facet",".pdf",sep = ""), width = 12, height = 14)



# #加载thin_chain
# thin_chain_H <- load_and_rename("thin_chain_H.RData", "thin_chain_H")
# thin_chain_UM <-load_and_rename("thin_chain_UM.RData", "thin_chain_UM")
# thin_chain_LM <-load_and_rename("thin_chain_LM.RData", "thin_chain_LM")

#2.1 Development and validation of a prediction tool for estimating RSV hospitalisation by birth month----

ICdata <-  read_excel("SCEQ_v4.2.xlsx", sheet = "-") %>% dplyr::select(Country,Income)
ValidationData <- read_excel("Validation.xlsx", sheet = "under1")[,-8] %>% left_join(ICdata) 
ValidationData$byBirthmonth.label =rep(1:12, length.out = nrow(ValidationData))
ValidationData_wide <- ValidationData[,-c(7,9)] %>% pivot_wider(names_from = Month,values_from = Seasonality) 

#### Validation of a prediction tool for estimating RSV hospitalisation by birth month using validation dataset
# model-predicted results
set.seed(2111)
ValidBirthResSummary <- genValidBirthRes(ValidationData_wide) 

# Cumulative proportion of RSV hospitalisations in the first year of life by birth month
ValidBirthRes12 <- ValidBirthResSummary %>% filter(byAge.label=="12") %>% 
  left_join(ValidationData,by=c("Study_ID","byBirthmonth.label","Country","Site"))
ValidBirthRes12 <- ValidBirthRes12 %>% mutate(Study_Label=paste(Study_ID," (",Site,",",Country,")",sep = ""))
ValidBirthRes12 <- ValidBirthRes12 %>% group_by(Study_ID) %>% mutate(total_RSV=sum(Seasonality))
ValidBirthRes12$IC <- factor(ValidBirthRes12$IC,levels =c("H","UM","LM","L"),labels=c("HICs","UMICs","LMICs/LICs","LICs"))

# Mark the study with the largest sample size in each group
ValidBirthRes12 <- do.call(rbind,by(ValidBirthRes12,ValidBirthRes12$IC,function(x){
  x$marking=ifelse(x$total_RSV==max(x$total_RSV),"1","0"); return(x)}))

# Visual comparison between the model-predictions and observed cumulative proportion of RSV hospitalisations in the first year of life by birth month
ggsave(ggplot(data=ValidBirthRes12, aes(x=byBirthmonth.label, y=p_cumsum.est, ymin=p_cumsum.lci, ymax=p_cumsum.uci)) +
         geom_line(colour="#CC0000",linetype="dashed",size=0.2)+
         geom_ribbon(alpha=0.5,fill="#d69999")+
         geom_line(colour="#3399FF",aes(x=byBirthmonth.label, y=Birthmonth_cumulative_12m),size=0.3)+
         labs(x="Birth month",y="Proportion")+
         scale_x_continuous(breaks = seq(1,12,1),labels = c("Jan","Feb","Mar","Apr","May","Jun",
                                                            "Jul","Aug","Sep","Oct","Nov","Dec"))+
         theme(legend.position = "right",
               text = element_text(size = 12),
               panel.background = element_blank(),
               axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
               axis.line.y = element_line(linetype=1,color="black",linewidth=0.25))+
         facet_wrap(~ IC+Study_Label,ncol = 4),  
       filename = paste("results/","Validplot",".pdf",sep = ""), width = 14, height = 16)

# Difference between model-predicted and observed peak birth months
byBirthmonth.label.real <- ValidBirthRes12  %>%
  group_by(Study_ID) %>% arrange(desc(Birthmonth_cumulative_12m),by_group=TRUE) %>% 
  slice(1) %>% dplyr::select(byBirthmonth.label)
byBirthmonth.label.predict <- ValidBirthRes12  %>%
  group_by(Study_ID) %>% arrange(desc(p_cumsum.est),by_group=TRUE) %>% 
  slice(1) %>% dplyr::select(byBirthmonth.label)
BirthMonth.diff <- left_join(byBirthmonth.label.real,byBirthmonth.label.predict ,by="Study_ID")
BirthMonth.diff$diff <- BirthMonth.diff$byBirthmonth.label.x-BirthMonth.diff$byBirthmonth.label.y
table(BirthMonth.diff$diff)

#2.2 correlation between the model-predicted and observed results by birth month----

#### Comparison between model-predicted and observed cumulative proportion of RSV hospitalisations in the first year of life 
# by birth month by country income level
by(ValidBirthRes12,ValidBirthRes12$Income,genCorPlots)


# selected sites with the highest number of RSV hospitalisations in the first year of life by income level for demonstration
ValidBirthRes12_most <- ValidBirthRes12 %>% group_by(Income) %>% filter(total_RSV==max(total_RSV))
ValidBirthRes12_most$Label <- paste(ValidBirthRes12_most$Site,", ",ValidBirthRes12_most$Country,sep = "")
by(ValidBirthRes12_most,ValidBirthRes12_most$Income,genCorPlots_add)


#### correlation between the model-predicted and the observed results of RSV hospitalisations by birth month for different months of age
# sample size more than 1000
ValidationData_byage <- read_excel("Validation.xlsx", sheet = "byage") #observed results
ValidBirthRes_byage <- ValidBirthResSummary %>% inner_join(ValidationData_byage,by=c("Study_ID","byBirthmonth.label","byAge.label"))

ValidBirthRes_byage <- do.call(rbind, by(ValidBirthRes_byage, 
                                         ValidBirthRes_byage[, c("Study_ID", "byBirthmonth.label")], 
                                         function(x) {x$p_cumsum = cumsum(x$p_combined)
                                           return(x)
                                         }))%>%  as.data.frame()

# Pearson’s correlation coefficients and 95% confidence intervals (CI)
# Proportion of RSV hospitalisations by birth month for different months of age
CorTest <- do.call(rbind, by(ValidBirthRes_byage, ValidBirthRes_byage[,c("IC","byAge.label")], function(data) {
  genCorTest(cor_data = data, value="p_combined") %>% as.data.frame() %>% mutate(group="Age")
}))

# Cumulative proportion of RSV hospitalisations by birth month for different months of age
CorTest_cum <- do.call(rbind, by(ValidBirthRes_byage, ValidBirthRes_byage[,c("IC","byAge.label")], function(data) {
  genCorTest(cor_data = data, value="p_cumsum") %>% as.data.frame() %>% mutate(group="cumAge")
}))

CorTest2 <- rbind(CorTest,CorTest_cum)
CorTest2$group <- factor(CorTest2$group,levels = c("Age","cumAge"),labels = c("Age","Cumulative age"))
CorTest2$Income <- factor(CorTest2$Income ,levels =c("H","UM","LM","L"),labels=c("HICs","UMICs","LMICs/LICs","LICs"))

ggsave(ggplot(data = CorTest2, aes(x = byAge.label, y = correlation, color=group)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(0.8), width = 0.8) +
  geom_point(position = position_dodge(0.8), size = 2) +
  labs(x="Age",y="Correlation coefficients (R)")+
  scale_x_continuous(breaks = seq(1,12,1),labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                                     "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
  ylim(0,1)+
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 10),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8),
    axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
    axis.line.y = element_line(linetype=1,color="black",linewidth=0.25)
  ) +
  facet_wrap(~Income,nrow = 3),
  filename = paste("results/CorTest",".pdf",sep = ""), width = 6, height = 8)



#3.1 Predicting RSV hospitalisation by birth month----

# Nationwide data only
SeasonData <- read_excel("Seasonality.xlsx", sheet = "Sheet1")
Country_location <- read_excel("Country_location.xlsx", sheet = "Sheet1")
ValidationData_nationwide <- ValidationData_wide %>% filter(Site=="Nationwide") %>% dplyr::select(Country,Income,Jan:Dec)
NationalSeasonData <- SeasonData %>% filter(Region=="nationwide") %>% dplyr::select(Country,Income,Jan:Dec) %>% rbind(ValidationData_nationwide)

# Predicting RSV hospitalisation by birth month
set.seed(3111)
BirthResSummary <- genCountryBirthRes(NationalSeasonData)
BirthResSummary <- BirthResSummary %>% left_join(Country_location,by="Country",multiple="any")
BirthResSummary$lat_group <- NA
for (i in 1:nrow(BirthResSummary)) {
  if (!is.na(BirthResSummary[i, "lat"])) {
    if (BirthResSummary[i, "lat"] > 35) {
      BirthResSummary[i, "lat_group"] <- "Temperate_region"
    } else if (BirthResSummary[i, "lat"] < 35 & BirthResSummary[i, "lat"] > 23.5) {
      BirthResSummary[i, "lat_group"] <- "Subtropical_region"
    } else if (BirthResSummary[i, "lat"] < 23.5 & BirthResSummary[i, "lat"] > -23.5) {
      BirthResSummary[i, "lat_group"] <- "Tropical_region"
    } else if (BirthResSummary[i, "lat"] < -23.5& BirthResSummary[i, "lat"] > -35) {
      BirthResSummary[i, "lat_group"] <- "Subtropical_region"
    }else if (BirthResSummary[i, "lat"] < -35) {
      BirthResSummary[i, "lat_group"] <- "Temperate_region"
    }
  }
}
BirthResSummary$IC <- factor(BirthResSummary$IC,levels =c("H","UM","LM","L"),labels=c("HICs","UMICs","LMICs","LICs"))
BirthResSummary$IC_group <- NA
BirthResSummary$IC_group <- ifelse(BirthResSummary$IC %in% c("LMICs","LICs"), "LMICs-LICs",as.character(BirthResSummary$IC))

# Describing the seasonality of nationwide data (observed)
Validation_national <- ValidationData %>% filter(Site=="Nationwide") %>% dplyr::select(Country,Month,Seasonality,byBirthmonth.label)
NationalSeasonData_long <-SeasonData %>% filter(Region=="nationwide") %>% dplyr::select(Country,Jan:Dec) %>% 
  pivot_longer(names_to = "Month",cols = Jan:Dec ,values_to = "Seasonality")
NationalSeasonData_long$byBirthmonth.label =rep(1:12, length.out = nrow(NationalSeasonData_long))
NationalSeasonData_long <- rbind(NationalSeasonData_long,Validation_national)
NationalSeasonData_long  <- do.call(rbind,by(NationalSeasonData_long , NationalSeasonData_long$Country, 
                                             function(x) {x$proportion = x$Seasonality/sum(x$Seasonality); return(x)}))
BirthResSummary <- BirthResSummary %>% left_join(NationalSeasonData_long,by=c("Country","byBirthmonth.label") )

#Peak months of age for different birth months for each country
BirthResSummary <- do.call(rbind,by(BirthResSummary, BirthResSummary[,c("Country","byBirthmonth.label")],
                                    function(x) {
                                      max_index = which.max(x$p_combined.est)
                                      x$peakAgeMonth = 0
                                      x$peakAgeMonth[max_index] = 1
                                      return(x)
                                    })) 

# Birth months with the highest RSV hospitalisation risk at the second month of life
peak2 <- BirthResSummary %>% filter(byAge.label==2)
num_birth <- peak2 %>% group_by(lat_group,Country) %>% summarise(num_birth=sum(peakAgeMonth))
summary(num_birth)

#Birth months with the highest RSV hospitalisation risk beyond the first six months of age,
peak6 <- BirthResSummary %>% filter(byAge.label>6)
num_birth6 <- peak6 %>% group_by(lat_group,Country) %>% summarise(num_birth=sum(peakAgeMonth))
summary(num_birth6)


# Predicted cumulative proportion of RSV hospitalisations in the first six month of life and first year of life by birth month and 
# proportion of RSV hospitalisations by calendar month, stratified by latitudinal group

by(BirthResSummary,BirthResSummary$lat_group,genPredictPlots_lat)

# Model-predicted month-by-month age distribution of RSV hospitalisations by birth month for different countries, stratified by latitudinal group
by(BirthResSummary,BirthResSummary$lat_group,genAgeBirthPlot_lat)

#legend
BirthResSummary_T <- BirthResSummary %>% filter(lat_group=="Temperate_region")
BirthResSummary_T$Month <- factor(BirthResSummary_T$Month,levels = c("Jan","Feb","Mar","Apr","May","Jun",
                                                           "Jul","Aug","Sep","Oct","Nov","Dec"))
mycol12 <- c("#0b0b0b","#4c9568", "#7fb961",
             "#b0d45d","#76afda","#abddff",
             "#b20000","#a14462","#f06152",
             "#e8743c","#ffc556","#ffe788")
plot1 <- ggplot() +
  geom_line(data=BirthResSummary_T,aes(x=byAge.label, y=p_combined.est,group=Month,colour=Month),size=1.2)+ 
  labs(x="Age",y="Proportion",colour="Birth month")+
  scale_x_continuous(breaks = seq(1,12,1),labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                                     "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
  scale_color_manual(values = mycol12) +
  theme(
    legend.position ="bottom",
    text = element_text(size = 25),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 60,vjust = 0.5,hjust = 0.5,size = 18),
    axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
    axis.line.y = element_line(linetype=1,color="black",linewidth=0.25),
    plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  guides(color = guide_legend(override.aes = list(size = 10),nrow = 1))+
  facet_wrap(~ Country,ncol = 6)

plot_legend <- get_legend(plot1)
plot_legend <- as_ggplot(plot_legend)
ggsave(plot_legend, filename = paste("results/","predict_AgeBirth","_","legend",".pdf",sep = ""), width = 24, height = 1)



#3.2 Risk ratios between the birth months with the highest and lowest cumulative RSV hospitalisation proportions----
# for the first 6 months of life and for the first year of life

# risk ratio for each country
set.seed(3112)
CountriesBirthResSummary612 <- do.call(rbind,by(NationalSeasonData,NationalSeasonData$Country,function(x){
  res <- genBirthRes(IC=x$Income,Seasonality_each=as.matrix(x[,c(3:14)]))
  res$Country <- x$Country
  return(res)
})) %>% filter(byAge.label==6 |byAge.label==12)
CountriesBirthResSummary612 <- CountriesBirthResSummary612 %>% left_join(BirthResSummary[,c(9:14)],by="Country",multiple="any") 
country.ratio <- do.call(rbind,by(CountriesBirthResSummary612,CountriesBirthResSummary612$Country,genBirthRes.Ratio.country ))
country.ratio.summary <- country.ratio %>% group_by(Country,byAge.label) %>%
  dplyr::summarise(ratio.est = median(ratio),
                   ratio.lci = quantile(ratio, 0.025),
                   ratio.uci = quantile(ratio, 0.975)
  )

country.ratio.summary$ratio_95CrI <- paste(round(country.ratio.summary$ratio.est,2)," (",round(country.ratio.summary$ratio.lci,2),"-",
                                          round(country.ratio.summary$ratio.uci,2),")",sep ="" )

# risk ratio by latitude group
# It is necessary to obtain risk ratios for all samples (1500 samples) in all countries and then combined by latitude.
RatioRes.summary.lat <- country.ratio %>% group_by(lat_group,byAge.label) %>%
  dplyr::summarise(ratio.est = median(ratio),
                   ratio.lci = quantile(ratio, 0.025),
                   ratio.uci = quantile(ratio, 0.975)
  )

RatioRes.summary.lat$ratio_95CrI <- paste(round(RatioRes.summary.lat$ratio.est,2)," (",round(RatioRes.summary.lat$ratio.lci,2),"-",
                                           round(RatioRes.summary.lat$ratio.uci,2),")",sep ="" )

#4. Descriptive analysis----
# observed data
#### Age distribution by years
AgeDistr_years <- read_excel("Validation.xlsx", sheet = "age_year")
AgeDistr_years <- do.call(rbind,by(AgeDistr_years, AgeDistr_years[,c("Study_ID","Year")], 
                                   function(x) {x$proportion = x$RSV_hospitalisations/sum(x$RSV_hospitalisations); return(x)}))

AgeDistr_years <- AgeDistr_years %>% left_join(ValidBirthRes12[,c("Study_ID","Study_Label")],by="Study_ID",multiple="any")

mycol15 <- c("#0b0b0b","#42465c","#356d67","#4c9568",
             "#7fb961","#b0d45d","#76afda",
             "#abddff","#b20000","#a14462",
             "#f06152","#eb998b","#e8743c",
             "#ffc556","#ffe788")
ggsave(ggplot(data=AgeDistr_years, aes(x=Age, y=proportion, color=factor(Year))) + 
         geom_line()+
         labs(x="Age",y="Proportion")+
         scale_x_continuous(breaks = seq(1,12,1),labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                                            "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
         scale_color_manual(values = mycol15) +
         theme(legend.position = "right",
               legend.title = element_blank(),
               legend.box = "vertical",
               text = element_text(size = 12),
               panel.background = element_blank(),
               axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8),
               axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
               axis.line.y = element_line(linetype=1,color="black",linewidth=0.25))+
         facet_wrap(~Study_Label,ncol = 2),
       filename = paste("results/years_age",".pdf",sep = ""), width = 8, height = 9 )

#### Summary of peak month of age with highest RSV hospitalisation proportion
AgeDistr_12 <- df.IP %>% filter(N_subAge==12)
AgeDistr_12 <- do.call(rbind,by(AgeDistr_12, AgeDistr_12$Study_ID, function(x) {x$proportion = x$RSV_count/sum(x$RSV_count); return(x)}))
peakMonth12 <- AgeDistr_12  %>%
  group_by(Study_ID) %>% arrange(desc(proportion),by_group=TRUE) %>% 
  slice(1) 
table(peakMonth12$Age_labeller)

#### Word map
World_ICdata <- read_excel("SCEQ_v4.2.xlsx", sheet = "World_IC") %>% dplyr::select(region=Country, Income) #Country income levels in 2019
World_ICdata$Income <- factor(World_ICdata$Income,levels =c("HICs","UMICs","LMICs","LICs"),labels=c("HICs","UMICs","LMICs","LICs"))
world_map <- map_data("world") %>% left_join(World_ICdata) %>% filter(region!="Antarctica")

ggsave(ggplot() +
         geom_map(
           data = world_map, map = world_map,
           aes(x=long, y=lat, map_id = region,fill=Income),
           color = "black",  linewidth = 0.1,alpha=0.7) +
         scale_fill_manual(values = c("#0392cf","#ffad60","#fdf498","#7bc043"),na.value = "#EEEEEE",na.translate=FALSE)+
         geom_point(data = LocationInfo,aes(x=long, y=lat),color="red",alpha = 0.7,size=3)+
         geom_hline(yintercept = c(23.5,-23.5), lty = 3, color = "darkgrey") +
         geom_hline(yintercept = c(35,-35), lty = 2, color = "darkgrey") +
         geom_hline(yintercept = 0, lty = 1, color = "darkgrey") +
         geom_text(aes(y = c(35,23.5,0,-23.5,-35), label =c("35°","23.5°","0°","-23.5°","-35°")), 
                   x = -Inf, hjust = -0.1, vjust = -1, size = 4,color = "darkgrey") +
         theme_void() +
         theme(legend.position = "bottom"),
       filename = paste("results/WorldMap_sites",".pdf",sep = ""), width = 12, height = 5)

#5. sensetivity analysis (age distribution by income level)----

#### Only studies containing PCR were included
df.IP.PCR <- df %>% filter(Setting=="Hospital"& s_analysis==1)

set.seed(1117)
paramSummary_PCR <- do.call(rbind, by(df.IP.PCR, df.IP.PCR$Income, function(data) {
  genRes2(inputdata = data, Group = data$Income[1], n.iteration = n.iteration, n.burnin = n.burnin, n.thin = n.thin,analysis="PCR")
}))

paramSummary_PCR$Group <- factor(paramSummary_PCR$Group,levels =c("H","UM","LM","L"),labels=c("HICs","UMICs","LMICs/LICs","LICs"))

genAgePlots(inputdata = paramSummary_PCR,type ="Income",analysis="PCR") 

# Exclude studies with total quality scores less than 6
Quality <- read_excel("SCEQ_v4.2.xlsx", sheet = "Quality Assessment") %>% dplyr::select(Study_ID,Score_total)
df.IP.High <- df.IP %>% left_join(Quality) %>% filter(Score_total>=6)

set.seed(1118)
paramSummary_High <- do.call(rbind, by(df.IP.High, df.IP.High$Income, function(data) {
  genRes2(inputdata = data, Group = data$Income[1], n.iteration = n.iteration, n.burnin = n.burnin, n.thin = n.thin,analysis="High")
}))

paramSummary_High$Group <- factor(paramSummary_High$Group,levels =c("H","UM","LM","L"),labels=c("HICs","UMICs","LMICs/LICs","LICs"))

genAgePlots(inputdata = paramSummary_High,type ="Income",analysis="High") 




