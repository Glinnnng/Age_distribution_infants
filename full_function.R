
#1. Estimating RSV hospitalisation by chronological month of age----

#### Calculate the likelihood for each study
genLikelihood_each <- function(each, param) {
  each <- cbind(
    each,
    data.frame(
      param = param,
      group = cut(1:12, c(0,each$Age_end), right = TRUE)
    ) %>% group_by(group) %>%
      dplyr::summarise(param_sum = sum(param))
  ) 
  return(dmultinom(each$RSV_count, prob = each$param_sum, log = TRUE))
} 

#### Calculate the sum of the likelihood of all studies
genLikelihood <- function(inputdata, id = "Study_ID", param) {
  return(sum(unlist(by(inputdata, inputdata[id], genLikelihood_each, param = param))))
}

#### Random walks on the parameters from a Dirichlet distribution parameterised based on the number of RSV cases of the data
genProposal <- function(inputdata = inputdata, param, step_dirichlet = step_dirichlet) {
  n <- sum(inputdata$RSV_count)
  alpha_dirichlet <- param * n
  alpha_dirichlet_update <- alpha_dirichlet + runif(12, -(n+12)*step_dirichlet, (n+12) * step_dirichlet)
  alpha_dirichlet_update <- ifelse(alpha_dirichlet_update<1, 1, alpha_dirichlet_update)
  return(as.vector(rdirichlet(1, alpha = alpha_dirichlet_update))) 
} 

#### Metropolis–Hasting’s algorithm for estimating model parameters
genMH <- function(log_likelihood_function, initial_params, n_iterations, inputdata, step_dirichlet) {
  chain <- matrix(NA, nrow = n_iterations, ncol = 12)
  chain[1, ] <- initial_params
  current_log_likelihood <- log_likelihood_function(param = initial_params, inputdata = inputdata)
  pb <- progress_bar$new(total = n_iterations, clear = TRUE, format = "  [:bar] :percent :etas")
  pb$tick()
  for (i in 2:n_iterations) {
    proposal <- genProposal(inputdata = inputdata, param = chain[i-1,], step_dirichlet = step_dirichlet)
    proposal_log_likelihood <- log_likelihood_function(param = proposal, inputdata = inputdata)
    acceptance_ratio <- exp(proposal_log_likelihood - current_log_likelihood) 
    if (runif(1) < acceptance_ratio) { 
      chain[i, ] <- proposal
      current_log_likelihood <- proposal_log_likelihood
    } else {
      chain[i, ] <- chain[i - 1, ]
    }
    pb$tick()
  }
  return(chain)
} 
#combine 3 chains
genMH_C <- function(log_likelihood_function = genLikelihood, initial_params, n_iterations=6000, inputdata, step_dirichlet = 0.005) {
  chain_1 <- genMH(initial_params = initial_params[[1]],log_likelihood_function = log_likelihood_function, n_iterations = n_iterations, 
                   inputdata = inputdata, step_dirichlet = step_dirichlet)
  chain_1 <- cbind(chain_1,matrix(1,nrow = n_iterations,ncol = 1))
  colnames(chain_1) <- c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","chain_index")
  chain_2 <- genMH(initial_params = initial_params[[2]],log_likelihood_function = log_likelihood_function, n_iterations = n_iterations, 
                   inputdata = inputdata, step_dirichlet = step_dirichlet)
  chain_2 <- cbind(chain_2,matrix(2,nrow = n_iterations,ncol = 1))
  colnames(chain_2) <- c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","chain_index")
  chain_3 <- genMH(initial_params = initial_params[[3]],log_likelihood_function = log_likelihood_function, n_iterations = n_iterations, 
                   inputdata = inputdata, step_dirichlet = step_dirichlet)
  chain_3 <- cbind(chain_3,matrix(3,nrow = n_iterations,ncol = 1))
  colnames(chain_3) <- c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","chain_index")
  result <- as.data.frame(rbind(chain_1,chain_2,chain_3))
  result$index <- 1:n_iterations
  return(result)
} 

#### burn in and thin for per chain (remaining 500 samples per chain) and calculate point estimate and 95% credible interval
# Three chains with different initial values were used: 
# in chain 1, all parameters were set to 1/12; 
# in chain 2, the first six parameters (corresponding to 0-<6 months) were each set to 0.1, and the second six parameters (corresponding to 6-<12 months) were set to 0.4/6; 
# in chain 3, the first six parameters were set to 0.4/6, and the second six parameters were set to 0.1.
genRes <- function(inputdata, n.iteration, n.burnin, n.thin, Group){
  chain_each <- genMH_C(inputdata=inputdata,initial_params = list(rep(1/12,12),c(rep(0.1, 6), rep(0.4/6,6)),c( rep(0.4/6,6), rep(0.1, 6))))
  chain_each$chain_index <- as.factor(chain_each$chain_index)
  chain_each_long <- chain_each %>% pivot_longer(cols = p1:p12,names_to = "parameter", 
                                       names_prefix = "p")
  chain_each_long$parameter <- as.numeric(chain_each_long$parameter)
  
  #Plot trace plots
  # IndexPlot <- ggplot(data = chain_each_long)+
  #   geom_point(aes(x=index, y=value,color=chain_index), size= 0.7, alpha=0.3)+
  #   theme_bw()+
  #   facet_wrap(~ parameter)
  # ggsave(IndexPlot,filename = paste("results/IndexPlots/","IndexPlot","_", Group ,".tiff",sep = ""), width = 10, height = 5) 
  
  thin_random <- sample(1:n.thin, size = 3, replace = TRUE) 
  thin_index <- data.frame(
    chain_index = as.factor(rep(1:3, each = (n.iteration - n.burnin)/n.thin)), 
    index_new = 1:((n.iteration - n.burnin)/n.thin),
    index_seq = 1:((n.iteration - n.burnin)/n.thin*3),
    index = c(
      seq(from  = n.burnin+thin_random[1], to = n.iteration, by = n.thin),
      seq(from  = n.burnin+thin_random[2], to = n.iteration, by = n.thin),
      seq(from  = n.burnin+thin_random[3], to = n.iteration, by = n.thin))
  ) 
  thin_chain <- left_join(thin_index, chain_each_long)
  thin_chain<- do.call(rbind,by(thin_chain, thin_chain$index_seq, function(x) {x$value_cumsum = cumsum(x$value); return(x)}))
  save(thin_chain,file = paste("thin_chain","_", Group ,".RData",sep = ""))
  paramSummary <- thin_chain %>% group_by(parameter) %>% dplyr::summarise(param.est=median(value),
                                                                              param.lci= quantile(value,0.025),
                                                                              param.uci= quantile(value,0.975),
                                                                              param_cumsum.est=median(value_cumsum),
                                                                              param_cumsum.lci= quantile(value_cumsum,0.025),
                                                                              param_cumsum.uci= quantile(value_cumsum,0.975)) %>% mutate(Group=Group)

  #Determine whether at least one study within the subgroup has 12 age groups so that the model can converge.
  if(sum(inputdata$N_subAge==12)>0){
  return(paramSummary)
}else{
  return(NULL)
  }  
} 

#### Plots to show age distribution
genAgePlots <- function(inputdata,type,analysis){
  ggsave(ggplot(data=inputdata, aes(x=parameter, y=param.est, ymin=param.lci, ymax=param.uci, fill=Group, linetype=Group)) + 
           geom_line()+
           geom_ribbon(alpha=0.2)+
           labs(x="Age",y="Proportion",fill=type,linetype=type)+
           scale_x_continuous(breaks = seq(1,12,1),labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                                              "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
           theme(legend.position = "bottom",
                 legend.justification = "center",
                 text = element_text(size = 20),
                 panel.background = element_blank(),
                 axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
                 axis.line.y = element_line(linetype=1,color="black",linewidth=0.25)),
           #guides(fill = guide_legend(nrow = 2),linetype = guide_legend(nrow = 2)),
         filename = paste("results/","AgeDistribution","_",type,"_",analysis,".pdf",sep = ""), width = 14, height = 7)
} 

# load_and_rename <- function(file_name, new_name) {
#   load(file_name)
#   new_name <-  thin_chain
#   rm(thin_chain)
#   return(new_name)
# }

#2. Development and validation of a prediction tool for estimating RSV hospitalisation by birth month----

#### Development a prediction tool for estimating RSV hospitalisation by birth month
#prediction tool calculating the proportion and cumulative proportion of RSV hospitalisation by birth month and month of age
#Risks and cumulative risks for different birth months in a single sample
genBirthEach <- function(byAge, byMonth) {
  res.matrix <- cbind(
    expand.grid(byAge.value = byAge, byMonth.value = byMonth),
    expand.grid(byAge.label = 1:12, byMonth.label = 1:12)
  )
  res.matrix$p_combined <- res.matrix$byAge.value *res.matrix$byMonth.value
  res.matrix$byBirthmonth.label <- ifelse(res.matrix$byMonth.label-res.matrix$byAge.label+1 > 0,
                                          res.matrix$byMonth.label-res.matrix$byAge.label+1, 
                                          res.matrix$byMonth.label-res.matrix$byAge.label+1+12)
  res.matrix <- res.matrix[order(res.matrix$byAge.label, decreasing = FALSE),]
  res.matrix <- do.call(rbind,by(res.matrix, res.matrix$byBirthmonth.label, function(x) {x$p_cumsum = cumsum(x$p_combined); return(x)}))
  return(res.matrix)
} 
# Risks and cumulative risks for different birth months for all the 1500 samples

genBirthRes <- function(IC,Seasonality_each){
  if(IC=="H"){age_input <- thin_chain_H
  }else{if(IC=="UM"){age_input <- thin_chain_UM
  }else{age_input <- thin_chain_LM}
  }
  
  age_input$parameter <- as.numeric(age_input$parameter)
  age_input <- age_input %>% arrange(index_seq,parameter)
  seasonality_input <- cbind(rdirichlet((n.iteration - n.burnin)/n.thin*3,Seasonality_each+1),
                             matrix(seq(1,(n.iteration - n.burnin)/n.thin*3,1), nrow = (n.iteration - n.burnin)/n.thin*3, ncol=1))
  colnames(seasonality_input) <- c("q1","q2","q3","q4","q5","q6","q7","q8","q9","q10","q11","q12","index_seq")
  seasonality_input <- as.data.frame(seasonality_input) %>% pivot_longer(cols = q1:q12, names_to = "parameter",
                                                                         names_prefix = "q",values_to = "value2") %>%  
    mutate(parameter=as.numeric(parameter)) %>% arrange(index_seq,parameter)
  
  
  input_combined <- left_join(age_input,seasonality_input)
  genBirthmonthMC <- function(input_combined_each){
    res <- genBirthEach(byAge = input_combined_each$value, byMonth = input_combined_each$value2)
    res$index_seq <- input_combined_each$index_seq[1]
    return(res)
  } 
  res.matrix.MC <- do.call(rbind, by(input_combined, input_combined$index_seq, genBirthmonthMC))
  } 
# Summary for the median and 95% credible interval
genBirthRes.summary <- function(IC,Seasonality_each){
  res.matrix.MC <- genBirthRes(IC=IC,Seasonality_each=Seasonality_each)
Res.summary <- res.matrix.MC %>% group_by(byBirthmonth.label, byAge.label) %>% 
  dplyr::summarise(p_combined.est = median(p_combined),
                   p_combined.lci = quantile(p_combined, 0.025),
                   p_combined.uci = quantile(p_combined, 0.975),
                   p_cumsum.est = median(p_cumsum),
                   p_cumsum.lci = quantile(p_cumsum, 0.025),
                   p_cumsum.uci = quantile(p_cumsum, 0.975)
                   )
Res.summary$IC <- IC
return(Res.summary)
} 

#### Validation of a prediction tool for estimating RSV hospitalisation by birth month using validation dataset.
genValidBirthRes <- function(seasonality_each){
  genBirthResEach <- function(seasonality){
    res <- genBirthRes.summary(IC=seasonality$Income, Seasonality_each=as.matrix(seasonality[,c(6:17)]))
    res$Country <- seasonality$Country
    res$Site <-  seasonality$Site
    res$Study_ID <-  seasonality$Study_ID
    return(res)
  }
  ValidRes <- do.call(rbind,by(seasonality_each, seasonality_each$Study_ID, genBirthResEach))
  return(ValidRes)
} 

####  correlation between the model-predicted and observed results
#  linear regression analysis and Pearson's correlation coefficients between the model-predicted and observed results 
# across different study sites stratified by country income level
genCorPlots <- function(cor_data){
  observed_cum12P <- cor_data$Birthmonth_cumulative_12m
  predicted_cum12P <- cor_data$p_cumsum.est
  CorPlot <- ggplot(data=cor_data,aes(x=observed_cum12P, y=predicted_cum12P))+
    geom_point(size=4,alpha=0.6,aes(color=cor_data$marking))+
    scale_color_manual(values = c("1"="#8A2BE2","0"="#b5d5e6"))+
    geom_smooth(method = "lm", color = "#b74c4c", fill = "#d69999")+ 
    labs(x="actual cumulative proportion",
         y="model-predicted cumulative proportion")+
    theme(
      legend.position ="none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      text = element_text(size = 25),
      axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
      axis.line.y = element_line(linetype=1,color="black",linewidth=0.25))+
    stat_cor(method = "pearson", hjust = 0.1, vjust = -1,label.x.npc = 0.1, label.y.npc = 0.8,size=8)+
    xlim(min(observed_cum12P),max(observed_cum12P)+0.001)+
    ylim(min(observed_cum12P),max(observed_cum12P)+0.001)
  ggsave(CorPlot ,filename = paste("results/","CorPlot","_",cor_data$Income[1],".pdf",sep = ""), width = 7, height = 7)
  
}

genCorPlots_add <- function(cor_data){
  ggsave(ggplot(data=cor_data, aes(x=byBirthmonth.label, y=p_cumsum.est, ymin=p_cumsum.lci, ymax=p_cumsum.uci)) +
           geom_line(colour="#CC0000",linetype="dashed",size=0.2)+
           geom_ribbon(alpha=0.5,fill="#d69999")+
           geom_line(colour="#3399FF",aes(x=byBirthmonth.label, y=Birthmonth_cumulative_12m),size=0.3)+
           labs(x="Birth month",y="Proportion")+
           scale_x_continuous(breaks = seq(1,12,1),labels = c("Jan","Feb","Mar","Apr","May","Jun",
                                                              "Jul","Aug","Sep","Oct","Nov","Dec"))+
           theme(legend.position = "right",
                 text = element_text(size = 25),
                 panel.background = element_blank(),
                 axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
                 axis.line.y = element_line(linetype=1,color="black",linewidth=0.25),
                 axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
           ylim(0,0.15)+
           facet_wrap(~ Country),  
         filename = paste("results/","Corplot_add","_",cor_data$Income[1],".pdf",sep = ""), width = 7, height =7)
  
}

#### Pearson’s correlation coefficients and 95% confidence intervals (CI)
# between the model-predicted and the observed proportion and cumulative proportion of RSV hospitalisations
# by birth month for different months of age
genCorTest <-function(cor_data, value) {
  
  if(value=="p_combined"){
    observed_p <- cor_data$p_combined
    predicted_p <- cor_data$p_combined.est
  }else{
    observed_p <- cor_data$p_cumsum
    predicted_p <- cor_data$p_cumsum.est
  }
  
  cor_df <- data.frame(byAge.label = NA, Income = NA, correlation = NA, ci_lower = NA, ci_upper = NA)
  row.names(cor_df) <- NULL
  
  
  cor_df$byAge.label <- cor_data$byAge.label[1]
  cor_df$correlation <- cor.test(observed_p, predicted_p)$estimate
  cor_df$ci_lower <-  cor.test(observed_p, predicted_p)$conf.int[1]
  cor_df$ci_upper <-  cor.test(observed_p, predicted_p)$conf.int[2]
  cor_df$Income <- cor_data$IC[1]
  return(cor_df)
}

#3. Predicting RSV hospitalisation by birth month----

genCountryBirthRes <- function(seasonality_each){
  genBirthResEach <- function(seasonality){
    res <- genBirthRes.summary(IC=seasonality$Income, Seasonality_each=as.matrix(seasonality[,c(3:14)]))
    res$Country <- seasonality$Country
    return(res)
  }
  CountryRes <- do.call(rbind,by(seasonality_each, seasonality_each$Country, genBirthResEach))
  return(CountryRes)
} 

# Predicted cumulative proportion of RSV hospitalisations in the first six month of life and first year of life by birth month
# and proportion of RSV hospitalisations by calendar month
genPredictPlots_lat <- function(predict_data){
  predict_data6 <- predict_data %>% filter(byAge.label == 6)
  predict_data12 <- predict_data %>% filter(byAge.label == 12)
  SeasonData <- predict_data %>% group_by(Country, byBirthmonth.label) %>% slice(1)
  ggsave(ggplot() +
           geom_line(data=predict_data6,aes(x=byBirthmonth.label, y=p_cumsum.est),colour="#008000",linetype="dashed",size=0.3)+ #橙色
           geom_ribbon(data=predict_data6,aes(x=byBirthmonth.label,ymin=p_cumsum.lci, ymax=p_cumsum.uci),alpha=0.4,fill="#99cc99")+
           geom_line(data=predict_data12,aes(x=byBirthmonth.label, y=p_cumsum.est),colour="#CC0000",linetype="dashed",size=0.3)+ #绿色
           geom_ribbon(data=predict_data12,aes(x=byBirthmonth.label,ymin=p_cumsum.lci, ymax=p_cumsum.uci),alpha=0.4,fill="#d69999")+
           geom_line(data = SeasonData,aes(x=byBirthmonth.label,y=proportion),colour="black",size=0.2)+
           labs(x="Month",y="Proportion")+
           scale_x_continuous(breaks = seq(1,12,1),labels = c("Jan","Feb","Mar","Apr","May","Jun",
                                                              "Jul","Aug","Sep","Oct","Nov","Dec"))+
           ylim(0,0.4)+
           theme(
             text = element_text(size = 27),
             panel.background = element_blank(),
             axis.text.x = element_text(angle = 60,vjust = 0.5,hjust = 0.5,size = 18),
             axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
             axis.line.y = element_line(linetype=1,color="black",linewidth=0.25))+
           facet_wrap(~ Country,ncol = 6),  
         filename = paste("results/","predict","_",predict_data$lat_group[1],".pdf",sep = ""), width = 20, height = ceiling(n_distinct(predict_data$Country) / 6)*3.8)
} 

# Model-predicted month-by-month age distribution of RSV hospitalisations by birth month for different countries

genAgeBirthPlot_lat <- function(predict_data){
  
  predict_data$Month <- factor(predict_data$Month,levels = c("Jan","Feb","Mar","Apr","May","Jun",
                                                             "Jul","Aug","Sep","Oct","Nov","Dec"))
  mycol12 <- c("#0b0b0b","#4c9568", "#7fb961",
               "#b0d45d","#76afda","#abddff",
               "#b20000","#a14462","#f06152",
               "#e8743c","#ffc556","#ffe788")
  ggsave(ggplot() +
           geom_line(data=predict_data,aes(x=byAge.label, y=p_combined.est,group=Month,colour=Month),size=0.5)+ 
           labs(x="Age",y="Proportion",colour="Birth month")+
           scale_x_continuous(breaks = seq(1,12,1),labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                                              "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
           scale_color_manual(values = mycol12, guide = guide_legend(nrow = 1)) +
           theme(
             legend.position ="none",
             text = element_text(size = 27),
             panel.background = element_blank(),
             axis.text.x = element_text(angle = 60,vjust = 0.5,hjust = 0.5,size = 18),
             axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
             axis.line.y = element_line(linetype=1,color="black",linewidth=0.25),
             plot.margin = unit(c(1, 1, 1, 1), "cm"))+
           facet_wrap(~ Country,ncol = 6),  
         filename = paste("results/","predict_AgeBirth","_",predict_data$lat_group[1],".pdf",sep = ""), width = 24, height = ceiling(n_distinct(predict_data$Country) / 6)*3.8)
} 

#### Risk ratios between the birth months with the highest and lowest cumulative RSV hospitalisation proportions for different month of age

# Risk ratio for each country
genBirthRes.Ratio.country <-function(x){
  res.matrix.max <- x %>% group_by(byAge.label,index_seq) %>% filter(p_cumsum==max(p_cumsum)) %>% ungroup()
  res.matrix.min <- x %>% group_by(byAge.label,index_seq) %>% filter(p_cumsum==min(p_cumsum)) %>% ungroup()
  res.matrix.ratio <-res.matrix.max %>%  left_join(res.matrix.min,by = c("byAge.label", "index_seq","Country", "IC", "lat", "long", "lat_group", "IC_group")) %>% 
    mutate(ratio=p_cumsum.x/p_cumsum.y) %>% ungroup()
  return(res.matrix.ratio )
}







