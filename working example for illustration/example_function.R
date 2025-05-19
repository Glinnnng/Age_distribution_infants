
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
  alpha <- 100 * each$param_sum +1
  likelihood <- extraDistr::ddirmnom(each$RSV_count, size = sum(each$RSV_count), alpha = alpha, log = TRUE) 
  return(likelihood)
}

#### Calculate the sum of the likelihood of all studies
genLikelihood <- function(inputdata, id = "Study_ID", param) {
  return(sum(unlist(by(inputdata, inputdata[id], genLikelihood_each, param = param))))
}

#### Random walks on the parameters from a Dirichlet distribution parameterised based on the number of RSV cases of the data
genProposal <- function(inputdata = inputdata, param, i, k) {
  n <- sum(inputdata$RSV_count)
  return(as.vector(gtools::rdirichlet(1, alpha = pmax(k * param, 1)))) 
}

#### Metropolis–Hasting’s algorithm for estimating model parameters
genMH <- function(log_likelihood_function, initial_params, n_iterations, inputdata, k_init) {
  chain <- matrix(NA, nrow = n_iterations, ncol = 12)
  chain[1, ] <- initial_params
  current_log_likelihood <- log_likelihood_function(param = initial_params, inputdata = inputdata)
  pb <- progress_bar$new(total = n_iterations, clear = TRUE, format = "  [:bar] :percent :etas")
  pb$tick()
  k <- k_init
  acceptance_count <- 0
  for (i in 2:n_iterations) {
    proposal <- genProposal(inputdata = inputdata, param = chain[i-1,], i = i, k = k)
    proposal_log_likelihood <- log_likelihood_function(param = proposal, inputdata = inputdata)
    acceptance_ratio <- exp(proposal_log_likelihood - current_log_likelihood)
    if (runif(1) < acceptance_ratio) {
      chain[i, ] <- proposal
      current_log_likelihood <- proposal_log_likelihood
      acceptance_count <- acceptance_count + 1
    } else {
      chain[i, ] <- chain[i - 1, ]
    }
    if (i %% 100 == 0) {
      current_acceptance <- acceptance_count / 100
      if (current_acceptance < 0.2) k <- k * 1.1
      if (current_acceptance > 0.5) k <- k * 0.9
      acceptance_count <- 0
      print(paste("Iter", i, "k =", round(k), "Acceptance =", round(current_acceptance, 2)))
    }
    pb$tick()
  }
  return(chain)
}

#combine 3 chains
genMH_C <- function(log_likelihood_function = genLikelihood, initial_params, n_iterations=600, inputdata, k_init = k_init) { 
  
  chain_1 <- genMH(initial_params = initial_params[[1]],log_likelihood_function = log_likelihood_function, n_iterations = n_iterations, 
                   inputdata = inputdata, k_init = k_init)
  chain_1 <- cbind(chain_1,matrix(1,nrow = n_iterations,ncol = 1))
  colnames(chain_1) <- c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","chain_index")
  chain_2 <- genMH(initial_params = initial_params[[2]],log_likelihood_function = log_likelihood_function, n_iterations = n_iterations, 
                   inputdata = inputdata, k_init = k_init)
  chain_2 <- cbind(chain_2,matrix(2,nrow = n_iterations,ncol = 1))
  colnames(chain_2) <- c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","chain_index")
  chain_3 <- genMH(initial_params = initial_params[[3]],log_likelihood_function = log_likelihood_function, n_iterations = n_iterations, 
                   inputdata = inputdata, k_init = k_init)
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
genRes <- function(inputdata, n.iteration, n.burnin, n.thin, Group, k_init ){
  chain_each <- genMH_C(inputdata=inputdata,initial_params = list(rep(1/12,12),c(rep(0.1, 6), rep(0.4/6,6)),c( rep(0.4/6,6), rep(0.1, 6))),
                        k_init=k_init )
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

