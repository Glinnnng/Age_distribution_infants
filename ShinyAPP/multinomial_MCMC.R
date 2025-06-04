n.iteration <- 6000
n.burnin <- 1000
n.thin <- 10

load("thin_chain.RData")
load("World_ICdata.RData")

genLikelihood_each <- function(each, param) {
  each <- cbind(
    each,
    data.frame(
      param = param,
      group = cut(1:12, c(0,each$Age_end), right = TRUE)
    ) %>% group_by(group) %>%
      dplyr::summarise(param_sum = sum(param))
  ) #给每一行数据加para_sum,及计算每个年龄组的概率
  return(dmultinom(each$RSV_count, prob = each$param_sum, log = TRUE))
} #计算每一个study的likelihood
genLikelihood <- function(inputdata, id = "Study_ID", param) {
  return(sum(unlist(by(inputdata, inputdata[id], genLikelihood_each, param = param))))
} #计算所有研究的likelihood之和
genProposal <- function(inputdata = inputdata, param, step_dirichlet = step_dirichlet) {
  n <- sum(inputdata$RSV_count)
  alpha_dirichlet <- param * n
  alpha_dirichlet_update <- alpha_dirichlet + runif(12, -(n+12)*step_dirichlet, (n+12) * step_dirichlet)
  alpha_dirichlet_update <- ifelse(alpha_dirichlet_update<1, 1, alpha_dirichlet_update)
  return(as.vector(rdirichlet(1, alpha = alpha_dirichlet_update))) #模拟从狄利克雷分布中抽样，狄利克雷分布为多项分布的概率分布，约束12个参数和为0-1，不是完全独立
} #计算下一步的方案（param应该是多少）
genMH <- function(log_likelihood_function, initial_params, n_iterations, inputdata, step_dirichlet) {
  chain <- matrix(NA, nrow = n_iterations, ncol = 12)
  chain[1, ] <- initial_params
  current_log_likelihood <- log_likelihood_function(param = initial_params, inputdata = inputdata)
  pb <- progress_bar$new(total = n_iterations, clear = TRUE, format = "  [:bar] :percent :etas")
  pb$tick()
  for (i in 2:n_iterations) {
    proposal <- genProposal(inputdata = inputdata, param = chain[i-1,], step_dirichlet = step_dirichlet)
    proposal_log_likelihood <- log_likelihood_function(param = proposal, inputdata = inputdata)
    acceptance_ratio <- exp(proposal_log_likelihood - current_log_likelihood) #当比值大于1时，一定往前走一步，当比值小于1时，以一定的概率往前走一步。
    if (runif(1) < acceptance_ratio) { 
      chain[i, ] <- proposal
      current_log_likelihood <- proposal_log_likelihood
    } else {
      chain[i, ] <- chain[i - 1, ]
    }
    pb$tick()
  }
  return(chain)
} #迭代6000次的parameter结果
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
} #合并三天链的结果


genRes <- function(chain_each, n.iteration, n.burnin, n.thin, IC){
  chain_each$chain_index <- as.factor(chain_each$chain_index)
  chain_each_long <- chain_each %>% pivot_longer(cols = p1:p12,names_to = "parameter", 
                                       names_prefix = "p")
  chain_each_long$parameter <- as.numeric(chain_each_long$parameter)
  IndexPlot <- ggplot(data = chain_each_long)+
    geom_point(aes(x=index, y=value,color=chain_index), size= 0.7, alpha=0.3)+
    theme_bw()+
    # theme(panel.grid.major = element_blank(), 
    #       panel.grid.minor = element_blank())
    facet_wrap(~ parameter)
  ggsave(IndexPlot,filename = paste("results/IndexPlots/","IndexPlot","_", IC ,".pdf",sep = ""), width = 10, height = 5) #绘制thin之前的parameter结果图
  
  thin_random <- sample(1:n.thin, size = 3, replace = TRUE) #十个里面随机抽一个，三条链抽三个
  thin_index <- data.frame(
    chain_index = as.factor(rep(1:3, each = (n.iteration - n.burnin)/n.thin)), #三条链，每条500个样
    index_new = 1:((n.iteration - n.burnin)/n.thin),
    index_seq = 1:((n.iteration - n.burnin)/n.thin*3),
    index = c(
      seq(from  = n.burnin+thin_random[1], to = n.iteration, by = n.thin),
      seq(from  = n.burnin+thin_random[2], to = n.iteration, by = n.thin),
      seq(from  = n.burnin+thin_random[3], to = n.iteration, by = n.thin))
  ) 
  return(left_join(thin_index, chain_each_long))
} #得到burnin和thin之后的数据框,app中没用上

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
} #单个样本不同出生月的风险和累积风险

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
                                                        names_prefix = "q",values_to = "value2") %>%  mutate(parameter=as.numeric(parameter)) %>% arrange(index_seq,parameter)
  

input_combined <- left_join(age_input,seasonality_input)
genBirthmonthMC <- function(input_combined_each){
  res <- genBirthEach(byAge = input_combined_each$value, byMonth = input_combined_each$value2)
  res$index_seq <- input_combined_each$index_seq[1]
  return(res)
} 
res.matrix.MC <- do.call(rbind, by(input_combined, input_combined$index_seq, genBirthmonthMC))
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
  } #1500个样本不同出生月的风险和累积风险的中位数和95%分位数

genCountryBirthRes <- function(seasonality_each){
  genBirthResEach <- function(seasonality){
    res <- genBirthRes(IC=seasonality$Income, Seasonality_each=as.matrix(seasonality[,c(3:14)]))
    res$Country <- seasonality$Country
    return(res)
  }
 CountryRes <- do.call(rbind,by(seasonality_each, seasonality_each$Country, genBirthResEach))
 return(CountryRes)
} #每个国家的144行结果（12个出生月*12个月份）,app中没用上
genBirthPlots <- function(BirthRes_each){
  BirthRes_each$byBirthmonth.abb <-factor(month.abb[BirthRes_each$byBirthmonth.label], levels = month.abb)
  RibbonPlot <- ggplot(data=BirthRes_each, aes(x=byAge.label)) + 
           geom_line(aes(y=p_combined.est),colour="#428bca",linetype="dashed",linewidth=0.2)+
           geom_ribbon(alpha=0.3,aes(ymin=p_combined.lci, ymax=p_combined.uci),fill="#428bca")+
           geom_line(aes(y=p_cumsum.est),colour="#5cb85c",linetype="dashed",linewidth=0.2)+
           geom_ribbon(alpha=0.3,aes(ymin=p_cumsum.lci, ymax=p_cumsum.uci),fill="#5cb85c")+        
           labs(x="Age",y="Proportion")+
           scale_x_continuous(breaks = seq(1,12,1),labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                                              "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
           theme(legend.position = "right",
                 text = element_text(size = 10),
                 panel.background = element_blank(),
                 axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8),
                 axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
                 axis.line.y = element_line(linetype=1,color="black",linewidth=0.25))+
           facet_wrap(~ byBirthmonth.abb)
  ggsave(RibbonPlot, filename = paste("results/RibbonPlots/",BirthRes_each$Country[1],"_","RP",".pdf",sep = ""), width = 10, height = 5)
  #热力图-瞬时
  HeatMaps <- ggplot(data = BirthRes_each, aes(x=as.factor(byAge.label), y=byBirthmonth.abb, fill=p_combined.est))+
    geom_tile()+
    scale_fill_gradient(low = "#f7fafc", high = "#005b96")+
    labs(x="Age",y="Birth Month",fill="Proportion")+
    # scale_y_discrete(labels=c("1"="Jan", "2"="Feb" ,"3"="Mar", "4"="Apr", "5"="May", "6"="Jun",
    #                           "7"="Jul" ,"8"="Aug" ,"9"="Sep", "10"="Oct","11"= "Nov","12"= "Dec"))+
    # geom_text(aes(label=p_combined.est),color = "black", size = 2)+
    scale_x_discrete(labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
    theme(panel.background = element_blank(),
          axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
  ggsave(HeatMaps, filename = paste("results/HeatMaps/",BirthRes_each$Country[1],"_","HM",".pdf",sep = ""), width = 7, height = 5)
  #热力图-累积
  HeatMaps_cum <- ggplot(data = BirthRes_each, aes(x=as.factor(byAge.label), y=byBirthmonth.abb, fill=p_cumsum.est))+
    geom_tile()+
    scale_fill_gradient(low = "#fffdfd", high = "#c0281d")+
    labs(x="Age",y="Birth Month",fill="Cumulative Proportion")+
    # scale_y_discrete(labels=c("1"="Jan", "2"="Feb" ,"3"="Mar", "4"="Apr", "5"="May", "6"="Jun",
    #                           "7"="Jul" ,"8"="Aug" ,"9"="Sep", "10"="Oct","11"= "Nov","12"= "Dec"))+
    # geom_text(aes(label=p_combined.est),color = "black", size = 2)+
    scale_x_discrete(labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
    theme(panel.background = element_blank(),
          axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
  ggsave(HeatMaps_cum, filename = paste("results/HeatMaps_cum/",BirthRes_each$Country[1],"_","HMc",".pdf",sep = ""), width = 7, height = 5)
 return(NULL)
  } #绘制出生月的可视化图,app中没用上

