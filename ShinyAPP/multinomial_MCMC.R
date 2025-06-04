n.iteration <- 6000
n.burnin <- 1000
n.thin <- 10

load("thin_chain.RData")
load("World_ICdata.RData")

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
  } 


