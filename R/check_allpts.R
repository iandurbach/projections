#####################################################################################
##### CHECK OUTPUT OF "INTERPOLATE_BETWEEN_ELIC_TIMES.R"
##### Checks correlations, marginal distributions, and sample trajectories
##### for generated FULL sequences 
#####################################################################################

check_allpts = function(y,nperiods)
{
  y = t(y)
  
  nsims = length(y[,1])
  
  yearseq = 1:sum(nperiods)
  
  y_means = c()
  y_actual = c()
  
  for(j in 1:length(nperiods))
  {
    
    n_to_fill = nperiods[j]
    
    start_col = ifelse(j==1,1,end_col+1)
    end_col = start_col + nperiods[j] - 1
    
    y_to_use = y[,start_col:end_col]
    
    y_means_j = apply(y_to_use,1,mean)
    period = rep(paste(yearseq[start_col],"-",yearseq[end_col]),length(y_means_j))
    y_means = rbind(y_means,cbind.data.frame(y_means_j,period))
    
    y_actual_j = stack(as.data.frame(y_to_use))[-2]
    period = rep(paste(yearseq[start_col],"-",yearseq[end_col]),length(y_actual_j))
    year = rep(yearseq[start_col]:yearseq[end_col],each=nsims)
    seq_id = rep(1:nsims,yearseq[end_col]-yearseq[start_col]+1)
    y_actual = rbind(y_actual,cbind.data.frame(y_actual_j,period,year,seq_id))
    
  }
  
  # histogram of means within periods
  colnames(y_means) = c("values","period")
  binwidth = (max(y_means$values)-min(y_means$values))/20
  y_means_hists = ggplot2::ggplot(y_means, ggplot2::aes(x=values,colour=period)) + 
    ggplot2::geom_freqpoly(binwidth=binwidth)
  
  # histogram of actual values within periods
  colnames(y_actual) = c("values","period","year","seq_id")
  binwidth = (max(y_actual$values)-min(y_actual$values))/20
  y_actual_hists = ggplot2::ggplot(y_actual, ggplot2::aes(x=values,colour=period)) + 
    ggplot2::geom_freqpoly(binwidth=binwidth)
  
  # few sample trajectories
  sampletraj = sample(1:nsims,min(20,nsims),replace=FALSE)
  y_red = subset(y_actual,seq_id %in% sampletraj)
  y_few_trajs = ggplot2::ggplot(y_red, ggplot2::aes(year,values,group=seq_id)) + ggplot2::geom_line()
  
  # all sample trajectories
  y_all_trajs = ggplot2::ggplot(y_actual, ggplot2::aes(year,values,group=seq_id)) + ggplot2::geom_line()
  
  myoutput = list(histograms_means=y_means_hists,histograms_points=y_actual_hists,
                  sample_trajectories=y_few_trajs,all_trajectories=y_all_trajs)
  return(myoutput)
}