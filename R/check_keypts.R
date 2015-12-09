#####################################################################################
##### CHECK OUTPUT OF "SIMULATE_AT_ELICITED_TIMES.R"
##### Checks correlations, marginal distributions, and sample trajectories
##### for generated sequences at key time-points
#####################################################################################

check_keypts <- function(y)
{
  nsims = length(y[,1])
  ntimes = length(y[1,])
  
  timelabels = 1:ntimes
  
  # correlation matrix
  y_cor = cor(y)
  
  # histograms
  seq_id = rep(1:nsims,ntimes)
  time <- rep(timelabels,each=nsims)
  y_new <- as.data.frame(cbind(seq_id,stack(as.data.frame(y))[,-2],time))
  colnames(y_new)[2] <- "values"
  binwidth_y <- (max(y_new$values)-min(y_new$values))/20
  y_new$time <- factor(y_new$time)
  y_hists <- ggplot2::ggplot(y_new, ggplot2::aes(x=values,colour=time)) + ggplot2::geom_freqpoly(binwidth=binwidth_y)
  
  # few sample trajectories
  ntraj = min(20,nsims)
  sampletraj = sample(1:nsims,ntraj,replace=FALSE)
  y_red = as.data.frame(y[sampletraj,])
  seq_id = rep(1:ntraj,ntimes)
  time = rep(timelabels,each=ntraj)
  y_new = as.data.frame(cbind(seq_id,stack(y_red)[,-2],time))
  colnames(y_new)[2] = "values"
  y_trajs = ggplot2::ggplot(y_new, ggplot2::aes(time,values,group=seq_id)) + ggplot2::geom_line() 
  
  myoutput = list(correlations=y_cor,histograms=y_hists,trajectories=y_trajs)
  
  return(myoutput)
}