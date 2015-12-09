#####################################################################################
##### CONSTRUCTS NICE PLOTS OF TRAJECTORIES WITH SUMMARY STATISTICS
##### Plots all simulated trajectories, with user-specified quantiles
#####################################################################################

plot_trajectories = function(X,xlabs=NULL,myprobs=c(0.025,0.10,0.5,0.90,0.975),xtitle="",
                             ytitle="",basesize=12){
  
  if(is.null(xlabs)==TRUE){xlabs = 1:length(X[,1])}
  
  X = as.data.frame(t(X))
  
  nquants <- length(myprobs)
  myquants <- function(x){return(quantile(x,probs=myprobs))}
  
  ntime <- length(X[1,])
  ntraj <- length(X[,1])
  
  qX <- as.data.frame(apply(X,2,myquants))
  
  pX <- cbind(rep(myprobs,times=ntime),
              rep(xlabs,each=nquants), stack(qX)[,1])
  colnames(pX) <- c("quantile","year","values")
  
  pX <- as.data.frame(pX)
  
  pX$quantile <- factor(pX$quantile,levels=myprobs,
                        labels=c("L95%","L80%","Med","U80%","U95%"))
  
  # all sample trajectories
  seq_id <- rep(1:ntraj,ntime)
  year <- rep(xlabs,each=ntraj)
  y_new <- as.data.frame(cbind(seq_id,stack(X)[,-2],year))
  colnames(y_new)[2] <- "values"
  
  Y = list(Xquants=pX,Xall=y_new)
  
  myfig = ggplot2::ggplot(Y[[2]], ggplot2::aes(year,values,group=seq_id)) + 
    ggplot2::geom_line(colour="black",alpha=0.2,size=0.3) +
    ggplot2::geom_line(data=Y[[1]],ggplot2::aes(year,values,group=quantile,colour=quantile)) +
    ggplot2::theme_gray(base_size=basesize) +
    ggplot2::scale_color_manual(values=c("red","blue","green","blue","red")) +
    ggplot2::theme(legend.position="none") + 
    ggplot2::scale_x_continuous(xtitle) + ggplot2::scale_y_continuous(ytitle)
  
  return(myfig)
  
}