#####################################################################################
##### SIMULATION OF VALUES AT ASSESSED TIME-POINTS
##### General procedure for generating possibly correlated sequences
##### at elicited time-points. 
##### Inter-temporal relationships can be specified in two different ways
##### 1) Correlation matrices showing correlations between any two time points
##### 2) Binary indicators of "mean reversion", where a "1" means that extreme
##### (high or low) values in that time period should be followed by less-extreme
##### values in the following period.
#####################################################################################

simulate_at_elic_times = function(nsims=1000,Sigma,x_min,x_mode1,x_mode2,x_max,mean_rev=0)
{
  n1 = length(x_min)
  n2 = length(x_mode1)  
  n3 = length(x_mode2)
  n4 = length(x_max)
  ns1 = length(Sigma[1,])
  ns2 = length(Sigma[,1])
  
  if(length(unique(c(n1,n2,n3,n4,ns1,ns2)))>1){
    stop("Check inputs: dimensions do not match") }
  
  n = n1
  
  x = MASS::mvrnorm(n=nsims, mu=rep(0,n), Sigma=Sigma)
  px <- pnorm(x)
  
  y <- matrix(0,nrow=nsims,ncol=n)
  for(j in 1:n)
  {
    y[,j] <- trapezoid::qtrapezoid(px[,j],min=x_min[j],mode1=x_mode1[j],
                        mode2=x_mode2[j],max=x_max[j]) 
  }
  
  if(sum(mean_rev)>0){
    
    if(length(mean_rev) != n){
      stop("Check inputs: length of mean_rev does not match length of other inputs") }
    
    if(mean_rev[n] == 1){
      warning("Last time period in mean_rev was set to 1: ignored!") }
    
    
    mean_rev_t = c(1:(n-1),0) * mean_rev
    mean_rev_t = mean_rev_t[mean_rev_t > 0]
    n_rev_t = length(mean_rev_t)
    
    ys = y
    
    for(j in 1:n)
    {
      ys[,j] = sort(y[,j],decreasing=T)
    }
    
    y_reorder = c()
    y_mr = c()
    y_cor = c()
    
    for(j in 1:n_rev_t)
    {
      
      j1 = mean_rev_t[j]
      j2 = j1 + 1
      
      ranks_ys = rank(-(abs(ys[,j1]-mean(ys[,j1]))))
      cutoffs = 1-ranks_ys/nsims
      mr = ifelse(runif(nsims)<cutoffs,1,0)
      
      ## create indices for mean reverting part
      mr_cor = correls[j1,j2]
      index_to_mr = mr * 1:nsims
      index_to_mr = index_to_mr[index_to_mr>0]
      targ_index_for_mr = seq(from=round(nsims/2-sum(mr)/2,0),by=1,length=sum(mr))
      
      x = MASS::mvrnorm(n=sum(mr), mu=c(0,0), 
                  Sigma=matrix(c(1,mr_cor,mr_cor,1),2,2))
      targ_cor_ranks = cbind(rank(-x[,1]),rank(-x[,2]))
      targ_cor_ranks = targ_cor_ranks[order(targ_cor_ranks[,1]),]
      targ_cor_index_for_mr = targ_index_for_mr[targ_cor_ranks[,2]]
      
      y_mr[[j]] = cbind(index_to_mr,targ_cor_index_for_mr)
      
      ## create indices for correlated part
      nmr_cor = correls[j1,j2]
      index_to_cor = (1:nsims)[-index_to_mr]
      targ_index_for_cor = (1:nsims)[-targ_index_for_mr]
      
      x = MASS::mvrnorm(n=(nsims-sum(mr)), mu=c(0,0), 
                  Sigma=matrix(c(1,nmr_cor,nmr_cor,1),2,2))
      targ_cor_ranks = cbind(rank(-x[,1]),rank(-x[,2]))
      targ_cor_ranks = targ_cor_ranks[order(targ_cor_ranks[,1]),]
      targ_cor_index_for_cor = targ_index_for_cor[targ_cor_ranks[,2]]
      
      y_cor[[j]] = cbind(index_to_cor,targ_cor_index_for_cor)
      
      y_reorder[[j]] = rbind(y_mr[[j]],y_cor[[j]])
      y_reorder[[j]] = y_reorder[[j]][order(y_reorder[[j]][,1]),]
      
    }
    
    y_mr = y
    for(j in 1:n_rev_t)
    {
      j1 = mean_rev_t[j]
      j2 = j1 + 1
      y_mr = y_mr[order(y_mr[,j]),]
      y_mr[,j2] = ys[y_reorder[[j]][,2],j2]
    }  
  }
  
  if(sum(mean_rev)>0){simdata = y_mr}else{simdata = y}

  return(simdata)
}