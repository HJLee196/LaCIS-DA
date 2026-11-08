ebh <- function(E, alpha=0.1){
  
  p <- length(E)
  E_ord <- order(E, decreasing = TRUE)
  E <- sort(E, decreasing = TRUE)
  comp <- E >= (p / alpha / (1:p))
  id <- suppressWarnings(max(which(comp>0)))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }
  return(list(rej = rej, cut = E_ord[id]))
}

bh <- function(E, alpha=0.1){
  
  p <- length(E)
  E_ord <- order(E, decreasing = FALSE)
  E <- sort(E, decreasing = FALSE)
  comp <- E <= (alpha * (1:p) / p)
  id <- suppressWarnings(max(which(comp>0)))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }
  return(list(rej = rej, cut = E_ord[id]))
}

multiknk = function(W_imp, alpha=0.1){
  m = dim(W_imp)[1] - 1
  
  W_max_ind = apply(W_imp,2,which.max)
  W_max = sapply(1:length(W_max_ind), function(x) W_imp[W_max_ind[x],x])
  W_med = sapply(1:length(W_max_ind), function(x) median(W_imp[-W_max_ind[x],x]))
  
  tau_stab <- (W_max - W_med)
  kap = (W_max_ind==1)
  
  ts = sort(c(0, abs(tau_stab)))
  ratio = sapply(ts, function(t) (1/m + 1/m*sum((tau_stab >= t)*(1-kap)))/max(1, 
                                                                              sum((tau_stab >= t)*(kap))))
  ok = which(ratio <= alpha)
  threshold <- ifelse(length(which(ts[ok]>0)) > 0, ts[ok][which(ts[ok]>0)[1]], Inf)
  
  selected <- which(tau_stab*kap >= threshold)
  return(list(rej = selected, cut = threshold))
}

## q value for BH
qbh <- function(p_val){
  p <- length(p_val)
  p_val_ord <- sort(p_val, decreasing = FALSE)
  q_val_ord <- p_val_ord*p/(1:p)
  
  q_val = rep(1,p)
  for(i in p:1){
    qi = min(q_val_ord[i],1)
    q_val[bh(p_val,qi)$rej] = qi
  }
  return(q_val)
}

# qvalue for ebh
qebh = function(W_imp, N = 100){
  p <- dim(W_imp)[2]
  q_val = rep(1,p)
  # W_imp = FinalResult[[i]][[2]]
  for(q in (seq(1,0,length.out = (N+1))[-(N+1)])){
    m = 10
    E <- matrix(NA, m, dim(W_imp)[2])
    for(w_it in 1:m){
      E[w_it,] <- W_imp[1,] - W_imp[(w_it+1),]
      tau <- knockoff.threshold(E[w_it,],fdr = q/2)
      E[w_it,] <- (E[w_it,] >= tau) / (1 + sum(E[w_it,] <= -tau))
    }
    
    E <- (dim(W_imp)[2]) * colMeans(E)
    
    q_val[ebh(E,q)$rej] = q
  }
  return(q_val)
}


# qvalue for multiple knockoffs
qmulti = function(W_imp, m = 10){
  p <- dim(W_imp)[2]
  
  W_max_ind = apply(W_imp,2,which.max)
  W_max = sapply(1:length(W_max_ind), function(x) W_imp[W_max_ind[x],x])
  W_med = sapply(1:length(W_max_ind), function(x) median(W_imp[-W_max_ind[x],x]))
  
  tau_stab <- (W_max - W_med)
  kap = (W_max_ind==1)
  
  ts = sort(c(0, abs(tau_stab)))
  ratio = sapply(ts, function(t) (1/m + 1/m*sum((tau_stab >= t)*(1-kap)))/max(1, 
                                                                              sum((tau_stab >= t)*(kap))))
  q_val_ord = unique(ratio)
  q_val_ord = sort(q_val_ord, decreasing = TRUE)
  q_val = rep(1,p)
  for(i in q_val_ord){
    ok = which(ratio <= i)
    threshold <- ifelse(length(which(ts[ok]>0)) > 0, ts[ok][which(ts[ok]>0)[1]], Inf)
    
    qi = min(i,1)
    q_val[which(tau_stab*kap >= threshold)] = qi
  }

  return(q_val)
}
