breslow <- function(ncc, n, no.boot=1000){
  ## ----------------------------------------------------------------------------
  ## Title: function_breslow.R
  ## ----------------------------------------------------------------------------
  ## Author: Jan Feifel
  ## ----------------------------------------------------------------------------
  ## Description: 
  ##       Derives the function of Equation (9) and its empirical variance
  ##              
  ## ----------------------------------------------------------------------------
  ## Required Packages: 
  ## ----------------------------------------------------------------------------
  ## Usage: 
  ##
  ##      data1: data file with structure similar to provided files
  ##      only data files with 2 covariates supported
  ##      
  ## ----------------------------------------------------------------------------
  ## Value: output (of type 'list'): 
  ##          times = unique event times 
  ##          cumHaz = estimator of cumulative hazard function
  ##          sigma.hat = variance estimator for (hat A - A) without star, (see 
  ##                      Eq. (4.6) Borgan et al. 1995) 
  ##          Boot.res = Eq (9) for each time
  ##          emp.var.boot = empirical variance of Eq (9) for each time
  ## ----------------------------------------------------------------------------
  
  # NCC-estimation (simple random sampling cancel out but for general)
  est.ncc <- coxph(Surv(Time, Fail)~ x.cov + x.cov2 + strata(Set) + offset(log(wt)), data = ncc)
  
  # Estimation helper function S_r^(gamma)=s_gamma
  s_0.ncc <- ncc[,(sum(exp(est.ncc$coefficients[1]* x.cov + est.ncc$coefficients[2]* x.cov2)* wt)), by=Set]$V1
  
  s_1.ncc <- cbind(ncc[,(sum(x.cov *exp(est.ncc$coefficients[1]* x.cov + est.ncc$coefficients[2]* x.cov2)* wt)), by=Set]$V1,
                   ncc[,(sum(x.cov2*exp(est.ncc$coefficients[1]* x.cov + est.ncc$coefficients[2]* x.cov2)* wt)), by=Set]$V1)
  
  ncc.S2 <- ncc[,multip:=exp(est.ncc$coefficients[1]*x.cov + est.ncc$coefficients[2]*x.cov2)* wt, by=.(Set,id)]
  
  ncc.S2 <- ncc.S2[,.( "ge11"=sum(x.cov^2* multip), 
                       "ge12"=sum(x.cov*x.cov2* multip), 
                       "ge21"=sum(x.cov*x.cov2* multip),  
                       "ge22"=sum(x.cov2^2* multip)),by=Set]
 # browser()
  s_2.ncc <- lapply(ncc.S2[unique(Set),  .N], function(i){
    d <- matrix(unlist(ncc.S2[Set==i, .(ge11, ge12, ge21, ge22)]), nrow = 2)
    colnames(d) <- c("x.cov", "x.cov2")
    return(d) })
  
  ncc.times <- unique(ncc$Time) # event times within ncc data set
  
  # NCC Breslow type estimator (see Section 2)
  A.ncc <- cumsum(1/s_0.ncc) 
  
  # components of variance estimator for (hat A - A) without star, (see Borgan et al. 1995)
  omega.NCC <- cumsum(1/s_0.ncc^2) # hat omega^2 from Theorem 1 
  B.hat <- apply(s_1.ncc/s_0.ncc^2,2,cumsum)
  E.ncc <- s_1.ncc/s_0.ncc
  
  n.sets <- ncc[unique(Set),.N] # Number of sampled sets
  
  colnames(E.ncc) <- paste("E",colnames(ncc[Fail==1,.(x.cov, x.cov2)]), sep=".")
  difZiMinE <- cbind(ncc[Fail==1,.(x.cov, x.cov2)], E.ncc)
  difZiMinE <- difZiMinE[, .(diff1=x.cov- E.x.cov, diff2=x.cov2-E.x.cov2)]
  
  V.ncc.helper <- data.table(E.ncc, Set=1:n.sets, "s0"=s_0.ncc)
  V.ncc.helper <- V.ncc.helper[,c("Ee11", "Ee12", "Ee21","Ee22"):=list(E.x.cov^2, 
                                                                        E.x.cov*E.x.cov2, 
                                                                        E.x.cov*E.x.cov2,  
                                                                        E.x.cov2^2)]
  V.ncc.helper <- merge(ncc.S2, V.ncc.helper, by = "Set")
  
  V.ncc.fast <- V.ncc.helper[, .( "Ve11"=ge11/s0-Ee11, "Ve12"=ge12/s0-Ee12, 
                                  "Ve21"=ge21/s0-Ee21, "Ve22"=ge22/s0-Ee22)]
  
  I.betahat <- matrix(colSums(V.ncc.fast),2,2) # observed information matrix evaluated at hat beta
  var.NCC <- sapply(1:dim(B.hat)[1], function(i){t(B.hat[i,]) %*% solve(I.betahat) %*% B.hat[i,]}) + 
                  omega.NCC # variance estimator for (hat A - A) without star, (see Eq. (4.6) Borgan et al. 1995) 
  
  #### Wild Bootstap with normal multiplier
  BootstrCB <- function(iter){
   
    Gi <- rnorm(n.sets,0,1) # normal distributed multiplier
    
    U_tau.star <- unlist(difZiMinE[,.(x.cov=sum(Gi*diff1), x.cov2=sum(Gi*diff2))]) # Eq (5) per bootstrap iteration and for hat beta 
    
    I_tau.star <- matrix(unlist(V.ncc.fast[,.("Ve11"=sum(Gi^2*Ve11),   
                                              "Ve12"=sum(Gi^2*Ve12), 
                                              "Ve21"=sum(Gi^2*Ve21),  
                                              "Ve22"=sum(Gi^2*Ve22))]), nrow = 2) # Eq (5) per bootstrap iteration and for hat beta 
    
    sqrtNbeta.starMinbeta.hat <- t(n*1/sqrt(n)*solve(I_tau.star)%*%U_tau.star) # Eq (8) per bootstrap iteration
    
    return(cumsum(Gi/s_0.ncc)*sqrt(n) - sqrtNbeta.starMinbeta.hat %*% t(B.hat) ) # Eq (9) per bootstrap iteration
  }
  
  Boot.res <- sapply(1:no.boot, FUN = BootstrCB)/sqrt(n) # Eq (9) for each time
  
  emp.var.boot <- apply(Boot.res, FUN =  var, MARGIN = 1) # empirical variance of (9) for each time
  
  return(list(times = ncc.times, 
              cumHaz = A.ncc,
              sigma.hat = var.NCC,
              Boot.res = Boot.res,
              emp.var.boot = emp.var.boot)
              )
}

