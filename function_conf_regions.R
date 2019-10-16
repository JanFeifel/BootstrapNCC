conf_regions <- function(bres, n, s1=0, s2=Inf, alpha = 0.05)
{
  ## ----------------------------------------------------------------------------
  ## Title: function_breslow.R
  ## ----------------------------------------------------------------------------
  ## Author: Jan Feifel
  ## ----------------------------------------------------------------------------
  ## Description: 
  ##       Generates depictable functions based on function_breslow incl. confi-
  ##       dence regions
  ##              
  ## ----------------------------------------------------------------------------
  ## Required Packages: 
  ## ----------------------------------------------------------------------------
  ## Usage: 
  ##
  ##      bres: output element of function_breslow.R
  ##      n: smaple size in entire cohort
  ##      s1: lower border of interval [s_1, s_2] (see Section 4)
  ##      s2: upper border of interval [s_1, s_2] (see Section 4)
  ##      alpha:   nominal level
  ##      
  ## ----------------------------------------------------------------------------
  ## Value: output (of type 'list'): 
  ##                    times = unique event times 
  ##                    cumHaz = estimator of cumulative hazard function
  ##                    CIp.lin = upper limit of linear transformed confidence interval
  ##                    CIm.lin = lower limit of linear transformed confidence interval
  ##                    CIp.log = upper limit of log-transformed confidence interval
  ##                    CIm.log = lower limit of log-transformed confidence interval
  ##                    CBp.lin = upper limit of linear transformed confidence band
  ##                    CBm.lin = upper limit of linear transformed confidence band
  ##                    CBp.EP = upper limit of log-transformed equal precision confidence band
  ##                    CBm.EP = lower limit of log-transformed equal precision confidence band
  ##                    CBp.HW = upper limit of log-transformed Hall-Wellner confidence band
  ##                    CBm.HW = lower limit of log-transformed Hall-Wellner confidence band
  ## ----------------------------------------------------------------------------
  
  # Restrict Eq (9) and its emprical variance to window [s1, s2]
  Boot.res.qu <- bres$Boot.res[bres$times<=s2 & bres$times>=s1,]
  emp.var.boot.qu <- bres$emp.var.boot[bres$times<=s2 & bres$times>=s1]
  
  
  ## calcualtion of c_{phi, 1-alpha}^{g_star}
  # when phi=phi1 and g_star=1 
  supreme.lin <- apply(Boot.res.qu, FUN =  function(x)max(abs(x)), MARGIN = 2)
  cb.lin <- quantile(supreme.lin, 1-alpha, na.rm = T)
  
  # when phi=phi2 and g_star = g1
  supreme.EP <- apply(1/sqrt(emp.var.boot.qu)*Boot.res.qu, FUN =  function(x)max(abs(x)), MARGIN = 2)
  cb.EP <- quantile(supreme.EP, 1-alpha, na.rm = T)
  
  # when phi=phi2 and g_star = g2
  supreme.HW <- apply(1/(1+n*emp.var.boot.qu)*Boot.res.qu, FUN =  function(x)max(abs(x)), MARGIN = 2)
  cb.HW <- quantile(x = supreme.HW, probs = 1-alpha, na.rm = T)
  
  # CB*** confidence band follows the equation in Section 4 with the respective functions
  # Ci** confidence interval relying on qunatile of normal distribution
  return(list(times =bres$times, 
              cumHaz =bres$cumHaz,
              CIp.lin = bres$cumHaz + qnorm(1-alpha/2) * sqrt(bres$sigma.hat),
              CIm.lin = bres$cumHaz - qnorm(1-alpha/2) * sqrt(bres$sigma.hat),
              CIp.log = bres$cumHaz * exp(+ qnorm(1-alpha/2)/(bres$cumHaz) * sqrt(bres$sigma.hat)),
              CIm.log = bres$cumHaz * exp(- qnorm(1-alpha/2)/(bres$cumHaz) * sqrt(bres$sigma.hat)),
              CBp.lin = bres$cumHaz + cb.lin,
              CBm.lin = bres$cumHaz - cb.lin,
              CBp.EP = bres$cumHaz * exp(+ cb.EP/(bres$cumHaz) * sqrt(bres$emp.var.boot)),
              CBm.EP = bres$cumHaz * exp(- cb.EP/(bres$cumHaz) * sqrt(bres$emp.var.boot)),
              CBp.HW = bres$cumHaz * exp(+ cb.HW/(bres$cumHaz) * (1+n*bres$emp.var.boot)), 
              CBm.HW = bres$cumHaz * exp(- cb.HW/(bres$cumHaz) * (1+n*bres$emp.var.boot))
  ))
}


