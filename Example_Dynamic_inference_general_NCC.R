#######################################################################################
#
#                                                                                     
#   Filename    :	      Example_Dynamic_inference_general_NCC.R    												  
#                                                                                     
#   Project     :       Biometrics article "Dynamic inference in general nested case-
#                       control designs"                                                             
#   Author      :       J. Feifel                                                                
#   Date        :       September 6, 2019
#   Purpose     :       Generate required functions for Figures 1 and 2
#																				  
#   R Version   :       R-3.6.1                                                                
#
#
#   Input source files  :  function_breslow.R, function_conf_regions.R, 
#   Input source files  :  data_SIR3_full.csv, data_SIR3_NCC.csv, data_SIR3_full_large.csv
#   Output data files   :  Web Figure 2.jpeg   Web figure 3.jpeg
#
#   Required R packages :  data.table  survival
#
#
########################################################################################

# Delete Work Space
rm(list = ls())

# packages required
library(data.table)
library(survival)

# include functions
source("function_breslow.R")
source("function_conf_regions.R")

# setting variables
no.boot <- 1000 # Number of realizations of hat W^star
alpha <- 0.05 # nominal level
s1 <- 9   # lower border of interval [s_1, s_2] (see Section 4)
s2 <- 60  # upper border of interval [s_1, s_2] (see Section 4)

# Read in data files
dat_full <- read.csv('data_SIR3_full.csv')
dat_ncc <- read.csv('data_SIR3_NCC.csv')
dat_full_large <- read.csv('data_SIR3_full_large.csv')

# data-preparation for breslow
dt <- data.table(dat_full)
dt[, x.cov := ifelse(sex=="M", 1,0)] # first covariate needs to be named x.cov
dt[, x.cov2 := age] # second covariate needs to be named x.cov2

ncc <- data.table(dat_ncc)
ncc[, x.cov := ifelse(sex=="M", 1,0)] # first covariate needs to be named x.cov
ncc[, x.cov2 := age] # second covariate needs to be named x.cov2

cox <- data.table(dat_full_large)
cox[, x.cov := ifelse(sex=="M", 1,0)] # first covariate needs to be named x.cov
cox[, x.cov2 := age] # second covariate needs to be named x.cov2

# Note: The original sample size (n in the manuscript) is required
# This is not necessarily computable from nested case-control data (NCC)
n <- dt[unique(id),.N] # calculate the origina cohort sample size

#############################
# NCC 
#############################
# coefficent estimation (simple random sampling weights cancel out but for general NCCs required)
# this is not required for the functions breslow or conf_regions
est.ncc <- coxph(Surv(Time, Fail)~ x.cov + x.cov2 + strata(Set) + offset(log(wt)), data = ncc)

bres.ncc <- breslow(ncc, n, no.boot=no.boot)

out.ncc <- conf_regions(bres.ncc, n = n, s1 = s1, s2 = s2)

###########################
## Cox 
###########################
# coefficent estimation (weights for simple random sampling not required in this step)
# Note: This command is not required for the later application of the functions breslow or conf_regions
## Cox using original data set
est.full <- coxph(Surv(Time, Fail)~ x.cov + x.cov2 , data = dt)
## Cox analysed in nested case-control style on large data set
est.cox <- coxph(Surv(Time, Fail)~ x.cov + x.cov2 + strata(Set) , data = cox)
# should be appriximately the same (simple check not required)
est.cox
est.full

bres.cox <- breslow(cox, n)

out.cox <- conf_regions(bres.cox, n = n, s1 = s1, s2 = s2)

###########################
## Plot 
###########################
jpeg(file = paste("Approximate Figure 2", "jpeg", sep = "."), height = 2500, width = 2500,  res = 400)
# cumulative hazard function
plot(out.ncc$times, out.ncc$cumHaz, type = 's', col='grey', lty=1, 
    ylab = 'Cumulative hazard for pneumonia', xlab = 'Time in intensive care unit (days)', xlim = c(10,60))
lines(out.cox$times, out.cox$cumHaz, type = 's', col='black', lty=1)
# log-transformed confidence intervals for classical Cox model
lines(out.cox$times, out.cox$CIp.log, type = 's', col='black', lty=3)
lines(out.cox$times, out.cox$CIm.log, type = 's', col='black', lty=3)
# Equal precision confidence bands for classical Cox model
lines(out.cox$times, out.cox$CBp.EP, type = 's', col='black', lty=2)
lines(out.cox$times, out.cox$CBm.EP, type = 's', col='black', lty=2)
# log-transformed confidence intervals for NCC data
lines(out.ncc$times, out.ncc$CIp.log, type = 's', col='grey', lty=3)
lines(out.ncc$times, out.ncc$CIm.log, type = 's', col='grey', lty=3)
# Equal precision confidence bands for NCC data
lines(out.ncc$times, out.ncc$CBp.EP, type = 's', col='grey', lty=2)
lines(out.ncc$times, out.ncc$CBm.EP, type = 's', col='grey', lty=2)

dev.off()