library(SLDAssay)
library(dplyr)
source("logL_fns.R")

#################
### screening ###
#################

###################################
### this function computes
###          1. mle of IUPM hat tau
###          2. maximum probability of an assay outcome over all values of tau
### for one sample
###################################
wrap.get.mle=function(pos,replicates,dilutions,monte=15000,conf.level=0.95,iupm=TRUE,na.rm=FALSE){
  if (sum(pos)==0){
    mle = 0
  } else if (sum(pos) == sum(replicates)){
    mle = 10^6
  } else {
    mle = exp(optim(par=10,likelihoodwrap,method="L-BFGS-B",pos=pos,replicates = replicates, dilutions = dilutions)$par)
  }
  
  maxprob = exp(-likelihood(pos,replicates,dilutions,mle))
  for (i in 1:length(replicates)){
    maxprob = maxprob * choose(replicates[i],pos[i])
  }
  result = c(mle,maxprob)
  return(result)
}
##################################
### This function performs the screening step, it returns
###        1. mle and fisher information for plausible outcomes for each sample
###        2. plausible combinations of outcomes from assay 1 and 2
##################################
findconfig = function(replicates1,dilutions1,replicates2,dilutions2,same.setting = FALSE,thres=0.001){
  temp = list(rep(NA,length(replicates1)))
  for (i in 1:length(replicates1)){
    temp[[i]] = seq(0,replicates1[i])
    names(temp)[i] = paste("level",i,sep="")
  }
  config1 = data.frame(expand.grid(temp))
  result1 = apply(config1,1,wrap.get.mle,replicates=replicates1,dilutions=dilutions1)
  config1$MLE = result1[1,]
  config1$maxprob = result1[2,]
  config1 = config1 %>% arrange(maxprob)
  index1 = sum(cumsum(config1$maxprob)<=thres)
  config1 = config1[-c(1:index1),]
  obsINFO1 = apply(config1,1,observedInfo,dilutions=dilutions1,replicates=replicates1)
  expINFO1 = apply(config1,1,expectedInfo,dilutions=dilutions1,replicates=replicates1)
  config1$info = expINFO1
  config1$obsinfo = obsINFO1
  
  if(same.setting == TRUE){
    config2 = config1
  }
  else{
    temp = list(rep(NA,length(replicates2)))
    for (i in 1:length(replicates2)){
      temp[[i]] = seq(0,replicates2[i])
      names(temp)[i] = paste("level",i,sep="")
    }
    config2 = data.frame(expand.grid(temp))
    result2 = apply(config2,1,wrap.get.mle,replicates=replicates2,dilutions=dilutions2)
    config2$MLE = result2[1,]
    config2$maxprob = result2[2,]
    config2 = config2 %>% arrange(maxprob)
    index2 = sum(cumsum(config2$maxprob)<=0.001)
    config2 = config2[-c(1:index2),]
    obsINFO2 = apply(config2,1,observedInfo,dilutions=dilutions2,replicates=replicates2)
    expINFO2 = apply(config2,1,expectedInfo,dilutions=dilutions2,replicates=replicates2)
    config2$info = expINFO2
    config2$obsinfo = obsINFO2
  }
  
  config = list()
  config[[1]] = config1
  config[[2]] = config2
  
  MLEgrid = expand.grid(config[[1]]$MLE,config[[2]]$MLE)
  MLEgrid$index1 = rep(seq(1:dim(config[[1]])[1]),dim(config[[2]])[1])
  MLEgrid$index2 = rep(seq(1:dim(config[[2]])[1]),each=dim(config[[1]])[1])
  config[[3]] = MLEgrid
  
  return(config)
}

##########################################
### This function finds extreme outcomes defined as
### those leading to a test statistic with a larger magnitude than the observed
### **Currently, only Wald test statistic is allowed**
##########################################
findExtreme = function(config,taudiff,pos1,pos2,testmethod,replicates1,replicates2,dilutions1,dilutions2){
  if (testmethod == "Wald"){
    reps1 = length(pos1)
    reps2 = length(pos2)
    
    index1 = apply(t(config[[1]][,1:reps1])==pos1,2,prod)==1
    tau1hat = config[[1]][index1,"MLE"]
    information1 = config[[1]][index1,"info"]
    
    index2 = apply(t(config[[2]][,1:reps2])==pos2,2,prod)==1
    tau2hat = config[[2]][index2,"MLE"]
    information2 = config[[2]][index2,"info"]
    
    currentStat = (tau1hat - tau2hat - taudiff)/sqrt(1/information1 + 1/information2)
    
    allinfo1 = config[[1]][config[[3]]$index1,"info"]
    allinfo2 = config[[2]][config[[3]]$index2,"info"]
    allstat = (config[[3]][,1]-config[[3]][,2]-taudiff)/sqrt(1/allinfo1+1/allinfo2)
    
    extremeConfig = config[[3]][abs(allstat)>=abs(currentStat),c("index1","index2")]
  } else{
    stop("Only Wald statistic is supported currently! Use testmethod == 'Wald'")
  }
  return(extremeConfig)
}

#######################################
### This function computes the total probability 
###      of extreme outcomes under given values of Delta and tau2
#######################################
calextremep = function(config,extremeConfig, taudiff, tau2,replicates1,replicates2,dilutions1,dilutions2){
  tau1 = taudiff+tau2
  pos1s = config[[1]][extremeConfig[,1],1:length(replicates1)]
  pos2s = config[[2]][extremeConfig[,2],1:length(replicates2)]
  prob = 1
  for (i in 1:length(replicates1)){
    poismean = tau1*dilutions1[i]/10^6
    prob = prob*dbinom(pos1s[,i],replicates1[i],1-exp(-poismean))
  }
  for (i in 1:length(replicates2)){
    poismean = tau2*dilutions2[i]/10^6
    prob = prob*dbinom(pos2s[,i],replicates2[i],1-exp(-poismean))
  }
  return(-sum(prob,na.rm=TRUE))
}

#####################################
### This function computes p-value using the M/BB/E approach
###      for testing a specific value of Delta
#####################################
mpvalue = function(config,extremeConfig,taudiff,replicates1,replicates2,dilutions1,dilutions2,taul=0,tauu=0,pvalmethod,pos1,pos2){
  if (pvalmethod=="BB"){
    ### taul and tauu are the lower and upper limit of a 99.5% CI for tau2
    pvalue = (-1)*(optim(par=10,fn=calextremep,config=config,extremeConfig=extremeConfig,taudiff=taudiff,
                    replicates1=replicates1,replicates2=replicates2,dilutions1=dilutions1,dilutions2=dilutions2,
                    method="L-BFGS-B",lower=taul,upper=tauu)$value)
    pvalue = pvalue + 0.005
  }
  else if (pvalmethod=="E"){
    CI2 = get.mle(pos2, replicates2, dilutions2, conf.level=0.99)$Exact_CI
    tau2hatnull = optim(par=10,fn=likenull,pos1 = pos1,replicates1 = replicates1,dilutions1 = dilutions1,
                 pos2 = pos2,replicates2 = replicates2, dilutions2 = dilutions2,taudiff = taudiff,
                 method="Brent",lower=CI2[1],upper=CI2[2])$par
    pvalue = calextremep(config,extremeConfig, taudiff, tau2hatnull,replicates1,replicates2,dilutions1,dilutions2)
    pvalue = pvalue*(-1)
  }
  else if (pvalmethod== "M"){
    pvalue = -optim(par=10,fn=calextremep,config=config,extremeConfig=extremeConfig,taudiff=taudiff,
                    replicates1=replicates1,replicates2=replicates2,dilutions1=dilutions1,dilutions2=dilutions2,
                    method="L-BFGS-B",lower=max(0,(-1)*taudiff),upper=min(10^6,10^6-taudiff))$value
    pvalue = pvalue
  }
  return(pvalue)
}

#######################################
### This function is a wrapper function to be passed into uniroot,
###      which will be used to find the exact CI
#######################################
onediffp = function(taudiff,config,pos1,pos2,replicates1,replicates2,dilutions1,dilutions2,pvalmethod,testmethod,thres=0.001){
  extremeConfig = findExtreme(config,taudiff, pos1=pos1, pos2=pos2,replicates1=replicates1,
                              replicates2=replicates2,dilutions1=dilutions1,dilutions2=dilutions2,
                              testmethod = testmethod)
  CI1 = get.mle(pos=pos1,dilutions=dilutions1,replicates=replicates1,conf.level=0.9975)$Exact_CI
  CI2 = get.mle(pos=pos2,dilutions=dilutions2,replicates=replicates2,conf.level=0.9975)$Exact_CI
  taul = max(CI2[1],CI1[1]-taudiff)
  tauu = min(CI2[2],CI1[2]-taudiff)
  p = mpvalue(config,extremeConfig,taudiff,replicates1,replicates2,dilutions1,dilutions2,taul,tauu,pvalmethod,pos1,pos2)
  return(p-0.05-2*thres)
}

#################################################################
### main function that computes MLE, exact and asymptotic CIs ###
#################################################################
twosample_mle = function(pos1,pos2,replicates1,replicates2,dilutions1,dilutions2,pvalmethod,testmethod){
  config = findconfig(replicates1,dilutions1,replicates2,dilutions2,same.setting=FALSE)
  tau1hat = wrap.get.mle(pos=pos1,replicates=replicates1,dilutions=dilutions1)[1] 
  tau2hat = wrap.get.mle(pos=pos2,replicates=replicates2,dilutions=dilutions2)[1]
  diffMLE = tau1hat - tau2hat
  info1 = expectedInfo(c(pos1,tau1hat),dilutions1,replicates1)
  info2 = expectedInfo(c(pos2,tau2hat),dilutions2,replicates2)
  se = sqrt(1/info1 + 1/info2)
  
  aCIlow = diffMLE - qnorm(0.975)*se
  aCIup = diffMLE + qnorm(0.975)*se
  CIlow = -10^6
  CIup = 10^6
  suppressWarnings({
    try({CIlow = uniroot(onediffp,interval=c(-100000,diffMLE),config=config, pos1=pos1,pos2=pos2,replicates1=replicates1,replicates2=replicates2,
                         dilutions1=dilutions1,dilutions2=dilutions2,pvalmethod=pvalmethod,testmethod=testmethod)$root})
    try({CIup = uniroot(onediffp,interval=c(diffMLE,100000),config=config, pos1=pos1,pos2=pos2,replicates1=replicates1,replicates2=replicates2,
                        dilutions1=dilutions1,dilutions2=dilutions2,pvalmethod=pvalmethod,testmethod=testmethod)$root})
  })
  result = vector("list",3)
  result[[1]] = diffMLE
  result[[2]] = c(CIlow,CIup)
  result[[3]] = c(aCIlow,aCIup)
  names(result) = c("delta","exactCI","asympCI")
  return(result)
}


########################
### a simple example ###
########################
replicates1 = rep(2,6)
replicates2 = rep(2,6)
dilutions1 = c(1e6,2e5,4e4,8e3,1600,320)
dilutions2 = dilutions1
pos1 = c(2,2,1,1,0,0)
pos2 = c(1,1,0,0,0,0)

### only Wald statistic is supported currently
twosample_mle(pos1,pos2,replicates1,replicates2,dilutions1,dilutions2,
              pvalmethod = "E",testmethod = "Wald")

