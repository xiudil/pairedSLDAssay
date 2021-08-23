library(SLDAssay)
library(dplyr)
source("logL_fns.R")

#################
### screening ###
#################
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
  result = c(mle,0,maxprob)
  return(result)
}

### screening step: find plausible configurations
findconfig = function(replicates1,dilutions1,replicates2,dilutions2,same.setting = FALSE,thres=0.001){
  temp = list(rep(NA,length(replicates1)))
  for (i in 1:length(replicates1)){
    temp[[i]] = seq(0,replicates1[i])
    names(temp)[i] = paste("level",i,sep="")
  }
  config1 = data.frame(expand.grid(temp))
  result1 = apply(config1,1,wrap.get.mle,replicates=replicates1,dilutions=dilutions1)
  config1$MLE = result1[1,]
  config1$maxprob = result1[3,]
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
    config2$maxprob = result2[3,]
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


### find extreme configurations
### this one depends on taudiff
findExtreme = function(config,taudiff,pos1,pos2,testmethod,replicates1,replicates2,dilutions1,dilutions2){
  if (testmethod=="naive"){
    reps1 = length(pos1)
    reps2 = length(pos2)
    tau1hat = config[[1]][apply(t(config[[1]][,1:reps1])==pos1,2,prod)==1,reps1+1]
    tau2hat = config[[2]][apply(t(config[[2]][,1:reps2])==pos2,2,prod)==1,reps2+1]
    currentStat = tau1hat-tau2hat-taudiff
    
    allstat = config[[3]][,1]-config[[3]][,2]-taudiff
    extremeConfig = config[[3]][abs(allstat)>=abs(currentStat),c("index1","index2")]
  }
  else if (testmethod == "Wald"){
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
  }
  
  ### this is still slow, and tau2hatnull often fall below 0
  ### change it to log scale and then transform it back
  else if (testmethod=="LRT"){
    LR = matrix(0,dim(config[[3]])[1],3)
    LR[,1] = rep(seq(1:dim(config[[1]])[1]),dim(config[[2]])[1])
    LR[,2] = rep(seq(1:dim(config[[2]])[1]),each=dim(config[[1]])[1])
    
    system.time({for (l in 1:dim(config[[3]])[1]){
        #print(l)
        i = LR[l,1]
        j = LR[l,2]
        tau2hat = config[[2]][j,"MLE"]
        tau1hat = config[[1]][i,"MLE"]
        info2 = config[[2]][j,"info"]
        info1 = config[[1]][i,"info"]
        
        tau2hatnull = exp(log(tau2hat)-info1/(info1+info2)/tau2hat*(taudiff-(tau1hat-tau2hat)))
        lnull = likenull(pos1=config[[1]][i,1:length(replicates1)],replicates1,dilutions1,
                         pos2=config[[2]][j,1:length(replicates2)],replicates2,dilutions2,tau2hatnull,taudiff)
        lfull = likelihood(config[[1]][i,1:length(replicates1)], replicates1, dilutions1, tau1hat)+
                likelihood(config[[2]][j,1:length(replicates2)], replicates2, dilutions2, tau2hat)
        LR[l,3] = as.numeric(lnull-lfull)
    }})
    index1 = which(apply(t(config[[1]][,1:length(pos1)])==pos1,2,prod)==1)
    index2 = which(apply(t(config[[2]][,1:length(pos2)])==pos2,2,prod)==1)
    current = LR[which(LR[,1]==index1 & LR[,2]==index2),3]
    extremeConfig = LR[LR[,3]>=current,c(1,2)]
  }
  return(extremeConfig)
}

### probability under a particular null
calextremep = function(config,extremeConfig, taudiff, tau2,replicates1,replicates2,dilutions1,dilutions2){
  tau1=taudiff+tau2
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

### compute p-value using the M/BB/E approach
mpvalue = function(config,extremeConfig,taudiff,replicates1,replicates2,dilutions1,dilutions2,taul=0,tauu=0,pvalmethod,pos1,pos2){
  if (pvalmethod=="BB"){
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

### this is a wrapper for uniroot
onediffp = function(taudiff,config,pos1,pos2,replicates1,replicates2,dilutions1,dilutions2,pvalmethod,testmethod,thres=0.001){
  extremeConfig = findExtreme(config,taudiff, pos1=pos1, pos2=pos2,replicates1=replicates1,
                              replicates2=replicates2,dilutions1=dilutions1,dilutions2=dilutions2,
                              testmethod = testmethod)
  CI1 = get.mle(pos=pos1,dilutions=dilutions1,replicates=replicates1,conf.level=0.9975)$Exact_CI
  CI2 = get.mle(pos=pos2,dilutions=dilutions2,replicates=replicates2,conf.level=0.9975)$Exact_CI
  taul = max(CI2[1],CI1[1]-taudiff)
  tauu = min(CI2[2],CI1[2]-taudiff)
  p = mpvalue(config,extremeConfig,taudiff,replicates1,replicates2,dilutions1,dilutions2,taul,tauu,pvalmethod,pos1,pos2)
  #return((p-0.05)^2)
  return(p-0.05-2*thres)
}


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
  try({CIlow = uniroot(onediffp,interval=c(-100000,diffMLE),config=config, pos1=pos1,pos2=pos2,replicates1=replicates1,replicates2=replicates2,
                       dilutions1=dilutions1,dilutions2=dilutions2,pvalmethod=pvalmethod,testmethod=testmethod)$root})
  try({CIup = uniroot(onediffp,interval=c(diffMLE,100000),config=config, pos1=pos1,pos2=pos2,replicates1=replicates1,replicates2=replicates2,
                      dilutions1=dilutions1,dilutions2=dilutions2,pvalmethod=pvalmethod,testmethod=testmethod)$root})
  result = vector("list",4)
  result[[1]] = diffMLE
  result[[2]] = c(CIlow,CIup)
  result[[3]] = c(aCIlow,aCIup)
  result[[4]] = c(pos1,pos2)
  names(result) = c("delta","exactCI","asympCI","rawAssayResult")
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

twosample_mle(pos1,pos2,replicates1,replicates2,dilutions1,dilutions2,"E","Wald")

