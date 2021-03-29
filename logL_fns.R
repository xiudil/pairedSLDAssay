### prepare the likelihood functions
### negative log-likelihood one sample
likelihood = function(pos, replicates, dilutions, tau){
  loglike = 0
  for (i in 1:length(pos)){
    if ((1-exp(-tau*dilutions[i]/10^6))!=0){
      loglike = loglike-(pos[i]*log(1-exp(-tau*dilutions[i]/10^6))
                         +(replicates[i]-pos[i])*(-tau*dilutions[i]/10^6))
    }
    else {
      loglike = loglike-(pos[i]*(-exp(-tau*dilutions[i]/10^6))
                         +(replicates[i]-pos[i])*(-tau*dilutions[i]/10^6))
    }
  }
  
  return(loglike)
}

likelihoodwrap = function(pos, replicates, dilutions, logtau){
  return(likelihood(pos, replicates, dilutions, exp(logtau)))
}

### negative log-likelihood two sample (tau2 and taudiff)
likenull = function(pos1,replicates1,dilutions1,pos2,replicates2,dilutions2,tau2,taudiff){
  loglikenull = likelihood(pos1,replicates1,dilutions1,tau2+taudiff)+
    likelihood(pos2,replicates2,dilutions2,tau2)
  return(as.numeric(loglikenull))
}

### observed and expected Fisher information
observedInfo = function(postau,dilutions,replicates){
  pos=postau[1:length(replicates)]
  tau=postau[length(replicates)+1]
  if (sum(pos)>0){
    info = sum(pos*dilutions^2/10^12*exp(-dilutions*tau/10^6)/(1-exp(-dilutions*tau/10^6))^2)
  }
  else{
    info = sum(replicates*dilutions^2/10^12*exp(-dilutions*tau/10^6)/(1-exp(-dilutions*tau/10^6)))
  }
  return(info)
}

expectedInfo = function(postau,dilutions,replicates){
  pos=postau[1:length(replicates)]
  tau=postau[length(replicates)+1]
  info = sum(replicates*dilutions^2/10^12*exp(-dilutions*tau/10^6)/(1-exp(-dilutions*tau/10^6)))
  return(info)
}