
#############################################
#  Calculate Bayes factor, AIC, BIC and DIC #
#############################################

library(MCMCpack)
library(MASS)
library(geoR)
library(mnormt)

#Loglik of inv-chisq distribution. If we use the default function "dinvchisq" in R, the outcome might be -Inf
log_dinvchisq=function(x,df,scale){
  return(df/2*log(df/2)-lgamma(df/2)+df/2*log(scale)-((df/2)+1)*log(x)-(df*scale)/(2*x))
}

#Given the regression parameters of each state, calculate the estimates of bike sharing demands for all time points
cal_mu=function(beta_matrix,x_matrix,num_state,num_period,num_covariate){
  mu_matrix=matrix(NA,num_period,num_state)
  for(i in 1:num_state){
    betatemp=beta_matrix[,i]
    mu_matrix[,i]=x_matrix%*%matrix(betatemp,num_covariate,1)
  }
  return(mu_matrix)
}

#Given the linear regression prediction of y in all states (which is defined as mu), calculate the likelihood of y_t
y_Ytheta=function(y_sequence,num_covariate,mu_matrix,num_state,num_period,sigma2_array){
  liklihood=matrix(NA,num_period,num_state)
  for(i in 1:num_period){
    ytemp=y_sequence[i]
    for(j in 1:num_state){
      liklihood[i,j]=dnorm(ytemp,mu_matrix[i,j],sqrt(sigma2_array[j]),log=T)
    }
    liklihood[i,]=liklihood[i,]-max(liklihood[i,])
    liklihood[i,]=exp(liklihood[i,])
  }
  return(liklihood)#Note that it is not the true log-likelihood. We divide a constant from them for calculating some probabilities via Bayes rule
}

#Given the linear regression prediction of y in all states (which is defined as mu), calculate the log-likelihood of y_t
logy_Ytheta2=function(y_sequence,num_covariate,mu_matrix,num_state,num_period,sigma2_array){
  liklihood=matrix(NA,num_period,num_state)
  for(i in 1:num_period){
    ytemp=y_sequence[i]
    for(j in 1:num_state){
      liklihood[i,j]=dnorm(ytemp,mu_matrix[i,j],sqrt(sigma2_array[j]),log=T)
    }
  }
  return(liklihood)
}

#Recursively calculate P(s_t|Y_t,P,beta,sigma2)
s_YthetaP=function(p_array,y_liklihood,
                   num_period,num_state,num_covariate){
  #p_array is the state self-transition probability.
  p_array2=p_array[-length(p_array)]
  s_store=matrix(0,num_period,num_state)
  s_store[1,1]=1
  for(i in 1:(num_period-1)){
    #Calculate P(s_t|Y_(t-1))
    temp1=s_store[i,]
    s_temp1=temp1*p_array
    temp1=temp1[-length(temp1)]
    s_temp2=c(0,temp1*(1-p_array2))
    s_temp=s_temp1+s_temp2#s_t|Y_(t-1)
    #Calculate P(s_t|Y_t)
    y_temp=y_liklihood[i+1,]
    sy_temp=s_temp*y_temp
    if(sum(sy_temp)!=0){
      sy_temp=exp(log(sy_temp)-log(sum(sy_temp)))
    }else{#If sum(sy_temp)=0, the normalization term is 0. Hence, we define P(s_t|Y_t) as 1/num_state.
      #Reference: X. Pang, B. Friedman, A. D. Martin, and K. M. Quinn, "Endogenous jurisprudential regimes," Political Analysis, vol. 20, no. 4, pp. 417-436, 2012
      sy_temp=rep(1/num_state,num_state)
    }
    s_store[i+1,]=sy_temp
  }
  return(s_store)
}

#Calculate P(s_t|Y_n, S^(t+1), theta, P), and sample S
simu_S=function(p_array,y_liklihood,num_period,num_state,num_covariate){
  s_sample=c(num_state)
  s_left=s_YthetaP(p_array,y_liklihood,num_period,num_state,num_covariate)
  #Backward sampling
  for(i in (num_period-1):1){
    #The reachable states
    s2=s_sample[1]
    if(s2==1){
      s_sample=c(1,s_sample)
      next
    }
    s1=s2-1#s_t=s_(t+1)-1
    #Right term of the Equ. (5) in reference: S. Chib, "Estimation and comparison of multiple change-point models," Journal of Econometrics, vol. 86, no. 2, pp. 221-241, 1998.
    #Note that it is possible that p equals to 0, which means there is only one sample in that state and the state is redundant
    p1=p_array[s1]
    p2=1-p1
    #Left term of the Equ. (5)
    s_left1=s_left[i,s1]
    s_left2=s_left[i,s2]
    #Calculate P(s_t|Y_n, S^(t+1), theta, P)
    prob=c(s_left1*p1,s_left2*p2)
    if((prob[1]==0&&prob[2]==0)){#When s_left1=0 and s_left1=0, randomly choose the next state
      r=runif(1)
      u=runif(1)
      if(r>u){
        s_sample=c(s1,s_sample)
      }else{
        s_sample=c(s2,s_sample)
      }
    }else{#Sample s
      prob=exp(log(prob)-log(sum(prob)))
      r=runif(1)
      if(r<prob[1]){
        s_sample=c(s1,s_sample)
      }else{
        s_sample=c(s2,s_sample)
      }
    }
  }
  return(s_sample)
}

#Given state sequence S, sample the state transition probability P
simu_P=function(s_sample,num_state,num_period,a0,b0){
  #The prior distribution of p is Beta(a0,b0)
  p_array=c()
  if(num_state>1){
    for(i in 1:(num_state-1)){
      nii=length(s_sample[which(s_sample==i)])-1
      x=rbeta(1,a0+nii,b0+1)
      p_array=c(p_array,x)
    }
  }
  p_array=c(p_array,1)
  return(p_array)
}

#Calculate the posterior log-likelihood of beta and sigma^2
calcu_pi_loglik=function(s_sample,y_sequence,x_matrix,num_state,num_period,num_covariate,beta_matrix,sigma2_array,
                         k0_prior,v0_prior,sigma02_prior){
  pi_loglik=0
  for(i in 1:num_state){
    if(i%in%s_sample){
      y_sequencetemp=y_sequence[which(s_sample==i)]
      y_sequencetemp=matrix(y_sequencetemp,length(y_sequencetemp),1)
      x_matrixtemp=x_matrix[which(s_sample==i),]
      if(!is.matrix(x_matrixtemp)){#When there is only one observation in state i
        x_matrixtemp=matrix(x_matrixtemp,1,num_covariate)
      }
      temp=t(x_matrixtemp)%*%x_matrixtemp+k0_prior*diag(1,num_covariate)
      temp=solve(temp)
      #The matrix should be symmetric and positive definite
      init=21
      signal=T
      while(signal){
        init=init-1
        temp=round(temp,init)#To ensure that the matrix would be symmetric
        if(!F%in%as.vector(t(temp)==temp)){
          signal=F
          while(min(eigen(temp)$values)<0){#When the matrix is not positive definite, we add a small value to its diagonal elements
            temp=temp+diag(10^(-init),num_covariate)
          }
        }
      }
      betatemp=temp%*%t(x_matrixtemp)%*%y_sequencetemp
      n=length(s_sample[which(s_sample==i)])
      vn=v0_prior+n
      mutemp=x_matrixtemp%*%matrix(betatemp,num_covariate,1)
      residual=y_sequencetemp-mutemp
      #When the demands are zero (missing data or the stations in the region are unavailable), the residuals might be zeros and the log-likelihood of sigma^2 would be -Inf 
      sn=t(residual)%*%residual+k0_prior*t(betatemp)%*%betatemp+v0_prior*sigma02_prior
      sigma2temp=temp*sigma2_array[i]
      #The matrix should be symmetric and positive definite
      init=21
      signal=T
      while(signal){
        init=init-1
        sigma2temp=round(sigma2temp,init)#To ensure that the matrix would be symmetric
        if(!F%in%as.vector(t(sigma2temp)==sigma2temp)){
          signal=F
          while(min(eigen(sigma2temp)$values)<0){#When the matrix is not positive definite, we add a small value to its diagonal elements
            sigma2temp=sigma2temp+diag(10^(-init),num_covariate)
          }
        }
      }
      #Calculate the posterior log-likelihood of beta and sigma^2
      pi_loglik=pi_loglik+log_dinvchisq(sigma2_array[i],df=vn,scale=sn/vn)+dmnorm(beta_matrix[,i],c(betatemp),sigma2temp,log=T)
    }else{
      #When state i is unreachable, calculate the prior log-likelihood of beta and sigma^2
      pi_loglik=pi_loglik+log_dinvchisq(sigma2_array[i],df=v0_prior,scale=sigma02_prior)+dmnorm(beta_matrix[,i],rep(0,num_covariate),diag(sigma2_array[i]/k0_prior,num_covariate),log=T)
    }
  }
  return(pi_loglik)
}


######################################################################################################################################
# Calculate the five terms of Equ (16) in the reference respectively, which are used for calculating Bayes factor, AIC, BIC and DIC
# Reference: S. Chib, "Estimation and comparison of multiple change-point models," Journal of Econometrics, vol. 86, no. 2, pp. 221-241, 1998.
######################################################################################################################################

#The first term: Calculate log P(Y|theta*,p*)
log_y=function(y_liklihood,y_logliklihood,p_array,
               num_period,num_covariate,num_state){
  #We need to calculate log P(y_t|Y_(t-1),theta*,P*)
  s_store=s_YthetaP(p_array,y_liklihood,num_period,num_state,num_covariate)#P(s_t|Y_t,P,beta,sigma2)
  p_array2=p_array[-length(p_array)]
  loglik=0
  temp5=c()
  for(i in 1:num_period){
    #P(s_t|Y_{t-1},P,beta,sigma2)
    if(i==1){
      s_temp=s_store[1,]
    }else{
      temp1=s_store[i-1,]
      s_temp1=temp1*p_array
      temp1=temp1[-length(temp1)]
      s_temp2=c(0,temp1*(1-p_array2))
      s_temp=s_temp1+s_temp2#s_t|Y_(t-1)
    }
    #P(y_t|Y_(t-1),theta*,P*,s_t)=P(y_t|s_t)
    temp=y_logliklihood[i,]#Note that it should be the true log-likelihood of y
    maxtemp=max(temp)
    temp=temp-maxtemp
    temp=sum(exp(temp)*s_temp)
    if(temp==0){#A strong outlier
      temp=-1000000
    }else{
      temp=log(temp)+maxtemp
    }
    temp5=c(temp5,temp)
    loglik=loglik+temp#temp is log P(y_t|Y_(t-1),theta*,P*,s_t), loglik is log(f(Y_n|theta*,P*))
  }
  return(loglik)
}

#The second term: Calculate log pi(beta*,theta2)
logpi_beta_theta2=function(beta_matrix,sigma2_array,num_state,num_covariate,
                           k0_prior,v0_prior,sigma02_prior){
  logpi_beta=0
  logpi_theta2=0
  for(i in 1:num_state){
    logpi_theta2=logpi_theta2+log_dinvchisq(sigma2_array[i],v0_prior,sigma02_prior)
    logpi_beta=logpi_beta+dmnorm(beta_matrix[,i],rep(0,num_covariate),diag(sigma2_array[i]/k0_prior,num_covariate),log=T)
  }
  return(logpi_beta+logpi_theta2)
}

#The third term: Calculate log pi(P*)
logpi_p=function(p_array,num_state,a0,b0){
  if(a0==1&&b0==1){#Beta(1,1) is U(0,1), which is the uniform distribution
    logpi3=(num_state-1)*dbeta(0.5,a0,b0)
  }else{
    logpi3=0
    if(num_state>1){
      for(i in 1:(num_state-1)){
        logpi3=logpi3+dbeta(p_array[i],a0,b0)
      }
    }
  }
  return(logpi3)
}

#The forth term: Calculate log pi(beta*,sigma2*|Y_n) from the posterior probabilities of beta and sigma^2 during MCMC simulation
logpi_beta_sigma=function(pi_loglik_store){
  maxtemp=max(pi_loglik_store)
  pi_loglik_store=pi_loglik_store-maxtemp
  posterior_p=log(mean(exp(pi_loglik_store)))+maxtemp
  return(posterior_p)
}

#The fifth term: Calculate log P(P*|Y,theta*), which is approximated using additional MCMC draws
log_posterior_p=function(y_liklihood,p_array,s_sample,num_state,num_period,num_covariate,a0,b0){
  num_sampling=5000+1000
  num_burnin=1000
  log_p=rep(0,num_sampling)
  #Start simulation
  for(i in 1:num_sampling){
    #Sample P and S
    p_arraytemp=simu_P(s_sample,num_state,num_period,a0,b0)
    s_sample=simu_S(p_arraytemp,y_liklihood,num_period,num_state,num_covariate)
    #Calculate pi(P|S)
    if(num_state>1){
      for(j in 1:(num_state-1)){
        nii=length(s_sample[which(s_sample==j)])-1
        log_p[i]=log_p[i]+dbeta(p_array[j],a0+nii,b0+1,log=T)#log [pi(p_11|S_(n,j))*pi(p_22|S_(n,j))*...*pi(p_KK|S_(n,j))]
      }
    }
  }
  log_p=log_p[-c(1:num_burnin)]
  log_p=log(mean(exp(log_p)))
  return(log_p)
}

#Calculate log marginal likelihood based on Equ (16)
log_marginal_lik=function(y_sequence,y_liklihood,p_array,mu_matrix,beta_matrix,sigma2_array,pi_loglik_store,s_sample,
                          num_period,num_covariate,num_state,k0_prior,v0_prior,sigma02_prior,a0,b0){
  y_logliklihood=logy_Ytheta2(y_sequence,num_covariate,mu_matrix,num_state,num_period,sigma2_array)
  logpi1=log_y(y_liklihood,y_logliklihood,p_array,num_period,num_covariate,num_state)
  logpi2=logpi_beta_theta2(beta_matrix,sigma2_array,num_state,num_covariate,k0_prior,v0_prior,sigma02_prior)
  logpi3=logpi_p(p_array,num_state,a0,b0)
  logpi4=logpi_beta_sigma(pi_loglik_store)
  logpi5=log_posterior_p(y_liklihood,p_array,s_sample,num_state,num_period,num_covariate,a0,b0)
  cat("loglik1:",logpi1,"loglik2:",logpi2,"loglik3:",logpi3,"loglik4:",logpi4,"loglik5:",logpi5,"\n")
  paralist=list(loglik_bar=logpi1,margin=logpi1+logpi2+logpi3-logpi4-logpi5)
  return(paralist)
}

################################################################
# Calculate bar(log-likelihood) in the expression of DIC
# Reference: D. J. Spiegelhalter, N. G. Best, B. P. Carlin, and A. Van Der Linde, "Bayesian measures of model complexity and fit," Journal of the Royal Statistical Society: Series b (Statistical Methodology), vol. 64, no. 4, pp. 583-639, 2002.
################################################################

#For each MCMC draws theta=(beta, sigma^2) and p, calculate log-likelihood of y
bar_loglik_y=function(y_sequence,x_matrix,
                      beta_store,sigma2_store,p_store,
                      num_period,num_covariate,num_state){
  num_sampling=length(beta_store[,1])
  loglik_storestore=rep(0,num_sampling)
  for(k in 1:num_sampling){
    beta_matrix=matrix(beta_store[k,],num_covariate,num_state)
    sigma2_array=sigma2_store[k,]
    p_array=p_store[k,]
    mu_matrix=cal_mu(beta_matrix,x_matrix,num_state,num_period,num_covariate)
    y_liklihood=y_Ytheta(y_sequence,num_covariate,mu_matrix,num_state,num_period,sigma2_array)
    y_logliklihood=logy_Ytheta2(y_sequence,num_covariate,mu_matrix,num_state,num_period,sigma2_array)
    loglik_storestore[k]=log_y(y_liklihood,y_logliklihood,p_array,
                               num_period,num_covariate,num_state)
  }
  loglik_storestore=mean(loglik_storestore)
  return(loglik_storestore)
}


##########################################################
#  Calculate the model comparison criteria of our models #
##########################################################

#Input data
osname="D:/change-point/robust_test_day/"#The absolute path that stores the input data
setwd(osname)
regiondata=read.csv("region_info.csv")
#The data frame that includes the information of regions (the regions of bike sharing stations).
#Alternatively, we can use regiondata=data.frame(region_label=c(1:num_region)), where num_region is the total number of regions
num_region=length(regiondata[,1])
data=read.csv("region_generation_final_daily.csv")
data=data[,-1]
weatherstr=c("aqi","prec","fog","rain","snow","press","rh","temp","wind","day")
data=data[which(data$isweekend==0),]#Remove weekend and holidays
data=data[which(data$isholiday==0),]
num_period=length(data[,1])

#Hyperparameters
a0=1
b0=1
k0_prior=1
v0_prior=1
sigma02_prior=1

#Read the files that store the information of demand change detection outcomes
osname_dailypara=paste(osname,"parastore_covariate/",sep="")#The absolute path that stores these files
osname_dailypara2=paste(osname,"parastore_covariate_unify/",sep="")#The absolute path that stores the new beta and sigma^2
setwd(osname)
storedata=read.csv("summary_null_model.csv",as.is=T)
storedata=storedata[,-1]
num_sampling=5000#The model comparison criteria are calculated based on the last "num_sampling" MCMC iterations
#Start calculation
for(regionindex in 1:num_region){
  #Input data
  num_state=storedata[regionindex,"num_change"]+1
  chosen_station=regiondata[regionindex,"region_label"]+1
  temp=data[,chosen_station]
  y_sequence=as.numeric(temp,num_period,1)
  datatemp=data.frame(data[,weatherstr])
  datatemp$intercept=1
  num_covariate=length(datatemp[1,])
  x_matrix=as.matrix(datatemp,num_period,num_covariate)
  #Read the MCMC draws of model parameters
  setwd(osname_dailypara2)
  beta_store=read.csv(paste("beta",as.character(as.integer(chosen_station)),".csv",sep=""),as.is=T)
  sigma2_store=read.csv(paste("sigma2",as.character(as.integer(chosen_station)),".csv",sep=""),as.is=T)
  setwd(osname_dailypara)
  s_store=read.csv(paste("s_store",as.character(as.integer(chosen_station)),".csv",sep=""),as.is=T)
  beta_store=as.matrix(beta_store[,-1])
  sigma2_store=as.matrix(sigma2_store[,-1])
  s_store=as.matrix(s_store[,-1])
  #Select the last "num_sampling" MCMC iterations
  num_samplingtemp=length(s_store[,1])
  beta_store=as.matrix(beta_store[(num_samplingtemp-num_sampling+1):num_samplingtemp,])
  sigma2_store=as.matrix(sigma2_store[(num_samplingtemp-num_sampling+1):num_samplingtemp,])
  s_store=as.matrix(s_store[(num_samplingtemp-num_sampling+1):num_samplingtemp,])
  #Rename the indexes of states. For example, the reachable states {1,3,4,6} would be renamed as {1,2,3,4}
  old_state=unique(as.integer(s_store[1,]))
  new_state=c(1:length(old_state))
  for(i in 1:length(s_store[1,])){
    for(j in 1:length(old_state)){
      s_store[which(s_store[,i]==old_state[j]),i]=new_state[j]
    }
  }
  #Posterior mean of regression parameters
  temp=as.numeric(matrix(1/num_sampling,1,num_sampling)%*%beta_store)
  beta_matrix=matrix(temp,num_covariate,num_state)
  sigma2_array=as.numeric(matrix(1/num_sampling,1,num_sampling)%*%sigma2_store)
  mu_matrix=cal_mu(beta_matrix,x_matrix,num_state,num_period,num_covariate)
  #Posterior estimates of s
  s_sample=c()
  for(i in 1:num_period){
    temp=as.integer(names(sort(table(s_store[,i]),decreasing=T))[1])
    s_sample=c(s_sample,temp)
  }
  #Posterior mean of p, which is obtained using additional MCMC simulation
  p_store=matrix(NA,num_sampling,num_state)
  for(i in 1:num_sampling){
    s_sampletemp=as.numeric(s_store[i,])
    p_store[i,]=simu_P(s_sampletemp,num_state,num_period,a0,b0)
  }
  p_array=matrix(1/num_sampling,1,num_sampling)%*%p_store
  p_array[num_state]=1
  #Calculate some log likelihoods
  y_liklihood=y_Ytheta(y_sequence,num_covariate,mu_matrix,num_state,num_period,sigma2_array)
  pi_loglik_store=rep(NA,num_sampling)
  for(i in 1:num_sampling){
    beta_matrixtemp=matrix(beta_store[i,],num_covariate,num_state)
    sigma2_arraytemp=sigma2_store[i,]
    pi_loglik_store[i]=calcu_pi_loglik(s_sample,y_sequence,x_matrix,num_state,num_period,num_covariate,
                                       beta_matrixtemp,sigma2_arraytemp,
                                       k0_prior,v0_prior,sigma02_prior)
  }
  #Calculate Bayes factor and the log-likelihood of y
  paralist=log_marginal_lik(y_sequence,y_liklihood,p_array,mu_matrix,beta_matrix,sigma2_array,pi_loglik_store,s_sample,
                            num_period,num_covariate,num_state,k0_prior,v0_prior,sigma02_prior,a0,b0)
  bayes_factor=paralist$margin
  storedata[regionindex,"bayes_factor"]=bayes_factor#Actually, it is log marginal likelihoods.
  #Calculate DIC
  barloglik=bar_loglik_y(y_sequence,x_matrix,beta_store,sigma2_store,p_store,
                         num_period,num_covariate,num_state)
  logpi1=paralist$loglik_bar
  storedata[regionindex,"loglik"]=logpi1#Theoretically, logpi1>barloglik
  barD=-2*barloglik
  Dbar=-2*logpi1
  pD=barD-Dbar#It should be a positive value
  DIC=barD+pD
  storedata[regionindex,"DIC"]=DIC
  cat(regionindex,bayes_factor,DIC,"\n")
  #Calculate AIC and BIC
  num_total_covariate=(num_covariate+1)*length(old_state)+length(old_state)-1#Total number of parameters in our model
  storedata[regionindex,"AIC"]=-2*logpi1+2*num_total_covariate
  storedata[regionindex,"BIC"]=-2*logpi1+num_total_covariate*log(num_period)
  cat(regionindex,storedata[regionindex,"AIC"],storedata[regionindex,"BIC"],"\n")
}
#Save outcomes
setwd(osname)
write.csv(storedata,"summary_covariate_unify.csv")
print("finish!")
#Print the average model comparion criteria per region
for(i in c("loglik","bayes_factor","DIC","AIC","BIC")){
  cat(i,mean(storedata[,i]),"\n")
}
