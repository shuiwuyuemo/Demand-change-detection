
#################################################################################################################################
#    Resample state-specific beta and sigma^2 to unify the independent variables and dependent variables in different models    #
#################################################################################################################################

library(MCMCpack)
library(MASS)
library(geoR)
library(truncnorm)

#Sample the regression parameters in each state
simu_beta_sigma=function(s_sample,y_sequence,x_matrix,num_state,num_period,num_covariate,
                         k0_prior,v0_prior,sigma02_prior){
  #N(0,sigma^2/k0_prior) is the prior distribution of beta
  #inv-chisq(v0_prior,sigma02_prior) is the prior distribution of sigma2
  beta_matrix=matrix(NA,num_covariate,num_state)
  sigma2_array=rep(NA,num_state)
  for(i in 1:num_state){
    if(i%in%s_sample){
      y_sequencetemp=y_sequence[which(s_sample==i)]
      y_sequencetemp=matrix(y_sequencetemp,length(y_sequencetemp),1)
      x_matrixtemp=x_matrix[which(s_sample==i),]
      if(!is.matrix(x_matrixtemp)){#when there is only one observation in state i
        x_matrixtemp=matrix(x_matrixtemp,1,num_covariate)
      }
      temp=solve(t(x_matrixtemp)%*%x_matrixtemp+k0_prior*diag(1,num_covariate))
      betatemp=temp%*%t(x_matrixtemp)%*%y_sequencetemp
      #Sample sigma^2
      n=length(s_sample[which(s_sample==i)])
      vn=v0_prior+n
      mutemp=x_matrixtemp%*%matrix(betatemp,num_covariate,1)
      residual=y_sequencetemp-mutemp
      sn=t(residual)%*%residual+k0_prior*t(betatemp)%*%betatemp+v0_prior*sigma02_prior
      sigma2_array[i]=rinvchisq(1,df=vn,scale=sn/vn)
      #Sample beta
      sigma2temp=temp*sigma2_array[i]
      betasample=mvrnorm(1,betatemp,sigma2temp)
      beta_matrix[,i]=betasample
    }else{
      #When state i is unreachable, sample the state from its prior distribution
      sigma2_array[i]=rinvchisq(1,df=v0_prior,scale=sigma02_prior)
      beta_matrix[,i]=mvrnorm(1,rep(0,num_covariate),diag(sigma2_array[i]/k0_prior,num_covariate))
    }
  }
  paralist=list(beta_matrix=beta_matrix,sigma2_array=sigma2_array)
  return(paralist)
}

###################Start resampling beta and sigma^2

#Input data
osname="D:/change-point/robust_test_day/"#The absolute path that stores the input data
osname_dailypara="D:/change-point/robust_test_day/parastore_covariate"#The absolute path that stores the model parameters
osname_dailypara2="D:/change-point/robust_test_day/parastore_covariate_unify"#The absolute path that stores the new beta and sigma^2
setwd(osname)
regiondata=read.csv("region_info.csv")
num_region=length(regiondata[,1])
data=read.csv("region_generation_final_daily.csv")
data=data[,-1]
weatherstr=c("aqi","prec","fog","rain","snow","press","rh","temp","wind","day")#Unified independent variables
data=data[which(data$isweekend==0),]
data=data[which(data$isholiday==0),]
num_period=length(data[,1])

#Hyperparameters
k0_prior=1#N(0,sigma^2/k0_prior) is the prior distribution of beta
v0_prior=1#inv-chisq(v0_prior,sigma02_prior) is the prior distribution of sigma2
sigma02_prior=1

for(regionindex in 1:num_region){
  #Input data
  temp=data[,regionindex]
  y_sequence=matrix(temp,num_period,1)
  datatemp=data.frame(data[,weatherstr])
  datatemp$intercept=1
  num_covariate=length(datatemp[1,])
  x_matrix=as.matrix(datatemp,num_period,num_covariate)
  #s sequence
  setwd(osname_dailypara)
  s_store=read.csv(paste("s_store",as.character(regionindex),".csv",sep=""))
  s_store=s_store[,-1]
  num_state=length(unique(as.integer(s_store[1,])))
  
  #Change the indexes of states. For example, we change state {1,3,4,6} into {1,2,3,4}
  old_state=unique(as.integer(s_store[1,]))
  new_state=c(1:length(old_state))
  for(i in 1:length(s_store[1,])){
    for(j in 1:length(old_state)){
      s_store[which(s_store[,i]==old_state[j]),i]=new_state[j]
    }
  }
  
  #Resample beta and sigma^2
  num_sampling=length(s_store[,1])
  beta_store=matrix(NA,num_sampling,num_state*num_covariate)
  sigma2_store=matrix(NA,num_sampling,num_state)
  #Start MCMC simulation
  for(i in 1:num_sampling){
    s_sample=as.integer(s_store[i,])
    paralist=simu_beta_sigma(s_sample,y_sequence,x_matrix,num_state,num_period,num_covariate,
                             k0_prior,v0_prior,sigma02_prior)
    beta_matrix=paralist$beta_matrix
    sigma2_array=paralist$sigma2_array
    beta_store[i,]=as.vector(beta_matrix)
    sigma2_store[i,]=sigma2_array
    if(i%%10000==0){
      cat("iteration",i,"\n")
    }
  }
  
  #Save beta and sigma^2
  setwd(osname_dailypara2)
  beta_storedata=data.frame(beta_store)
  write.csv(beta_storedata,paste("beta",as.character(regionindex),".csv",sep="",collapse=""))
  sigma2_storedata=data.frame(sigma2_store)
  write.csv(sigma2_storedata,paste("sigma2",as.character(regionindex),".csv",sep="",collapse=""))
  print(regionindex)
}
print("finish!!")
