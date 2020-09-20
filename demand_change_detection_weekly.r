
##################################################################
#       demand change detection with weekly-integrated data      #
##################################################################

library(MCMCpack)
library(MASS)
library(geoR)
library(truncnorm)

#Given the regression parameters of each state, calculate the estimates of bike sharing demands for all time points
cal_mu=function(beta_matrix,x_matrix,num_state,num_period,num_covariate){
  mu_matrix=matrix(NA,num_period,num_state)
  for(i in 1:num_state){
    betatemp=beta_matrix[,i]
    mu_matrix[,i]=x_matrix%*%matrix(betatemp,num_covariate,1)
  }
  return(mu_matrix)
}

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

#Sample state s. Note that we only need to sample s_t if s_t!=s_{t-1} or s_t!=s_{t+1}
simu_S=function(s_sample,y_sequence,mu_matrix,sigma2_array,num_period,a0,b0){
  #Beta(a0,b0) is the prior distribution of state transition probabilities
  #1. When s_1!=s_2, sample s_1
  if(s_sample[1]!=s_sample[2]){
    #calculate P(y_1|Y_{-1},theta_?)
    state1=s_sample[1]
    ytemp=y_sequence[1]
    mutemp=mu_matrix[1,state1]
    sigmatemp=sqrt(sigma2_array[state1])
    ytemp1=dnorm(ytemp,mutemp,sigmatemp,log=T)
    state2=s_sample[2]
    mutemp=mu_matrix[1,state2]
    sigmatemp=sqrt(sigma2_array[state2])
    ytemp2=sum(dnorm(ytemp,mutemp,sigmatemp,log=T))
    maxtemp=max(ytemp1,ytemp2)#To avoid infinite values, we minus the larger value, which does not influence the final outcome
    ytemp1=ytemp1-maxtemp
    ytemp2=ytemp2-maxtemp
    #Calculate P(s_1|s_2,s_3,...) and sample s_1
    p1=a0*b0/(a0+b0)^2*exp(ytemp1)
    n_s2s2=length(s_sample[which(s_sample==state2)])-1
    p2=b0/(a0+b0)*(n_s2s2+a0)/(n_s2s2+b0+a0)*exp(ytemp2)
    p1=p1/(p1+p2)
    u=runif(1)
    if(u<p1){
      s_sample[1]=state1
    }else{
      s_sample[1]=state2
    }
  }
  #2. When s_(t-1)=i and s_(t+1)=i+1, sample s_t
  for(i in 2:(num_period-1)){
    if(s_sample[i-1]==s_sample[i+1]){
      next
    }
    s_sample[i]=0
    #Calculate P(y_t|Y_(-t),theta_?)
    state1=s_sample[i-1]
    ytemp=y_sequence[i]
    mutemp=mu_matrix[i,state1]
    sigmatemp=sqrt(sigma2_array[state1])
    ytemp1=sum(dnorm(ytemp,mutemp,sigmatemp,log=T))
    state2=s_sample[i+1]
    mutemp=mu_matrix[i,state2]
    sigmatemp=sqrt(sigma2_array[state2])
    ytemp2=sum(dnorm(ytemp,mutemp,sigmatemp,log=T))
    maxtemp=max(ytemp1,ytemp2)#To avoid infinite values, we minus the larger value, which does not influence the final outcome
    ytemp1=ytemp1-maxtemp
    ytemp2=ytemp2-maxtemp
    #Calculate p(s_t|s_(t-1)) and p(s_(t+1)|s_t)
    n11=length(s_sample[which(s_sample==state1)])-1
    p_font1=(n11+a0)/(n11+b0+a0)
    p_back1=b0/(n11+1+b0+a0)
    p_font2=b0/(n11+b0+a0)
    n22=length(s_sample[which(s_sample==state2)])-1
    p_back2=(n22+a0)/(n22+b0+a0)
    #Calculate P(s_t|s_(-t)) and sample s_t
    p1=p_font1*p_back1*exp(ytemp1)
    p2=p_font2*p_back2*exp(ytemp2)
    p1=p1/(p1+p2)
    u=runif(1)
    if(u<p1){
      s_sample[i]=state1
    }else{
      s_sample[i]=state2
    }
  }
  #3. When s_T!=s_(T-1), sample s_T
  if(s_sample[num_period-1]!=s_sample[num_period]){
    state1=s_sample[num_period-1]
    n_s1s1=length(s_sample[which(s_sample==state1)])-1
    #Calculate P(y_T|Y_(-T),theta_?)
    ytemp=y_sequence[num_period]
    mutemp=mu_matrix[num_period,state1]
    sigmatemp=sqrt(sigma2_array[state1])
    ytemp1=sum(dnorm(ytemp,mutemp,sigmatemp,log=T))
    state2=s_sample[num_period]
    mutemp=mu_matrix[num_period,state2]
    sigmatemp=sqrt(sigma2_array[state2])
    ytemp2=sum(dnorm(ytemp,mutemp,sigmatemp,log=T))
    maxtemp=max(ytemp1,ytemp2)#To avoid infinite values, we minus the larger value, which does not influence the final outcome
    ytemp1=ytemp1-maxtemp
    ytemp2=ytemp2-maxtemp
    #Calculate P(s_T|s_(-T)) and sample s_T
    p1=(n_s1s1+a0)/(n_s1s1+b0+a0)*exp(ytemp1)
    p2=b0/(n_s1s1+b0+a0)*exp(ytemp2)
    p1=p1/(p1+p2)
    u=runif(1)
    if(u<p1){
      s_sample[num_period]=state1
    }else{
      s_sample[num_period]=state2
    }
  }
  return(s_sample)
}


####################################
#      Numerical experiments       #
#  We build M_{basic} in this code #
####################################

osname="D:/change-point/robust_test_weekly/"#The absolute path that stores the input data
setwd(osname)
regiondata=read.csv("region_info.csv")
#The data frame that includes the information of regions (the regions of bike sharing stations).
#Alternatively, we can use regiondata=data.frame(region_label=c(1:num_region)), where num_region is the total number of regions
num_region=length(regiondata[,1])
data=read.csv("region_generation_final_weekly.csv")
data=data[,-1]
weatherstr=c("aqi","prec","fog","rain","snow","press","rh","temp","wind")
num_period=length(data[,1])

#Hyperparameters
a0=1#Beta(a0,b0) is the prior distribution of state transition probabilities
b0=1
k0_prior=1#N(0,sigma^2/k0_prior) is the prior distribution of beta
v0_prior=1#inv-chisq(v0_prior,sigma02_prior) is the prior distribution of sigma2
sigma02_prior=1

#The information of demand changes which have been detected by our model.
if(T){#When we start demand change detection
  #Initialization
  storedata=data.frame(regiondata$region_label)#The data frame that includes: num_sampling, num_change, compute_time
  storedata$loc_change=NA
  storedata$num_change=NA
  storedata$num_try=0
  locstoredata=data.frame(matrix(NA,num_region,num_period))#The data frame that stores the final posterior estimates of the timing of demand changes
  locstoredata0=locstoredata#The data frame that stores the estimates of the timing of demand changes in the first stage of our implementation procedure
  num_converge=0
}else{#After we detect the demand changes of several regions, and ready to detect the demand changes of the rest regions
  #The information of the regions with detected demand changes
  storedata=read.csv("summary_covariate.csv",as.is=T)
  locstoredata=read.csv("loc_change_covariate.csv",as.is=T)
  locstoredata0=read.csv("loc_change0_covariate.csv",as.is=T)
  storedata=storedata[,-1]
  locstoredata=locstoredata[,-1]
  locstoredata0=locstoredata0[,-1]
  num_converge=length(na.omit(storedata$num_change))
}
names(locstoredata)=c(1:num_period)
names(locstoredata0)=c(1:num_period)
num_circle=0
while(num_converge!=num_region){
  num_circle=num_circle+1
  for(regionindex in 1:num_region){
    if(!is.na(storedata[regionindex,"num_change"])){
      next
    }
    storedata[regionindex,"num_try"]=storedata[regionindex,"num_try"]+1
    chosen_station=regiondata[regionindex,"region_label"]+1
    #input data
    y_sequence=matrix(data[,chosen_station],num_period,1)
    datatemp=data.frame(data[,weatherstr])
    datatemp$intercept=1
    num_covariate=length(datatemp[1,])
    x_matrix=as.matrix(datatemp,num_period,num_covariate)
    num_state=20#Initial number of states
    
    ###########The first stage of our implementation procedure: deriving the candicate demand changes###########
    #The number of iteration
    max_sampling=200000#The upper limit for the number of iteration
    num_sampling=10000#Examine the reachable states per "num_sampling" iterations
    num_chain=4#The number of MCMC simulation chains
    #Store the demand changes detected by these MCMC simulation chains
    chain_s_store=matrix(NA,num_chain,num_sampling*num_period)#Store the state sequence
    chain_state_store=matrix(NA,num_chain,num_sampling*num_state)#Store the indexes of states, since some states will become unreachable
    chain_loc_change_store=matrix(NA,num_chain,num_sampling*num_state)#Store the timing of demand changes
    last_s_store=matrix(NA,num_chain,num_period)#Store the state sequence of the last MCMC draw. We use it to start the next "num_sampling" simulations
    ###Start simulation###
    first_iter=rep(T,num_chain)#Indicate whether an MCMC chain has started simulation
    start_time=proc.time()
    signaltemp=F
    for(chain in 1:num_chain){
      stop_criteria=0#Whether the MCMC simulation roughly converges
      current_sample=0#Number of iterations
      while(stop_criteria<1){
        if(first_iter[chain]){#When the MCMC chain has not started yet, we apply an equidistant state sequence
          s_sample=rep(1,num_period)
          start1=floor(num_period/num_state)
          for(i in 2:num_state){
            start2=floor(i*num_period/num_state)
            s_sample[start1:start2]=i
            start1=start2+1
          }
          first_iter[chain]=F
        }else{
          #When the MCMC chain has started, use the state sequence of the last iteration as the initial state sequence
          s_sample=last_s_store[chain,]
        }
        #Store the next "num_sampling" MCMC draws of state sequences
        s_store=matrix(NA,num_sampling,num_period)
        state_store=matrix(NA,num_sampling,num_state)
        loc_change_store=matrix(NA,num_sampling,num_state)
        #Start MCMC sampling
        for(i in 1:num_sampling){
          paralist=simu_beta_sigma(s_sample,y_sequence,x_matrix,num_state,num_period,num_covariate,
                                   k0_prior,v0_prior,sigma02_prior)
          beta_matrix=paralist$beta_matrix
          sigma2_array=paralist$sigma2_array
          mu_matrix=cal_mu(beta_matrix,x_matrix,num_state,num_period,num_covariate)
          s_sample=simu_S(s_sample,y_sequence,mu_matrix,sigma2_array,num_period,a0,b0)
          #Store the MCMC draws
          s_store[i,]=s_sample
          stateset=unique(s_sample)
          for(j in 1:length(stateset)){
            state_store[i,j]=stateset[j]
            loc_change_store[i,j]=length(s_sample[which(s_sample<=stateset[j])])+1#Note that these include T+1
          }
          #Print per 1000 iterations
          if(i%%1000==0){
            cat("chain",chain,", iteration",i,"\n")
          }
        }
        last_s_store[chain,]=s_sample#Store the state sequence of the last iteration
        current_sample=current_sample+num_sampling
        
        #######Check whether the number of demand changes remains constant for "2*num_sampling" iterations#######
        if(current_sample>max_sampling){#If the number of iterations exceeds the upper limit, skip this region
          signaltemp=T
          print("The number of iterations exceeds the threshold")
          break
        }else{
          if(length(unique(s_store[1,]))!=length(unique(s_store[num_sampling,]))){#The number of demand changes is not constant
            cat("The number of demand changes before:",length(unique(s_store[1,])),"; The number of demand changes after:",length(unique(s_store[num_sampling,])),"\n")
            stop_criteria=0
          }else{#The number of demand changes is constant
            final_num_sampling=current_sample
            is_converge=T
            cat("Converge after",current_sample,"iterations","\n")
            stop_criteria=stop_criteria+0.5#If the number of demand changes is constant for "2*num_sampling" iterations, the algorithm roughly converge
          }
        }
      }
      #Store the MCMC draws of this MCMC chain
      chain_s_store[chain,]=as.vector(s_store)
      chain_state_store[chain,]=as.vector(state_store)
      chain_loc_change_store[chain,]=as.vector(loc_change_store)
    }
    if(signaltemp){next}#The number of iterations exceeds the upper limit, skip this region
    duration_time=(proc.time()-start_time)/3600
    cat("Duration:",duration_time[3],"\n")
    
    ###########Derive candidate demand changes###########
    mode_precent=0.5#We need to select the time points whose detected probabilities are over "mode_precent"
    thres=num_sampling*num_chain*mode_precent
    dist_thres=2#The degree of dispersion. When the two candicate demand changes are within dist_thres time points, they would be treated as the same demand change
    ###Start calculation###
    chain_loc_change_store2=as.vector(chain_loc_change_store)
    chain_loc_change_store2=na.omit(chain_loc_change_store2)
    chain_loc_change_store2=chain_loc_change_store2[which(chain_loc_change_store2!=num_period+1)]
    #Select the time points whose detected probabilities are over "mode_precent". They are parts of the candicate demand changes
    temp=table(chain_loc_change_store2)
    temp=sort(temp,decreasing=T)
    mode_select=as.numeric(names(temp[which(temp>thres)]))
    mode_select=as.integer(mode_select)
    mode_select=sort(mode_select)
    #Select the time points which are the modes of the detected probability distribution
    density_temp=density(chain_loc_change_store2)
    peak=c()
    valleyleft=c(density_temp$x[1])
    valleyright=c()
    findleft=T
    findcenter=F
    findright=F
    for(j in 2:(length(density_temp$x)-1)){
      if(density_temp$y[j-1]<density_temp$y[j]&&density_temp$y[j]>=density_temp$y[j+1]&&findleft){
        peak=c(peak,density_temp$x[j])
        findleft=F
        findcenter=T
      }else if(density_temp$y[j-1]>density_temp$y[j]&&density_temp$y[j]<=density_temp$y[j+1]&&findcenter){
        valleyright=c(valleyright,density_temp$x[j])
        findright=T
        findcenter=F
      }else if(density_temp$y[j-1]>=density_temp$y[j]&&density_temp$y[j]<density_temp$y[j+1]&&findright){
        valleyleft=c(valleyleft,density_temp$x[j])
        findright=F
        findleft=T
      }
    }
    valleyright=c(valleyright,max(density_temp$x))
    #Add these modes as the candicate demand changes if the detection probabilities within "dist_thres" from them exceed "mode_precent"
    for(j in 1:length(peak)){
      temp=chain_loc_change_store2[which(chain_loc_change_store2>=peak[j]-dist_thres)]
      temp=temp[which(temp<=peak[j]+dist_thres)]
      #Whether there are candidate demand changes within "dist_thres" from them
      temp2=mode_select[which(mode_select>valleyleft[j])]
      temp2=temp2[which(temp2<valleyright[j])]
      temp3=mode_select[which(mode_select>peak[j]-dist_thres)]
      temp3=temp3[which(temp3<peak[j]+dist_thres)]
      if(length(temp)>thres&&length(temp2)==0&&length(temp3)==0){
        mode_select=c(mode_select,as.integer(peak[j]))
      }
    }
    mode_select=unique(mode_select)
    mode_select=sort(mode_select)
    #Save the candidate demand changes
    storedata[regionindex,"num_change0"]=length(mode_select)
    storedata[regionindex,"loc_change0"]=paste(as.character(mode_select),sep="",collapse=",")
    locstoredata0[regionindex,]=0
    locstoredata0[regionindex,mode_select]=1
    
    #######The second stage of our implementation procedure: deriving the final demand changes###########
    num_state=length(mode_select)+1
    #Number of MCMC iterations
    max_sampling=1000000#The upper limits of the number of iterations before the number of demand changes remains constant
    max_sampling2=2000000#The upper limits of the number of iterations before the algorithm converges
    final_num_sampling=0#The number of iterations after the burnin period
    current_sample=0#The number of iterations have been processed
    num_sampling=10000#Examine the algorithm convergence per "num_sampling" iterations
    thres_sampling=20000#The lower limits of the number of iterations
    num_chain=2#Number of MCMC chains
    #Store the MCMC draws in these MCMC simulation chains
    chain_s_store=matrix(NA,num_chain,num_sampling*num_period)#State sequence
    chain_beta_store=matrix(NA,num_chain,num_sampling*num_state*num_covariate)
    chain_sigma2_store=matrix(NA,num_chain,num_sampling*num_state)
    chain_state_store=matrix(NA,num_chain,num_sampling*num_state)#The indexes of states
    chain_loc_change_store=matrix(NA,num_chain,num_sampling*num_state)#Timing of demand changes
    chain_num_change_store=matrix(NA,num_chain,2)#Store the numbers of demand changes in the first and the last iterations of the current "num_sampling" iterations
    last_s_store=matrix(NA,num_chain,num_period)#The state sequence in the last iteration
    #For convergence inference
    stop_criteria=T#Whether the algorithm converges
    num_stable=F#Whether the number of demand changes is stable: all MCMC chains have the same numbers of demand changes, and the number of demand changes remains constant in the current "num_sampling" iterations.
    num_temp=0#Half the number of iterations after all chains have the same numbers of demand changes
    first_iter=rep(T,num_chain)#Whether the MCMC chain needs to be initialized
    final_num_change=-1#The final number of demand changes
    R_thres=1.2#When the corrected potential scale reductions of all parameters are less than "R_thres", the algorithm converges
    start_time=proc.time()
    is_converge=F#Whether the algorithm converges
    start_time=proc.time()
    ###Start simulation###
    while(stop_criteria){
      if(current_sample>max_sampling2){
        stop_criteria=F
        print("The number of iterations exceeds the threshold")
        break
      }
      #After the number of demand changes is stable, we use the last half of MCMC iterations to judge whether the algorithm converges
      if(num_stable){#The number of demand changes is stable
        #The matrices that store the model parameters
        lentemp=as.integer(num_temp+num_sampling/2)
        chain_s_store_temp=chain_s_store
        chain_s_store=matrix(NA,num_chain,lentemp*num_period)
        chain_beta_store_temp=chain_beta_store
        chain_beta_store=matrix(NA,num_chain,lentemp*num_state*num_covariate)
        chain_sigma2_store_temp=chain_sigma2_store
        chain_sigma2_store=matrix(NA,num_chain,lentemp*num_state)
        chain_state_store_temp=chain_state_store
        chain_state_store=matrix(NA,num_chain,lentemp*num_state)
        chain_loc_change_store_temp=chain_loc_change_store
        chain_loc_change_store=matrix(NA,num_chain,lentemp*num_state)
      }else{
        #the number of demand changes is NOT stable
        chain_s_store=matrix(NA,num_chain,num_sampling*num_period)
        chain_beta_store=matrix(NA,num_chain,num_sampling*num_state*num_covariate)
        chain_sigma2_store=matrix(NA,num_chain,num_sampling*num_state)
        chain_state_store=matrix(NA,num_chain,num_sampling*num_state)
        chain_loc_change_store=matrix(NA,num_chain,num_sampling*num_state)
      }
      
      ###MCMC simulation###
      for(chain in 1:num_chain){
        if(first_iter[chain]){#Initialize the state sequence using the candidate demand changes
          s_sample=rep(1,num_period)
          if(num_state>2){
            for(j in 1:(num_state-2)){
              s_sample[mode_select[j]:(mode_select[j+1]-1)]=j+1
            }
          }
          if(num_state>1){
            s_sample[mode_select[num_state-1]:num_period]=num_state
          }
          first_iter[chain]=F
        }else{
          #Use the state sequence in the last iteration
          s_sample=last_s_store[chain,]
        }
        #Store the MCMC draws in the next "num_sampling" iterations
        s_store=matrix(NA,num_sampling,num_period)
        state_store=matrix(NA,num_sampling,num_state)
        loc_change_store=matrix(NA,num_sampling,num_state)
        beta_store=matrix(NA,num_sampling,num_state*num_covariate)
        sigma2_store=matrix(NA,num_sampling,num_state)
        for(i in 1:num_sampling){
          #MCMC sampling
          paralist=simu_beta_sigma(s_sample,y_sequence,x_matrix,num_state,num_period,num_covariate,
                                   k0_prior,v0_prior,sigma02_prior)
          beta_matrix=paralist$beta_matrix
          sigma2_array=paralist$sigma2_array
          mu_matrix=cal_mu(beta_matrix,x_matrix,num_state,num_period,num_covariate)
          s_sample=simu_S(s_sample,y_sequence,mu_matrix,sigma2_array,num_period,a0,b0)
          #Store
          s_store[i,]=s_sample
          stateset=unique(s_sample)
          for(j in 1:length(stateset)){
            state_store[i,j]=stateset[j]
            loc_change_store[i,j]=length(s_sample[which(s_sample<=stateset[j])])+1#Note that these include T+1
          }
          beta_store[i,]=as.vector(beta_matrix)
          sigma2_store[i,]=sigma2_array
          #Print per 1000 iteration
          if(i%%1000==0){
            cat("chain",chain,", iteration",i,"\n")
          }
        }
        #Store the number of demand changes for convergence inference
        temp=unique(na.omit(state_store[1,]))
        chain_num_change_store[chain,1]=length(temp)-1
        temp=unique(na.omit(state_store[num_sampling,]))
        chain_num_change_store[chain,2]=length(temp)-1
        #Store the simulation outcome of this chain
        if(!num_stable){#When the number of demand changes is not stable, abandon the MCMC draws in previous iterations
          chain_s_store[chain,]=as.vector(s_store)
          chain_state_store[chain,]=as.vector(state_store)
          chain_loc_change_store[chain,]=as.vector(loc_change_store)
          chain_beta_store[chain,]=as.vector(beta_store)
          chain_sigma2_store[chain,]=as.vector(sigma2_store)
        }else{#When the number of demand changes is stable, abandon the first "num_sampling/2" iterations and enlarge the storage matrices
          temp=chain_s_store_temp[chain,]
          temp=matrix(temp,num_temp,num_period)
          temp=as.matrix(temp[-c(1:as.integer(num_sampling/2)),])
          temp=rbind(temp,s_store)
          temp=as.vector(temp)
          chain_s_store[chain,]=temp#state sequence
          temp=chain_beta_store_temp[chain,]
          temp=matrix(temp,num_temp,num_state*num_covariate)
          temp=as.matrix(temp[-c(1:as.integer(num_sampling/2)),])
          temp=rbind(temp,beta_store)
          temp=as.vector(temp)
          chain_beta_store[chain,]=temp#beta
          temp=chain_sigma2_store_temp[chain,]
          temp=matrix(temp,num_temp,num_state)
          temp=as.matrix(temp[-c(1:as.integer(num_sampling/2)),])
          temp=rbind(temp,sigma2_store)
          temp=as.vector(temp)
          chain_sigma2_store[chain,]=temp#sigma^2
          temp=chain_state_store_temp[chain,]
          temp=matrix(temp,num_temp,num_state)
          temp=as.matrix(temp[-c(1:as.integer(num_sampling/2)),])
          temp=rbind(temp,state_store)
          temp=as.vector(temp)
          chain_state_store[chain,]=temp#The indexes of states
          temp=chain_loc_change_store_temp[chain,]
          temp=matrix(temp,num_temp,num_state)
          temp=as.matrix(temp[-c(1:as.integer(num_sampling/2)),])
          temp=rbind(temp,loc_change_store)
          temp=as.vector(temp)
          chain_loc_change_store[chain,]=temp#Timing of demand changes
        }
        last_s_store[chain,]=s_sample#The state sequence of the last iteration
      }
      current_sample=current_sample+num_sampling
      
      ###Convergence inference###
      # [1] A. Gelman and D. Rubin, "Inference from iterative simulation using multiple sequences," Statistical Science, vol. 7, no. 4, pp.45-472, 1992
      # [2] S. P. Brooks and A. Gelman, "General methods for monitoring convergence of iterative simulations," Journal of Computational and Graphical Statistics, vol. 7, no. 4, pp. 434-455, 1998.
      ########################################
      if(current_sample>=thres_sampling){#The number of iterations exceeds the lower limit of iterations
        if(max(chain_num_change_store[1:num_chain,])!=min(chain_num_change_store[1:num_chain,])){#The number of demand changes is not stable
          num_stable=F
          num_temp=0
          cat("Iteration",current_sample,": The number of demand changes is not stable","\n")
          print(chain_num_change_store)
          if(current_sample>max_sampling){#The number of iterations exceeds the upper limit of iterations
            stop_criteria=F
            print("The number of iterations exceeds the threshold")
            break
          }
          next
        }
        final_num_change=min(chain_num_change_store[1:num_chain,])
        #Convergence inference
        if(final_num_change>0){
          #1. Test the timing of demand changes
          meantemp=matrix(0,num_chain,final_num_change)
          W=rep(0,final_num_change)
          s2temp=matrix(0,num_chain,final_num_change)
          R_hat=rep(0,final_num_change)
          for(chain in 1:num_chain){
            temp=matrix(chain_loc_change_store[chain,],num_temp,num_state)
            temp=as.matrix(temp[,1:final_num_change])
            for(j in 1:final_num_change){
              meantemp[chain,j]=mean(temp[,j])
              s2temp[chain,j]=var(temp[,j])
              W[j]=W[j]+s2temp[chain,j]
            }
          }
          W=W/num_chain
          for(j in 1:final_num_change){
            Bn=var(meantemp[,j])
            sigma2_plus=(num_temp-1)/num_temp*W[j]+Bn
            V_hat2=sigma2_plus+Bn/num_chain
            var_V_hat=((num_temp-1)/num_temp)^2/num_chain*var(s2temp[,j])+
              ((num_chain+1)/num_chain)^2*2/(num_chain-1)*Bn^2+
              2*(num_chain+1)*(num_temp-1)/(num_chain^2)/num_temp*
              (cov(s2temp[,j],meantemp[,j]^2)-2*mean(meantemp[,j])*cov(s2temp[,j],meantemp[,j]))
            df=2*V_hat2/var_V_hat
            #R_hat[j]=V_hat2/W[j]*df/(df-2)
            R_hat[j]=V_hat2/W[j]*(df+3)/(df+1)
            if(V_hat2==0||W[j]==0||is.na(R_hat[j])){
              R_hat[j]=1
            }
          }
          R_hat=sqrt(R_hat)
          if(max(R_hat)>R_thres){
            num_stable=F
            num_temp=0
            cat("Iteration",current_sample,": The timing of demand changes does not converge","\n")
            print(R_hat)
            next
          }else{
            num_stable=T
            num_temp=as.integer(length(chain_s_store[1,])/num_period)#Number of iterations used for convergence inference (or after the burnin period)
          }
        }else{
          num_stable=T
          num_temp=as.integer(length(chain_s_store[1,])/num_period)#Number of iterations used for convergence inference (or after the burnin period)
        }
        #2. state-specific beta, (final_num_change+1)*num_covariate in total
        meantemp=matrix(0,num_chain,(final_num_change+1)*num_covariate)
        W=rep(0,(final_num_change+1)*num_covariate)
        s2temp=matrix(0,num_chain,(final_num_change+1)*num_covariate)
        R_hat=rep(0,(final_num_change+1)*num_covariate)
        for(chain in 1:num_chain){
          #state
          statetemp=matrix(chain_state_store[chain,],num_temp,num_state)
          statetemp=statetemp[num_temp,c(1:(final_num_change+1))]
          #beta
          temp=matrix(chain_beta_store[chain,],num_temp,num_state*num_covariate)
          indextemp=c()
          for(j in statetemp){
            indextemp=c(indextemp,c(1:num_covariate)+(j-1)*num_covariate)
          }
          temp=temp[,indextemp]
          for(j in 1:((final_num_change+1)*num_covariate)){
            meantemp[chain,j]=mean(temp[,j])
            s2temp[chain,j]=var(temp[,j])
            W[j]=W[j]+s2temp[chain,j]
          }
        }
        W=W/num_chain
        for(j in 1:((final_num_change+1)*num_covariate)){
          Bn=var(meantemp[,j])
          sigma2_plus=(num_temp-1)/num_temp*W[j]+Bn
          V_hat2=sigma2_plus+Bn/num_chain
          var_V_hat=((num_temp-1)/num_temp)^2/num_chain*var(s2temp[,j])+
            ((num_chain+1)/num_chain)^2*2/(num_chain-1)*Bn^2+
            2*(num_chain+1)*(num_temp-1)/(num_chain^2)/num_temp*
            (cov(s2temp[,j],meantemp[,j]^2)-2*mean(meantemp[,j])*cov(s2temp[,j],meantemp[,j]))
          df=2*V_hat2/var_V_hat
          #R_hat[j]=V_hat2/W[j]*df/(df-2)
          R_hat[j]=V_hat2/W[j]*(df+3)/(df+1)
          if(V_hat2==0||W[j]==0||is.na(R_hat[j])){
            R_hat[j]=1
          }
        }
        R_hat=sqrt(R_hat)
        if(max(R_hat)>R_thres){
          cat("Iteration",current_sample,": beta does not converge","\n")
          print(R_hat)
          next
        }
        #3. state-specific sigma, final_num_change+1 in total
        meantemp=matrix(0,num_chain,final_num_change+1)
        W=rep(0,final_num_change+1)
        s2temp=matrix(0,num_chain,final_num_change+1)
        R_hat=rep(0,final_num_change+1)
        for(chain in 1:num_chain){
          #state
          statetemp=matrix(chain_state_store[chain,],num_temp,num_state)
          statetemp=statetemp[num_temp,c(1:(final_num_change+1))]
          #sigma
          temp=matrix(chain_sigma2_store[chain,],num_temp,num_state)
          temp=as.matrix(temp[,statetemp])
          temp=sqrt(temp)
          for(j in 1:(final_num_change+1)){
            meantemp[chain,j]=mean(temp[,j])
            s2temp[chain,j]=var(temp[,j])
            W[j]=W[j]+s2temp[chain,j]
          }
        }
        W=W/num_chain
        for(j in 1:(final_num_change+1)){
          Bn=var(meantemp[,j])
          sigma2_plus=(num_temp-1)/num_temp*W[j]+Bn
          V_hat2=sigma2_plus+Bn/num_chain
          var_V_hat=((num_temp-1)/num_temp)^2/num_chain*var(s2temp[,j])+
            ((num_chain+1)/num_chain)^2*2/(num_chain-1)*Bn^2+
            2*(num_chain+1)*(num_temp-1)/(num_chain^2)/num_temp*
            (cov(s2temp[,j],meantemp[,j]^2)-2*mean(meantemp[,j])*cov(s2temp[,j],meantemp[,j]))
          df=2*V_hat2/var_V_hat
          #R_hat[j]=V_hat2/W[j]*df/(df-2)
          R_hat[j]=V_hat2/W[j]*(df+3)/(df+1)
          if(V_hat2==0||W[j]==0||is.na(R_hat[j])){
            R_hat[j]=1
          }
        }
        R_hat=sqrt(R_hat)
        if(max(R_hat)>R_thres){
          cat("Iteration",current_sample,": sigma does not converge","\n")
          print(R_hat)
          next
        }
        #If all parameters converge, the algorithm converges
        final_num_sampling=current_sample
        is_converge=T
        cat("Iteration",current_sample,", the MCMC algorithm converges","\n")
        stop_criteria=F
      }
    }
    duration_time=(proc.time()-start_time)/3600
    cat("Duration:",duration_time[3],"\n")
    
    ###Store the MCMC draws of the last MCMC chain if the algorithm converges###
    if(is_converge){#The algorithm converges
      num_converge=num_converge+1
      #The MCMC draws of the last MCMC chain
      num_sampling_temp=num_temp
      s_store=matrix(chain_s_store[chain,],num_sampling_temp,num_period)
      state_store=matrix(chain_state_store[chain,],num_sampling_temp,num_state)
      loc_change_store=matrix(chain_loc_change_store[chain,],num_sampling_temp,num_state)
      beta_store=matrix(chain_beta_store[chain,],num_sampling_temp,num_state*num_covariate)
      sigma2_store=matrix(chain_sigma2_store[chain,],num_sampling_temp,num_state)
      final_num_change=length(unique(na.omit(s_store[num_sampling_temp,])))-1
      loc_change_store=as.matrix(loc_change_store[,c(1:final_num_change)])
      #Store the information of demand changes
      storedata[regionindex,"num_change"]=final_num_change
      storedata[regionindex,"compute_time"]=duration_time[3]
      locstoredata[regionindex,]=0
      storedata[regionindex,"loc_change"]=""
      if(final_num_change>0){
        #The posterior modes of the timing of demand changes
        changeloc=c()
        for(i in 1:final_num_change){
          temp=loc_change_store[,i]
          temp=table(temp)
          temp=sort(temp,decreasing=T)
          changeloc=c(changeloc,as.integer(names(temp)[1]))
        }
        locstoredata[regionindex,changeloc]=1
        storedata[regionindex,"loc_change"]=paste(as.character(changeloc),sep="",collapse=",")
      }
      #Save
      setwd(osname)
      write.csv(storedata,"summary_covariate.csv")
      write.csv(locstoredata,"loc_change_covariate.csv")
      write.csv(locstoredata0,"loc_change0_covariate.csv")
      #Store the MCMC draws of the last MCMC chain
      setwd(paste(osname,"parastore_covariate",sep=""))#The absolute path that stores the model parameters
      loc_changedata=data.frame(loc_change_store[,c(1:final_num_change)])
      write.csv(loc_changedata,paste("loc_change",as.character(chosen_station),".csv",sep="",collapse=""))
      s_sampledata=data.frame(s_store)
      write.csv(s_sampledata,paste("s_store",as.character(chosen_station),".csv",sep="",collapse=""))
      beta_storedata=data.frame(beta_store)
      statetemp=state_store[num_sampling_temp,c(1:(final_num_change+1))]
      indextemp=c()
      for(j in statetemp){
        indextemp=c(indextemp,c(1:num_covariate)+(j-1)*num_covariate)
      }
      beta_storedata=beta_storedata[,indextemp]
      write.csv(beta_storedata,paste("beta",as.character(chosen_station),".csv",sep="",collapse=""))
      sigma2_storedata=data.frame(sigma2_store)
      sigma2_storedata=sigma2_storedata[,statetemp]
      write.csv(sigma2_storedata,paste("sigma2",as.character(chosen_station),".csv",sep="",collapse=""))
      
      ###Visualization###
      setwd(paste(osname,"figurestore_covariate",sep=""))#The absolute path that stores the figures
      #1. Plot: the timing of demand changes, the sequence of true bike sharing demands (or y), and the sequence of bike sharing demands predicted by the beta of each state
      png(paste("MCMC",as.character(chosen_station),".png",sep="",collapse=""),height=500,width=1500,res=100)
      par(mfrow=c(1,2))
      #The timing of demand changes
      ts.plot(s_sample,main=paste(as.character(i),"iterations",sep=" "))
      print(paste(as.character(i),"iterations",sep=" "))
      #the sequence of true bike sharing demands and the bike sharing demands predicted by the beta in the last MCMC iteration
      strtemp=c()
      plot(c(1:num_period),y_sequence,main="True y and predicted y",type="l",col=1,xlim=c(1,num_period+60),xlab="period")
      for(j in 1:length(statetemp)){
        lines(c(1:num_period),mu_matrix[,statetemp[j]],col=j+1)
        strtemp=c(strtemp,paste("Regime",as.character((statetemp[j]))))
      }
      lines(c(1:num_period),y_sequence,col=1)
      legend("topright",legend=c("True y",strtemp),col=c(1:(length(statetemp)+1)),lty=1,bg="white")
      dev.off()
      #2. Plot: the timing of demand changes, the sequence of true bike sharing demands (or y), and the sequence of independent variables
      png(paste("TS",as.character(chosen_station),".png",sep="",collapse=""),height=500,width=1500,res=100)
      par(mfrow=c(1,2))
      ts.plot(y_sequence)
      for(j in changeloc){
        abline(v=j,col="blue")
      }
      ymax=c()
      for(i in 1:(num_covariate-1)){
        temp=na.omit(x_matrix[,i]/y_sequence)
        temp=temp[temp!=Inf]
        ymax=c(ymax,quantile(temp,0.95))
      }
      ymax=max(ymax)
      ts.plot(x_matrix[,1]/y_sequence,col=1,ylim=c(0,ymax),xlim=c(1,350))
      if(num_covariate>1){
        for(i in 2:(num_covariate-1)){
          lines(x_matrix[,i]/y_sequence,col=i)
        }
      }
      for(j in changeloc){
        abline(v=j,col="blue")
      }
      legend("topright",legend=c(as.character(chosen_station),as.character(weatherstr)),col=c(1:(1+num_covariate)),lty=1,bg="white")
      dev.off()
    }
  }
}
print("finish!!")
