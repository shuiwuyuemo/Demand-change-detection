
############################################################
# Clusters the data segments divided by the demand changes #
############################################################

library(ggplot2)
library(apcluster)

#Absolute paths
osname0="D:/change-point/robust_test_day/"#The absolute path stores input data
osname=paste(osname0,"regime_cluster/",sep="")#The absolute path stores clustering results
osname0="E:/research/bicycle-share/change_point/revision/回复1"
osname="E:/research/bicycle-share/change_point/regime_cluster"
covariate_choicetemp="adaptive"#The model. covariate: M_{basic}, covariate_trend: M_{trend}, covariate_piece: M_{piecewise}

###Calculate the state-specific regression prediction of bike sharing demands###
setwd(osname0)
data=read.csv("region_generation_final_daily.csv")
data=data[which(data$isweekend==0),]
data=data[which(data$isholiday==0),]
data$intercept=1
daydata=data$day
data=data[,-1]
num_period=length(data[,1])
num_region=60#Number of region
covariate=c("aqi","prec","fog","rain","snow","press","rh","temp","wind","day","intercept")
num_covariate=length(covariate)
datastore=data.frame(region=NA,regime=NA,hastrip=NA,startdate=NA,enddate=NA,num_day=NA,average_trip=NA)
for(i in covariate){
  datastore[,i]=NA
}
datastore=datastore[-1,]
count=1
for(region in 1:num_region){
  setwd(paste(osname0,"parastore_",covariate_choicetemp,"_unify",sep=""))
  betatemp=read.csv(paste("beta",as.character(region),".csv",sep=""))
  betatemp=betatemp[,-1]
  setwd(paste(osname0,"parastore_",covariate_choicetemp,sep=""))
  loctemp=read.csv(paste("loc_change",as.character(region),".csv",sep=""))
  loctemp=data.frame(loctemp[,-1])
  if(loctemp[1,1]!=num_period+1){#When there is no demand change, 
    startloc=c()
    endloc=c()
    for(i in 1:(length(loctemp[1,]))){
      startloc=c(startloc,max(loctemp[,i]))
      endloc=c(endloc,max(loctemp[,i]))
    }
    startloc=c(1,startloc)
    endloc=c(endloc,1501)
    #Input data
    temp=data[,region]
    datatemp=data.frame(data[,covariate])
    datatemp$y=temp
    num_sampling=length(betatemp[,1])
    covariate_estimation=as.numeric(matrix(1/num_sampling,1,num_sampling)%*%as.matrix(betatemp))
    covariate_estimation=matrix(covariate_estimation,num_covariate,length(startloc))
    for(i in 1:length(startloc)){
      #The time period of state i
      datatemptemp=datatemp[startloc[i]:(endloc[i]-1),]
      #Store the information of state i
      datastore[count,covariate]=as.numeric(covariate_estimation[,i])
      datastore[count,"region"]=region
      datastore[count,"regime"]=i
      datastore[count,"startdate"]=startloc[i]
      datastore[count,"enddate"]=endloc[i]-1
      datastore[count,"startdate_temp"]=startloc[i]
      datastore[count,"enddate_temp"]=endloc[i]-1
      datastore[count,"num_day"]=datastore[count,"enddate"]-datastore[count,"startdate"]+1
      datastore[count,"average_trip"]=mean(datatemptemp$y)
      if(mean(datatemptemp$y)<0.1){#When there is no bike sharing demand in this data segment
        datastore[count,"hastrip"]=0 
      }else{
        datastore[count,"hastrip"]=1
      }
      count=count+1
    }
  }else{
    setwd(paste(osname0,"parastore_",covariate_choicetemp,"_unify",sep=""))
    betatemp=read.csv(paste("beta",as.character(region),".csv",sep=""))
    betatemp=betatemp[,-1]
    #Input data
    temp=data[,region]
    datatemp=data.frame(data[,covariate])
    datatemp$y=temp
    num_sampling=length(betatemp[,1])
    #Store the information of state 1
    datastore[count,covariate]=as.numeric(matrix(1/num_sampling,1,num_sampling)%*%as.matrix(betatemp))
    datastore[count,"region"]=region
    datastore[count,"regime"]=1
    datastore[count,"startdate"]=1
    datastore[count,"enddate"]=length(daydata)
    datastore[count,"startdate_temp"]=1
    datastore[count,"enddate_temp"]=length(daydata)
    datastore[count,"num_day"]=datastore[count,"enddate"]-datastore[count,"startdate"]+1
    datastore[count,"average_trip"]=mean(datatemptemp$y)
    if(mean(datatemptemp$y)<0.1){#When there is no bike sharing demand in this data segment
      datastore[count,"hastrip"]=0 
    }else{
      datastore[count,"hastrip"]=1
    }
    count=count+1
  }
  print(region)
}
datastore[is.na(datastore)]=0
setwd(osname)
write.csv(datastore,paste("mu_col_",covariate_choicetemp,".csv",sep=""))

###Calculate the dissimilarity matrix (asymmetrical) based on Canberra distance###
#Dissimilarity function
Canberradist_dist=function(seq1,seq2){
  temp=which(!(seq1==0&seq2==0))
  if(length(temp)>0){
    seq1new=seq1[temp]
    seq2new=seq2[temp]
    L2dist=sum(abs(seq1new-seq2new)/(abs(seq1new)+abs(seq2new)))/length(seq1)
  }else{
    L2dist=0
  }
  return(L2dist)
}
#Input data
setwd(osname0)
data=read.csv("region_generation_final_daily.csv")
data=data[which(data$isweekend==0),]
data=data[which(data$isholiday==0),]
data=data[,-1]
data$intercept=1
setwd(osname)
datastore=read.csv(paste("mu_col_",covariate_choicetemp,".csv",sep=""))
datastore=datastore[,-1]
covariate=c("aqi","prec","fog","rain","snow","press","rh","temp","wind","day","intercept")
num_period=1501#T
num_covariate=length(covariate)
num_regime=length(datastore[,1])
distmatrix=matrix(NA,num_regime,num_regime) 
for(i in 1:num_regime){
  regiontemp=datastore[i,"region"]
  coef_temp=matrix(as.numeric(datastore[i,covariate]),num_covariate,1)
  count=1
  for(region in 1:60){
    loctemp=datastore[which(datastore$region==region),]
    #Input data
    y_true=data[,region]
    x_temp=data.frame(data[,covariate])
    y_pred=as.numeric(as.matrix(x_temp)%*%coef_temp)
    if(length(loctemp[,1])>1){#Have demand changes
      for(j in 1:length(loctemp[,1])){
        #The true bike sharing demands and the predicted bike sharing demands
        y_truetemp=y_true[loctemp[j,"startdate_temp"]:(loctemp[j,"enddate_temp"])]
        y_predtemp=y_pred[loctemp[j,"startdate_temp"]:(loctemp[j,"enddate_temp"])]
        y_predtemp=round(y_predtemp,6)
        distmatrix[i,count]=Canberradist_dist(y_truetemp,y_predtemp)#Calculate the dissimilarity. Index i refers to the model of state i, and index j refers to the bike sharing demands in data segment j
        count=count+1
      }
    }else{
      #The true bike sharing demands and the predicted bike sharing demands
      y_truetemp=y_true
      y_predtemp=y_pred
      y_predtemp=round(y_predtemp,6)
      distmatrix[i,count]=Canberradist_dist(y_truetemp,y_predtemp)#Calculate the dissimilarity. Index i refers to the model of state i, and index j refers to the bike sharing demands in data segment j
      count=count+1
    }
  }
  print(i)
}
setwd(osname)
#Rename the columns of the data frame
for(i in 1:length(datastore[,1])){
  datastore[i,"colname"]=paste(as.character(datastore[i,"region"]),"_",as.character(datastore[i,"regime"]),sep="")
}
#Save
write.csv(datastore,paste("mu_col_new_",covariate_choicetemp,".csv",sep=""))
distmatrix=data.frame(distmatrix)
names(distmatrix)=datastore$colname
write.csv(distmatrix,paste("regime_dist_new_",covariate_choicetemp,".csv",sep=""))

###Calculate the symmetric similarity matrix and cluster the data segments###
setwd(osname)
data=read.csv(paste("regime_dist_new_",covariate_choicetemp,".csv",sep=""))
data=data[,-1]
print(mean(diag(as.matrix(data))))
data[is.na(data)]=0#Change the NA in the dissimilarity matrix into 0
coldata=read.csv(paste("mu_col_new_",covariate_choicetemp,".csv",sep=""))
coldata=coldata[,-1]
colname=as.character(coldata$colname)
names(data)=colname
#Derive the symmetric similarity matrix based on Equ (19)
for(i in 1:length(data[,1])){
  data[i,]=data[i,]-data[i,i]
}
for(i in 1:(length(data[,1])-1)){
  for(j in (i+1):length(data[,1])){
    temp=max(data[i,j],data[j,i])
    data[i,j]=temp
    data[j,i]=temp
  }
}
data=-as.matrix(data)
#Apply affinity propagation method to cluster the data segments
data[is.na(data)]=0
cat(quantile(data,0.1),quantile(data,0.25),"\n")
data=as.SparseSimilarityMatrix(data)
ap_clust=apcluster(s=data)
clusterlist=ap_clust@clusters
print(length(clusterlist))#Print the number of clusters
coldata$cluster=NA
for(i in 1:length(clusterlist)){
  indextemp=clusterlist[[i]]
  indextemp=as.numeric(indextemp)
  coldata[indextemp,"cluster"]=i
}
print(table(coldata$cluster))
setwd(osname)
write.csv(coldata,paste("mu_col_new_",covariate_choicetemp,".csv",sep=""))

###Summarize the statistics of the bike sharing demands in each cluster###
#Including: 
#1. The number of data segments; 
#2. the summation, median, minimum and maximum of the duration, daily demands and growth trends for the data segments in each cluster.
setwd(osname0)
tripdata=read.csv("region_generation_final_daily.csv")
tripdata=tripdata[which(tripdata$isweekend==0),]
tripdata=tripdata[which(tripdata$isholiday==0),]
daydata=tripdata$day
tripdata=tripdata[,-1]
data=coldata
clusterset=unique(data$cluster)
clusterset=sort(clusterset)
datastore=data.frame(cluster=clusterset)
for(cluster in clusterset){
  numdaytemp=data[which(data$cluster==cluster),"num_day"]
  trendtemp=data[which(data$cluster==cluster),"day"]
  #Collect the demand records of the data segments in a cluster
  triptemp=c()
  datatemp=data[which(data$cluster==cluster),]
  for(i in 1:length(datatemp[,1])){
    startdate=datatemp[i,"startdate"]
    enddate=datatemp[i,"enddate"]
    temp=tripdata[startdate:enddate,]
    region=datatemp[i,"region"]
    triptemp=c(triptemp,temp[,region])
  }
  #Calculate the statistics
  datastore[cluster,"num_regime"]=length(datatemp[,1])
  datastore[cluster,"numday_sum"]=sum(numdaytemp)
  datastore[cluster,"numday_q50"]=median(numdaytemp)
  datastore[cluster,"numday_min"]=min(numdaytemp)
  datastore[cluster,"numday_max"]=max(numdaytemp)
  datastore[cluster,"trip_q50"]=median(triptemp)
  datastore[cluster,"trip_min"]=min(triptemp)
  datastore[cluster,"trip_max"]=max(triptemp)
  datastore[cluster,"trend_q50"]=median(trendtemp)*365
  datastore[cluster,"trend_min"]=min(trendtemp)*365
  datastore[cluster,"trend_max"]=max(trendtemp)*365
}
#Sort the data by the order of median bike sharing demands, and find the demand changes in category 1
cluster_sort=datastore[order(datastore$trip_q50),"cluster"]
datastore=datastore[cluster_sort,]
#Sort the data by the order of median growth trends, and find the demand changes in category 2
datastore=datastore[order(datastore$cluster),]
cluster_sort=datastore[order(datastore$trend_q50,decreasing=T),"cluster"]
datastore=datastore[cluster_sort,]
setwd(osname)
write.csv(datastore,paste("regime_summary_",covariate_choicetemp,".csv",sep=""))

###Visualize the clusters to identify the regimes with zero demand or high growth trends
setwd(osname0)
data=read.csv("region_generation_final_daily.csv")
data=data[which(data$isweekend==0),]
data=data[which(data$isholiday==0),]
data=data[,-1]
names(data)[1:60]=as.character(c(1:60))
data=data[,c(as.character(c(1:60)),"day")]
setwd(osname)
#Sort the clusters by the order of median bike sharing demands
datastore=read.csv(paste("regime_summary_",covariate_choicetemp,".csv",sep=""))
datastore=datastore[order(datastore$cluster),]
cluster_sort=datastore[order(datastore$trip_q50),"cluster"]
#Sort the data segments in each cluster by the order of start dates and end dates
indexdata=read.csv(paste("mu_col_new_",covariate_choicetemp,".csv",sep=""))
indexdata$order=NA
count=1
for(cluster_id in cluster_sort){
  datatemp=indexdata[which(indexdata$cluster==cluster_id),]
  datatemp=datatemp[order(datatemp$startdate,datatemp$enddate),]
  for(i in 1:length(datatemp[,1])){
    indexdata[which(indexdata$region==datatemp[i,"region"]&indexdata$regime==datatemp[i,"regime"]),"order"]=count
    count=count+1
  }
}
#Arrange daily bike sharing demands
num_period=length(data[,1])
num_regime=length(indexdata[,1])
data$index=c(1:num_period)
datastoretemp=matrix(NA,num_period,num_regime)
for(i in 1:num_regime){
  indextemp=indexdata[i,"order"]
  regiontemp=indexdata[i,"region"]
  datatemp=data[indexdata[i,"startdate"]:indexdata[i,"enddate"],]
  datastoretemp[c(min(datatemp$index):max(datatemp$index)),indextemp]=datatemp[,regiontemp]
}
colname=names(data)
datastore=data.frame(label1=rep(1:num_period,num_regime),label2=rep(1:num_regime,each=num_period),value=as.vector(datastoretemp))
#The lines that divide the clusters
lentemp=c()
temp=0
num_cluster=length(unique(indexdata$cluster))
for(i in cluster_sort){
  temp=temp+length(indexdata[which(indexdata$cluster==i),1])
  lentemp=c(lentemp,temp)
}
linedata=data.frame(x=rep(c(0.5,num_period+0.5),num_cluster),y=rep(lentemp+0.5,each=2),group=as.factor(rep(c(1:num_cluster),each=2)))
lentemp=c(0.5,lentemp,num_regime+0.5)
lentemp2=c()
for(i in 1:num_cluster){
  lentemp2=c(lentemp2,(lentemp[i]+lentemp[i+1])/2)
}
lentemp2=lentemp2+0.5
#Plot
pdf(paste("segment_cluster_",covariate_choicetemp,".pdf",sep=""),width=7.2,height=8)
ggplot()+
  geom_raster(data=datastore,aes(x=label1,y=label2,fill=value))+
  geom_line(data=linedata,col="black",linetype="dashed",aes(x=x,y=y,group=group))+
  geom_hline(yintercept=num_regime+0.5)+
  geom_vline(xintercept=num_period+0.5)+
  scale_fill_gradientn(colours=c("#F2E3E3","#A13030","grey20","#616130","#949449","#AFAF61"),
                       breaks=c(0,1000,2000,3000,4000,5000),
                       label=as.character(c(0,1000,2000,3000,4000,5000)),
                       na.value="white",limit=c(0,5000),space="Lab")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),
        legend.key=element_rect(fill="white"),legend.background=element_rect(fill="white",color="black"),
        axis.line=element_line(colour="black"),axis.text=element_text(colour="black"),#axis.text.y=element_blank(),
        axis.text.x=element_text(angle=5),
        legend.text=element_text(size=11),text=element_text(size=14),legend.title=element_text(size=12),
        axis.ticks.y=element_blank(),#axis.ticks.x=element_blank(),
        legend.key.size=unit(0.5,"cm"),legend.key.height=unit(0.5,"cm"),
        strip.background=element_rect(fill="white"),panel.spacing.y=unit(0.6,"cm"))+
  guides(fill=guide_colorbar(title="Daily\ndemands\n",title.hjust=0.5,title.vjust=0.5),
         color=guide_legend(title=element_blank()))+
  scale_x_continuous(expand=c(0,0),breaks=c(184,549,914,1280,1645,2010),
                     labels=c("2014-01-01","2015-01-01","2016-01-01","2017-01-01","2018-01-01","2019-01-01"))+
  scale_y_continuous(expand=c(0,0),breaks=lentemp2,labels=as.character(c(1:length(unique(indexdata$cluster)))))+
  xlab("Date")+ylab("The cluster of data segments")
dev.off()
#Manually examine whether the demand changes in categories 1 and 2 can be identified
#Find the clusters that correspond to categories 1 and 2
cat("category 1: cluster",cluster_sort[1],"\n")
cat("category 2: cluster",cluster_sort[c(7,8,11)],"\n")

###Classify the demand changes into 3 categories. Collect the timing of demand changes in category 3. This is used for investigating their temporal similarity to events in our paper########
setwd(paste(osname0,"/loc_change_data",sep=""))
locdata=read.csv(paste("loc_change_",covariate_choicetemp,".csv",sep=""))
locdata=locdata[,-1]
locdata2=locdata
locdata3=locdata
setwd(osname)
#The demand changes in category 1
data=read.csv(paste("mu_col_new_",covariate_choicetemp,".csv",sep=""))
num_regime=length(data[,1])
zero_cluster=13#Prespecify the index of the cluster in category 1 based on regime_summary.csv
data$zero_cluster=0
for(i in 1:length(data[,1])){
  if(data[i,"cluster"]==zero_cluster){
    data[i,"zero_cluster"]=1#The demand changes transits from the data segments with zero demands
    if(i>1&&data[i-1,"region"]==data[i,"region"]){#The demand changes transits to the data segments with zero demands
      data[i-1,"zero_cluster"]=1
    }
  }
}
data=data[which(data$zero_cluster==1),]
for(i in 1:length(data[,1])){
  region=data[i,"region"]
  regime=data[i,"regime"]
  loctemp=locdata[region,]
  loctemp=data.frame(loctemp)
  loctemp=names(loctemp)[which(loctemp[1,]==1)]
  if(regime<=length(loctemp)){
    locdata2[region,loctemp[regime]]=0
  }
}
sum1=sum(locdata)
sum2=sum(locdata2)
num_category1=sum1-sum2
cat("Number of demand changes in category 1:",num_category1,"\n")
#The demand changes in category 2
data=read.csv(paste("mu_col_new_",covariate_choicetemp,".csv",sep=""))
rapidgrowth_cluster=c(8,7,1)
data=data[which(data$cluster%in%rapidgrowth_cluster),] 
for(i in 1:length(data[,1])){
  region=data[i,"region"]
  regime=data[i,"regime"]
  loctemp=locdata[region,]
  loctemp=data.frame(loctemp)
  loctemp=names(loctemp)[which(loctemp[1,]==1)]
  if(regime<=length(loctemp)){
    locdata2[region,loctemp[regime]]=0
    locdata3[region,loctemp[regime]]=0
  }
}
num_region=60
num_category3=sum(locdata2)
num_category2=num_regime-num_category1-num_category3-num_region
cat("Number of demand changes in category 2:",num_category2,"\n")
cat("Number of demand changes in category 3:",num_category3,"\n")
setwd(paste(osname0,"/loc_change_data",sep=""))
write.csv(locdata2,paste("loc_change_",covariate_choicetemp,"_category3.csv",sep=""))#Category 3
write.csv(locdata3,paste("loc_change_",covariate_choicetemp,"_normal.csv",sep=""))#Category 1 and category 3
