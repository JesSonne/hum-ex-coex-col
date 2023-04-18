library(gdata)
library(reshape2)
library(vegan)
library(MASS)
library(sf)
library(rgeos)

z_score=function(x,pop){abs(x-mean(pop,na.rm=T))/sd(pop,na.rm=T)}

eu_dist=function(x1,x2,y1,y2){
  aa=abs(x1-x2)
  bb=abs(y1-y2)
  cc=sqrt(aa^2+bb^2)
  return(cc)
  
}

z_to_prob=function(z){pnorm(z,lower.tail = T)-pnorm(z,lower.tail = F)}

range01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
col_model=function(net_nam,n_reps=500,network_folder="...network_folder",buff=10){
  if(buff==10){load("data/source_pool_list_10km.RData")
    source_pool_list=source_pool_list_10km}
  if(buff==100){load("data/source_pool_list_100km.RData")
    source_pool_list=source_pool_list_100km
  }
  
  net_id=as.numeric(strsplit(net_nam,"_")[[1]][[1]])
  
  res_ex_coex=data.frame(net=net_nam,ex=NA,co_ex=NA,co_exr=NA,col=NA)
  
  #loading network
  Network = read.csv(paste0(network_folder,"/",net_nam),row.names = 1,h=T,sep=";")
  colnames(Network)=gsub("\\."," ", colnames(Network))
  Network=Network[,complete.cases(t(Network))]
  Network=Network[rowSums(Network)>0,]
  
  #hummingbird species in the network
  hum_sp=colnames(Network)
  
  #community clim
  com_cor=dat[which(dat$Network.ID.in.database ==net_id),c("Longitude","Latitude")]
  
  #removing temperate migrants from Central American networks
  if(paste(net_nam) %in% C_Am_nets){
    if(T %in% (as.character(hum_sp) %in% migra)){
      ppr=which(hum_sp %in% migra)
      Network=Network[,-ppr]
      hum_sp=hum_sp[-ppr]
      
      ppr1=which(rowSums( Network)==0)
      if(length(ppr1>0)){ Network=Network[-ppr1,]}
      
    }}
  
  
  # determining source pool
  sp_names=strsplit(names(source_pool_list),"_")
  sp_names=as.numeric(unlist(lapply(sp_names,"[[",1)))
  
  ft=which(sp_names==net_id)
  sp=source_pool_list[[ft]]
  
  #determining species-level colonization rate
  col_rate=rep(NA,length(sp))
  for(y in 1:length(sp)){
    foc_sp=as.character(sp[y])
    sp_shp=names(sp_clim)[which(names(sp_clim)==foc_sp)]
    
    pos=which(names(sp_clim)==sp_shp)
    cont=sp_clim[[pos]]
    
    
    if(nrow(cont)>10){
      cc_PCA=cont
      
      
      fut_clim=data.frame(extract(fut_ras,com_cor))[,-1]
      colnames(fut_clim)=c("bio1","bio4","bio12","bio15")
      pred <- predict(pca, newdata=fut_clim)[,c(1:2)] # predicting future data using the contemporary climates' principal componants 
      
      
      obs_PCA1_fut=(pred[1]-mean(cc_PCA[,1]))/sd(cc_PCA[,1])
      obs_PCA2_fut=(pred[2]-mean(cc_PCA[,2]))/sd(cc_PCA[,2])
      
      cc_PCA1=decostand(cc_PCA[,1],"standardize")
      cc_PCA2=decostand(cc_PCA[,2],"standardize")
      
      fut_dist=sqrt(z_score(obs_PCA1_fut,cc_PCA1)^2+z_score(obs_PCA2_fut,cc_PCA2)^2)
      
      col_rate[y]=1-z_to_prob(fut_dist)
      
    }else {col_rate[y]=0}
    
  } # end y
  
  names(col_rate)=sp
  
  #run simulation 
  res_col_net=rep(NA,n_reps)
  res_col_reg=rep(NA,n_reps)
  for(f in 1:n_reps){
    
    col_sp=rep(0,length(sp))
    for(q in 1:length(sp)){
      col_sp[q]=sample(c(1,0),1,prob =c(col_rate[q],1-col_rate[q]))
      
    }
    res_col_net[f]=sum(col_sp)/ncol(Network)
    res_col_reg[f]=sum(col_sp)/length(col_sp)
    
  }
  
  res_col_total=data.frame(net=net_nam,col_net=NA,col_reg=NA)
  res_col_total$col_net[1]=mean(res_col_net)
  res_col_total$col_reg[1]=mean(res_col_reg)
  
  
  return(res_col_total)
  #print(i)
  
}
