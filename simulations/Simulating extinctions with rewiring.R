library(gdata)
library(reshape2)
library(vegan)
library(MASS)
library(sf)

z_score=function(x,pop){abs(x-mean(pop,na.rm=T))/sd(pop,na.rm=T)}

eu_dist=function(x1,x2,y1,y2){
  aa=abs(x1-x2)
  bb=abs(y1-y2)
  cc=sqrt(aa^2+bb^2)
  return(cc)
  
}

z_to_prob=function(z){pnorm(z,lower.tail = T)-pnorm(z,lower.tail = F)}

range01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}

coex_model_rewireing=function(net_nam,n_reps=500,network_folder="...network_folder",r_flex=0.5){
  
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
  if(paste(net_id) %in% C_Am_nets){
    if(T %in% (as.character(hum_sp) %in% migra)){
      ppr=which(hum_sp %in% migra)
      Network=Network[,-ppr]
      hum_sp=hum_sp[-ppr]
      
      ppr1=which(rowSums( Network)==0)
      if(length(ppr1>0)){ Network=Network[-ppr1,]}
      
    }}
  
  
  # calculating climate-driven extinction risk
  survivability=rep(NA,length(hum_sp))
  clim_impact=rep(NA,length(hum_sp))
  for(y in 1:length(hum_sp)){
    foc_sp=as.character(hum_sp[y])
    sp_shp=which(names(sp_clim)==foc_sp)
    
    cont=sp_clim[[sp_shp]]
    
    if(nrow(cont)>10){
      cc_PCA=cont
      
      obs_PCA1=(extract(clim_PCA1_reg,com_cor)[,2]-mean(cc_PCA[,1]))/sd(cc_PCA[,1])
      obs_PCA2=(extract(clim_PCA2_reg,com_cor)[,2]-mean(cc_PCA[,2]))/sd(cc_PCA[,2])
      
      
      fut_clim=data.frame(extract(fut_ras,com_cor))[,-1]
      colnames(fut_clim)=c("bio1","bio4","bio12","bio15")
      pred <- predict(pca, newdata=fut_clim)[,c(1:2)] # predicting future data using the contemporary climates' principal componants 
      
      
      obs_PCA1_fut=(pred[1]-mean(cc_PCA[,1]))/sd(cc_PCA[,1])
      obs_PCA2_fut=(pred[2]-mean(cc_PCA[,2]))/sd(cc_PCA[,2])
      
      cc_PCA1=decostand(cc_PCA[,1],"standardize")
      cc_PCA2=decostand(cc_PCA[,2],"standardize")
      
      
      obs_dist=sqrt(z_score(obs_PCA1,cc_PCA1)^2+z_score(obs_PCA2,cc_PCA2)^2)
      fut_dist=sqrt(z_score(obs_PCA1_fut,cc_PCA1)^2+z_score(obs_PCA2_fut,cc_PCA2)^2)
      
      clim_dist=(fut_dist-obs_dist)
      
      clim_impact[y]=z_to_prob(clim_dist)*z_to_prob(fut_dist)
      survivability[y]=z_to_prob(fut_dist)   
      
      if(fut_dist-obs_dist<0){clim_impact[y]=0;survivability[y]=1}
      
    }else {clim_impact[y]=0;survivability[y]=1}
    
    
    
    
  } # end y
  
  
  res_coex=rep(NA,n_reps)
  clim_res=rep(NA,n_reps)
  
  for(x in 1:n_reps){
    
    ex=rep(0,length(hum_sp))
    for(u in 1:length(hum_sp)){
      ex[u]=sample(c(1,0),size = 1, prob = c(clim_impact[u],1-clim_impact[u]))
      
    } # end u
    
    net=Network
    orig_net=Network
    
    clim_res[x]=sum(ex)/length(hum_sp)# proportion of climate-driven extinctions
    
    #clim_res[x]=sum(ex)
    
    net[,ex==1]=0
    net1=orig_net
    #establish new network with rewireing
    net=rewire(ori=net1,novel=net,lev="high",flex=r_flex)
    
    
    #remove extinct species
    while(sum(ex)>0){
      
      
      che=1-(rowSums(net)/rowSums(net1))
      che[is.na(che)]=0
      #probability of removing plants
      p_probs=1-(rowSums(net)/rowSums(orig_net))
      p_probs=che
      
      
      p_ex=rep(0,length(p_probs))
      for(u in 1:length(p_probs)){
        p_ex[u]=sample(c(1,0),size = 1, prob = c(p_probs[u],1-p_probs[u]))
      } # end u
      
      net1=net
      net[p_ex==1,]=0
      #rewire low level
      net=rewire(ori=net1,novel=net,lev="low",flex=r_flex)
      
      che=1-(colSums(net)/colSums(net1))
      che[is.na(che)]=0
      
      #probability of removing hummingbird
      h_probs=1-(colSums(net)/colSums(orig_net))
      h_probs=che
      
      ex_h=rep(0,length(h_probs))
      for(u in 1:length(h_probs)){
        ex_h[u]=sample(c(1,0),size = 1, prob = c(h_probs[u],1-h_probs[u]))
      } # end u
      
      
      net1=net
      net[,ex_h==1]=0
      # rewire again for high level
      net=rewire(ori=net1,novel=net,lev="high",flex=r_flex)
      
      e_b=rep(1,length(ex_h))
      e_a=rep(1,length(ex_h))
      
      e_b[which(colSums(net1)>0)]=0
      e_a[which(colSums(net)>0)]=0
      
      
      ex=rep(0,length(ex_h))
      ex[which(rowSums(cbind(e_b,e_a))==1)]=1
      
    }
    
    
    
    res_coex[x]=length(which(colSums(net)==0))/length(hum_sp)   # proportion of climate-driven extinctions + coextinctions
    #res_coex[x]=length(which(colSums(net)==0))
    
    
    
  } # end x
  
  
  
  output=list()
  output[[1]]=clim_res
  output[[2]]=res_coex
  
  return(output)
  #print(i)
  
}