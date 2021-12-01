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

col_model=function(net_nam,n_reps=500,buff=10){
   if(buff==10){load("source_pool_list_10km.RData")
    source_pool_list=source_pool_list_10km}
  if(buff==100){load("source_pool_list_100km.RData")
    source_pool_list=source_pool_list_100km
    }
   
  path=paste0(getwd(),"/nets/",net_nam)
  
  res_col_total=data.frame(net=net_nam,ex=NA,co_ex=NA,co_exr=NA,col_net=NA,col_reg=NA)
  
  #loading network
  Network = read.csv(path,row.names = 1,h=T,sep=";")
  colnames(Network)=gsub("\\."," ", colnames(Network))
  Network=Network[,which(colnames(Network) %in% names(shp))] # removing species that are not in the polygons
  Network=Network[,complete.cases(t(Network))]
  Network=Network[rowSums(Network)>0,]
  
  
  
  #substitute names with 
  sub=subset(sp_nam,sp_nam$file==net_nam)
  sub=subset(sub,sub$name %in% names(shp))
  
  
  #community clim
  com_cor=dat[which(dat$file ==net_nam),c("Longitude","Latitude")]
  
  
  #removing migrants from temperate North America
  if(paste(net_nam) %in% cas){
    if(T %in% (as.character(sub$name) %in% migra)){
      ppr=which(sub$name %in% migra)
      Network=Network[,-ppr]
      sub=sub[-ppr,]
      
      ppr1=which(rowSums( Network)==0)
      if(length(ppr1>0)){ Network=Network[-ppr1,]}
      
    }}
  
  
  # determining source pool
  ft=which(names(source_pool_list)==net_nam)
  sp=source_pool_list[[ft]]
  
  #determining species-level colonization rate
  col_rate=rep(NA,length(sp))
  for(y in 1:length(sp)){
    foc_sp=as.character(sp[y])
    sp_shp=names(shp)[which(names(shp)==foc_sp)]
    
    pos=which(names(sp_clim)==sp_shp)
    cont=sp_clim[[pos]]

    
    if(nrow(cont)>10){
      cc_PCA=cont
      
      obs_PCA1_fut=(extract(fut_ras[[1]],com_cor)-mean(cc_PCA[,1],na.rm=T))/sd(cc_PCA[,1],na.rm=T)
      obs_PCA2_fut=(extract(fut_ras[[2]],com_cor)-mean(cc_PCA[,2],na.rm=T))/sd(cc_PCA[,2],na.rm=T)
      
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
  
  res_col_total$col_net[1]=mean(res_col_net)
  res_col_total$col_reg[1]=mean(res_col_reg)
  
  
  return(res_col_total)
  #print(i)
  
}
