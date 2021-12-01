rewire=function(
ori=orig_net, # original network
novel=net,  # netoork with lost interactions
lev="high", #"high": rewiring for hummingbirds; "low" revireing for plants
flex=0.5){  #fraction of interactions to rewire

if(lev=="high"){
  aa=rowSums(ori)   #
  bb=rowSums(novel) #
  
  dd=data.frame(sp=names(aa),int_ori=aa)
  pos=match(names(aa),names(bb))
  dd$int_ex=bb[pos]
  dd$int_lost=dd$int_ori-dd$int_ex
  dd$int_rew=round(dd$int_lost*flex) # rounding to whole interaction
  
  #identifying the rewireing species
  rew_sp=dd$sp[which(dd$int_ex>0&dd$int_rew>0)]
  
  #repeat for each rewireing species
  if(sum(dd$int_rew)>0&length(rew_sp)>0){
  for(q in 1:length(rew_sp)){
    ssp=rew_sp[q]
    sub=novel[which(rownames(novel)==ssp),] #
    sub_sp=colnames(sub)[which(sub>0)] #
    ii=dd$int_rew[which(dd$sp==ssp)]
    
    #create rewireing sample
    tab=table(sample(sub_sp,ii,replace=T))
    for(u in 1:length(tab)){
      pos_c=which(colnames(novel)==names(tab)[u]) #
      pos_r=which(rownames(novel)==ssp)           #
      novel[pos_r,pos_c]=novel[pos_r,pos_c]+tab[u]
      
    } # end u
    
  } # end q
  
  
}} # end if

  
  if(lev=="low"){
    aa=colSums(ori)   #
    bb=colSums(novel) #
    
    dd=data.frame(sp=names(aa),int_ori=aa)
    pos=match(names(aa),names(bb))
    dd$int_ex=bb[pos]
    dd$int_lost=dd$int_ori-dd$int_ex
    dd$int_rew=round(dd$int_lost*flex) # rounding to whole interaction
    
    #identifying the rewireing species
    rew_sp=dd$sp[which(dd$int_ex>0&dd$int_rew>0)]
    
    #repeat for each rewireing species
    if(sum(dd$int_rew)>0&length(rew_sp)>0){
      for(q in 1:length(rew_sp)){
        ssp=rew_sp[q]
        sub=novel[,which(colnames(novel)==ssp)] #
        names(sub)=rownames(novel)
        sub_sp=names(sub)[which(sub>0)] #
        ii=dd$int_rew[which(dd$sp==ssp)]
        
        #create rewireing sample
        tab=table(sample(sub_sp,ii,replace=T))
        for(u in 1:length(tab)){
          pos_r=which(rownames(novel)==names(tab)[u]) #
          pos_c=which(colnames(novel)==ssp)           #
          novel[pos_r,pos_c]=novel[pos_r,pos_c]+tab[u]
          
        } # end u
        
      } # end q
      
  
  
}} # end both if
  
  
return(novel) 
}  
