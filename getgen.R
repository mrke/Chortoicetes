getgen<-function(hatchdates, stagefreq1){
  n<-0 
  p<-0
  starve<-24 # hrs locusts can go without food
  for(m in 1:length(hatchdates)){ #loop through hatch dates and test for successful reproduction
    row.names(DATA) <- NULL 
    masses  <- rep(0,nrow(DATA))
    n<-n+1
    ############## START HERE ################# 
    nograss<-0
    start<-hatchdates[m]
    drymass<-exp(pars_grow[1]) # fitted initial mass of 1 mg (egg mass)
    for(i in (start+1):length(masses)){
      if(i==start+1){
        masses[start]<-exp(pars_grow[1])
      }
      ddrymass<-dmass(pars_grow,drymass, environ[i,6])
      masses[i]<-masses[i-1]+ddrymass/24.
      instar<-getinstarfrommass(masses[i])
      # increment frequency of stage in stagefreq data frame
      eval(parse(text=paste("stagefreq1[",i,",",instar,"col]<-stagefreq1[",i,",",instar,"col]+1", sep="")))
      drymass<-masses[i]
      if(nograss>0){ # reset 'nograss' if there's been some growth
        if(plant.pres[i]==1){
          nograss<-0
        }
      }
      if(plant.pres[i]==0){
        nograss<-nograss+1
      }
      if(nograss<starve){
        if(masses[i]>55){ # if development complete, record hatch date 
          reprodate<-i
          p<-p+1
          if(p==1){
            reprodates<-c(reprodate,reprodate+7*24,reprodate+14*24) # three pods
            #reprodates<-c(reprodate) # one pod
          }else{
            reprodates<-c(reprodates,c(reprodate,reprodate+7*24,reprodate+14*24)) # three pods
            #reprodates<-c(reprodates,reprodate) # one pod
          }
          break
        }
      }else{
        stagefreq1[i,Deathcol]<- stagefreq1[i,Deathcol] + 1
        break
      }
    } # end loop through hatch dates 
    
    devs1<-as.data.frame(cbind(as.data.frame(environ$dates),masses))
    colnames(devs1)<-c('dates','mass')
    devs<-subset(devs1,mass>0)
    if(n==1){
      gens<-devs
    }else{
      gens<-rbind(gens,devs)
    }
    
  }
  if(exists("reprodates")==FALSE){
    reprodates<-NULL
  }
  return(list(reprodates=reprodates,gens=gens, stagefreq = stagefreq1))
}
