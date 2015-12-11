gethatch<-
  function(ovirows, stagefreq1, DATA){
  n<-0 
  p<-0
  for(m in 1:nrow(ovirows)){ #loop through oviposition dates and test for successful development
    row.names(DATA) <- NULL 
    dev  <- rep(0,nrow(DATA))
    # set developmental thresholds
    # create vector for previous day photoperiod
    DATA$PHOTOPREV <- DATA$PHOTO
    DATA$PHOTOPREV[2:nrow(DATA)]<-DATA$PHOTO[1:(nrow(DATA)-1)] 
    
    egg <- egglay(DATA$PHOTO[ovirows[m,2]],DATA$PHOTOPREV[ovirows[m,2]], DATA$D5cm, DATA$D10cm, DATA$MOIST5cm, DATA$MOIST10cm)
    # make counter for total number of consecutive egg hatches
    egg_gen <- 1
    # make empty vector for development 
    dev  <- rep(0,nrow(DATA))
    # make empty vectore for diapause potential
    dp   <- rep(0,nrow(DATA))
    # make empty state variable
    state <- rep("empty",nrow(DATA))
    # make empty variable for diapause egg (TRUE\FALSE) 
    d_e <- rep(FALSE,nrow(DATA))
    if(ovirows[m,2]+recover*24<=length(plant.pres) & ((ovirows[m,1]==1 & sum(plant.pres[ovirows[m,2]:(ovirows[m,2]+recover*24)])>=24*recover) | (ovirows[m,1]==0  ))){ # grass present for at least time to first lay
      ############## START HERE ################# 
      n<-n+1
      # loop through each date and update egg development vector (dev)
      if(ovirows[m,1]==1){ # check if arriving for first time or if last clutch at end of grass growth
        start<-ovirows[m,2]+recover*24
      }else{
        start<-ovirows[m,2]
      }
      for(i in start:length(dev)){
        # if cold, increment cold hours 
        if(egg$temp[i]<cold_thresh){
          egg$cold_hours <- egg$cold_hours + 1  
        }
        # if dry, increment dry hours, if moist, reset dry hours
        if(egg$moist[i]<dry_thresh){
          egg$dry_hours <- egg$dry_hours + 1  
        }else{
          egg$dry_hours <- 0
        }  
        #if too long in cold conditions or in diapause, avert diapause potential 
        if(egg$cold_hours>cold_hours_thresh || egg$diapause_hours>diapause_hours_thresh){
          egg$diapause_pot <- FALSE
        }
        
        if(egg$diapause_pot){
          dp[i]<-1
        }  else{ 
          dp[i]<-0
        }
        
        # If egg has diapause potential and at 45% dev, then egg is in diapause
        if(dev[i-1]>0.45 && dev[i-1]<0.46 && egg$diapause_pot == TRUE){ 
          egg$diapause_hours  <- egg$diapause_hours + 1
          egg$in_diapause <- TRUE
        } else{ 
          egg$in_diapause <- FALSE
        }
        # if not quiescent (Q1:dry and 30% dev)
        if(egg$moist[i]<dry_thresh && dev[i-1]>0.30 && dev[i-1]<0.31){
          egg$Q1 <- TRUE
          egg$Q2 <- FALSE
          state[i]<- "Q1"
          stagefreq1[i,Q1col]<-stagefreq1[i, Q1col]+1
          # else is not Q2:not diapause, dry, 45% dev  
        }else if(!egg$in_diapause && egg$moist[i]<dry_thresh && dev[i-1]>0.45 && dev[i-1]<0.46){
          egg$Q1 <- FALSE 
          egg$Q2 <- TRUE
          state[i] <- "Q2"
          stagefreq1[i,Q2col]<-stagefreq1[i,Q2col]+1
        }else{
          egg$Q1 <- FALSE
          egg$Q2 <- FALSE
          if(egg$in_diapause){
          state[i]<- "diapause"
          stagefreq1[i,Diapausecol]<-stagefreq1[i,Diapausecol]+1
          }else{
           state[i]<- "non-diapause"
          }
        }
        # nor in diapause  (cold_hours <720, diapuase = 1, and 45% dev) 
        if(i>1){ # check that isn't first position in vector
        if(!(egg$Q1||egg$Q2||egg$in_diapause)){
          # increment development at a 0.08 increase rate per day at 32C or use temp correct factor
          egg$dev <- dev[i-1] + devrate(egg$temp[i]) 
          state[i]<-"dev"
          stagefreq1[i,Eggdevcol]<-stagefreq1[i,Eggdevcol]+1
          # else in diapause/quiescent, so do not increment 
        } else{  egg$dev <- dev[i-1]
        }  
        }
        
        # if development complete, record hatch date  
        if(egg$dev>0.99){
          egg_gen <- egg_gen + 1
          state[i]<-"fin"
          hatchdate<-i
          if(n==1){
            hatchdates<-hatchdate
          }else{
            hatchdates<-c(hatchdates,hatchdate)
          }
          break
          # if consecutive dry hours is more than 6 months, terminate egg and start again   
        }else if(egg$dry_hours>dry_hours_thresh){
          stagefreq1[i,Deathcol]<-stagefreq1[i,Deathcol]+1
          break
        }
        # update egg development vector and diapause egg vector
        dev[i]   <- egg$dev
        d_e[i]   <- egg$diapause_egg
      }
      
    } # end check if within recovery time
    if(ovirows[m,1]!=1){
      p<-p+1
    }
  } # end loop through ovip dates
  return(list(hatchdates=hatchdates, stagefreq=stagefreq1))
}
