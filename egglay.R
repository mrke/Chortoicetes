# if photo period (DATA%PHOTO) < 13 and decreasing during egg lay, then make diapause egg 
# (lay at 2.5cm and stop development at 45% under ideal conditions)
egglay <- function(photoperiod,previous_photoperiod, SHALtemp, DEEPtemp, SHALmoist, DEEPmoist){
  if(photoperiod<photo_thresh & photoperiod<=previous_photoperiod){
    #set diapause potential
    out.diapause_pot <-TRUE
    out.temp <- SHALtemp
    out.moist <- SHALmoist
    out.diapause_egg <-TRUE
  }else{
    out.diapause_pot <- FALSE
    out.temp <- DEEPtemp
    out.moist<- DEEPmoist
    out.diapause_egg <-FALSE
  }
  out.diapause_hours <- 0
  out.cold_hours     <- 0
  out.dev            <- 0
  out.dry_hours      <- 0
  out.in_diapause    <- FALSE
  
  list(moist = out.moist                  ,
       temp = out.temp                    ,
       diapause_pot = out.diapause_pot    ,
       in_diapause  = out.in_diapause     ,
       diapause_hours = out.diapause_hours,
       cold_hours   = out.cold_hours      ,
       dry_hours    = out.dry_hours       ,
       diapause_egg = out.diapause_egg    ,
       dev = out.dev)
}

# Temp dependent function for C. terminifera egg dev rate fitted from Gregg (1984) egg dev. rates 
devrate <- function(temp){
  if(temp<=32){
    y = -12334*(1/(temp+273)) + 38.027 
    return(exp(y)/24) # divide by 24 to get units in 1/h
  } else{
    return(devrate(32))
  }
}