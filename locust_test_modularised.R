maindir<-getwd()
source('microclim.R')
source('ectotherm.R')
survey_freq<-read.csv('survey_freq_98-09.csv')
locustdata<-read.csv('locustdata_90-09.csv')
locustdata$DATE_<-as.POSIXct(locustdata$DATE_,format="%Y-%m-%d")
ii<-5
for(ii in 1:50){
  
  longlat<-c(survey_freq[ii,2],survey_freq[ii,3]) # type a long/lat here in decimal degrees
  
  #get soil data
  source("../micro_australia/get.soil.R")
  soil.hydro<-get.soil()
  PE<-soil.hydro$PE
  BB<-soil.hydro$BB
  BD<-soil.hydro$BD
  KS<-soil.hydro$KS
  PE[1:9]<-CampNormTbl9_1[3,4] #air entry potential J/kg 
  KS[1:9]<-CampNormTbl9_1[3,6] #saturated conductivity, kg s/m3
  BB[1:9]<-CampNormTbl9_1[3,5] #soil 'b' parameter
  PE[10:13]<-CampNormTbl9_1[4,4] #air entry potential J/kg 
  KS[10:13]<-CampNormTbl9_1[4,6] #saturated conductivity, kg s/m3
  BB[10:13]<-CampNormTbl9_1[4,5] #soil 'b' parameter
  BulkDensity <- BD[seq(1,19,2)]*1000#rep(1360,10) # soil bulk density (kg/m3)
  
  # run microclimate model to get microclimate for ectotherm model and soil temps for predicting egg development and food availability
  micro<-micro_aust(loc = longlat, ystart = 1990, yfinish = 1990, PE = PE, BB = BB, BD = 
      BD, KS = KS, BulkDensity = BulkDensity)
  metout<-micro$metout
  shadmet<-micro$shadmet
  soil<-micro$soil
  shadsoil<-micro$shadsoil
  soilmoist<-micro$soilmoist
  shadmoist<-micro$shadmoist
  humid<-micro$humid
  shadhumid<-micro$shadhumid
  soilpot<-micro$soilpot
  shadpot<-micro$shadpot
  rainfall<-micro$RAINFALL
  MAXSHADES<-micro$MAXSHADES
  elev<-as.numeric(micro$ALTT)
  REFL<-as.numeric(micro$REFL)
  dim<-micro$dim
  
  RAINFALL<-rainfall
  
  # run ectotherm model to get body temperatures
  ecto<-ectotherm(lometry = 1, VTMAX = 43, VTMIN = 29, TBASK = 25, 
    TEMERGE = 15, ctmax = 47, ctmin = 1, TPREF = 37)
  
  environ<-as.data.frame(ecto$environ)
  tzone<-paste("Etc/GMT-",10,sep="") # doing it this way ignores daylight savings!
  dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
  dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
  
  environ<-cbind(dates,environ)
  soil<-cbind(dates,soil)
  metout<-cbind(dates,metout)
  shadsoil<-cbind(dates,shadsoil)
  shadmet<-cbind(dates,shadmet)
  
  
  dates2<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days") 
  dates2<-subset(dates2, format(dates2, "%m/%d")!= "02/29") # remove leap years
  rainfall<-as.data.frame(cbind(dates2,rainfall))
  colnames(rainfall)<-c("dates","rainfall")
}