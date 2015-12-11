maindir<-getwd()
library(NicheMapR)
library(plyr)
survey_freq<-read.csv('survey_freq_98-09.csv')
locustdata<-read.csv('locustdata_90-09.csv')
locustdata$DATE_<-as.POSIXct(locustdata$DATE_,format="%Y-%m-%d")
ystart<-1934
yfinish<-2014
nyears<-yfinish-ystart+1
ii<-19
for(ii in 1:50){
  
  longlat<-c(survey_freq[ii,2],survey_freq[ii,3]) # type a long/lat here in decimal degrees
  #longlat<-c(142,-12)
  loc<-longlat
  source("../micro_australia/get.soil.R")
  source("../micro_australia/micro_aust.R")
  DEP = c(0., 1.,  3.,  5,  10,  20.,  30.,  60.,  90.,  200.)  
  DEP = c(0., 2.5, 5, 7.5, 10,  20.,  30.,  60.,  90.,  200.)  
 DEP <- c(0., 2.5, 5,  10, 15, 20.,  30.,  60.,  90.,  200.) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
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
  micro<-micro_aust(loc = longlat, ystart = ystart, yfinish = yfinish, PE = PE, BB = BB, BD = 
      BD, KS = KS, BulkDensity = BulkDensity, dailywind = 0, DEP = DEP)
  metout<-as.data.frame(micro$metout)
  shadmet<-as.data.frame(micro$shadmet)
  soil<-as.data.frame(micro$soil)
  shadsoil<-as.data.frame(micro$shadsoil)
  soilmoist<-as.data.frame(micro$soilmoist)
  shadmoist<-as.data.frame(micro$shadmoist)
  humid<-as.data.frame(micro$humid)
  shadhumid<-as.data.frame(micro$shadhumid)
  soilpot<-as.data.frame(micro$soilpot)
  shadpot<-as.data.frame(micro$shadpot)
  rainfall<-as.data.frame(micro$RAINFALL)
  MAXSHADES<-as.data.frame(micro$MAXSHADES)
  elev<-as.numeric(micro$ALTT)
  REFL<-as.numeric(micro$REFL)
  dim<-micro$dim
  
  RAINFALL<-rainfall
  
  # run ectotherm model to get body temperatures
  source("c:/git/NicheMapR/R/ectotherm.R")
  ecto<-ectotherm(lometry = 1, VTMAX = 43, VTMIN = 29, TBASK = 25, 
    TEMERGE = 15, ctmax = 47, ctmin = 1, TPREF = 37)
  
  environ<-as.data.frame(ecto$environ)
  environ2<-subset(environ,TIME==12)
  photop<-spline(rep(seq(1,365),nyears)[1:(nrow(environ)/24)],environ2$DAYLENGTH,n=nrow(environ),xmin=1,xmax=365*nyears,method="periodic")$y
  environ$DAYLENGTH<-photop
  masbal<-as.data.frame(ecto$masbal)
  tzone<-paste("Etc/GMT-",10,sep="") # doing it this way ignores daylight savings!
  if(yfinish!=as.numeric(substr(Sys.time(),1,4))){
  dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
  }else{
  dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate(as.numeric(substr(Sys.time(),1,4)),as.numeric(substr(Sys.time(),6,7)),as.numeric(substr(Sys.time(),9,10)),tz=tzone)-3600*13, by="hours")
  }
  dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
  
  environ<-cbind(dates,environ)
  masbal<-cbind(dates,masbal)
  soil<-cbind(dates,soil)
  metout<-cbind(dates,metout)
  shadsoil<-cbind(dates,shadsoil)
  shadmet<-cbind(dates,shadmet)
  
  if(yfinish!=as.numeric(substr(Sys.time(),1,4))){
  dates2<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days") 
  }else{
  dates2<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate(as.numeric(substr(Sys.time(),1,4)),as.numeric(substr(Sys.time(),6,7)),as.numeric(substr(Sys.time(),9,10)),tz=tzone)-3600*13, by="days") 
  }
  dates2<-subset(dates2, format(dates2, "%m/%d")!= "02/29") # remove leap years
  rainfall<-as.data.frame(cbind(dates2,rainfall))
  colnames(rainfall)<-c("dates","rainfall")
  
  # get subset for current site 
  attach(locustdata)
  locustsub<-locustdata[(long == longlat[1])&(lat==longlat[2]),]
  detach(locustdata)
  locustsub$date<-as.character(locustsub$DATE_,format="%Y-%m-%d")
  
  root_shallow<-3#kk-1#2#5 # how shallow do the roots go? 2 to 10, corresopnding to 2.5, 5, 10, 15, 20, 30, 60, 90 and 200 cm
  root_deep<-6#kk#6 # how deep do the roots go? 2 to 10, corresopnding to 2.5, 5, 10, 15, 20, 30, 60, 90 and 200 cm
  growth_delay<-0# days after suitable soil moisture that new growth occurs
  wilting_thresh<-200*-1 # water potential for wilting point J/kg (divide by 100 to get bar)
  permanent_wilting_point<-1500*-1 # water potential for permanent wilting point (PWP) J/kg (divide by 100 to get bar)
  FoodWater<-82 # water content of fully hydrated veg
  
  source("plantgro.R")
  plantout<-plantgro(soilpot = soilpot, soilmoist = soilmoist, root_shallow = root_shallow, 
    root_deep = root_deep, growth_delay = growth_delay, wilting_thresh=wilting_thresh, permanent_wilting_point=permanent_wilting_point)

  plant.pres<-plantout$plant.pres
  plant.moist<-plantout$pct.water
  
#   plot(plant.pres,type='h')
#   
#   plantout2<-plantout
#   plantout2$date<-as.character(as.Date(plantout$date1,format="%Y-%m-%h"))  
# 
#   merge_results_e<-merge(plantout2,locustsub,by="date")
#   merge_results_e<-subset(merge_results_e,Enew!='NA')
#   pred<-aggregate(merge_results_e$moist.index,by=list(merge_results_e$date),FUN=max)
#   obser<-aggregate(merge_results_e$Enew,by=list(merge_results_e$date),FUN=max)
#   r_ephem<-round(cor(round(pred$x),obser$x,method="spearman"),4)
#   p_ephem<-round(as.numeric(cor.test(round(pred$x),obser$x)[3]),4)
#   plot(pred$x~jitter(obser$x),col='blue',ylab='predicted greenness',xlab='observed greenness',main=paste("ephemeral plants, roots ",DEP[root_deep]," cm",sep=""),xlim=c(0,11))
#   text(2,max(pred$x)-1,paste("r=",r_ephem,"p=",p_ephem))
#   
  ###################### run locust simulation ##################################
  environ_agg<-aggregate(environ[,2:23],by=list(format(environ$date,'%Y-%m-%d')),max)
  TCs<-cbind(rep(1:length(environ_agg$TC),24),rep(environ_agg$TC,24))
  TCs<-TCs[order(TCs[,1]),]
  oviposit<-plant.pres
  oviposit[TCs[,2]<29]<-0 # prevent oviposition occurring on days when too cold (max daily TC less than 29 degrees)
  recover<-1 # time locust needs to build up resources to lay first batch
  ovidates<-as.data.frame(cbind(oviposit[1:length(oviposit)],c(0,oviposit[1:(length(oviposit)-1)])))
  ovidates<-cbind(dates[1:length(dates)],ovidates)
  ovidates_start<-subset(ovidates, V1-V2==1)
  ovirows_start<-cbind(1,as.numeric(rownames(ovidates_start))) # rows at which grass growth started
  ovidates_finish<-subset(ovidates, V1-V2==-1)
  ovirows_finish<-cbind(0,as.numeric(rownames(ovidates_finish))) # rows at which grass growth started
  if(oviposit[length(oviposit)]==1){
    ovirows_finish<-rbind(ovirows_finish,cbind(0,length(oviposit)))
  }
  ovirows<-rbind(ovirows_start,ovirows_finish)
  relays<-ovirows_start[(ovirows_finish[,2]-ovirows_start[,2])/(7*24)>=2,]
  clutches<-(ovirows_finish[,2]-ovirows_start[,2])/(7*24)
  clutches<-floor(subset(clutches,clutches>=2))
  if(length(clutches)>0){
  relays2<-rep(0,sum(clutches))
  counter=0
  for(mm in 1:nrow(relays)){
    for(mn in 1:clutches[mm]){
      counter=counter+1
      relays2[counter]<-relays[mm,2]+7*24*mn
    }
  }
  ovirows_mid<-cbind(1,relays2)
  ovirows<-rbind(ovirows,ovirows_mid)
  }
  ovirows<-ovirows[order(-ovirows[,1],ovirows[,2]),]
  ########### chortoicetes oviposition model##########
  
  # to do, for now assume on first and last days of grass growth periods
  
  ############### chortoicetes egg model #############
  
  DATA1<-cbind(soil[,1:8],metout[,13:14],environ[,19],soilmoist[,4:6]*100)
  colnames(DATA1)<-c('DATE','JULDAY','TIME','D0cm','D2.5cm','D5cm', 'D10cm','D15cm','ZEN','SOLR','PHOTO','MOIST3cm','MOIST5cm','MOIST10cm')
  
  source("egglay.R")
  source("devrate.R")
  source("gethatch.R")
  source("dmass.R")
  source("dsurvival.R")
  source("getinstarfrommass.R")
  source("getgen.R")
  
  photo_thresh <- 13 # h, below this threshold, eggs are laid with diapause potential
  dry_thresh   <- 9  # %, below this threshold, the low soil moisture may trigger quiescence or desiccation
  dry_hours_thresh <- 365*24 # h,  above this threshold, egg dessicates 
  cold_thresh  <- 15 # C, below this threshold, total cold hours accumulate
  cold_hours_thresh <- 60*24 # h, above this threshold of cumulative cold hours, diapause potential is lost
  diapause_hours_thresh <- 7*7*24 # h,  above this threshold of cumulative diapause hours, diapause potential is lost 
  DATA<-DATA1
  pars_grow<-read.csv('growth_coeff.csv')[,2]
  pars_survive<-read.csv('survive_coeff.csv')[,2]
  
  # make matrix for frequency of different stages through time
  stagefreqnames<-c('Eggdev', 'Q1', 'Q2', 'Diapause', 'N1','N2','N3','N4','N5','Adult','Death')
  for(name in stagefreqnames){ # make variables to easily access matrix column by name
    eval(parse(text=paste(name,"col<-",which(stagefreqnames==name), sep = ""))) 
  }
  stagefreq<-rep(0,nrow(DATA1)*length(stagefreqnames))
  dim(stagefreq)<-c(nrow(DATA1),length(stagefreqnames))
  
  
  for(g in 1:15){ # cap at max of 15 generations
    if(g==1){
      hatchings<-gethatch(ovirows, stagefreq, DATA=DATA1)
      hatchdates<-hatchings$hatchdates 
      generations<-getgen(hatchdates, hatchings$stagefreq)
      reprodates<-generations$reprodates[generations$reprodates<nrow(generations$stagefreq)] # make sure not going beyond time series
      reprodates1<-reprodates
      gens<-generations$gens
    }else{
      if(length(reprodates)>0){
        hatchings<-gethatch(cbind(0,reprodates),generations$stagefreq, DATA=DATA1)
        hatchdates<-hatchings$hatchdates
        hatchdates<-hatchdates[!duplicated(hatchdates)]
        generations<-getgen(hatchdates, generations$stagefreq)
        reprodates<-generations$reprodates[generations$reprodates<nrow(generations$stagefreq)] # make sure not going beyond time series
        gens<-generations$gens
      }else{
        break
      }
    }
    if(length(reprodates)>0){
      if(g==1){
        allgens<-gens
      }else{
        if(sum(reprodates1)!=sum(reprodates)){
          allgens<-rbind(allgens,gens)
          reprodates1<-reprodates
        }else{
          break
        }
      }
      cat('gen ',g,'\n')
      plot(environ$TC~environ$dates,type='l',col='orange',ylim=c(0,60))
      points(plant.moist~environ$dates,type='l',col='green')
      points(allgens$mass~allgens$dates,type='p',col='blue',cex=0.1)
    }else{
      break
    }
  }
  
  success<-as.data.frame(subset(allgens,mass>55))
  
  stagefreqDF<-as.data.frame(generations$stagefreq)
  
  stagefreqDF<-cbind(dates,longlat[1],longlat[2],stagefreqDF)
  names(stagefreqDF)<-c('date','long','lat',stagefreqnames)
  
  stagefreqDF_agg<-aggregate(stagefreqDF[,2:14],by=list(format(dates,'%Y-%m-%d')),max)
  allgens_agg<-aggregate(allgens[,2],by=list(format(allgens$dates,'%Y-%m-%d')),max)
  colnames(stagefreqDF_agg)[1]<-'date'
  
  stagefreqDF_agg$nymphs<-rowSums(stagefreqDF_agg[,9:12]) # exclude 1st instars as they may have just hatched and died if no food
  egglaytimes<-plantout[ovirows[,2],]
  egglaytimes<-aggregate(egglaytimes[,2:5],by=list(format(egglaytimes$date1,'%Y-%m-%d')),max)
  plantout_agg<-aggregate(plantout[,2:5],by=list(format(plantout$date1,'%Y-%m-%d')),max)
  stagefreqDF_agg$SoilMoist<-plantout_agg$soilmoist
  stagefreqDF_agg$PlantWater<-plantout_agg$pct.water
  stagefreqDF_agg$PlantPres<-plantout_agg$plant.pres
  
  write.csv(stagefreqDF_agg,paste("results/stagefreq_site",ii,".csv",sep=""))
  write.table(cbind(ii,g),file = "results/gens.csv",append = TRUE, col.names = F, sep = ",")
  

  plotstagefreqDF_agg<-subset(stagefreqDF_agg,as.numeric(substr(stagefreqDF_agg$date,start = 1, stop = 4 ))>=1990)

  
  year_vals<-subset(environ,as.character(dates,format='%d/%m')=="01/01")
  year_vals<-subset(year_vals,as.character(year_vals$dates,format='%H')=="00") # get midnight
 
  
  plot(stagefreqDF_agg$Diapause~as.POSIXct(plotstagefreqDF_agg$date),ylab="",xlab="",type='l',ylim=c(0,24),col="NA")
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(egglaytimes$pct.water*0~as.POSIXct(egglaytimes$Group.1),type='p',col="black",pch=16)  
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year
  points(environ_agg$DAYLENGTH~dates2,type='l',col='blue')
  points(stagefreqDF_agg$Diapause~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  points(allgens$mass/10~allgens$dates,type='p',col='brown',cex=0.1)

  plot(stagefreqDF_agg$Diapause~as.POSIXct(plotstagefreqDF_agg$date),ylab="",xlab="",type='l',ylim=c(0,24),col="NA")
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(egglaytimes$pct.water*0~as.POSIXct(egglaytimes$Group.1),type='p',col="black",pch=16)  
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year
  points(stagefreqDF_agg$N1~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  points(environ_agg$DAYLENGTH~dates2,type='l',col='blue')

  plot(stagefreqDF_agg$Diapause~as.POSIXct(plotstagefreqDF_agg$date),ylab="",xlab="",type='l',ylim=c(0,24),col="NA")
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(egglaytimes$pct.water*0~as.POSIXct(egglaytimes$Group.1),type='p',col="black",pch=16)  
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year
  points(stagefreqDF_agg$N2~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  points(environ_agg$DAYLENGTH~dates2,type='l',col='blue')

  plot(stagefreqDF_agg$Diapause~as.POSIXct(plotstagefreqDF_agg$date),ylab="",xlab="",type='l',ylim=c(0,24),col="NA")
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(egglaytimes$pct.water*0~as.POSIXct(egglaytimes$Group.1),type='p',col="black",pch=16)  
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year
  points(stagefreqDF_agg$N3~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  points(environ_agg$DAYLENGTH~dates2,type='l',col='blue')
    
  plot(stagefreqDF_agg$Diapause~as.POSIXct(plotstagefreqDF_agg$date),ylab="",xlab="",type='l',ylim=c(0,24),col="NA")
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(egglaytimes$pct.water*0~as.POSIXct(egglaytimes$Group.1),type='p',col="black",pch=16)  
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year
  points(stagefreqDF_agg$N4~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  points(environ_agg$DAYLENGTH~dates2,type='l',col='blue')
  
  plot(stagefreqDF_agg$Diapause~as.POSIXct(plotstagefreqDF_agg$date),ylab="",xlab="",type='l',ylim=c(0,24),col="NA")
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(egglaytimes$pct.water*0~as.POSIXct(egglaytimes$Group.1),type='p',col="black",pch=16)  
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year
  points(stagefreqDF_agg$N5~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  points(environ_agg$DAYLENGTH~dates2,type='l',col='blue')
  
      plot(environ$TC~environ$dates,type='l',col='grey',ylab="",xlab="",ylim=c(0,50))
      abline(29,0,col='blue',lty=2)
      abline(43,0,col='red',lty=2)
      abline(37,0,col='orange',lty=2,lwd=2)
      points(environ$SHADE/10~environ$dates,type='l',col='dark green')
      points(environ$DEP/10~environ$dates,type='l',col='brown')
      
      points(allgens$mass~allgens$dates,type='p',col='blue',cex=0.1)
  
  plot(stagefreqDF_agg$Q1~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(stagefreqDF_agg$Q1~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year
   
   
  plot(stagefreqDF_agg$Q2~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(stagefreqDF_agg$Q2~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year
    
  filename<-paste("results/locust time series roots",DEP[root_deep]," cm regrow thresh ",growth_delay," days site ",ii,".pdf",sep="")
  pdf(filename,paper="A4r",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(3,1)) # set up for 12 plots in 2 columns
  par(oma = c(2,2,2,2) + 1) # margin spacing stuff
  par(mar = c(3,3,1,1) + 1) # margin spacing stuff 
  par(mgp = c(3,1,0) ) # margin spacing stuff   
  
  plotstagefreqDF_agg<-subset(stagefreqDF_agg,as.numeric(substr(stagefreqDF_agg$date,1,4))>1989)
  plot(plotstagefreqDF_agg$nymphs/max(plotstagefreqDF_agg$nymphs)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l',ylim=c(0,8),ylab="nymph density",xlab="")
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(locustsub$NDENS-0.5~locustsub$DATE_,type='h',col='orange',lwd=2)  
  points(plotstagefreqDF_agg$nymphs/max(plotstagefreqDF_agg$nymphs)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year

  
  plot(plotstagefreqDF_agg$Adult/max(plotstagefreqDF_agg$Adult)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l',ylim=c(0,8),ylab="adult density",xlab="")
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(locustsub$ADENS-0.5~locustsub$DATE_,type='h',col='orange',lwd=2)  
  points(plotstagefreqDF_agg$Adult/max(plotstagefreqDF_agg$Adult)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year
  
  plot((plotstagefreqDF_agg$Adult+plotstagefreqDF_agg$nymphs)/max(plotstagefreqDF_agg$Adult+plotstagefreqDF_agg$nymphs)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l',ylim=c(0,8),ylab="adult and nymph density",xlab="")
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(locustsub$ADENS + locustsub$NDENS-0.5~locustsub$DATE_,type='h',col='orange',lwd=2)  
  points((plotstagefreqDF_agg$Adult+plotstagefreqDF_agg$nymphs)/max(plotstagefreqDF_agg$Adult+plotstagefreqDF_agg$nymphs)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l')  
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year
  title(paste("lat/long ",longlat[2],",",longlat[1],sep=""))
  dev.off()
  
filename<-paste("results/locust time series 1934-2014 roots",DEP[root_deep]," cm regrow thresh ",growth_delay," days site ",ii,".pdf",sep="")
  pdf(filename,paper="A4r",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(3,1)) # set up for 12 plots in 2 columns
  par(oma = c(2,2,2,2) + 1) # margin spacing stuff
  par(mar = c(3,3,1,1) + 1) # margin spacing stuff 
  par(mgp = c(3,1,0) ) # margin spacing stuff   
  plotstagefreqDF_agg<-stagefreqDF_agg
  points(plotstagefreqDF_agg$nymphs/max(plotstagefreqDF_agg$nymphs)*5~as.POSIXct(plotstagefreqDF_agg$date),type='h',lwd=1,ylim=c(0,5),ylab="nymph density",xlab="",col='blue')
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*5*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  #points(locustsub$NDENS-0.5~locustsub$DATE_,type='h',col='orange',lwd=2)  
 # points(plotstagefreqDF_agg$nymphs/max(plotstagefreqDF_agg$nymphs)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  abline(v=year_vals$dates,col='grey',lty=2) # add lines to show beginning of each year
  
  obs_outbreaks<-read.csv("locust_outbreaks.csv")
  plot(plotstagefreqDF_agg$nymphs/max(plotstagefreqDF_agg$nymphs)*5~as.POSIXct(plotstagefreqDF_agg$date),type='l',ylim=c(0,5),ylab="nymph density",xlab="",col="NA")
  points(obs_outbreaks$outbreak~as.POSIXct(obs_outbreaks$date,format="%d/%m/%Y"),type='h',col='blue',lwd=2)
  points(plotstagefreqDF_agg$nymphs/max(plotstagefreqDF_agg$nymphs)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  
  obs<-aggregate(obs_outbreaks$outbreak,by=list(substr(as.Date(obs_outbreaks$date,format="%d/%m/%Y"),1,4)),max)
  pred<-aggregate(plotstagefreqDF_agg$nymphs,by=list(substr(as.Date(plotstagefreqDF_agg$date,format="%Y-%m-%d"),1,4)),max)
  plot(as.numeric(pred$Group.1),pred$x/max(pred$x)*5,type='h')
  points((as.numeric(obs$Group.1)+0.2),obs$x,type='h',col='red')
  obspred<-merge(obs,pred,by = "Group.1")
  #obspred$x.y=log10(obspred$x.y+1)
  obspred$x.y=obspred$x.y/max(obspred$x.y)*5
  obspred<-as.data.frame(obspred)
  plot(jitter(obspred$x.x)~obspred$x.y)
  cor(obspred$x.x,obspred$x.y)
  cor.test(obspred$x.x,obspred$x.y)

  plot(plotstagefreqDF_agg$Adult/max(plotstagefreqDF_agg$Adult)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l',ylim=c(0,8),ylab="adult density",xlab="")
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(locustsub$ADENS-0.5~locustsub$DATE_,type='h',col='orange',lwd=2)  
  points(plotstagefreqDF_agg$Adult/max(plotstagefreqDF_agg$Adult)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l')
  
  plot((plotstagefreqDF_agg$Adult+plotstagefreqDF_agg$nymphs)/max(plotstagefreqDF_agg$Adult+plotstagefreqDF_agg$nymphs)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l',ylim=c(0,8),ylab="adult and nymph density",xlab="")
  points(plantout_agg$pct.water/max(plantout_agg$pct.water)*6*plantout_agg$plant.pres~dates2,type='h',col="darkolivegreen2")  
  points(locustsub$ADENS + locustsub$NDENS-0.5~locustsub$DATE_,type='h',col='orange',lwd=2)  
  points((plotstagefreqDF_agg$Adult+plotstagefreqDF_agg$nymphs)/max(plotstagefreqDF_agg$Adult+plotstagefreqDF_agg$nymphs)*8~as.POSIXct(plotstagefreqDF_agg$date),type='l')  
  title(paste("lat/long ",longlat[2],",",longlat[1],sep=""))
  dev.off()

  filename<-paste("results/grass correl roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days site ",ii,".pdf",sep="")
  pdf(filename,paper="A4",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(2,1)) # set up for 12 plots in 2 columns
  par(oma = c(2,2,2,2) + 1) # margin spacing stuff
  par(mar = c(3,3,1,1) + 1) # margin spacing stuff 
  par(mgp = c(3,1,0) ) # margin spacing stuff 
  
  
  
  plantout2<-plantout
  plantout2$date<-as.character(as.Date(plantout$date1,format="%Y-%m-%h"))  

  merge_results_e<-merge(plantout2,locustsub,by="date")
  merge_results_e<-subset(merge_results_e,Enew!='NA')
  pred<-aggregate(merge_results_e$moist.index,by=list(merge_results_e$date),FUN=max)
  obser<-aggregate(merge_results_e$Enew,by=list(merge_results_e$date),FUN=max)
  r_ephem<-round(cor(round(pred$x),obser$x,method="spearman"),4)
  p_ephem<-round(as.numeric(cor.test(round(pred$x),obser$x)[3]),4)
  plot(pred$x~jitter(obser$x),col='blue',ylab='predicted greenness',xlab='observed greenness',main=paste("ephemeral plants, roots ",DEP[root_deep]," cm",sep=""),xlim=c(0,11))
  text(2,max(pred$x)-1,paste("r=",r_ephem,"p=",p_ephem))

  merge_results_p<-merge(plantout2,locustsub,by="date")
  merge_results_p<-subset(merge_results_p,PERENnew!='NA')
  pred<-aggregate(merge_results_p$moist.index,by=list(merge_results_p$date),FUN=max)
  obser<-aggregate(merge_results_p$PERENnew,by=list(merge_results_p$date),FUN=max)
  r_peren<-round(cor(round(pred$x),obser$x,method="spearman"),4)
  p_peren<-round(as.numeric(cor.test(round(pred$x),obser$x)[3]),4)
  plot(pred$x~jitter(obser$x),col='blue',ylab='predicted greenness',xlab='observed greenness',main=paste("perrenial plants, roots ",DEP[root_deep]," cm",sep=""),xlim=c(0,11))
  text(2,max(pred$x)-1,paste("r=",r_peren,"p=",p_peren))
    dev.off()

  
  filename<-paste("results/locust correl roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days site ",ii,".pdf",sep="")
  pdf(filename,paper="A4",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(3,1)) # set up for 12 plots in 2 columns
  par(oma = c(2,2,2,2) + 1) # margin spacing stuff
  par(mar = c(3,3,1,1) + 1) # margin spacing stuff 
  par(mgp = c(3,1,0) ) # margin spacing stuff 
  stagefreqDF_agg$date1<-stagefreqDF_agg$date
  stagefreqDF_agg$locusts_pred<-stagefreqDF_agg$nymphs+stagefreqDF_agg$Adult
  locustsub$locusts_obs<-locustsub$NDENS+locustsub$ADENS
  locustsub_nymphs<-ddply(locustsub,.(date),summarise,NDENS = max(NDENS))
  colnames(locustsub_nymphs)[1]<-'date1'
  locustsub_locusts<-ddply(locustsub,.(date),summarise,locusts_obs = max(locusts_obs))
  colnames(locustsub_locusts)[1]<-'date1'  
  locustsub_adults<-ddply(locustsub,.(date),summarise,ADENS = max(ADENS))
  colnames(locustsub_adults)[1]<-'date1'
  merge_results_n<-merge(stagefreqDF_agg,locustsub_nymphs,by="date1")
  merge_results_a<-merge(stagefreqDF_agg,locustsub_adults,by="date1")
  merge_results_l<-merge(stagefreqDF_agg,locustsub_locusts,by="date1")
  r_nymphs<-round(cor(merge_results_n$nymphs,merge_results_n$NDENS),3)
  r_adults<-round(cor(merge_results_a$Adult,merge_results_a$ADENS),3)
  r_locusts<-round(cor(merge_results_l$locusts_pred,merge_results_l$locusts_obs),3)
  p_nymphs<-round(as.numeric(cor.test(merge_results_n$nymphs,merge_results_n$NDENS)[3]),3)
  p_adults<-round(as.numeric(cor.test(merge_results_a$Adult,merge_results_a$ADENS)[3]),3)
  p_locusts<-round(as.numeric(cor.test(merge_results_l$locusts_pred,merge_results_l$locusts_obs)[3]),3)
  plot(merge_results_n$nymphs~jitter(merge_results_n$NDENS),col='blue',ylab='predicted nymphs',xlab='observed nymphs',main='nymphs')
  text(2,max(merge_results_n$nymphs)*.75,paste("r=",r_nymphs,"p=",p_nymphs))
  plot(merge_results_a$Adult~jitter(merge_results_a$ADENS),col='blue',ylab='predicted adults',xlab='observed adults',main='adults')
  text(2,max(merge_results_a$Adult)*.75,paste("r=",r_adults,"p=",p_adults))
  plot(merge_results_l$locusts_pred~jitter(merge_results_l$locusts_obs),col='blue',ylab='predicted locusts',xlab='observed locusts',main='adults')
  text(2,max(merge_results_l$locusts_pred)*.75,paste("r=",r_locusts,"p=",p_locusts))
  dev.off()

  all_e<-as.data.frame(cbind(DEP[root_deep],merge_results_e))
  all_p<-as.data.frame(cbind(DEP[root_deep],merge_results_p))
  all_n<-as.data.frame(cbind(DEP[root_deep],merge_results_n))
  all_a<-as.data.frame(cbind(DEP[root_deep],merge_results_a))
  all_l<-as.data.frame(cbind(DEP[root_deep],merge_results_l))
  all_r_e<-as.data.frame(cbind(ii,DEP[root_deep],r_ephem,p_ephem))
  all_r_p<-as.data.frame(cbind(ii,DEP[root_deep],r_peren,p_peren))
  all_r_n<-as.data.frame(cbind(ii,DEP[root_deep],r_nymphs,p_nymphs))
  all_r_a<-as.data.frame(cbind(ii,DEP[root_deep],r_adults,p_adults))
  all_r_l<-as.data.frame(cbind(ii,DEP[root_deep],r_locusts,p_locusts))

  write.table(all_e,"results/all_e.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_p,"results/all_p.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_n,"results/all_n.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_a,"results/all_a.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_l,"results/all_l.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_r_e,"results/all_r_e.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_r_p,"results/all_r_p.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_r_n,"results/all_r_n.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_r_a,"results/all_r_a.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_r_l,"results/all_r_l.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)

}