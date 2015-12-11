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
  loc<-longlat
  #get soil data
  source("../micro_australia/get.soil.R")
  source("../micro_australia/micro_aust.R")
  DEP = c(0., 1.,  3.,  5,  10,  20.,  30.,  60.,  90.,  200.)  
  soil.hydro<-get.soil()
  #soil.hydro<-read.csv("data.csv")
 # PE<-rep(0,19)
  #KS<-PE
  #BB<-PE
  #BD<-PE
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
      BD, KS = KS, BulkDensity = BulkDensity, dailywind = 0)
  #micro<-micro_global(loc = longlat, runmoist = 1, nyears = nyears, timeinterval = 365, PE = PE, BB = BB, BD = 
  #    BD, KS = KS)
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
  
  # get subset for current site 
  attach(locustdata)
  locustsub<-locustdata[(long == longlat[1])&(lat==longlat[2]),]
  detach(locustdata)
  locustsub$date<-as.character(locustsub$DATE_,format="%Y-%m-%d")
  
  root_shallow<-3#kk-1#2#5 # how shallow do the roots go? 2 to 10, corresopnding to 1, 3, 5, 10, 20, 30, 60, 90 and 200 cm
  root_deep<-6#kk#6 # how deep do the roots go? 2 to 10, corresopnding to 1, 3, 5, 10, 20, 30, 60, 90 and 200 cm
  growth_delay<-1# days after suitable soil moisture that new growth occurs
  wilting_thresh<-200*-1 # water potential for wilting point J/kg (divide by 100 to get bar)
  permanent_wilting_point<-1500*-1 # water potential for permanent wilting point (PWP) J/kg (divide by 100 to get bar)
  FoodWater<-82 # water content of fully hydrated veg
  
  source("plantgro.R")
  plantout<-plantgro(soilpot = soilpot, soilmoist = soilmoist, root_shallow = root_shallow, 
    root_deep = root_deep, growth_delay = growth_delay, wilting_thresh=wilting_thresh, permanent_wilting_point=permanent_wilting_point)
  #plantout<-subset(plantout, format(plantout$date1,"%H")=="12")

  plant.pres<-plantout$plant.pres
  plant.moist<-plantout$pct.water
 
  
 #plot(plantout$moist*plant.pres~plantout$date1,type='h',col='green')
 #points(locustsub$Enew~as.POSIXct(locustsub$date),type = 'h')
  
  ###################### run locust simulation ##################################
  
  recover<-1 # time locust needs to build up resources to lay first batch
  ovidates<-as.data.frame(cbind(plant.pres[1:(length(plant.pres)-1)],plant.pres[2:length(plant.pres)]))
  ovidates<-cbind(dates[2:length(dates)],ovidates)
  ovidates2<-subset(ovidates, V1-V2==-1)
  ovirows_start<-cbind(1,as.numeric(rownames(ovidates2))) # rows at which grass growth started
  ovidates2<-subset(ovidates, V2-V1==-1)
  ovirows_finish<-cbind(0,as.numeric(rownames(ovidates2))) # rows at which grass growth finish
  ovirows<-rbind(ovirows_start,ovirows_finish)
  ovirows<-rbind(c(0,recover),ovirows) # add first day of sim as an oviposition date
  
  ########### chortoicetes oviposition model##########
  
  # to do, for now assume on first and last days of grass growth periods
  
  ############### chortoicetes egg model #############
  
  DATA1<-cbind(soil[,1:8],metout[,13:14],environ[,19],soilmoist[,5:7]*100)
  colnames(DATA1)<-c('DATE','JULDAY','TIME','D0cm','D1cm','D3cm','D5cm','D10cm','ZEN','SOLR','PHOTO','MOIST3cm','MOIST5cm','MOIST10cm')
  
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
  stagefreqnames<-c('Eggdev', 'Q1', 'Q2', 'Diapause', 'N1','N2','N3','N4','N5','Adult','Death' )
  for(name in stagefreqnames){ # make variables to easily access matrix column by name
    eval(parse(text=paste(name,"col<-",which(stagefreqnames==name), sep = ""))) 
  }
  stagefreq<-rep(0,nrow(DATA1)*length(stagefreqnames))
  dim(stagefreq)<-c(nrow(DATA1),length(stagefreqnames))
  
  
  for(g in 1:15){
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
      #plot(environ$TC~environ$dates,type='l',col='orange',ylim=c(0,60))
      #points(plant.moist~environ$dates,type='l',col='green')
      #points(allgens$mass~allgens$dates,type='p',col='blue',cex=0.1)
    }else{
      break
    }
  }
  
  success<-as.data.frame(subset(allgens,mass>55))
  
  
  
  stagefreqDF<-as.data.frame(generations$stagefreq)
  
  stagefreqDF<-cbind(dates,longlat[1],longlat[2],stagefreqDF)
  names(stagefreqDF)<-c('date','long','lat',stagefreqnames)
  
  stagefreqDF_agg<-aggregate(stagefreqDF[,2:14],by=list(format(dates,'%Y-%m-%d')),max)
  colnames(stagefreqDF_agg)[1]<-'date'
  
  stagefreqDF_agg$nymphs<-rowSums(stagefreqDF_agg[,9:12]) # exclude 1st instars as they may have just hatched and died if no food
  
  filename<-paste("plant growth test output/time series ",DEP[root_deep]," cm regrow thresh ",growth_delay," days site ",ii,".pdf",sep="")
  pdf(filename,paper="A4r",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(3,1)) # set up for 12 plots in 2 columns
  par(oma = c(2,2,2,2) + 1) # margin spacing stuff
  par(mar = c(3,3,1,1) + 1) # margin spacing stuff 
  par(mgp = c(3,1,0) ) # margin spacing stuff   
  plot(stagefreqDF_agg$nymphs/max(stagefreqDF_agg$nymphs)*8~as.POSIXct(stagefreqDF_agg$date),type='l',ylim=c(0,8),ylab="nymph density",xlab="")
  points(plant.moist/max(plant.moist)*3*plant.pres~plantout$date1,type='h',col="darkolivegreen2")  
  points(locustsub$NDENS-0.5~locustsub$DATE_,type='h',col='orange',lwd=2)  
  points(stagefreqDF_agg$nymphs/max(stagefreqDF_agg$nymphs)*8~as.POSIXct(stagefreqDF_agg$date),type='l')
  


  
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
  plot(merge_results_n$nymphs~jitter(merge_results_n$NDENS),col='blue',ylab='predicted nymphs',xlab='observed nymphs',main='nymphs')
  text(min(merge_results_n$NDENS)+1*1.2,max(merge_results_n$nymphs)*.75,paste("r=",r_nymphs))

#   plot(stagefreqDF_agg$Adult/max(stagefreqDF_agg$Adult)*8~as.POSIXct(stagefreqDF_agg$date),type='l',ylim=c(0,8),ylab="adult density",xlab="")
#   points(plant.moist/max(plant.moist)*3*plant.pres~plantout$date1,type='h',col="darkolivegreen2")  
#   points(locustsub$ADENS-0.5~locustsub$DATE_,type='h',col='orange',lwd=2)  
#   points(stagefreqDF_agg$Adult/max(stagefreqDF_agg$Adult)*8~as.POSIXct(stagefreqDF_agg$date),type='l')
#   
#   plot((stagefreqDF_agg$Adult+stagefreqDF_agg$nymphs)/max(stagefreqDF_agg$Adult+stagefreqDF_agg$nymphs)*8~as.POSIXct(stagefreqDF_agg$date),type='l',ylim=c(0,8),ylab="adult and nymph density",xlab="")
#   points(plant.moist/max(plant.moist)*3*plant.pres~plantout$date1,type='h',col="darkolivegreen2")  
#   points(locustsub$ADENS + locustsub$NDENS-0.5~locustsub$DATE_,type='h',col='orange',lwd=2)  
#   points((stagefreqDF_agg$Adult+stagefreqDF_agg$nymphs)/max(stagefreqDF_agg$Adult+stagefreqDF_agg$nymphs)*8~as.POSIXct(stagefreqDF_agg$date),type='l')  
  title(paste("lat/long ",longlat[2],",",longlat[1],sep=""))
  dev.off()
  
}
  
  #plot plant growth metric against observed plant growth index
  filename<-paste("plant growth test output/perennial roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days site ",ii,"_1990_1999.pdf",sep="")
  pdf(filename,paper="A4",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(5,2)) # set up for 12 plots in 2 columns
  par(oma = c(2,2,2,2) + 0.1) # margin spacing stuff
  par(mar = c(3,3,1,1) + 0.1) # margin spacing stuff 
  par(mgp = c(3,1,0) ) # margin spacing stuff 
  
  for(yr in 1990:1999){
    plotplant.moist<-subset(plantout,as.numeric(format(plantout$date1, "%Y"))==yr)
    plot(plotplant.moist$moist~plotplant.moist$date1,type='l',col='dark green',main=yr,ylim=c(0,11))
    points(stagefreqDF_agg$nymphs~as.POSIXct(stagefreqDF_agg$date),type='h',col='orange')
    points(stagefreqDF_agg$Adult~as.POSIXct(stagefreqDF_agg$date),type='h',col='red',lty=2)  
    #points(locustsub$Enew~locustsub$DATE_,col='red',type='h',lwd=2)
    points(locustsub$PERENnew~locustsub$DATE_,col='blue',type='h',lwd=2)
    #points(rainfall$rainfall/10~rainfall$dates,lty=1,type='h',col='light blue')
    points(locustsub$NDENS~as.POSIXct(locustsub$DATE_), type='p',col='black',pch=16,cex=2)
    #points(locustsub$ADENS~as.POSIXct(locustsub$DATE_), type='p',col='grey',pch=16)  
    points(plotplant.moist$moist~plotplant.moist$date1,type='l',col='dark green',main=yr,ylim=c(0,11))
    
  }
  title(paste("perennial plants, roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days, lat/long ",longlat[2],",",longlat[1],sep=""),outer=TRUE)
  dev.off()
  
  #plot plant growth metric against observed plant growth index
  filename<-paste("plant growth test output/perennial roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days site ",ii,"_2000_2009.pdf",sep="")
  pdf(filename,paper="A4",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(5,2)) # set up for 12 plots in 2 columns
  par(oma = c(2,2,2,2) + 0.1) # margin spacing stuff
  par(mar = c(3,3,1,1) + 0.1) # margin spacing stuff 
  par(mgp = c(3,1,0) ) # margin spacing stuff 
  
  for(yr in 2000:2009){
    plotplant.moist<-subset(plantout,as.numeric(format(plantout$date1, "%Y"))==yr)
    plot(plotplant.moist$moist~plotplant.moist$date1,type='l',col='dark green',main=yr,ylim=c(0,11))
    points(stagefreqDF_agg$nymphs~as.POSIXct(stagefreqDF_agg$date),type='h',col='orange')
    points(stagefreqDF_agg$Adult~as.POSIXct(stagefreqDF_agg$date),type='h',col='red',lty=2)  
    #points(locustsub$Enew~locustsub$DATE_,col='red',type='h',lwd=2)
    points(locustsub$PERENnew~locustsub$DATE_,col='blue',type='h',lwd=2)
    #points(rainfall$rainfall/10~rainfall$dates,lty=1,type='h',col='light blue')
    points(locustsub$NDENS~as.POSIXct(locustsub$DATE_), type='p',col='black',pch=16,cex=2)
    points(locustsub$ADENS~as.POSIXct(locustsub$DATE_), type='p',col='grey',pch=16)  
    points(plotplant.moist$moist~plotplant.moist$date1,type='l',col='dark green',main=yr,ylim=c(0,11))
  }
  title(paste("perennial plants, roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days, lat/long ",longlat[2],",",longlat[1],sep=""),outer=TRUE)
  dev.off()
  
  filename<-paste("plant growth test output/ephemeral roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days site ",ii,".pdf",sep="")
  pdf(filename,paper="A4",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(6,2)) # set up for 12 plots in 2 columns
  par(oma = c(2,2,2,2) + 0.01) # margin spacing stuff
  par(mar = c(3,3,1,1) + 0.01) # margin spacing stuff 
  par(mgp = c(3,1,0) ) # margin spacing stuff 
  
  for(yr in 1998:2009){
    plotplant.moist<-subset(plantout,as.numeric(format(plantout$date1, "%Y"))==yr)
    plot(plotplant.moist$moist~plotplant.moist$date1,type='l',col='dark green',main=yr,ylim=c(0,11))
    points(stagefreqDF_agg$nymphs~as.POSIXct(stagefreqDF_agg$date),type='h',col='orange')
    points(stagefreqDF_agg$Adult~as.POSIXct(stagefreqDF_agg$date),type='h',col='red',lty=2)     
    points(locustsub$Enew~locustsub$DATE_,col='blue',type='h',lwd=1)
    #points(locustsub$PERENnew~locustsub$DATE_,col='blue',type='h',lwd=1)
    #points(rainfall$rainfall/10~rainfall$dates,lty=1,type='h',col='light blue')
    points(locustsub$NDENS~as.POSIXct(locustsub$DATE_), type='p',col='black',pch=16,cex=2)
    points(locustsub$ADENS~as.POSIXct(locustsub$DATE_), type='p',col='grey',pch=16)  
  }
  title(paste("ephemeral plants, roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days, lat/long ",longlat[2],",",longlat[1],sep=""),outer=TRUE)
  dev.off()
  
  filename<-paste("plant growth test output/correl roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days site ",ii,".pdf",sep="")
  pdf(filename,paper="A4",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(2,1)) # set up for 12 plots in 2 columns
  par(oma = c(2,2,2,2) + 1) # margin spacing stuff
  par(mar = c(3,3,1,1) + 1) # margin spacing stuff 
  par(mgp = c(3,1,0) ) # margin spacing stuff 
  
  
  
  plantout2<-plantout
  plantout2$moist<-plantout2$moist*plantout2$plant.pres
  plantout2$date<-as.character(as.Date(plantout$date1,format="%Y-%m-%h"))  

  merge_results_e<-merge(plantout2,locustsub,by="date")
  merge_results_e<-subset(merge_results_e,Enew!='NA')
  pred<-aggregate(merge_results_e$moist,by=list(merge_results_e$date),FUN=max)
  obser<-aggregate(merge_results_e$Enew,by=list(merge_results_e$date),FUN=max)
  r_ephem<-round(cor(round(pred$x),obser$x,method="spearman"),4)
  p_ephem<-round(as.numeric(cor.test(round(pred$x),obser$x)[3]),4)
  plot(pred$x~jitter(obser$x),col='blue',ylab='predicted greenness',xlab='observed greenness',main=paste("ephemeral plants, roots ",DEP[root_deep]," cm",sep=""),xlim=c(0,11))
  text(1,max(pred$x)-1,paste("r=",r_ephem,"p=",p_ephem))

  merge_results_p<-merge(plantout2,locustsub,by="date")
  merge_results_p<-subset(merge_results_p,PERENnew!='NA')
  pred<-aggregate(merge_results_p$moist,by=list(merge_results_p$date),FUN=max)
  obser<-aggregate(merge_results_p$PERENnew,by=list(merge_results_p$date),FUN=max)
  r_peren<-round(cor(round(pred$x),obser$x,method="spearman"),4)
  p_peren<-round(as.numeric(cor.test(round(pred$x),obser$x)[3]),4)
  plot(pred$x~jitter(obser$x),col='blue',ylab='predicted greenness',xlab='observed greenness',main=paste("perrenial plants, roots ",DEP[root_deep]," cm",sep=""),xlim=c(0,11))
  text(1,max(pred$x)-1,paste("r=",r_peren,"p=",p_peren))
  
    
  library(plyr)
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
  plot(merge_results_n$nymphs~jitter(merge_results_n$NDENS),col='blue',ylab='predicted nymphs',xlab='observed nymphs',main='nymphs')
  text(min(merge_results_n$NDENS)+1*1.2,max(merge_results_n$nymphs)*.75,paste("r=",r_nymphs))
  plot(merge_results_a$Adult~jitter(merge_results_a$ADENS),col='blue',ylab='predicted adults',xlab='observed adults',main='adults')
  text(min(merge_results_a$ADENS)+1*1.2,max(merge_results_a$Adult)*.75,paste("r=",r_adults))
  plot(merge_results_l$locusts_pred~jitter(merge_results_l$locusts_obs),col='blue',ylab='predicted locusts',xlab='observed locusts',main='adults')
  text(min(merge_results_l$locusts_obs)+1*1.2,max(merge_results_l$locusts_pred)*.75,paste("r=",r_locusts))
  dev.off()
  
  #if(ii==1 & kk==nodestart){
  all_e<-as.data.frame(cbind(DEP[root_deep],merge_results_e))
  all_p<-as.data.frame(cbind(DEP[root_deep],merge_results_p))
  all_l<-as.data.frame(cbind(DEP[root_deep],merge_results))
  all_r_e<-as.data.frame(cbind(ii,DEP[root_deep],r_ephem))
  all_r_p<-as.data.frame(cbind(ii,DEP[root_deep],r_peren))
  all_r_n<-as.data.frame(cbind(ii,DEP[root_deep],r_nymphs))
  all_r_a<-as.data.frame(cbind(ii,DEP[root_deep],r_adults))
  #}else{
  #all_e<-rbind(all_e,as.data.frame(cbind(DEP[root_deep],merge_results_e)))
  #all_p<-rbind(all_p,as.data.frame(cbind(DEP[root_deep],merge_results_p)))
  #  all_r_e<-rbind(all_r_e,as.data.frame(cbind(ii,DEP[root_deep],r_ephem)))
  #  all_r_p<-rbind(all_r_p,as.data.frame(cbind(ii,DEP[root_deep],r_peren)))
  #}  
  
  write.table(all_e,"plant growth test output/all_e.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_p,"plant growth test output/all_p.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_l,"plant growth test output/all_l.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_r_e,"plant growth test output/all_r_e.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_r_p,"plant growth test output/all_r_p.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_r_n,"plant growth test output/all_r_n.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(all_r_a,"plant growth test output/all_r_a.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  
  
  
}