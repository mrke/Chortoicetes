
spatial<-"c:/Australian Environment/" # place where climate input files are kept
maindir<-getwd()
survey_freq<-read.csv('survey_freq_98-09.csv')
locustdata<-read.csv('locustdata_90-09.csv')
locustdata$DATE_<-as.POSIXct(locustdata$DATE_,format="%Y-%m-%d")
ii<-1
for(ii in 1:50){
  
############## location and climatic data  ###################################
sitemethod <- 0 # 0=specified single site long/lat, 1=place name search using geodis (needs internet)
longlat<-c(survey_freq[ii,2],survey_freq[ii,3]) # type a long/lat here in decimal degrees
loc <- "Broken Hill, Australia" # type in a location here, used if option 1 is chosen above
timezone<-0 # if timezone=1 (needs internet), uses GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)
rungads<-1 # use the Global Aerosol Database?
dailywind<-1 # use daily windspeed database?
terrain<-0 # include terrain (slope, aspect, horizon angles) (1) or not (0)?
soildata<-1 # include soil data for Australia (1) or not (0)?
snowmodel<-0 # run snow version? (slower!)
ystart<-1998
yfinish<-2009
nyears<-yfinish-ystart+1# integer, number of years for which to run the microclimate model, only for AWAP data (!!max 10 years!!)

############# microclimate model parameters ################################
EC <- 0.0167238 # Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)
RUF <- 0.004 # Roughness height (m), , e.g. sand is 0.05, grass may be 2.0, current allowed range: 0.001 (snow) - 2.0 cm.
# Next for parameters are segmented velocity profiles due to bushes, rocks etc. on the surface, IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO!
Z01 <- 0. # Top (1st) segment roughness height(m)
Z02 <- 0. # 2nd segment roughness height(m)
ZH1 <- 0. # Top of (1st) segment, height above surface(m)
ZH2 <- 0. # 2nd segment, height above surface(m)
SLE <- 0.96 # Substrate longwave IR emissivity (decimal %), typically close to 1
ERR <- 2.0 # Integrator error for soil temperature calculations
DEP <- c(0., 1.,  3.,  5,  10,  20.,  30.,  60.,  90.,  200.) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
# 
library(dismo)
if(sitemethod==1){
longlat <- geocode(loc)[1, 3:4] # assumes first geocode match is correct
}
prevdir<-getwd()
setwd('x:')
cmd<-paste("R --no-save --args ",longlat[1]," ",longlat[2]," < extract.R",sep='')
system(cmd)
soilpro<-read.csv('data.csv')
FC<-(7.561+1.176*soilpro$clay-0.009843*soilpro$clay^2+0.2132*soilpro$silt)/100
PWP<-(-1.304+1.117*soilpro$clay-0.009309*soilpro$clay^2)/100
setwd(prevdir)
#write.table(cbind(isite,longlat[1],longlat[2],soilpro), file = "c:/git/micro_australia/soilprops.txt", sep = ",", col.names = F, qmethod = "double", append = T)
   
# pre-extracted
#soilpro<-read.csv("c:/git/ectotherm/sleepy/soilprops.txt",header=FALSE)
#colnames(soilpro)<-c('i','site','long','lat','desc','blkdens','clay','silt','sand')
#soilpro<-subset(soilpro,site==1)
#soilpro<-soilpro[,5:9]
   
soil_depths<-c(2.5,7.5,22.5,45,80,150)
plot(soilpro$clay~soil_depths,ylim=c(0,100),col='red',type='l')
points(soilpro$sand~soil_depths,ylim=c(0,100),col='orange',type='l')
points(soilpro$silt~soil_depths,ylim=c(0,100),col='grey',type='l')
title(main=loc)
legend("topleft", inset=.05,
       legend=round(soilpro[1,3:5],1),bty="n", 
       horiz=TRUE, bg=NULL, cex=0.8)

DEP2<-rep(0,18)
j<-1
for(i in 1:length(DEP2)){
  if(i%%2==0){
    DEP2[i]<-DEP2[i-1]+(DEP[j]-DEP2[i-1])/2
  }else{
    DEP2[i]<-DEP[j]
    j<-j+1
  }
}
DEP2<-as.data.frame(floor(DEP2))
colnames(DEP2)<-"DEPTH"
 
# par(mfrow=c(2,2))
# par(mar = c(2,2,1,2) + 0.1) # margin spacing stuff 
# par(mar = c(5,5,5,5) + 0.1) # margin spacing stuff 

CampNormTbl9_1<-read.csv('../micro_australia/CampNormTbl9_1.csv')
dclay<-0.001 #mm
dsilt<-0.026 #mm
dsand<-1.05 #mm
a<-(soilpro$clay/100)*log(dclay) + (soilpro$sand/100)*log(dsand) + (soilpro$silt/100)*log(dsilt)
b.1<-(((soilpro$clay/100)*log(dclay)^2+(soilpro$sand/100)*log(dsand)^2+(soilpro$silt/100)*log(dsilt)^2)-a^2)^(1/2)
dg<-exp(a)
sigma_g<-exp(b.1)
PES<-(0.5*dg^(-1/2))*-1
b<--2*PES+0.2*sigma_g
PE<-PES*(soilpro$blkdens/1.3)^(0.67*b)
KS<-0.004*(1.3/soilpro$blkdens)^(1.3*b)*exp(-6.9*soilpro$clay/100-3.7*soilpro$silt/100)
BD<-soilpro$blkdens
   

plot(KS~soil_depths,xlim=c(-1,200),ylim=c(0.000017,0.0058))
KS_spline <-spline(soil_depths,KS,n=201,xmin=0,xmax=200,method='natural')
points(KS_spline$y,col='red',type='l')
KS_spline<-as.data.frame(cbind(KS_spline$x,KS_spline$y))
colnames(KS_spline)<-c('DEPTH','VALUE')
KS<-merge(DEP2,KS_spline)
KS<-c(KS[1,2],KS[,2])
KS[KS<0.000017]<-0.000017
   
plot(PE~soil_depths,xlim=c(-1,200),ylim=c(-15,0))
PE_spline <-spline(soil_depths,PE,n=201,xmin=0,xmax=200,method='natural')
points(PE_spline$y,col='red',type='l')
PE_spline<-as.data.frame(cbind(PE_spline$x,PE_spline$y))
colnames(PE_spline)<-c('DEPTH','VALUE')
PE<-merge(DEP2,PE_spline)
PE<-c(-1*PE[1,2],-1*PE[,2])
   
plot(b~soil_depths,xlim=c(-1,200),ylim=c(2,24))
b_spline <-spline(soil_depths,b,n=201,xmin=0,xmax=200,method='natural')
points(b_spline$y,col='red',type='l')
b_spline<-as.data.frame(cbind(b_spline$x,b_spline$y))
colnames(b_spline)<-c('DEPTH','VALUE')
b<-merge(DEP2,b_spline)
BB<-c(b[1,2],b[,2])
   
plot(BD~soil_depths,xlim=c(-1,200),ylim=c(1,1.6))
BD_spline <-spline(soil_depths,BD,n=201,xmin=0,xmax=200,method='natural')
points(BD_spline$y,col='red',type='l')
BD_spline<-as.data.frame(cbind(BD_spline$x,BD_spline$y))
colnames(BD_spline)<-c('DEPTH','VALUE')
BD<-merge(DEP2,BD_spline)
BD<-c(BD[1,2],BD[,2])
  
Thcond <- 2.5 # soil minerals thermal conductivity (W/mC)
Density <- 2560. # soil minerals density (kg/m3)
SpecHeat <- 870. # soil minerals specific heat (J/kg-K)
BulkDensity <- 1300 # soil bulk density (kg/m3)
cap<-1 # organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)
SatWater <- 0.26 # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
Clay <- 20 # clay content for matric potential calculations (%)
SoilMoist <- 0 # fractional soil moisture (decimal %)
rainmult<-1 # rain multiplier for surface soil moisture (use to induce runoff), proportion
runmoist<-1 # run soil moisture model (0=no, 1=yes)?
SoilMoist_Init<-c(0.1,0.12,0.15,0.3,0.4,0.4,0.4,0.4,0.4,0.4) # initial soil water content, m3/m3
evenrain<-0 # spread daily rainfall evenly across 24hrs (1) or one event at midnight (0)
maxpool<-10000#6 # max depth for water pooling on the surface, mm (to account for runoff)
L<-c(0,0,rep(4,9),1.8,0.95,0.85,0.8,0.4,0.366,0,0)*10000
L<-c(0,0,8.18990859,7.991299442,7.796891252,7.420411664,7.059944542,6.385001059,5.768074989,4.816673431,4.0121088,1.833554792,0.946862989,0.635260544,0.804575,0.43525621,0.366052856,0,0)*10000
LAI<-0.1 # leaf area index, used to partition traspiration/evaporation from PET
PE[1:9]<-CampNormTbl9_1[3,4] #air entry potential J/kg 
KS[1:9]<-CampNormTbl9_1[3,6] #saturated conductivity, kg s/m3
BB[1:9]<-CampNormTbl9_1[3,5] #soil 'b' parameter
PE[10:13]<-CampNormTbl9_1[4,4] #air entry potential J/kg 
KS[10:13]<-CampNormTbl9_1[4,6] #saturated conductivity, kg s/m3
BB[10:13]<-CampNormTbl9_1[4,5] #soil 'b' parameter
REFL<-0.15 # soil reflectance (decimal %)
slope<-0. # slope (degrees, range 0-90)
aspect<-180. # aspect (degrees, 0 = North, range 0-360)
hori<-rep(0,24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
PCTWET<-40 # percentage of surface area acting as a free water surface (%)
CMH2O <- 1. # precipitable cm H2O in air column, 0.1 = VERY DRY; 1.0 = MOIST AIR CONDITIONS; 2.0 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)  
TIMAXS <- c(1.0, 1.0, 0.0, 0.0)   # Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise            											
TIMINS <- c(0.0, 0.0, 1.0, 1.0)   # Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon
minshade<-0. # minimum available shade (%)
maxshade<-70. # maximum available shade (%)
runshade<-1. # run the model twice, once for each shade level (1) or just for the first shade level (0)?
manualshade<-1 # if using soildata, which includes shade, this will override the data from the database and force max shade to be the number specified above
Usrhyt <- 3# local height (cm) at which air temperature, relative humidity and wind speed calculatinos will be made 
rainwet<-1.5 # mm rain that causes soil to become 90% wet
snowtemp<-1.5 # temperature at which precipitation falls as snow (used for snow model)
snowdens<-0.4 # snow density (mg/m3)
snowmelt<-1. # proportion of calculated snowmelt that doesn't refreeze
undercatch<-1. # undercatch multipier for converting rainfall to snow
rainmelt<-0.016 # paramter in equation that melts snow with rainfall as a function of air temp
write_input<-1 # write csv files of final input to working directory? 1=yes, 0=no.
warm<-0 # uniform warming of air temperature input to simulate climate change
loop<-0 # if doing multiple years, this shifts the starting year by the integer value

# run the model
niche<-list(L=L,LAI=LAI,SoilMoist_Init=SoilMoist_Init,evenrain=evenrain,runmoist=runmoist,maxpool=maxpool,PE=PE,KS=KS,BB=BB,BD=BD,loop=loop,warm=warm,rainwet=rainwet,manualshade=manualshade,dailywind=dailywind,terrain=terrain,soildata=soildata,loc=loc,ystart=ystart,yfinish=yfinish,nyears=nyears,RUF=RUF,SLE=SLE,ERR=ERR,DEP=DEP,Thcond=Thcond,Density=Density,SpecHeat=SpecHeat,BulkDensity=BulkDensity,Clay=Clay,SatWater=SatWater,SoilMoist=SoilMoist,CMH2O=CMH2O,TIMAXS=TIMAXS,TIMINS=TIMINS,minshade=minshade,maxshade=maxshade,Usrhyt=Usrhyt,REFL=REFL,slope=slope,aspect=aspect,hori=hori,rungads=rungads,cap=cap,write_input=write_input,spatial=spatial,snowmodel=snowmodel,snowtemp=snowtemp,snowdens=snowdens,snowmelt=snowmelt,undercatch=undercatch,rainmelt=rainmelt,rainmult=rainmult,runshade=runshade)
setwd('../micro_australia/')
source('NicheMapR_Setup_micro.R')
nicheout<-NicheMapR(niche)
setwd(maindir)

# get output
metout<-as.data.frame(nicheout$metout[1:(365*24*nyears),]) # above ground microclimatic conditions, min shade
shadmet<-as.data.frame(nicheout$shadmet[1:(365*24*nyears),]) # above ground microclimatic conditions, max shade
soil<-as.data.frame(nicheout$soil[1:(365*24*nyears),]) # soil temperatures, minimum shade
shadsoil<-as.data.frame(nicheout$shadsoil[1:(365*24*nyears),]) # soil temperatures, maximum shade
soilmoist<-as.data.frame(nicheout$soilmoist[1:(365*24*nyears),]) # soil water content, minimum shade
shadmoist<-as.data.frame(nicheout$shadmoist[1:(365*24*nyears),]) # soil water content, maximum shade
humid<-as.data.frame(nicheout$humid[1:(365*24*nyears),]) # soil humidity, minimum shade
shadhumid<-as.data.frame(nicheout$shadhumid[1:(365*24*nyears),]) # soil humidity, maximum shade
soilpot<-as.data.frame(nicheout$soilpot[1:(365*24*nyears),]) # soil water potential, minimum shade
shadpot<-as.data.frame(nicheout$shadpot[1:(365*24*nyears),]) # soil water potential, maximum shade
rainfall<-as.data.frame(nicheout$RAINFALL)
MAXSHADES<-as.data.frame(nicheout$MAXSHADES)
elev<-as.numeric(nicheout$ALTT)
REFL<-as.numeric(nicheout$REFL)
longlat<-as.matrix(nicheout$longlat)
fieldcap<-as.numeric(nicheout$fieldcap)
wilting<-as.numeric(nicheout$wilting)
ectoin<-c(elev,REFL,longlat,fieldcap,wilting,ystart,yfinish)

tzone<-paste("Etc/GMT-",10,sep="") # doing it this way ignores daylight savings!
dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
dates2<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days") 
dates2<-subset(dates2, format(dates2, "%m/%d")!= "02/29") # remove leap years
rainfall<-as.data.frame(cbind(dates2,rainfall))
colnames(rainfall)<-c("dates","rainfall")

# get subset for current site 
attach(locustdata)
locustsub<-locustdata[(long == longlat[1])&(lat==longlat[2]),]
detach(locustdata)
locustsub$date<-as.character(locustsub$DATE_,format="%Y-%m-%d")

#locustsub$Enew[is.na(locustsub$PERENnew)==FALSE & is.na(locustsub$Enew)==TRUE] <- 0 # assume that if perennials are observed and no data for ephemerals, then no ephemerals 
#locustsub$PERENnew[is.na(locustsub$Enew)==FALSE & is.na(locustsub$PERENnew)==TRUE] <- 0 # assume that if ephemerals are observed and no data for perennials, then no perennials
  
nodestart<-3
nodefinish<-9
for(kk in nodestart:nodefinish){
  ##################### parameters:
  root_deep<-kk#6 # how deep do the roots go? 2 to 10, corresopnding to 1, 3, 5, 10, 20, 30, 60, 90 and 200 cm
  root_shallow<-kk-1#2#5 # how shallow do the roots go? 2 to 10, corresopnding to 1, 3, 5, 10, 20, 30, 60, 90 and 200 cm
  growth_delay<-1 # days after suitable soil moisture that new growth occurs
  wilting_thresh<-200*-1 # water potential for wilting point J/kg (divide by 100 to get bar)
  permanent_wilting_point<-1500*-1 # water potential for permanent wilting point (PWP) J/kg (divide by 100 to get bar)
  FoodWater<-82 # water content of fully hydrated veg
  #################################
  
  grassgrowths<-as.data.frame(soilpot)
  soilmoist2b<-as.data.frame(soilmoist)
  soilmoist2b<-subset(soilmoist2b,TIME==720)
  grassgrowths<-subset(grassgrowths,TIME==720)
  grassgrowths<-grassgrowths[,((root_shallow+3):(3+root_deep))] # get range of depths to take max of
  grassgrowths<-apply(grassgrowths, 1, mean)
  
  grow<-grassgrowths
  grow[grow>permanent_wilting_point]<-1 # find times when below the PWP
  grow[grow<=permanent_wilting_point]<-0 # find times when above the PWP
  counter<-0
  grow2<-grow*0
  for(j in 1:length(grow)){ # accumulate runs of days above the permanent wilting point (growth possible)
    if(j==1){
      if(grow[j]==1){
        counter<-counter+1
      }
      grow2[j]<-counter
    }else{
      if(grow[j-1]>0 & grow[j]==1){
        counter<-counter+1
      }else{
        counter<-0
      }
      grow2[j]<-counter
    }
  }
  grow3<-grow2
  grow3[grow3<growth_delay]<-0 # apply growth delay specified by the user for time required for plats to come back after PWP has been hit
  grow3[grow3>0]<-1 # make vector of 0 and 1 where 1 means plants could have come back from drought
  
  soilmoist2b<-soilmoist2b[,((root_shallow+3):(3+root_deep))] # get range of depths to take max of
  soilmoist2b<-apply(soilmoist2b, 1, mean)
  
  grassgrowths<-as.data.frame(cbind(grassgrowths,soilmoist2b))
  colnames(grassgrowths)<-c('pot','moist')
  grassgrowths$pot[grassgrowths$pot>wilting_thresh]<-FoodWater # assume plants start wilting at about 2 bar, but above this they are at max water content
  grassgrowths$moist<-grassgrowths$moist*100 # convert to percent
  potmult<-grassgrowths$pot
  potmult[potmult!=82]<-0
  potmult[potmult!=0]<-1  
  wilting<-subset(grassgrowths,pot==FoodWater) # find soil moisture range corresponding to values above the wilting point
  wilting<-min(wilting$moist) # get the min soil moisture at which plants aren't wilting
  grassgrowths<-grassgrowths$moist
  grassgrowths[grassgrowths>wilting]<-FoodWater # now have vector of either max plant water content or soil moisture content - need to convert the latter into a smooth decline to zero from max value
  minmoist<-0
  grassgrowths[grassgrowths<FoodWater]<-(grassgrowths[grassgrowths<FoodWater]-minmoist)/(wilting-minmoist)*FoodWater # for just the values less than max water content, make them equal to the ratio of moisture level and wilting point moisture and then multiplied by food water content so have moisture content ranging from max level down to min moisture possible
  grassgrowths<-grassgrowths/100*grow3 # now convert to proportion and cut out times below PWP (including regrowth penalty)
  grasstsdms<-grassgrowths
  grassmoist<-grassgrowths
  grassmoist<-as.data.frame(cbind(dates2,grassgrowths))
  colnames(grassmoist)<-c('date1','moist')
  grassmoist$date1<-dates2
  grassmoist$moist<-grassmoist$moist/max(grassmoist$moist)*11 # put in units scaling from 0-11
  # next four lines spread the values out more evenly over 11 categories
  minval<-min(grassmoist$moist[grassmoist$moist!=0])
  grassmoist$moist[grassmoist$moist==0]<-minval
  grassmoist$moist<-grassmoist$moist-minval
  grassmoist$moist<-grassmoist$moist/max(grassmoist$moist)*11 # put in units scaling from 0-11

  #plot plant growth metric against observed plant growth index
  filename<-paste("plant growth test output/perennial roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days site ",ii,".pdf",sep="")
  pdf(filename,paper="A4",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(6,2)) # set up for 12 plots in 2 columns
  par(oma = c(2,2,2,2) + 0.1) # margin spacing stuff
  par(mar = c(3,3,1,1) + 0.1) # margin spacing stuff 
  par(mgp = c(3,1,0) ) # margin spacing stuff 
  
  for(yr in 1998:2009){
  plotgrassmoist<-subset(grassmoist,as.numeric(format(grassmoist$date1, "%Y"))==yr)
  plot(plotgrassmoist$moist~plotgrassmoist$date1,type='l',col='dark green',main=yr,ylim=c(0,11))
  #points(locustsub$Enew~locustsub$DATE_,col='red',type='h',lwd=2)
  points(locustsub$PERENnew~locustsub$DATE_,col='blue',type='h',lwd=1)
  points(rainfall$rainfall/10~rainfall$dates,lty=1,type='h',col='light blue')
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
  plotgrassmoist<-subset(grassmoist,as.numeric(format(grassmoist$date1, "%Y"))==yr)
  plot(plotgrassmoist$moist~plotgrassmoist$date1,type='l',col='dark green',main=yr,ylim=c(0,11))
  points(locustsub$Enew~locustsub$DATE_,col='blue',type='h',lwd=1)
  #points(locustsub$PERENnew~locustsub$DATE_,col='blue',type='h',lwd=1)
  points(rainfall$rainfall/10~rainfall$dates,lty=1,type='h',col='light blue')
  }
  title(paste("ephemeral plants, roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days, lat/long ",longlat[2],",",longlat[1],sep=""),outer=TRUE)
  dev.off()
  
  filename<-paste("plant growth test output/correl roots ",DEP[root_deep]," cm regrow thresh ",growth_delay," days site ",ii,".pdf",sep="")
  pdf(filename,paper="A4",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(2,1)) # set up for 12 plots in 2 columns
  par(oma = c(2,2,2,2) + 1) # margin spacing stuff
  par(mar = c(3,3,1,1) + 1) # margin spacing stuff 
  par(mgp = c(3,1,0) ) # margin spacing stuff 
  
  
  grassmoist$date<-as.character(grassmoist$date1)
  merge_results_p<-merge(grassmoist,locustsub,by="date")
  merge_results_p<-subset(merge_results_p,PERENnew!='NA')
  pred<-aggregate(merge_results_p$moist,by=list(merge_results_p$date),FUN=max)
  obser<-aggregate(merge_results_p$PERENnew,by=list(merge_results_p$date),FUN=max)
  r_peren<-round(cor(round(pred$x),obser$x,method="spearman"),2)
  plot(pred$x~jitter(obser$x),col='blue',ylab='predicted greenness',xlab='observed greenness',xlim=c(0,11),ylim=c(0,11),main=paste("perennial plants, roots ",DEP[root_deep]," cm",sep=""))
  text(3,8,paste("r=",r_peren))
  
  merge_results_e<-merge(grassmoist,locustsub,by="date")
  merge_results_e<-subset(merge_results_e,Enew!='NA')
  pred<-aggregate(merge_results_e$moist,by=list(merge_results_e$date),FUN=max)
  obser<-aggregate(merge_results_e$Enew,by=list(merge_results_e$date),FUN=max)
  r_ephem<-round(cor(round(pred$x),obser$x,method="spearman"),2)
  plot(pred$x~jitter(obser$x),col='blue',ylab='predicted greenness',xlab='observed greenness',xlim=c(0,11),ylim=c(0,11),main=paste("ephemeral plants, roots ",DEP[root_deep]," cm",sep=""))
  text(3,8,paste("r=",r_ephem))
  dev.off()

if(ii==1 & kk==nodestart){
  all_e<-as.data.frame(cbind(DEP[root_deep],merge_results_e))
  all_p<-as.data.frame(cbind(DEP[root_deep],merge_results_p))
  all_r_e<-as.data.frame(cbind(ii,DEP[root_deep],r_ephem))
  all_r_p<-as.data.frame(cbind(ii,DEP[root_deep],r_peren))
}else{
  all_e<-rbind(all_e,as.data.frame(cbind(DEP[root_deep],merge_results_e)))
  all_p<-rbind(all_p,as.data.frame(cbind(DEP[root_deep],merge_results_p)))
  all_r_e<-rbind(all_r_e,as.data.frame(cbind(ii,DEP[root_deep],r_ephem)))
  all_r_p<-rbind(all_r_p,as.data.frame(cbind(ii,DEP[root_deep],r_peren)))
}  
  
} # end loop through grass model params



} # end loop through sites




# #plot for earlier years, with just text code for grass condition for now, per year
# for(yr in 1990:1997){
#   plotgrassmoist<-subset(grassmoist,as.numeric(format(grassmoist$date, "%Y"))==yr)
#   plot(plotgrassmoist$moist~plotgrassmoist$date,type='l',col='dark green',main=yr,ylim=c(0,11))
#   text(xlocustsub$DATE_,y=jitter(rep(5,nrow(locustsub)),amount=2),labels=dataset_ephem$OLDVEG,cex=0.5)
#   #text(x=dataset_ephem$DATE_,y=10,labels=dataset_ephem$OLDVEG,cex=0.5)
#   points(rainfall$rainfall/10~rainfall$dates,lty=2,type='h',col='light blue')
# }

# # now per month for a closer look
# mons<-c("01","02","03","04","05","06","07","08","09","10","11","12")
# #plot plant growth
# for(yr in 1990:1997){
#   for(mon in 1:12){
#     plotgrassmoist<-subset(grassmoist,format(grassmoist$date, "%Y-%m")==paste(yr,"-",mons[mon],sep=""))
#     plot(plotgrassmoist$moist~plotgrassmoist$date,type='l',col='dark green',main=yr,ylim=c(0,11))
#     text(x=dataset_ephem$DATE_,y=jitter(rep(5,nrow(dataset_ephem)),amount=2),labels=dataset_ephem$OLDVEG,cex=0.5)
#     #text(x=dataset_ephem$DATE_,y=10,labels=dataset_ephem$OLDVEG,cex=0.5)
#     points(rainfall$rainfall/10~rainfall$dates,lty=2,type='h',col='light blue')
#   }
# }






























write.csv(metout,'metout.csv')
write.csv(shadmet,'shadmet.csv')
write.csv(soil,'soil.csv')
write.csv(shadsoil,'shadsoil.csv')
write.csv(soilmoist,'soilmoist.csv')
write.csv(soilpot,'soilpot.csv')
write.csv(humid,'humid.csv')
write.csv(shadmoist,'shadmoist.csv')
write.csv(shadhumid,'shadhumid.csv')
write.csv(shadpot,'shadpot.csv')
write.csv(rainfall,'rainfall.csv')
write.csv(ectoin,'ectoin.csv')
write.csv(DEP,'DEP.csv')
write.csv(MAXSHADES,'MAXSHADES.csv')

tzone<-paste("Etc/GMT-",10,sep="") # doing it this way ignores daylight savings!
dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
soil<-cbind(dates,soil)
metout<-cbind(dates,metout)
shadsoil<-cbind(dates,shadsoil)
shadmet<-cbind(dates,shadmet)
soilmoist<-cbind(dates,soilmoist)
shadmoist<-cbind(dates,shadmoist)
humid<-cbind(dates,humid)
shadhumid<-cbind(dates,shadhumid)
soilpot<-cbind(dates,soilpot)
shadpot<-cbind(dates,shadpot)

dates2<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days") 
dates2<-subset(dates2, format(dates2, "%m/%d")!= "02/29") # remove leap years
rainfall<-as.data.frame(cbind(dates2,rainfall))
colnames(rainfall)<-c("dates","rainfall")

dstart<-as.POSIXct(as.Date('01/01/2000', "%d/%m/%Y"))-3600*11
dfinish<-as.POSIXct(as.Date('31/12/2000', "%d/%m/%Y"))-3600*10
plotsoilmoist<-subset(soilmoist,  soilmoist$dates > dstart & soilmoist$dates < dfinish )
plothumid<-subset(humid,  humid$dates > dstart & humid$dates < dfinish )
plotsoilpot<-subset(soilpot,  soilpot$dates > dstart & soilpot$dates < dfinish )
plotsoil<-subset(soil,  soil$dates > dstart & soil$dates < dfinish )

#plot(plotsoilmoist$dates, plotsoilmoist[,4]*100,type='l',col = "red",lty=1,ylim = c(0,50),ylab='moisture (% vol)',xlab='date')
plot(plotsoilmoist$dates, plotsoilmoist[,5]*100,type='l',col = 3,lty=1,ylim = c(0,50),ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,6]*100,type='l',col = 4,lty=1,ylim = c(0,50),ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,7]*100,type='l',col = 5,lty=1,ylim = c(0,50),ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,8]*100,type='l',col = 6,lty=1,ylim = c(0,50),ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,9]*100,type='l',col = 7,lty=1,ylim = c(0,50),ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,10]*100,type='l',col = 8,lty=1,ylim = c(0,50),ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,11]*100,type='l',col = 9,lty=1,ylim = c(0,50),ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,12]*100,type='l',col = 10,lty=1,ylim = c(0,50),ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,13]*100,type='l',col = 11,lty=1,ylim = c(0,50),ylab='moisture (% vol)',xlab='date')
points(rainfall$rainfall~rainfall$dates,type='h',col='dark blue')
# 
# plot(plothumid$dates, plothumid[,4]*100,type='l',col = "red",lty=1,ylim = c(0,100),ylab='relative humdity (%)',xlab='date')
# points(plothumid$dates, plothumid[,5]*100,type='l',col = 3,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plothumid$dates, plothumid[,6]*100,type='l',col = 4,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plothumid$dates, plothumid[,7]*100,type='l',col = 5,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plothumid$dates, plothumid[,8]*100,type='l',col = 6,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plothumid$dates, plothumid[,9]*100,type='l',col = 7,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plothumid$dates, plothumid[,10]*100,type='l',col = 8,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plothumid$dates, plothumid[,11]*100,type='l',col = 9,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plothumid$dates, plothumid[,12]*100,type='l',col = 10,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plothumid$dates, plothumid[,13]*100,type='l',col = 11,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# 
plot(plotsoilpot$dates, plotsoilpot[,4],type='l',col = "red",lty=1,ylim = c(-5000,0),ylab='water potential (J/kg)',xlab='date')
plot(plotsoilpot$dates, plotsoilpot[,5],type='l',col = 3,lty=1,ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,6],type='l',col = 4,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,7],type='l',col = 5,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,8],type='l',col = 6,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,9],type='l',col = 7,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,10],type='l',col = 8,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,11],type='l',col = 9,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,12],type='l',col = 10,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,13],type='l',col = 11,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# 
# plot(plotsoil$dates, plotsoil[,4],type='l',col = "red",lty=1,ylim = c(-10,80),ylab='temperature (C)',xlab='date')
# points(plotsoil$dates, plotsoil[,5],type='l',col = 3,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plotsoil$dates, plotsoil[,6],type='l',col = 4,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plotsoil$dates, plotsoil[,7],type='l',col = 5,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plotsoil$dates, plotsoil[,8],type='l',col = 6,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plotsoil$dates, plotsoil[,9],type='l',col = 7,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plotsoil$dates, plotsoil[,10],type='l',col = 8,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plotsoil$dates, plotsoil[,11],type='l',col = 9,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plotsoil$dates, plotsoil[,12],type='l',col = 10,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
# points(plotsoil$dates, plotsoil[,13],type='l',col = 11,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')