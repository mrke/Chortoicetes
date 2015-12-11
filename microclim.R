
microclim<-function(longlat=c(-89.40123,43.07305),sitemethod=0,ystart=1990,yfinish=1990,timezone=0,EC=0.0167238,
  rainfrac=0.5,densfun=c(0,0),writecsv=0,tides=matrix(data = 0., nrow = 24*dim, ncol = 3),shore=0,
  L=c(0,0,8.18990859,7.991299442,7.796891252,7.420411664,7.059944542,6.385001059,5.768074989,4.816673431,4.0121088,1.833554792,0.946862989,0.635260544,0.804575,0.43525621,0.366052856,0,0)*10000,
  LAI=0.1,evenrain=0,runmoist=0,maxpool=10000,PE=rep(1.1,19),KS=rep(0.0037,19),BB=rep(4.5,19),
  BD=rep(1.3,19),loc="Madison Wisconsin",timeinterval=12,nyears=1,RUF=0.004,SLE=0.95,ERR=1.5,
  DEP=c(0., 2.5,  5.,  10.,  15.,  20.,  30.,  50.,  100.,  200.),
  Thcond=2.5,Density=2560,SpecHeat=870,BulkDensity=1300,Clay=20,CMH2O=1.,
  TIMAXS=c(1.0, 1.0, 0.0, 0.0),TIMINS=c(0.0, 0.0, 1.0, 1.0),minshade=0,maxshade=90,Usrhyt=1,
  REFL=0.15,slope=0,aspect=0,hori=rep(0,24),rungads=1,cap=1,write_input=0,spatial="c:/Australian Environment/",
  snowmodel=0,snowtemp=1.5,snowdens=0.375,snowmelt=0.9,undercatch=1,rainmult=1,
  rainmelt=0.0125,runshade=1,mac=0,PCTWET=0,
  SoilMoist_Init=c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3),soiltype=4, ...){
######################### model modes ###########################################################
mac<-0 # choose mac (1) or pc (0) 
writecsv<-0 # make Fortran code write output as csv files
write_input<-0 # write csv files of final input to working directory? 1=yes, 0=no.
runshade<-1 # run the model twice, once for each shade level (1) or just for the first shade level (0)?
runmoist<-1 # run soil moisture model (0=no, 1=yes)?
snowmodel<-0 # run the snow model (0=no, 1=yes)? - note that this runs slower
shore<-0 # include tide effects (if 0, an empty matrix of tide effects is created)
rungads<-1 # use the Global Aerosol Database?
#########################################################################################################

############## location and climatic data  ###################################
spatial<-"c:/Australian Environment/" # place where climate input files are kept
sitemethod <- 0 # 0=specified single site long/lat, 1=place name search using geodis (needs internet)
#longlat<-c(survey_freq[ii,2],survey_freq[ii,3]) # type a long/lat here in decimal degrees
loc <- "Broken Hill, Australia" # type in a location here, used if option 1 is chosen above
timezone<-0 # if timezone=1 (needs internet), uses GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)
dailywind<-1 # use daily windspeed database?
terrain<-0 # include terrain (slope, aspect, horizon angles) (1) or not (0)?
soildata<-1 # include soil data for Australia (1) or not (0)?
#ystart <- 1990# start year for weather generator calibration dataset or AWAP database
#yfinish <- 2009# end year for weather generator calibration dataset
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
  
Thcond <- rep(2.5,10) # soil minerals thermal conductivity (W/mC)
Density <- rep(2560.,10) # soil minerals density (kg/m3)
SpecHeat <- rep(870.,10) # soil minerals specific heat (J/kg-K)
BulkDensity <- BD[seq(1,19,2)]*1000#rep(1360,10) # soil bulk density (kg/m3)
cap<-1 # organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)
SatWater <- rep(0.26,10) # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
Clay <- rep(22,10) # clay content for matric potential calculations (%)
SoilMoist <- 0.2 # fractional soil moisture (decimal %)
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
snowdens<-0.325 # snow density (mg/m3)
densfun<-c(0.001369,0.1095) # slope and intercept of linear model of snow density as a function of day of year - if it is c(0,0) then fixed density used
snowmelt<-0.9 # proportion of calculated snowmelt that doesn't refreeze
undercatch<-1.0 # undercatch multipier for converting rainfall to snow
rainmelt<-0.0125#85 # paramter in equation that melts snow with rainfall as a function of air temp, start with 0.0125
write_input<-0 # write csv files of final input to working directory? 1=yes, 0=no.
warm<-0 # uniform warming of air temperature input to simulate climate change
loop<-0 # if doing multiple years, this shifts the starting year by the integer value

# run the model
niche<-list(sitemethod=sitemethod,L=L,LAI=LAI,SoilMoist_Init=SoilMoist_Init,evenrain=evenrain,runmoist=runmoist,maxpool=maxpool,PE=PE,KS=KS,BB=BB,BD=BD,loop=loop,warm=warm,rainwet=rainwet,manualshade=manualshade,dailywind=dailywind,terrain=terrain,soildata=soildata,loc=loc,ystart=ystart,yfinish=yfinish,nyears=nyears,RUF=RUF,SLE=SLE,ERR=ERR,DEP=DEP,Thcond=Thcond,Density=Density,SpecHeat=SpecHeat,BulkDensity=BulkDensity,Clay=Clay,SatWater=SatWater,SoilMoist=SoilMoist,CMH2O=CMH2O,TIMAXS=TIMAXS,TIMINS=TIMINS,minshade=minshade,maxshade=maxshade,Usrhyt=Usrhyt,REFL=REFL,slope=slope,aspect=aspect,hori=hori,rungads=rungads,cap=cap,write_input=write_input,spatial=spatial,snowmodel=snowmodel,snowtemp=snowtemp,snowdens=snowdens,snowmelt=snowmelt,undercatch=undercatch,rainmelt=rainmelt,rainmult=rainmult,runshade=runshade)
setwd('../micro_australia/')
source('NicheMapR_Setup_micro.R')
nicheout<-NicheMapR(sitemethod=sitemethod,L=L,LAI=LAI,SoilMoist_Init=SoilMoist_Init,evenrain=evenrain,runmoist=runmoist,maxpool=maxpool,PE=PE,KS=KS,BB=BB,BD=BD,loop=loop,warm=warm,rainwet=rainwet,manualshade=manualshade,dailywind=dailywind,terrain=terrain,soildata=soildata,loc=loc,ystart=ystart,yfinish=yfinish,nyears=nyears,RUF=RUF,SLE=SLE,ERR=ERR,DEP=DEP,Thcond=Thcond,Density=Density,SpecHeat=SpecHeat,BulkDensity=BulkDensity,Clay=Clay,SatWater=SatWater,SoilMoist=SoilMoist,CMH2O=CMH2O,TIMAXS=TIMAXS,TIMINS=TIMINS,minshade=minshade,maxshade=maxshade,Usrhyt=Usrhyt,REFL=REFL,slope=slope,aspect=aspect,hori=hori,rungads=rungads,cap=cap,write_input=write_input,spatial=spatial,snowmodel=snowmodel,snowtemp=snowtemp,snowdens=snowdens,snowmelt=snowmelt,undercatch=undercatch,rainmelt=rainmelt,rainmult=rainmult,runshade=runshade)
setwd(maindir)

# get output
metout<-as.matrix(nicheout$metout) # above ground microclimatic conditions, min shade
shadmet<-as.matrix(nicheout$shadmet) # above ground microclimatic conditions, max shade
soil<-as.matrix(nicheout$soil) # soil temperatures, minimum shade
shadsoil<-as.matrix(nicheout$shadsoil) # soil temperatures, maximum shade
soilmoist<-as.matrix(nicheout$soilmoist) # soil water content, minimum shade
shadmoist<-as.matrix(nicheout$shadmoist) # soil water content, maximum shade
humid<-as.matrix(nicheout$humid) # soil humidity, minimum shade
shadhumid<-as.matrix(nicheout$shadhumid) # soil humidity, maximum shade
soilpot<-as.matrix(nicheout$soilpot) # soil water potential, minimum shade
shadpot<-as.matrix(nicheout$shadpot) # soil water potential, maximum shade
rainfall<-as.matrix(nicheout$RAINFALL)
MAXSHADES<-as.matrix(nicheout$MAXSHADES)
elev<-as.numeric(nicheout$ALTT)
REFL<-as.numeric(nicheout$REFL)
longlat<-as.matrix(nicheout$longlat)
fieldcap<-as.numeric(nicheout$fieldcap)
wilting<-as.numeric(nicheout$wilting)
ectoin<-c(elev,REFL,longlat,fieldcap,wilting,ystart,yfinish)

# write.csv(metout,'microclimate/metout.csv')
# write.csv(shadmet,'microclimate/shadmet.csv')
# write.csv(soil,'microclimate/soil.csv')
# write.csv(shadsoil,'microclimate/shadsoil.csv')
# write.csv(soilmoist,'microclimate/soilmoist.csv')
# write.csv(soilpot,'microclimate/soilpot.csv')
# write.csv(humid,'microclimate/humid.csv')
# write.csv(shadmoist,'microclimate/shadmoist.csv')
# write.csv(shadhumid,'microclimate/shadhumid.csv')
# write.csv(shadpot,'microclimate/shadpot.csv')
# write.csv(rainfall,'microclimate/rainfall.csv')
# write.csv(ectoin,'microclimate/ectoin.csv')
# write.csv(DEP,'microclimate/DEP.csv')
# write.csv(MAXSHADES,'microclimate/MAXSHADES.csv')

return(list(metout=metout,shadmet=shadmet,soil=soil,shadsoil=shadsoil,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,MAXSHADES=MAXSHADES,elev=elev,REFL=REFL,fieldcap=fieldcap,wilting=wilting,ectoin=ectoin,rainfall=rainfall))
}