## make animation of plague locust through time

require(foreign) # read/write dbf
library(ggplot2) # Plotting maps.
library(maps)    # Map data.
library(oz)      # Map data for Australia.
library(scales)  # Functions: alpha() comma()
library(ggmap)   # Google maps.
library(zoom)
library(rgeos)   # for gContains(...)


# read australia map data
ds <- read.csv(file.path("APLC data/", "ozdata.csv"))
# ds  <- subset(ds, border == "coast") # dont keep state lines
mapdata<-ds[,c(2,3,6)] # long, lat, state
mapdata$ausmap <- rep(TRUE,nrow(mapdata)) # ausmap
mapdata$NDENS <- rep(-1,nrow(mapdata)) # nymph density
mapdata$ADENS <- rep(NA,nrow(mapdata)) # adult density
mapdata$DATE_ <- rep(NA,nrow(mapdata)) # sampling date
mapdata$OLDVEG     <- rep(NA,nrow(mapdata))  # old veg condition index (does not separate between perenial and ephem)
mapdata$E     <- rep(NA,nrow(mapdata))  # ephemeral plant greeness
mapdata$PEREN <- rep(NA,nrow(mapdata)) # perennial plant greeness
mapdata<-subset(mapdata, state != "TAS")

# Load locust survey data
surveylist <- c("survey90-92.dbf", "survey93-95.dbf","surveys96-98.dbf", "survey98-01.dbf", "survey02-05.dbf" , "survey05-09.dbf")
datalist <- list()
datalength <- 0
for(i in 1:length(surveylist)){
  dat <- read.dbf(file.path("APLC data/",surveylist[i]))
  if(i==1){
    datalist[length(datalist)+1]<-list(data.frame(dat$X_COORD, dat$Y_COORD, 
                                                  rep(NA,nrow(dat)),rep(FALSE,nrow(dat)),
                                                  dat$NDENS, dat$ADENS, dat$DATE_, 
                                                  dat$CONDITION, rep(NA,nrow(dat)),rep(NA,nrow(dat)),rep(NA,nrow(dat)),rep(NA,nrow(dat))))
 
  }
  if(i==2){
    datalist[length(datalist)+1]<-list(data.frame(dat$X_COORD, dat$Y_COORD, 
                                                  rep(NA,nrow(dat)),rep(FALSE,nrow(dat)),
                                                  dat$NDENS, dat$ADENS, dat$DATE_, 
                                                  dat$CONDITION, rep(NA,nrow(dat)),rep(NA,nrow(dat)),rep(NA,nrow(dat)),rep(NA,nrow(dat))))
    
  }
  if(i==3){
    datalist[length(datalist)+1]<-list(data.frame(dat$X_COORD, dat$Y_COORD, 
                                                  rep(NA,nrow(dat)),rep(FALSE,nrow(dat)),
                                                  dat$NDENS, dat$ADENS, dat$DATE_, 
                                                  dat$CONDITION, rep(NA,nrow(dat)),rep(NA,nrow(dat)),rep(NA,nrow(dat)),rep(NA,nrow(dat))))
    
    
  }
  if(i==4){
    datalist[length(datalist)+1]<-list(data.frame(dat$X_COORD, dat$Y_COORD, 
                                                  rep(NA,nrow(dat)),rep(FALSE,nrow(dat)),
                                                  dat$NDENS, dat$ADENS, dat$DATE_, 
                                                  rep(NA,nrow(dat)),E=dat$E, PEREN=dat$P,
                                                  Enew=rep(NA,nrow(dat)),PERENnew=rep(NA,nrow(dat))))
#     word1 <- unique(datalist[[i]]$OLDVEG)
#     word2 <- unique(datalist[[i]]$E)
#     word3 <- unique(datalist[[i]]$PEREN)
#     words<-unique(word1,word2,word3)
    # code grass labels as numeric
  words<-c('NA', 'GS', 'LG', 'G', 'DO', 'DGB','GT','D','VD','DGB-LG','GT-LG','D_LG')
  key  <-c(NA,10,8,11,9,7,5,2,1,6,4,3)
    
    # Replace vegetation code with numerical approximation
    for(j in 1:nrow(dat)){
      Eno<-which(words == datalist[[i]]$E[j])
      Pno<-which(words == datalist[[i]]$PEREN[j])
      if(length(Eno)==0){# catch integer(0) 
        datalist[[i]]$Enew[j]<-NA
      }else{
        datalist[[i]]$Enew[j]<-key[Eno]
      }
      if(length(Pno)==0){# catch integer(0) 
        datalist[[i]]$PERENnew[j]<-NA
      }else{
        datalist[[i]]$PERENnew[j]<-key[Pno]
      }
    }
  }
  if(i==5){
    datalist[length(datalist)+1]<-list(data.frame(dat$X_COORD, dat$Y_COORD, 
                                                  rep(NA,nrow(dat)),rep(FALSE,nrow(dat)),
                                                  dat$NDENS, dat$ADENS, dat$DATE_, 
                                                  rep(NA,nrow(dat)),E=dat$E, PEREN=dat$P,
                                                  Enew=rep(NA,nrow(dat)),PERENnew=rep(NA,nrow(dat))))
  #         word1 <- unique(datalist[[i]]$OLDVEG)
  #         word2 <- unique(datalist[[i]]$E)
  #         word3 <- unique(datalist[[i]]$PEREN)
  #         words<-unique(word1,word2,word3)
  
  # code grass labels as numeric
  words<-c('NA', 'GS', 'LG', 'G', 'DO', 'DGB','GT','D','VD','DGB-LG','GT-LG','D_LG')
  key  <-c(NA,10,8,11,9,7,5,2,1,6,4,3)

  # Replace vegetation code with numerical approximation
  for(j in 1:nrow(dat)){
    Eno<-which(words == datalist[[i]]$E[j])
    Pno<-which(words == datalist[[i]]$PEREN[j])
    if(length(Eno)==0){# catch integer(0) 
      datalist[[i]]$Enew[j]<-NA
    }else{
      datalist[[i]]$Enew[j]<-key[Eno]
    }
    if(length(Pno)==0){# catch integer(0) 
      datalist[[i]]$PERENnew[j]<-NA
    }else{
      datalist[[i]]$PERENnew[j]<-key[Pno]
    }
  }
  }
  if(i==6){
    datalist[length(datalist)+1]<-list(data.frame(dat$POINT_X, dat$POINT_Y, 
                                                  rep(NA,nrow(dat)),rep(FALSE,nrow(dat)),
                                                  dat$NDENS, dat$ADENS, dat$DATE_, 
                                                  rep(FALSE,nrow(dat)),E=dat$E, PEREN=dat$PEREN,
                                                  Enew=rep(NA,nrow(dat)),PERENnew=rep(NA,nrow(dat))))
    # code grass labels as numeric
  words<-c('NA', 'GS', 'LG', 'G', 'DO', 'DGB','GT','D','VD','DGB-LG','GT-LG','D_LG')
  key  <-c(NA,10,8,11,9,7,5,2,1,6,4,3)
    
    # Replace vegetation code with numerical approximation
    for(j in 1:nrow(dat)){
      Eno<-which(words == datalist[[i]]$E[j])
      Pno<-which(words == datalist[[i]]$PEREN[j])
      if(length(Eno)==0){# catch integer(0) 
        datalist[[i]]$Enew[j]<-NA
      }else{
        datalist[[i]]$Enew[j]<-key[Eno]
      }
      if(length(Pno)==0){# catch integer(0) 
        datalist[[i]]$PERENnew[j]<-NA
      }else{
        datalist[[i]]$PERENnew[j]<-key[Pno]
      }
    }
    
  }
  names(datalist[[i]])<-c("long",   "lat",    "state",  "ausmap", "NDENS",  "ADENS",  "DATE_",  
                          "OLDVEG", "E",      "PEREN", "Enew", "PERENnew" )
  datalength<-datalength + nrow(dat) # sum of all data rows
}

# Show sample of coded veg data for "survey98-01.dbf", "survey02-05.dbf" , "survey05-09.dbf"
# E (ephemeral veg) and PEREN (perenial veg) are coded by APLC while Enew and PERENnew are 
# numerical conversions specified by Mike
head(datalist[[4]])
head(datalist[[5]])
head(datalist[[6]])

# # stitch together data
# locustdata<-rep(NA,datalength*ncol(mapdata))
# dim(locustdata)<-c(datalength,ncol(mapdata))
# start<-1
# for(i in 1:length(datalist)){
#   end<-start+nrow(datalist[[i]])
#   locustdata[start:(end-1),]<-datalist[[i]]
#   start <- end
# }


locustdata<-datalist[[1]]
for(i in 2:6){
locustdata<-rbind(locustdata,datalist[[i]])
}

locustdata<-as.data.frame(locustdata)
locustdata$DATE_<-as.Date(locustdata$DATE_, ,origin="1970-01-01")

# get rid of anomolous data long/lat = 1 = -54
locustdata<-subset(locustdata, lat!=-54)
locustdata<-subset(locustdata, long!=1)


# round long lat measurements to nearest 0.5 degrees 
locustdata1 <- as.data.frame(locustdata)
locustdata1$lat <- round(locustdata$lat*2,0)/2 # round to nearest 0.5
locustdata1$long <- round(locustdata$long*2,0)/2
surveyfreq<-as.data.frame(table(locustdata1$long,locustdata1$lat))
names(surveyfreq)<-c('long','lat', 'freq')
surveyfreq$lat<-as.numeric(as.character(surveyfreq$lat))
surveyfreq$long<-as.numeric(as.character(surveyfreq$long))
surveyfreq$freq<-as.numeric(as.character(surveyfreq$freq))

# plot survey distribution all data
p <- ggplot() + coord_fixed()
basemap <- p + geom_polygon(data=mapdata, aes(x=long,y=lat,group=state, fill = 'red'), colour = "white")
p <- basemap + geom_point(data=surveyfreq[which(surveyfreq$freq>0),],
                          aes(x=long,
                              y=lat,
                              size = freq))
p+ggtitle("1990-2009")

# order long/lats by survey freq
attach(surveyfreq)
sforder <-surveyfreq[order(-freq),]
detach(surveyfreq)
write.csv(sforder,'survey_freq_90-09.csv')
# print top surveyed coordinates
head(sforder)

# get subset for most sampled longlat 
attach(locustdata1)
locustsub<-locustdata1[which((long == sforder[1,1])&(lat==sforder[1,2])),]
detach(locustdata1)

plot(locustsub$DATE_, locustsub$NDENS, type='h')









locustdata98 <- as.data.frame(locustdata)
locustdata98$lat <- round(locustdata$lat*2,0)/2 # round to nearest 0.5
locustdata98$long <- round(locustdata$long*2,0)/2
locustdata98<-subset(locustdata98,DATE_ > as.Date("1998-01-01", ,origin="1970-01-01")) # just 1998 onwards, when veg estimates are cleaner
surveyfreq98<-as.data.frame(table(locustdata98$long,locustdata98$lat))
names(surveyfreq98)<-c('long','lat', 'freq')
surveyfreq98$lat<-as.numeric(as.character(surveyfreq98$lat))
surveyfreq98$long<-as.numeric(as.character(surveyfreq98$long))
surveyfreq98$freq<-as.numeric(as.character(surveyfreq98$freq))



# plot survey distribution 1998-2009 data
p2 <- ggplot() + coord_fixed()
basemap <- p2 + geom_polygon(data=mapdata, aes(x=long,y=lat,group=state, fill = 'red'), colour = "white")
p2 <- basemap + geom_point(data=surveyfreq98[which(surveyfreq98$freq>0),],
                          aes(x=long,
                              y=lat,
                              size = freq))
p2+ggtitle("1998-2009")

# order long/lats by survey freq
attach(surveyfreq98)
sforder98 <-surveyfreq98[order(-freq),]
detach(surveyfreq98)
write.csv(sforder98,'survey_freq_98-09.csv')
# print top surveyed coordinates
head(sforder98)

# get subset for most sampled longlat 
attach(locustdata98)
locustsub98<-locustdata98[which((long == sforder98[1,1])&(lat==sforder98[1,2])),]
detach(locustdata98)

plot(locustsub98$DATE_, locustsub98$NDENS, type='h')
points(locustsub98$DATE_, locustsub98$Enew/3, col='green',type='p')


