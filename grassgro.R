plantgro<-function(soilpot=soilpot,soilmoist=soilmoist, root_shallow=4, root_deep=8, growth_delay=1, wilting_thresh=-200, permanent_wilting_point=-1500, FoodWater=82){
  
  meanpot<-as.data.frame(soilpot)[,((root_shallow+3):(3+root_deep))] # get range of soil water potential depths to take mean of
  meanpot<-apply(meanpot, 1, mean) # get average soil water potential across chosen depth range
  
  grow<-meanpot
  grow[grow>permanent_wilting_point]<-1 # find times above the PWP (growth possible)
  grow[grow<=permanent_wilting_point]<-0 # find times when below the PWP (plant dead)
  counter<-0
  grow2<-grow*0 # create empty vector
  for(j in 1:length(grow)){ # accumulate runs of days above the permanent wilting point (growth possible)
    if(j==1){ # check, if first hour, whether we're starting with a growth day or not
      if(grow[j]==1){
        counter<-counter+1 # if so, increment the counter
      }
      grow2[j]<-counter
    }else{ # otherwised, onyl increment the counter if growth on both present hour and the hour just past
      if(grow[j-1]>0 & grow[j]==1){
        counter<-counter+1
      }else{
        counter<-0
      }
      grow2[j]<-counter
    }
  }
  grow3<-grow2
  grow3[grow3<growth_delay]<-0 # apply growth delay specified by the user for time required for plants to come back after PWP has been hit
  grow3[grow3>0]<-1 # make vector of 0 and 1 where 1 means plants could have come back from drought
  
  meanmoist<-as.data.frame(soilmoist)[,((root_shallow+3):(3+root_deep))]  # get range of soil water moisture depths to take mean of
  meanmoist<-apply(meanmoist, 1, mean)
  
  mean.moist.pot<-as.data.frame(cbind(meanpot,meanmoist))
  colnames(mean.moist.pot)<-c('pot','moist')
  mean.moist.pot$pot[mean.moist.pot$pot>wilting_thresh]<-FoodWater # assume plants start wilting at about 2 bar, but above this they are at max water content
  mean.moist.pot$moist<-mean.moist.pot$moist*100 # convert to percent
  potmult<-mean.moist.pot$pot
  potmult[potmult!=82]<-0
  potmult[potmult!=0]<-1  
  wilting<-subset(mean.moist.pot,pot==FoodWater) # find soil moisture range corresponding to values above the wilting point
  wilting<-min(wilting$moist) # get the min soil moisture at which plants aren't wilting
  pct.water<-mean.moist.pot$moist
  pct.water[pct.water>wilting]<-FoodWater # now have vector of either max plant water content or soil moisture content - need to convert the latter into a smooth decline to zero from max value
  minmoist<-0
  pct.water[pct.water<FoodWater]<-(pct.water[pct.water<FoodWater]-minmoist)/(wilting-minmoist)*FoodWater # for just the values less than max water content, make them equal to the ratio of moisture level and wilting point moisture and then multiplied by food water content so have moisture content ranging from max level down to min moisture possible
  pct.water<-pct.water/100*grow3 # now convert to proportion and cut out times below PWP (including regrowth penalty)
  plantmoist<-as.data.frame(cbind(dates,mean.moist.pot))
  colnames(plantmoist)<-c('date1','moist')
  plantmoist$date1<-dates
  plantmoist$moist<-plantmoist$moist/max(plantmoist$moist)*11 # put in units scaling from 0-11
  # next four lines spread the values out more evenly over 11 categories
  minval<-min(plantmoist$moist[plantmoist$moist!=0])
  plantmoist$moist[plantmoist$moist==0]<-minval
  plantmoist$moist<-plantmoist$moist-minval
  plantmoist$moist<-plantmoist$moist/max(plantmoist$moist)*11 # put in units scaling from 0-11
  plantmoist2<-plantmoist
  plantmoist2$moist[plantmoist2$moist>0]<-1

  plant.pres<-plantmoist2$moist
  plant.moist<-subset(plantmoist, format(plantmoist$date1,"%H")=="12")
 
  return(list(plant.pres=plant.pres,plant.moist=plant.moist))  
}