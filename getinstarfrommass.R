getinstarfrommass<-function(mass){
  if(mass<2.88){
    instar<-"N1"
  }else if(mass<5.56){
    instar<-"N2"
  }else if(mass<11.71){
    instar<-"N3"
  }else if(mass<25.53){
    instar<-"N4"
  }else if(mass<55.){
    instar<-"N5"
  }else if(mass>=55.){
    instar<-"Adult"
  }
  return(instar)
}
