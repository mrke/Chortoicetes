dmass <- function(pars_grow,mass, temp){ #units of mg dry weight/day!
  if(temp>39){
    temp<-39
  }
  if(temp<25.9 & temp>20){
    temp<-25.9
  }
  # change in mass, as a function of fitted parameters, mass, and temp C
  dm<- mass*(pars_grow[2]+pars_grow[3]*temp**1+
               pars_grow[4]*temp**2+
               pars_grow[5]*temp**3+
               pars_grow[6]*temp**4+
               pars_grow[7]*temp**5)
  if(temp<=20.0){
    dm<-0
  }
  return(dm)
}
