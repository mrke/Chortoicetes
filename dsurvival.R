dsurvival<-function(pars_survive,survival, temp){ # units suvival prob/day
  # change in proportion surviving as a function of fitted pars, proportion surviving, and temp C  
  temp[temp>39]<-39
  temp[temp<29.5]<-29.5
  dy <- -survival*(      pars_survive[2]*temp**1+
                           pars_survive[3]*temp**2+
                           pars_survive[4]*temp**3+
                           pars_survive[5]*temp**4)
  
  return(dy)}