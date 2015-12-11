# Temp dependent function for C. terminifera egg dev rate fitted from Gregg (1984) egg dev. rates 
devrate <- function(temp){
  if(temp<=32){
    y = -12334*(1/(temp+273)) + 38.027 
    return(exp(y)/24) # divide by 24 to get units in 1/h
  } else{
    return(devrate(32))
  }
}