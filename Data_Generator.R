generate_data <- function(N=200, dist="bern", parV = 0.5){
  if(dist=="bern"){
    return(rbinom(N,1, parV))
  }
  if(dist=="norm"){
    return(rnorm(n=N,mean=parV[1],sd=parV[2]))
  }
  if(dist=="Lickert"){
    return(round(rbeta(N, 5, parV)*5)+1)
  }
  if(dist=="expon"){
    return(rexp(n = N, rate = parV))
  }
  if(dist=="beta"){
    return(rbeta(N, 5, parV))
  }
}

superior <- function(par1=c(0,0), par2=FALSE, dist="bern"){
  if(dist=="bern"){
    if(par1[1]<par1[2]){
      return(1)
    }
    if(par1[1]>par1[2]){
      return(0)
    }
    if(par1[1]==par1[2]){
      return(2)
    }
  }
  if(dist=="norm"){
    if(par1[1]<par1[2]){
      return(1)
    }
    if(par1[1]>par1[2]){
      return(0)
    }
    if(par1[1]==par1[2]){
      return(2)
    }
  }
  if(dist=="Lickert" || dist=="beta"){
    exp1 <- 5/(5+par1[1])
    exp2 <- 5/(5+par1[2])
    if(exp1<exp2){
      return(1)
    }
    if(exp1>exp2){
      return(0)
    }
    if(exp1==exp2){
      return(2)
    }
  }
  if(dist=="expon"){
    if(par1[1]>par1[2]){
      return(1)
    }
    if(par1[1]<par1[2]){
      return(0)
    }
    if(par1[1]==par1[2]){
      return(2)
    }
  }
}