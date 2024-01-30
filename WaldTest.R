wald.test.binary <- function(x0,x1, measure="sd"){
  
  p0 <- mean(x0); p1 <- mean(x1)
  q0 <- 1-p0 ; q1 <- 1-p1
  sdx0 <- p0*q0; sdx1 <- p1*(1-p1)
  n0 <- length(x0); n1 <- length(x1)
  
  if(measure=="sd"){
    if((sdx0==0 && sdx1==0)){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    } else {
      Z <- (p0-p1)/sqrt(sdx0/n0+sdx1/n1)
    }
  } else if(measure=="rr") {
    if(p0==1 || (p0==0 && p1*q1==0)){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    }
    Z <- (q1/q0-1)/sqrt(p0/(n0*q0^3)+p1*q1/n1)
  } else if(measure=="or"){
    if((q0*p1==1) || (p0==1 && p1==1) || (p0==0 && p1==0)){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    }
    else{
      Z <- (p0*q1/(q0*p1)-1)/sqrt(p0/(n0*q0^3)+p1/(n1*q1^3))
    }
  } else if(measure=="lrr"){
    if(p0==1 || p1==1 || (p0==0 && p1==0) ){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    } else {
      Z <- log(q1/q0)/ sqrt(p0/(n0*q0)+p1/(n1*q1))
    }
  } else if(measure=="lor"){
    if(p0==0 || q0==0 || p1==0 || q1==0){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    } else {
      Z <- log(q1/q0)/ sqrt(p0/(n0*q0)+p1/(n1*q1))
    }
  }
  return(2*pnorm(-1*abs(Z)))
}

wald.test <- function(x0,x1, measure="sd"){
  p0 <- mean(x0); p1 <- mean(x1)
  sdx0 <- sd(x0); sdx1 <- sd(x1)
  n0 <- length(x0); n1 <- length(x1)
  
  
  if((sdx0==0 && sdx1==0)){
    if(p0==p1){Z = 0}
    if(p0<p1){Z = -Inf}
    if(p0>p1){Z = Inf}
  } else {
    if(measure=="sd"){
      Z <- (p0-p1)/sqrt(sdx0/n0+sdx1/n1)
    }
  }
  return(2*pnorm(-1*abs(Z)))
}