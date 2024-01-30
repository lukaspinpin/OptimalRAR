library(foreach)
library(doParallel)
library(doRNG)
library(Rfast)

simulation <- function(N=200, dist="bern", K=1, par1 = c(0,0),  par2 = NULL, burnin=3, ar=NULL, nsim= 10^4, method="ER", measure="sd"){
  
  #Generate Result Data Frame
  colnames <-  c()
  for(i in 1:ncol(par1)){
    colnames <- c(colnames, paste("p1", i, sep = "", collapse = NULL))
  }
  if(!is.null(par2)){
    for(i in 1:ncol(par2)){
      colnames <- c(colnames, paste("p2", i, sep = "", collapse = NULL))
    }
  }
  for(j in 1:(K+1)){
    colnames <- c(colnames, paste( "n", j-1, sep = "", collapse = NULL))
    colnames <- c(colnames, paste("Var(n", j-1, ")",sep = "", collapse = NULL))
  }
  colnames <- c(colnames, "EMR", "Var_EMR", "Z_Type1", "BM_Type1", "Z_Power", "BM_Power", "Per_Superior", "Var_Superior")
  if(!is.null(ar)){
    colnames <- c(colnames, ar)
    colnames <- c(colnames, paste("BIAS(", ar, ")",sep = "", collapse = NULL))
  }
  result <- data.frame(matrix(ncol = length(colnames), nrow = max(nrow(par1),nrow(par2))))
  colnames(result) <- colnames
  
  cl <- parallel::makeCluster(4)  #increase according to the limitations of your cluster if you run it on the cluster
  doParallel::registerDoParallel(cl)
  
  #2-Armed Or Multi-armed
  ######  2  Armed  ##########################################################
  if(K==1){
    ########## Theoretical Variance Calculation ############################
    if(dist=="bern"){
      mean <- par1
      fun <- function(x){x*(1-x)}
      variance <-  data.frame(lapply(par1,fun))
      theta <- par1[,2]*(1-par1[,1]) + 0.5*(par1[,1]*par1[,2]+(1-par1[,1])*(1-par1[,2]))
      psi <- 0.5*theta +0.25
    }
    if(dist=="norm"){
      mean <- par1
      fun <- function(x){x^2}
      variance <-  data.frame(lapply(par2,fun))
      xn <- (mean[,1] - mean[,2])/sqrt(variance[,1]+variance[,2])
      theta <-  1 - dnorm(xn)
      psi <- 0.5*theta +0.25
    }
    if(dist=="expon"){
      fun <- function(x){1/x}
      mean <- data.frame(lapply(par1,fun))
      fun <- function(x){1/(x^2)}
      variance <-  data.frame(lapply(par1,fun))
      theta <-  par1[,1]/(par1[,1]+par1[,2])
      psi <- 0.5*theta +0.25
    }
    if(dist=="beta"){
      fun <- function(x){5/(x+5)}
      mean <- data.frame(lapply(par1,fun))
      fun <- function(x){5*x/((5+x)^2+(5+x+1))}
      variance <-  data.frame(lapply(par1,fun))
      theta <- c(0.5,0.5631895, 0.6317465, 0.704089, 0.7777655, 0.813794,0.8486965,
                 0.881401, 0.91132, 0.937633, 0.959517, 0.976518, 0.988456, 0.991853,
                  0.994524, 0.995621, 0.997357, 0.998557, 0.998991, 0.999326, 0.999584)
      psi <- 0.5*theta +0.25
      psi <- c(rep(0.5, length(theta)),theta[-1])
    }
    if(dist=="Lickert"){
      ## Variance needed to be estimated through simulations in extra Code
      mean <-  read.csv("Lickert_Mean.csv", header = TRUE)
      mean <- mean[,2:3]
      variance <-  read.csv("Lickert_Variance.csv", header = TRUE)
      variance <- variance[,2:3]
      theta <- c(0.5, 0.544, 0.593, 0.6446757, 0.696642, 0.7469708, 0.796875, 0.818026,
                  0.840608, 0.852698, 0.8653635, 0.878488, 0.892124, 0.906136, 0.92029,
                  0.93435, 0.9479362, 0.9606685, 0.9720625, 0.9817345, 0.9893505,
                  0.9947455, 0.998026, 0.999626)
      psi <- 0.5*theta +0.25
      psi <- c(rep(0.5, length(psi)),psi[-1])
    }
    ############## ER ######################################################
    if(method == "ER" ){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.ER <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%dorng% 
        
        result_part <-result_part_emp
        
        source('ER.R')
        p <- as.numeric(par1[i,])
        start <- length(p)
        result_part[1,1:start] <- p
        if(is.null(par2)){
          sim.ER <- sim_ER(N=N, dist=dist, K=1, par1 = p,  burnin= burnin, nsim=nsim, measure=measure)
        } else{
          sig <- as.numeric(par2[i,])
          start <- start+ length(sig)
          result_part[1,1:start] <- c(p, sig)
          sim.ER <- sim_ER(N=N, dist=dist, K=1, par1 = p, burnin= burnin, par2=sig, nsim=nsim)
        }
        
        
        result_part[1,(start+1)] <- mean(sim.ER[,5]); result_part[1,(start+2)] <- var(sim.ER[,5]) 
        result_part[1,(start+3)] <- mean(sim.ER[,6]); result_part[1,(start+4)] <- var(sim.ER[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.ER[,1])
        result_part[1,(start+6)] <- var(sim.ER[,1])
        
        reject_Z <- sum(sim.ER[,2] < 0.05)
        reject_BM <- sum(sim.ER[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.ER[,4])
        result_part[1,(start+12)] <- var(sim.ER[,4])
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.ER
    }
    ############## Sequential ######################################################
    if(method == "SMLE_R_minF"){  ### change
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.SMLE_R <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%dorng% 
        source('SMLE.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        start <- length(p)
        result_part[1,1:start] <- p
        sim.SMLE_R <- sim_SMLE(N=N, dist=dist, K=1, par1 = p,  nsim=nsim, measure=measure, ar=ar, burnin=burnin)
        result_part[1,(start+1)] <- mean(sim.SMLE_R[,5]); result_part[1,(start+2)] <- var(sim.SMLE_R[,5]) 
        result_part[1,(start+3)] <- mean(sim.SMLE_R[,6]); result_part[1,(start+4)] <- var(sim.SMLE_R[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.SMLE_R[,1])
        result_part[1,(start+6)] <- var(sim.SMLE_R[,1])
        
        reject_Z <-  sum(sim.SMLE_R[,2] < 0.05)
        reject_BM <- sum(sim.SMLE_R[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.SMLE_R[,4])
        result_part[1,(start+12)] <- var(sim.SMLE_R[,4])
        
        result_part[1,(start+13)] <- 1- sqrt(p[1])/(sqrt(p[1]) + sqrt(p[2]))

        result_part[1,(start+14)] <- mean(sim.SMLE_R[,7]) - result_part[1,(start+13)]
        result_part
      }
      result[,1:(length(colnames))] <- result.SMLE_R 
    }
    if(method == "SMLE_Neyman"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.SMLE_Neyman <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%dorng% 
        source('SMLE.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.SMLE_Neyman <- sim_SMLE(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.SMLE_Neyman <- sim_SMLE(N=N, dist=dist, K=1, par1 = p,  measure=measure, nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- p
        }
        result_part[1,(start+1)] <- mean(sim.SMLE_Neyman[,5]); result_part[1,(start+2)] <- var(sim.SMLE_Neyman[,5]) 
        result_part[1,(start+3)] <- mean(sim.SMLE_Neyman[,6]); result_part[1,(start+4)] <- var(sim.SMLE_Neyman[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.SMLE_Neyman[,1])
        result_part[1,(start+6)] <- var(sim.SMLE_Neyman[,1])
        
        reject_Z <-  sum(sim.SMLE_Neyman[,2] < 0.05)
        reject_BM <- sum(sim.SMLE_Neyman[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.SMLE_Neyman[,4])
        result_part[1,(start+12)] <- var(sim.SMLE_Neyman[,4])
        
        result_part[1,(start+13)] <- 1 - var[1]/(var[1]+var[2])
        result_part[1,(start+14)] <- mean(sim.SMLE_Neyman[,7]) - result_part[1,(start+13)]
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.SMLE_Neyman 
    }
    if( method == "ERADE"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.ERADE <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD", "pseudorank")) %dorng% { #%dorng% 
        source('ERADE.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        mu <- as.numeric(mean[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.ERADE <- sim_ERADE(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.ERADE <- sim_ERADE(N=N, dist=dist, K=1, par1 = p,   measure= measure, nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- p
        }
        result_part[1,(start+1)] <- mean(sim.ERADE[,5]); result_part[1,(start+2)] <- var(sim.ERADE[,5]) 
        result_part[1,(start+3)] <- mean(sim.ERADE[,6]); result_part[1,(start+4)] <- var(sim.ERADE[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.ERADE[,1])
        result_part[1,(start+6)] <- var(sim.ERADE[,1])
        
        reject_Z <-  sum(sim.ERADE[,2] < 0.05)
        reject_BM <- sum(sim.ERADE[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.ERADE[,4])
        result_part[1,(start+12)] <- var(sim.ERADE[,4])
        
        if(ar=="R_minF"){
          result_part[1,(start+13)] <- 1- sqrt(p[1])/(sqrt(p[1]) + sqrt(p[2]))
        } 
        if(ar=="Neyman"){
          result_part[1,(start+13)] <- 1 - var[1]/(var[1]+var[2])
        }
        result_part[1,(start+14)] <- mean(sim.ERADE[,7]) - result_part[1,(start+13)]
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.ERADE 
    }
    if( method == "DBCD"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.DBCD <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD", "pseudorank")) %dorng% { #%dorng% 
        source('DBCD.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        mu <- as.numeric(mean[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.DBCD <- sim_DBCD(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.DBCD <- sim_DBCD(N=N, dist=dist, K=1, par1 = p,  measure= measure, nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- p
        }
        result_part[1,(start+1)] <- mean(sim.DBCD[,5]); result_part[1,(start+2)] <- var(sim.DBCD[,5]) 
        result_part[1,(start+3)] <- mean(sim.DBCD[,6]); result_part[1,(start+4)] <- var(sim.DBCD[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.DBCD[,1])
        result_part[1,(start+6)] <- var(sim.DBCD[,1])
        
        reject_Z <-  sum(sim.DBCD[,2] < 0.05)
        reject_BM <- sum(sim.DBCD[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.DBCD[,4])
        result_part[1,(start+12)] <- var(sim.DBCD[,4])
        
        if(ar=="R_minF"){
          result_part[1,(start+13)] <- 1- sqrt(p[1])/(sqrt(p[1]) + sqrt(p[2]))
        } 
        if(ar=="Neyman" ){
          result_part[1,(start+13)] <- 1 - var[1]/(var[1]+var[2])
        }
        result_part[1,(start+14)] <- mean(sim.DBCD[,7]) - result_part[1,(start+13)]
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.DBCD 
    }
    
  } else {
    print("Unfortunouatly, the multi-armed case has not been implemented yet")
  }
  
  #Print and save Results
  if(method=="ERADE" || method=="DBCD"){
    filename <-  paste(method, ar, N, dist,"All.csv",sep="_")
  } else {
    filename <-  paste(method, N, dist,"All.csv",sep="_")
  }
  
  write.csv(round(result,4), filename, row.names=TRUE)
  print(round(result,4))
  print(round(N*(1-result[7])))
  stopCluster(cl)
}

# Generate Parameter Data Frame
# Bernoulli
par1_Ber <- data.frame(0.941,0.991)


#Scenarios
N <- c(420)
burn <- round(N/12)
Distributions <- c("bern")
param_list <- expand.grid(N = N, Distributions = Distributions)


for(i in 1:dim(param_list)[1]){
  nsim = 5*10^4 # or 10^4!
  if(param_list[i,]$Distributions=="bern") {
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, burnin=N/2, nsim=nsim, method = "ER", measure = "lrr") #10^4
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "SMLE_R_minF", ar="R_minF", measure = "lrr")
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "SMLE_Neyman", ar="Neyman", measure = "lrr") #10^4
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "DBCD", ar="R_minF", measure = "lrr") #10^4
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "DBCD", ar="Neyman", measure = "lrr") #10^4
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "ERADE", ar="R_minF", measure = "lrr") #10^4
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "ERADE", ar="Neyman", measure = "lrr") #10^4
  }
}