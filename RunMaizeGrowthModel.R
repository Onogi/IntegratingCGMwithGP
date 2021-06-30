#Execute either row to import the function
source("MaizeGrowthModel.R")
Rcpp::sourceCpp("MaizeGrowthModel.cpp")

#Function to fit the model to observations (GW) with MSE
FitMGM <- function(Parameter, Temp, SR, Population, GW){
  
  TLN <- Parameter[1]
  AM <- Parameter[2]
  SRE <- Parameter[3]
  MTU <- Parameter[4]
  
  Result <- MGM(Temp, SR, Population, TLN, AM, SRE, MTU)
  mean((GW - Result)^2)
}

#Optimize parameters using PSO
library(pso)

#Define parameter ranges
Range <- matrix(c(6, 23, 700, 800, 1.5, 1.7, 1050, 1250), nr=2)

#Repeat 10 times and choose the best result
Estimate <- NULL
MSE <- 1e+6
for(set in 1:10){
  Result <- psoptim(rep(NA, 4), FitMGM, 
                    Temp=Temp, SR=SR, Population=Population, GW=GW,
                    lower=Range[1, ], upper=Range[2, ])
  if(Result$value < MSE){
    Estimate <- Result$par
    MSE <- Result$value
  }
}

