#Execute either row to import the function
source("DVRmodel.R")
Rcpp::sourceCpp("DVRmodel.cpp")

#Function to fit the model to observations (DTH) with MSE
FitDVR <- function(Parameter, Temp, Photo, DTH){
  
  G <- Parameter[1]
  Alpha <- Parameter[2]
  Beta <- Parameter[3]
  Po <- Parameter[4]
  To <- Parameter[5]
  
  Result <- DVRmodel(Temp, Photo, G, Alpha, Beta, Po, To)
  mean((DTH - Result)^2)
}

#Optimize parameters using PSO
library(pso)

#Define parameter ranges
Range <- matrix(c(20, 120, 1e-5, 10, 1e-5, 10, 0, 24, 8, 42), nr=2)

#Repeat 100 times and choose the best result
Estimate <- NULL
MSE <- 1e+6
for(set in 1:10){
  Result <- psoptim(rep(NA, 5), FitDVR, 
                    Temp=Temp, Photo=Photo, DTH=DTH,
                    lower=Range[1, ], upper=Range[2, ])
  if(Result$value < MSE){
    Estimate <- Result$par
    MSE <- Result$value
  }
}