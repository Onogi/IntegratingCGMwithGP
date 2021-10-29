#Simulate data sets#############################################################
#Here daily mean temperature and solar radiation at Tsukuba (Japan) 
#measured for 2015-2020 are used:

##Tsukuba2015.csv
##Tsukuba2016.csv
##Tsukuba2017.csv
##Tsukuba2018.csv
##Tsukuba2019.csv
##Tsukuba2020.csv

##These files were downloaded from the web site of Japan Meteorological Agency
##https://www.data.jma.go.jp/obd/stats/etrn/index.php

##Each file contains observations from 1st May to 30th Nov.
##1st May is the sowing date

##The six years are treated as "environments" hereafter.
Ne <- 6#Number of environments

##Read files
Temperature <- NULL
Solarradiation <- NULL
for(i in 1:Ne){
  Meteoinfo <- read.csv(paste("Tsukuba",c(2015:2020)[i],".csv",sep=""), skip=5, header=F)
  Temperature <- cbind(Temperature, Meteoinfo[,2])
  Solarradiation <- cbind(Solarradiation, Meteoinfo[,4])
}
any(is.na(Temperature))
FALSE
any(is.na(Solarradiation))
FALSE

#Population (plant density, plants/m2) at each environment
Population <- rep(6, Ne)

#Maximum number of days
Md <- nrow(Temperature)

#As genome data, rice SNP data in Onogi et al. (2016, TAG) are used
##The columns correspond with lines (176 lines)
##The rows correspond with SNPs (162 SNPs)
Geno <- as.matrix(read.csv("Geno.csv", row.names = 1))
Nl <- ncol(Geno)
Nm <- nrow(Geno)

#Create parameter values for each genotype (line)
##Ranges of parameter values
Range <- matrix(c(6, 23, 700, 800, 1.5, 1.7, 1050, 1250), nr=2)
colnames(Range) <- c("TLN", "AM", "SRE", "MTU")
Np <- 4#Number of parameters

##Create parameters
QTL <- matrix(sample(1:Nm, 10*Np, replace = T), nc = Np)#10 QTLs for each parameter
Beta <- matrix(rnorm(10*Np), nc = Np)#QTL effects are drawn
Parameter <- matrix(0, nr = Nl, nc = Np)
for(para in 1:Np){
  U <- t(Geno)[, QTL[, para]] %*% Beta[, para, drop = F]
  U <- U + rnorm(Nl, 0, sqrt(var(U)*0.25))#Add noises
  v <- range(U)
  U <- U * (Range[2, para] - Range[1, para])/(v[2] - v[1])#Adjust ranges
  U <- U + (Range[1, para] - min(U))
  Parameter[, para] <- U
}
colnames(Parameter) <- colnames(QTL) <- colnames(Range)

#Create phenotypes (grain weight)
##Use cpp version here
Rcpp::sourceCpp("MaizeGrowthModel.cpp")
GW <- matrix(NA, nr = Ne, nc = Nl)
for(line in 1:Nl){
  GW[, line] <- MGM(Temperature, Solarradiation, Population, 
                    Parameter[line, 1],
                    Parameter[line, 2],
                    Parameter[line, 3],
                    Parameter[line, 4])
}
GW <- GW + rnorm(Nl * Ne, 0, sqrt(var(GW) * 0.5))#Add noises


#Optimize parameters in the maize growth model##################################
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

#Optimize for all lines
Estimate.All <- matrix(NA, nrow = Nl, ncol = Np)
colnames(Estimate.All) <- c("TLN", "AM",  "SRE", "MTU")
for(line in 1:Nl){
  
  #Repeat 10 times and choose the best result
  Estimate <- NULL
  MSE <- 1e+6
  for(set in 1:10){
    Result <- psoptim(rep(NA, 4), FitMGM, 
                      Temp = Temperature, SR = Solarradiation,
                      Population = Population,
                      GW = GW[, line],
                      lower = Range[1, ], upper = Range[2, ])
    if(Result$value < MSE){
      Estimate <- Result$par
      MSE <- Result$value
    }
  }
  
  Estimate.All[line, ] <- Estimate
}

#Compare with the true values
##For parameter TLN
plot(Parameter[, "TLN"], Estimate.All[, "TLN"]);abline(0,1)
sqrt(mean((Parameter[, "TLN"]-Estimate.All[, "TLN"])^2))

##For parameter AM
plot(Parameter[, "AM"], Estimate.All[, "AM"]);abline(0,1)
sqrt(mean((Parameter[, "AM"]-Estimate.All[, "AM"])^2))

##For parameter SRE
plot(Parameter[, "SRE"], Estimate.All[, "SRE"]);abline(0,1)
sqrt(mean((Parameter[, "SRE"]-Estimate.All[, "SRE"])^2))

##For parameter MTU
plot(Parameter[, "MTU"], Estimate.All[, "MTU"]);abline(0,1)
sqrt(mean((Parameter[, "MTU"]-Estimate.All[, "MTU"])^2))


#The resulting estimates (Estimate.All) can be used for genetic analyses
#Estimate heritability of parameter TLN
library(rrBLUP)#Can be installed from CRAN
K <- A.mat(t(Geno))
rrBLUP.TLN <- mixed.solve(y = Estimate.All[, "TLN"], K = K)
rrBLUP.TLN$Vu/(rrBLUP.TLN$Vu + rrBLUP.TLN$Ve)

#GWAS using mixed model
pheno <- data.frame(ID=paste("V", 1:Nl, sep=""), Estimate.All)
geno <- data.frame(ID=1:nrow(Geno), Chr=0, Pos=0, Geno)
rrBLUP.GWAS <- GWAS(pheno = pheno, geno = geno, K = K, plot = FALSE)
for(para in 1:4) {
  plot(rrBLUP.GWAS[[3 + para]])
  for(i in 1:10) abline(v = QTL[i, para], lty = 2, col = 2)#positions of true QTLs
}

#GWAS using Bayesian regression
library(VIGoR)#Can be installed from CRAN
ETA <- list(list(model = "BayesC", X = t(Geno)))
for(para in 1:4){
  VIGoR.GWAS <- vigor(Estimate.All[, para], ETA)
  plot(VIGoR.GWAS$ETA[[1]]$Rho)
  for(i in 1:10) abline(v = QTL[i, para], lty = 2, col = 4)
}


#Optimize parameters using GenomeBasedModel#####################################
#Create input objects of GenomeBasedModel
##Create Y
Y <- GW

##Create Input
Input <- rbind(matrix(as.numeric(Temperature), nr = Md * Ne, nc = 176),
               matrix(as.numeric(Solarradiation), nr = Md * Ne, nc = 176))

##Create Freevec
Freevec <- c(Nl, Ne, Md, Population)

##Reference values
Referencevalues <- apply(Range, 2, mean)

##Add line IDs to Y and Geno
Y <- rbind(1:Nl, Y)
X <- rbind(1:Nl, Geno)

##Missing values (not included in this simulation)
Missing <- 999999

#Run GenomeBasedModel
Rcpp::sourceCpp("MaizeGrowthModel_GBM.cpp")
library(GenomeBasedModel)
Methodcode <- c(1, 2, 7, 4)
#1: Bayesian lasso (BL)
#2: Extended Bayesian lasso (EBL)
#7: BayesB
#4: BayesC

start<-proc.time()[3]
Result <- GenomeBasedModel(Input, Freevec, Y, Missing, Np, X,
                           Methodcode, Referencevalues, MGM_GBM,
                           PassMatrix = TRUE)
end<-proc.time()[3]
(end-start)/3600
2.668261
##Took ~2.6 hours

#Compare with the true values
##For parameter TLN
Posterior.mean <- apply(Result$Para[[1]], 2, mean)
plot(Parameter[, "TLN"], Posterior.mean);abline(0,1)
sqrt(mean((Parameter[, "TLN"]-Posterior.mean)^2))

##For parameter AM
Posterior.mean <- apply(Result$Para[[2]], 2, mean)
plot(Parameter[, "AM"], Posterior.mean);abline(0,1)
sqrt(mean((Parameter[, "AM"]-Posterior.mean)^2))

##For parameter SRE
Posterior.mean <- apply(Result$Para[[3]], 2, mean)
plot(Parameter[, "SRE"], Posterior.mean);abline(0,1)
sqrt(mean((Parameter[, "SRE"]-Posterior.mean)^2))

##For parameter MTU
Posterior.mean <- apply(Result$Para[[4]], 2, mean)
plot(Parameter[, "MTU"], Posterior.mean);abline(0,1)
sqrt(mean((Parameter[, "MTU"]-Posterior.mean)^2))
#=>Accurate estimation is not easy both for PSO and GenomeBasedModel
#But GenomeBasedModel will return more accurate estimates for all parameters.


#Marker effects on the 1st and 2nd parameters
#that were estimated with Bayesian lasso and extended Bayesian lasso
plot(Result$Genome[[1]]$Beta)
for(i in 1:10) abline(v = QTL[i, 1], lty = 2, col = 2)#positions of true QTLs
plot(Result$Genome[[2]]$Beta)
for(i in 1:10) abline(v = QTL[i, 2], lty = 2, col = 2)

#Inclusion probabilities on the 3rd and 4th parameters
plot(Result$Genome[[3]]$Rho)
for(i in 1:10) abline(v = QTL[i, 3], lty = 2, col = 2)
plot(Result$Genome[[4]]$Rho)
for(i in 1:10) abline(v = QTL[i, 3], lty = 2, col = 2)

#=>Shrinkage magnitures of these regression methods can be changed by
#modifying hyperparameter values.

#For example,
Hyperpara <- list (list(Phi = 1, Omega = 0.1), #BL
                   list(Phi = 1, Omega = 0.1, Psi = 0.1, Theta = 0.1), #EBL
                   list(Nu=5, S2=1, Kappa=0.05), #BayesB
                   list(Nu=5, S2=1, Kappa=0.05)) #BayesC
#default values are shown in the vignette
#For regression methods and their hyperparameters,
#see the pdf manual of my another package VIGoR (https://github.com/Onogi/VIGoR)

#Result <- GenomeBasedModel(Input, Freevec, Y, Missing, Np, X,
#                            Methodcode, Referencevalues, MGM_GBM,
#                            Hyperpara = Hyperpara,
#                            PassMatrix = TRUE)







