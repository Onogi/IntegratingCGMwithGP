#Import sample data#############################################################
#All data was used in Onogi et al. (2016) TAG.

#Emergence and heading dates
##The columns correspond with environments (9 environments)
##The rows correspond with lines (176 lines)
Emergence <- as.matrix(read.csv("EmergenceDate.csv", row.names = 1))
Heading <- as.matrix(read.csv("HeadingDate.csv", row.names = 1))

#Photoperiod and temperature at each environment
##The columns correspond with environments (9 environments)
##The rows correspond with days (197 days)
Photoperiod <- as.matrix(read.csv("Photoperiod.csv", row.names = 1))
Temperature <- as.matrix(read.csv("Temperature.csv", row.names = 1))

#SNP genotypes (-1, 0, 1)
##The columns correspond with lines (176 lines)
##The rows correspond with SNPs (162 SNPs)
Geno <- as.matrix(read.csv("Geno.csv", row.names = 1))

#Number of lines
Nl <- nrow(Emergence)

#Number of environments
Ne <- ncol(Emergence)

#Number of days
Md <- nrow(Temperature)

#Days to heading
DaysToHeading <- Heading - Emergence + 1


#DVR model######################################################################
#Use cpp version here
Rcpp::sourceCpp("DVRmodel.cpp")


#Optimize parameters in the DVR model###########################################
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

#Optimize for all lines
Estimate.All <- matrix(NA, nrow = Nl, ncol = 5)
colnames(Estimate.All) <- c("G", "Alpha", "Beta", "Po", "To")
for(line in 1:Nl){
  
  #Extract temperature for this line
  Temp <- NULL
  for(env in 1:Ne){
    if(!is.na(Emergence[line, env])){
      w <- Temperature[Emergence[line, env]:Md, env]
      if(length(w) < Md) w <- c(w, rep(0, Md - length(w)))
      Temp <- cbind(Temp, w)
    }
  }
  
  #Extract photoperiod for this line
  Photo <- NULL
  for(env in 1:Ne){
    if(!is.na(Emergence[line, env])){
      w <- Photoperiod[Emergence[line, env]:Md, env]
      if(length(w) < Md) w <- c(w, rep(0, Md - length(w)))
      Photo <- cbind(Photo, w)
    }
  }
  
  #Extract days to heading (NA removed)
  DTH <- DaysToHeading[line, ]
  DTH <- DTH[!is.na(DTH)]
  
  #Repeat 10 times and choose the best result
  Estimate <- NULL
  MSE <- 1e+6
  for(set in 1:10){
    Result <- psoptim(rep(NA, 5), FitDVR, 
                      Temp = Temp, Photo = Photo, DTH = DTH,
                      lower = Range[1, ], upper = Range[2, ])
    if(Result$value < MSE){
      Estimate <- Result$par
      MSE <- Result$value
    }
  }
  
  Estimate.All[line, ] <- Estimate
}

#The resulting estimates (Estimate.All) can be used for genetic analyses
#Estimate heritability of parameter G
library(rrBLUP)#Can be installed from CRAN
K <- A.mat(t(Geno))
rrBLUP.G <- mixed.solve(y = Estimate.All[, "G"], K = K)
rrBLUP.G$Vu/(rrBLUP.G$Vu + rrBLUP.G$Ve)

#GWAS using mixed model
pheno <- data.frame(ID=paste("V", 1:Nl, sep=""), Estimate.All)
geno <- data.frame(ID=1:nrow(Geno), Chr=0, Pos=0, Geno)
rrBLUP.GWAS <- GWAS(pheno = pheno, geno = geno, K = K, plot = FALSE)
for(para in 1:5) plot(rrBLUP.GWAS[[3 + para]])

#GWAS using Bayesian regression
library(VIGoR)#Can be installed from CRAN
ETA <- list(list(model = "BayesC", X = t(Geno)))
for(para in 1:5){
  VIGoR.GWAS <- vigor(Estimate.All[, para], ETA)
  plot(VIGoR.GWAS$ETA[[1]]$Rho)
}


#Optimize parameters using GenomeBasedModel#####################################
#create input objects of GenomeBasedModel
##observed days to heading (make all emergence days 1)
Y <- t(DaysToHeading)
Missing <- 999999
Y[is.na(Y)] <- Missing

##Extract temperature for all lines and create Input
##Each column of Input correspond with each line
##The number of rows of Input is Ne * Md
Input <- NULL
x <- NULL
for(line in 1:Nl){
  v <- NULL
  for(env in 1:Ne){
    if(is.na(Emergence[line,env])){
      w <- Temperature[1:Md, env]#assume to be sowed at the first day when missing
    }else{
      w <- Temperature[Emergence[line, env]:Md, env]
    }
    if(length(w) < Md) w <- c(w, rep(0, Md - length(w)))
    v <- c(v, w)
  }
  x <- cbind(x, v)
}
Input <- rbind(Input, x)

##Extract photoperiod for all lines and add to Input
x <- NULL
for(line in 1:Nl){
  v <- NULL
  for(env in 1:Ne){
    if(is.na(Emergence[line, env])){
      w <- Photoperiod[1:Md, env]#assume to be sowed at the first day when missing
    }else{
      w <- Photoperiod[Emergence[line, env]:Md, env]
    }
    if(length(w) < Md) w <- c(w, rep(0, Md - length(w)))
    v <- c(v, w)
  }
  x <- cbind(x, v)
}
Input <- rbind(Input, x)
rm(v, w, x)
##Now the number of rows of Input is 2 * Ne * Md
##The upper half includes temperature and the lower half includes photoperiod

##Create Freevec
Freevec <- c(Nl, Ne, Md)

##Number of parameters in the DVR model
Np <- 5

##Reference values of the parameters (taken from Onogi et al. 2016)
Referencevalues <- c(4, 1.61, 1.61, 10, 30)

##Add line IDs to Y and Geno
Y <- rbind(1:Nl, Y)
X <- rbind(1:Nl, Geno)

##Transformation
##The firt 3 parameters are log-transformed.
##The last 2 parameters are not transformed.
Transformation <- c(rep("log",3),"nt","nt")

#Run GenomeBasedModel
Rcpp::sourceCpp("DVRmodel_GBM.cpp")
library(GenomeBasedModel)
Methodcode <- c(1, 2, 4, 7, 8)
#1: Bayesian lasso
#2: Extended Bayesian lasso
#4: BayesC
#7: BayesB
#8: GBLUP
#=>Various methods are used here because this is an example.
#But variable selection methods (EBL, BayesB, or BayesC) will be suitable for this data.

start<-proc.time()[3]
Result <- GenomeBasedModel(Input, Freevec, Y, Missing, Np, X,
                           Methodcode, Referencevalues, DVRmodel_GBM,
                           Transformation = Transformation,
                           PassMatrix = TRUE)
end<-proc.time()[3]
(end-start)/3600
2.812911
##Took ~3 hours
##Now the implementation of GBLUP is slow and this would be a bottleneck.

#Compare with the results of PSO
##For parameter G
Originalscale <- exp(Result$Para[[1]])#Back to original scale
Posterior.mean <- apply(Originalscale, 2, mean)
plot(Posterior.mean, Estimate.All[, "G"]);abline(0,1)

##For parameter Alpha
Originalscale <- exp(Result$Para[[2]])
Posterior.mean <- apply(Originalscale, 2, mean)
plot(Posterior.mean, Estimate.All[, "Alpha"]);abline(0,1)

##For parameter Beta
Originalscale <- exp(Result$Para[[3]])
Posterior.mean <- apply(Originalscale, 2, mean)
plot(Posterior.mean, Estimate.All[, "Beta"]);abline(0,1)

##For parameter Po
Posterior.mean <- apply(Result$Para[[4]], 2, mean)
plot(Posterior.mean, Estimate.All[, "Po"]);abline(0,1)

##For parameter To
Posterior.mean <- apply(Result$Para[[5]], 2, mean)
plot(Posterior.mean, Estimate.All[, "To"]);abline(0,1)

#Marker effects on the 1st and 2nd parameters 
#that are estimated with Bayesian lasso and extended Bayesian lasso
#Note that these effects are for the log-transformed parameters.
plot(Result$Genome[[1]]$Beta)
plot(Result$Genome[[2]]$Beta)

#Inclusion probabilities on the 3rd and 4th parameters
#that are estimated with BayesC and BayesB
plot(Result$Genome[[3]]$Rho)
plot(Result$Genome[[4]]$Rho)

#Total additive effects on the 5th parameter
#that are estimated with GBLUP
plot(Result$Genome[[5]]$U)


#Hyperparametes of regression methods can be modified.
#For example,
Hyperpara <- list (list(Phi = 1, Omega = 0.1), #BL
                   list(Phi = 1, Omega = 0.1, Psi = 0.1, Theta = 0.1), #EBL
                   list(Nu=5, S2=1, Kappa=0.05), #BayesC
                   list(Nu=5, S2=1, Kappa=0.05), #BayesB
                   list(Nu=5, S2=1)) #BLUP
#default values are shown in the vignette
#For regression methods and their hyperparameters,
#see the pdf manual of my another package VIGoR (https://github.com/Onogi/VIGoR)

#Result <- GenomeBasedModel(Input, Freevec, Y, Missing, Np, X,
#                           Methodcode, Referencevalues, DVRmodel_GBM,
#                           Transformation = Transformation,
#                           Hyperpara = Hyperpara,
#                           PassMatrix = TRUE)

