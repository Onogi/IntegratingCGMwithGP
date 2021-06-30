Rcpp::sourceCpp("MaizeGrowthModel_GBM.cpp")
library(GenomeBasedModel)
Methodcode <- c(1, 2, 4, 8)
#1: Bayesian lasso
#2: Extended Bayesian lasso
#4: BayesC
#8: GBLUP

Result <- GenomeBasedModel(Input, Freevec, Y, Missing, Np, Geno,
                           Methodcode, Referencevalues, MGM_GBM,
                           PassMatrix = TRUE)

#Posterior means of the 2nd parameter
apply(Result$Para[[2]], 2, mean)

#Marker effects on the 1st and 2nd parameters 
Result$Genome[[1]]$Beta
Result$Genome[[2]]$Beta

#Marker inclusion probabilities for the 3rd parameter
Result$Genome[[3]]$Rho

#Total additive effects on the 4th parameter
Result$Genome[[4]]$U


