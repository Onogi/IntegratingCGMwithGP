Rcpp::sourceCpp("DVRmodel_GBM.cpp")
library(GenomeBasedModel)
Methodcode <- c(1, 2, 4, 7, 8)
#1: Bayesian lasso
#2: Extended Bayesian lasso
#4: BayesC
#7: BayesB
#8: GBLUP

Result <- GenomeBasedModel(Input, Freevec, Y, Missing, Np, Geno,
                           Methodcode, Referencevalues, DVRmodel_GBM,
                           Transformation = c(rep("log",3),"nt","nt"),
                           PassMatrix = TRUE)

#Posterior means of the 1st parameter
apply(Result$Para[[1]], 2, mean)

#Marker effects on the 1st and 2nd parameters 
#that are estimated with Bayesian lasso and extended Bayesian lasso
Result$Genome[[1]]$Beta
Result$Genome[[2]]$Beta

#Marker inclusion probabilities for the 3rd and 4th parameters
#that are estimated with BayesC and BayesB
Result$Genome[[3]]$Rho
Result$Genome[[4]]$Rho

#Total additive effects on the 5th parameter
#that are estimated with GBLUP
Result$Genome[[5]]$U
