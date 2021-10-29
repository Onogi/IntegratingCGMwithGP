# IntegratingCGMwithGP

These are the R and Rcpp scripts developed for illustration of crop growth models (CGMs) and genomic prediction-assisted CGMs in "Integration of crop growth models and genomic prediction" by Akio Onogi, a chapter of "Genomic prediction of complex traits" editted by Nour Ahmadi and Jerome Bartholomé.  
  
The files and Boxes in the chapter correspond with each other as follows.  
  
DVRmodel.R                 : Box 1.  
DVRmodel.cpp               : Box 2.  
RunDVRmodel.R              : Box 5.  
DVRmodel_GBM.cpp           : Box 7.  
RunDVRmodel_GBM.R          : Box 9.  
  
MaizeGrowthModel.R         : Box 3.  
MaizeGrowthModel.cpp       : Box 4.  
RunMaizeGrowthModel.R      : Box 6.  
MaizeGrowthModel_GBM.cpp   : Box 8.  
RunMaizeGrowthModel_GBM.cpp: Box 10.  

Analysis examples using these scripts are put in AnalysisExample_DVRmodel and AnalysisExample_MaizeGrowthModel.
R scritps in the folders (Example_DVRmodel.R and Example_MaizeGrowthModel.R) are illustrations how to optimize CGMs, connect with genome data, and use GenomeBasedModel. The results are stored in Example_DVRmodel.RData and Example_MaizeGrowthModel.RData, respectively.
