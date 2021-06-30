MGM <- function(Temp, SR, Population, TLN, AM, SRE, MTU){

  Ne <- ncol(Temp)#number of environments
  Md <- nrow(Temp)#maximum day
  stopifnot(ncol(Temp) == ncol(SR))
  stopifnot(ncol(Temp) == length(Population))
  stopifnot(nrow(Temp) == nrow(SR))
  
  #Variables determined beforehand
  LNM <- 3.53 + 0.46 * TLN
  Tb.pheno <- 8
  Tb.grain <- 0
  SRE.reduction <- 0.75
  TU.emergence <- 87
  TU.silking <- 67
  TU.radiation <- 500
  
  #Functions
  tu <- function(Temp, d, Tb) {Temp[d] - Tb}
  ln <- function(TU) {2.5 * exp(TU * 0.00225)}
  pla <- function(LN, LNM, AM, Width) {
    Steps <- seq(2.5, LN, by=Width)
    Height <- AM * exp(-0.0344 * (Steps-LNM)^2 + 0.000731 * (Steps-LNM)^3)
    sum(Height * Width)
  }
  fas <- function(TU) {min(0.00161 * exp(0.00328 * TU), 1)}
  lai <- function(PLA, FAS, Population) {Population * (PLA - PLA * FAS) / 10000}
  ph <- function(SR, SRE, LAI) {SR * SRE * (1 - exp(-0.4 * LAI))}
  hi <- function(DFSm3) {min(DFSm3 * 0.015, 0.5)}
  
  #Objects returned
  GW.maturity <- numeric(Ne)

  #Run for each environment
  for(env in 1:Ne){
    
    #State variables
    TU <- TU.emergence * -1
    TU.forsilking <- 0
    TU.grain <- 0
    DFS <- -3
    LN <- 0
    AR <- 0
    PLA <- 0
    FAS <- 0
    LAI <- 0
    BM <- numeric(Md)
    GW <- numeric(Md)
    
    #Growth stage
    Stage <- 1
    #1: vagetative phase
    #2: reproductive phase (after silking)
    
    for(d in 1:Md){
      
      if(Stage == 1){
        
        TU.d <- tu(Temp[, env], d, Tb.pheno)
        TU <- TU + TU.d
        if(TU >= 0){
          
          if(LN < TLN){
            LN <- ln(TU)
            if(LN + 2 < TLN) 
              PLA <- pla(LN + 2, LNM, AM, 0.5)
            else 
              PLA <- pla(TLN, LNM, AM, 0.5)
          } else {
            TU.forsilking <- TU.forsilking + TU.d
          }
          
          FAS <- fas(TU)
          LAI <- lai(PLA, FAS, Population[env])
          BM[d] <- BM[d-1] + ph(SR[d, env], SRE, LAI)
          if(TU.forsilking >= TU.silking) Stage <- 2
        }
      } else {
        
        TU <- TU + tu(Temp[, env], d, Tb.pheno)
        TU.grain <- TU.grain + tu(Temp[, env], d, Tb.grain)
        FAS <- fas(TU)
        LAI <- lai(PLA, FAS, Population[env])
        
        if(TU.grain > TU.radiation)
          SRE.grain <- SRE * SRE.reduction
        else
          SRE.grain <- SRE
        
        BM[d] <- BM[d-1] + ph(SR[d, env], SRE.grain, LAI)
        DFS <- DFS + 1
        if(DFS > 0) GW[d] <- BM[d] * hi(DFS)
        if(TU.grain > MTU) break
      }
    }
    GW.maturity[env] <- GW[d]
  }
  GW.maturity
}
