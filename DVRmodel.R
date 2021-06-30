DVRmodel <- function(Temp, Photo, G, Alpha, Beta, Po, To){
  
  Ne <- ncol(Temp)#number of environments
  Md <- nrow(Temp)#maximum day
  stopifnot(ncol(Temp)==ncol(Photo))
  stopifnot(nrow(Temp)==nrow(Photo))
  
  #Variables determined beforehand
  Pb <- 0.0
  Pc <- 24.0
  Tb <- 8.0
  Tc <- 42.0
  Popb <- Po - Pb
  Pcpo <- Pc - Po
  Pratio <- Pcpo / Popb
  Totb <- To - Tb
  Tcto <- Tc - To
  Tratio <- Tcto / Totb
  Dvs1 <- 0.145*G + 0.005*G*G
  Dvs2 <- 0.345*G + 0.005*G*G
  
  #Objects returned
  DTH <- numeric(Ne)
  
  #Run for each environment
  for(env in 1:Ne){
    
    #Growth stage
    Photosensitive <- 0
    Dvs <- 0
    
    for(d in 1:Md){
      
      if(Photosensitive){
        
        Pd <- Photo[d, env]
        if(Pd<Po) 
          gPd <- 1.0 
        else 
          gPd <- (Pd - Pb)/Popb * ((Pc - Pd)/Pcpo)^Pratio
        
        Td <- Temp[d, env]
        if(Td<Tb|Td>Tc) 
          fTd <- 0.0 
        else 
          fTd <- (Td - Tb)/Totb * ((Tc - Td)/Tcto)^Tratio
        
        Dvr <- fTd^Alpha * gPd^Beta
      } else {
        
        Td <- Temp[d, env]
        if(Td<Tb|Td>Tc) 
          fTd <- 0.0 
        else 
          fTd <- (Td - Tb)/Totb * ((Tc - Td)/Tcto)^Tratio
        
        Dvr <- fTd^Alpha
      }
      Dvs <- Dvs + Dvr
      
      if(Dvs >= Dvs1 & Dvs <= Dvs2) 
        Photosensitive <- 1 
      else 
        Photosensitive <- 0
      
      if(Dvs > G) {Nothead <- 0; break} else {Nothead <- 1}
    }
    if(Nothead) DTH[env] <- Md + 1 else DTH[env] <- d
  }
  DTH
}
