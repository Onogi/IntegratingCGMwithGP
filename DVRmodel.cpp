#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
NumericVector DVRmodel(NumericMatrix Temp, NumericMatrix Photo,
                       double G, double Alpha, double Beta, double Po, double To){
  
  double Dvs, Dvr, Dvs1, Dvs2;
  double Pb=0.0, Pc=24.0, Popb, Pcpo, Pratio, Pd, gPd;
  double Tb=8.0, Tc=42.0, Totb, Tcto, Tratio, Td, fTd;
  int env, d;
  int Ne = Temp.ncol(), Md = Temp.nrow();
  int Photosensitive;
  
  Popb = Po - Pb;
  Pcpo = Pc - Po;
  Pratio = Pcpo / Popb;
  Totb = To - Tb;
  Tcto = Tc - To;
  Tratio = Tcto / Totb;
  Dvs1 = 0.145*G + 0.005*G*G;
  Dvs2 = 0.345*G + 0.005*G*G;

  NumericVector DTH(Ne);
  
  for(env = 0; env < Ne; ++env){
    
    for(d = 0, Dvs = 0.0, Photosensitive = 0; d < Md; ++d){
      
      if(Photosensitive){
        
        Pd = Photo(d, env);
        if(Pd < Po) 
          gPd = 1.0; 
        else 
          gPd = (Pd - Pb) / Popb * pow((Pc - Pd) / Pcpo, Pratio);
        
        Td = Temp(d, env);
        if(Td < Tb||Td > Tc) 
          fTd = 0.0; 
        else 
          fTd = (Td - Tb) / Totb * pow((Tc - Td) / Tcto, Tratio);
        
        Dvr = pow(fTd, Alpha) * pow(gPd, Beta);
      } else {
        
        Td = Temp(d, env);
        if(Td < Tb||Td > Tc) 
          fTd = 0.0; 
        else 
          fTd = (Td - Tb) / Totb * pow((Tc - Td) / Tcto, Tratio);
        
        Dvr = pow(fTd, Alpha);
      }
      Dvs += Dvr;
      
      if(Dvs >= Dvs1 && Dvs <= Dvs2) 
        Photosensitive = 1; 
      else 
        Photosensitive = 0;
      
      if(Dvs > G) break;
    }
    DTH(env) = d + 1;
  }
  return DTH;
}
