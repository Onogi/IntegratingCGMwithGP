#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
NumericMatrix DVRmodel_GBM(NumericMatrix Input, NumericVector Freevec,
                           NumericMatrix Parameter){
  
  double G, Alpha, Beta, Po, To;
  double Dvs, Dvr, Dvs1, Dvs2;
  double Pb=0.0, Pc=24.0, Popb, Pcpo, Pratio, Pd, gPd;
  double Tb=8.0, Tc=42.0, Totb, Tcto, Tratio, Td, fTd;
  int line, env, d, envMd;
  int Nl, Ne, Md, NeMd;
  int Photosensitive;

  Nl = (int)Freevec(0);//number of lines
  Ne = (int)Freevec(1);//number of environments
  Md = (int)Freevec(2);//maximum day
  NeMd = Ne * Md;//used to extract photoperiod from Input
  
  //Object returned
  NumericMatrix DTH(Ne, Nl);
  
  //Run for each line
  for(line=0; line<Nl; ++line){
    
    //Parameters of current line
    G = exp(Parameter(0, line));
    Alpha = exp(Parameter(1, line));
    Beta = exp(Parameter(2, line));
    Po = Parameter(3, line);
    To = Parameter(4, line);
    
    Popb = Po - Pb;
    Pcpo = Pc - Po;
    Pratio = Pcpo / Popb;
    Totb = To - Tb;
    Tcto = Tc - To;
    Tratio = Tcto / Totb;
    Dvs1 = 0.145*G + 0.005*G*G;
    Dvs2 = 0.345*G + 0.005*G*G;
    
    //Run for each enironment
    for(env = 0; env < Ne; ++env){
      
      envMd = env * Md;//used to extract values from Input
      
      for(d = 0, Dvs = 0.0, Photosensitive = 0; d < Md; ++d){
        
        if(Photosensitive){
          
          Pd = Input(NeMd + envMd + d, line);
          if(Pd < Po) 
            gPd = 1.0; 
          else 
            gPd = (Pd - Pb) / Popb * pow((Pc - Pd) / Pcpo, Pratio);
          
          Td = Input(envMd + d, line);
          if(Td < Tb||Td > Tc) 
            fTd = 0.0; 
          else 
            fTd = (Td - Tb) / Totb * pow((Tc - Td) / Tcto, Tratio);
          
          Dvr = pow(fTd, Alpha) * pow(gPd, Beta);
        } else {
          
          Td = Input(envMd + d, line);
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
      DTH(env, line) = d + 1;
    }
  }
  return DTH;
}
