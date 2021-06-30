#include <Rcpp.h>
using namespace Rcpp;
//Functions
double tu(double Temp, double Tb) {
  return Temp - Tb;
}
double ln(double TU) {
  return 2.5 * exp(TU * 0.00225);
}
double ar(double LN, double LNM, double AM) {
  return AM * exp(-0.0344 * pow(LN - LNM, 2.0) + 
                  0.000731 * pow(LN - LNM, 3.0));
}
double pla(double LN, double LNM, double AM, double Width) {
  double Sum = 0.0, l = 2.5;
  
  do {
    Sum += Width * ar(l, LNM, AM);    
    l += Width;
  } while (l <= LN);
  
  return Sum;
}
double fas(double TU) {
  double temp = 0.00161 * exp(0.00328 * TU);
  return temp < 1.0 ? temp: 1.0;
}
double lai(double PLA, double FAS, double Population) {
  return Population * (PLA - PLA * FAS) / 10000;
}
double ph(double SR, double SRE, double LAI) {
  return SR * SRE * (1 - exp(-0.4 * LAI));
}
double hi(double DFSm3) {
  double temp = DFSm3 * 0.015;
  return temp < 0.5 ? temp: 0.5;
}
//[[Rcpp::export]]
NumericMatrix MGM_GBM(NumericMatrix Input, NumericVector Freevec, 
                      NumericMatrix Parameter){
  
  int line, env, d, DFS, Stage, envMd;
  double TLN, AM, SRE, MTU;
  double LNM, Tb_pheno, Tb_grain, SRE_reduction, TU_emergence; 
  double TU_silking, TU_radiation;
  double TU, TU_d, TU_forsilking, TU_grain, SRE_grain;
  double LN, PLA, FAS, LAI;
  int Nl, Ne, Md, NeMd;

  Nl = (int)Freevec(0);//number of lines
  Ne = (int)Freevec(1);//number of environments
  Md = (int)Freevec(2);//maximum day
  NeMd = Ne * Md;//used to extract solar radiation from Input
  
  NumericVector BM(Md), GW(Md), Population(Ne);
  
  //Extract populations from Freevec
  for(env=0; env<Ne; ++env) Population(env) = Freevec(3 + env);
  
  //Object returned;
  NumericMatrix GW_maturity(Ne, Nl);
  
  //Variables determined beforehand
  Tb_pheno = 8.0;
  Tb_grain = 0.0;
  SRE_reduction = 0.75;
  TU_emergence = 87.0;
  TU_silking = 67.0;
  TU_radiation = 500.0;
  
  //Run for each line
  for(line=0; line<Nl; ++line){
    
    //Parameters of current line
    TLN = Parameter(0, line);
    AM = Parameter(1, line);
    SRE = Parameter(2, line);
    MTU = Parameter(3, line);
    LNM = 3.53 + 0.46 * TLN;
    
    //Run for each environment
    for(env=0; env<Ne; ++env){
      
      envMd = env * Md;//used to extract values from Input
      
      //State variables
      TU = TU_emergence * -1.0;
      TU_forsilking = 0.0;
      TU_grain = 0.0;
      DFS = -3;
      LN = 0.0;
      PLA = 0.0;
      FAS = 0.0;
      LAI = 0.0;
      for(d=0; d<Md; ++d){
        BM(d) = 0.0; GW(d) = 0.0;
      }
      
      //Growth stage
      Stage = 1;
      
      //Run
      for(d=0; d<Md; ++d){
        
        if(Stage == 1){
          
          TU_d = tu(Input(envMd + d, line), Tb_pheno);
          TU += TU_d;
          if(TU >= 0.0){
            
            if(LN < TLN){
              LN = ln(TU);
              if(LN + 2.0 < TLN)
                PLA = pla(LN + 2.0, LNM, AM, 0.5);
              else
                PLA = pla(TLN, LNM, AM, 0.5);
            } else {
              TU_forsilking += TU_d;
            }
            
            FAS = fas(TU);
            LAI = lai(PLA, FAS, Population(env));
            BM(d) = BM(d-1) + ph(Input(NeMd + envMd + d, line), SRE, LAI);
            if(TU_forsilking >= TU_silking) Stage = 2;
          }
        }else{
          
          TU += tu(Input(envMd + d, line), Tb_pheno);
          TU_grain += tu(Input(envMd + d, line), Tb_grain);
          FAS = fas(TU);
          LAI = lai(PLA, FAS, Population(env));
          
          if(TU_grain > TU_radiation){
            SRE_grain = SRE * SRE_reduction;
          } else {
            SRE_grain = SRE;
          }
          BM(d) = BM(d-1) + ph(Input(NeMd + envMd + d, line), SRE_grain, LAI);
          
          DFS = DFS + 1;
          if(DFS > 0) GW(d) = BM(d) * hi(DFS);
          if(TU_grain > MTU) break;
        }
      }
      if(d == Md) 
        GW_maturity(env, line) = GW(d-1); 
      else 
        GW_maturity(env, line) = GW(d);
    }
  }
  return GW_maturity;
}
