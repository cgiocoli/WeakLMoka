#include "cMZhao.h"

const int n=128,ni=512;
const size_t limit = n;
std:: vector<double> Zdci(n),Zlzi(n),Zzi(n);
std:: vector<double> vii(ni),lmii(ni);
const double mmin=1e3,mmax=1e20;


/*
  Initialise vectors for MAH
*/
void iniTables(cosmology *cosmology,double z0){
  double lz0 = log10(1+z0);
  fill_linear(Zlzi,n,0.+lz0,1.+lz0);
  for(int i=0;i<n;i++){
    Zzi[i]=-1+pow(10,Zlzi[i]);
    Zdci[i] = cosmology->deltaC(Zzi[i])/
      (cosmology->growthFactor(1./(1+Zzi[i]))/cosmology->growthFactor(1.)/(1.+Zzi[i]));
  }
  fill_linear(lmii,ni,log10(mmin),log10(mmax));
  for(int i=0;i<ni;i++){
    // vii[i] = cdm->variance(pow(10.,lmii[i]));
  }
}


double getCZhao(cosmology *cosmology, double mass, double redshift){
  return 0;
}







