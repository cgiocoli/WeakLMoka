#include "cMZhao.h"

const int n=512;
std:: vector<double> Zdci(n),Zlzi(n),Zzi(n);
  std:: vector<double> Zlm,Zls;
const double mmin=1e3,mmax=1e20;

double getS(double mi){
  double lmi = log10(mi);
  double lsi = getY(Zlm,Zls,lmi);
  return pow(10,lsi);
}

/*
  Initialise vectors for MAH
*/
void iniTables(cosmology *cosmology, std:: vector<double> lm, std:: vector<double> ls){
  fill_linear(Zlzi,n,0.,2.);
  for(int i=0;i<n;i++){
    Zzi[i]=-1+pow(10,Zlzi[i]);
    Zdci[i] = cosmology->deltaC(Zzi[i])/
      (cosmology->growthFactor(1./(1+Zzi[i]))/cosmology->growthFactor(1.)/(1.+Zzi[i]));
  }
  for(int i=0;i<lm.size();i++){
    Zlm.push_back(lm[i]);
    Zls.push_back(ls[i]);
  }
}

double getCZhao(cosmology *cosmology, double mass, double redshift,std:: string vir){
  if(Zlm.size()<1){
    std:: cout << " masses and variances not properly initialized! " << std:: endl;
    std:: cout << " check this out! ... or lead the correct file for this " << std:: endl;
    std:: cout << " .... I will STOP here!!! " << std:: endl;
    exit(1);
  }
  double f = 0.04;
  double alphaf;
  if(vir=="200"){
    alphaf = 1.356/pow(f,0.65)*exp(-2.0*f*f*f);
  }else{
    alphaf = 0.815/pow(f,0.707)*exp(-2.0*f*f*f);
  }
  double wf = sqrt(2*log(alphaf+1));
  double dcf = cosmology->deltaC(redshift)/
    (cosmology->growthFactor(1./(1+redshift))/cosmology->growthFactor(1.)/(1.+redshift)) + 
    wf*sqrt(getS(f*mass) - getS(mass));
  double zf = -1 + pow(10.,getY(Zdci,Zlzi,dcf));
  double t0 = cosmology->time(redshift);
  double tf = cosmology->time(zf);
  if(vir=="200"){
    // equation (23) Giocoli+13
    return 4.*pow(  1 + pow(t0/3.2/tf,8.),1./13.);
  }else{
    // standard
    return 4.*pow(  1 + pow(t0/3.75/tf,8.4),1./8.);
    // increased by \sigma_8
    // return 4.78*pow(  1 + pow(t0/3.75/tf,10.41),1./8.);
  }
}







