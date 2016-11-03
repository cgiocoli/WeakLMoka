#include "power2D.h"

void powerl(std::valarray<float> a,std::valarray<float> b,
	int nx, int ny, 
	double boxlx, double boxly,
	double *ll, double *Pl,int nl){
  // go in the fourir space doing the zero padding
  int zerosize = 2;

  // size of the new map in x and y directions, factor by which each size is increased
  int Nnx=int(zerosize*nx);
  int Nny=int(zerosize*ny);
  double Nboxlx = boxlx*zerosize;
  double Nboxly = boxly*zerosize;
  
  std:: valarray<float> Na,Nb;
  Na.resize( Nnx*Nny );
  Nb.resize( Nnx*Nny );

  // assume locate in a rectangular map and build up the new one
  for( int i=0; i<Nnx; i++ ) for( int j=0; j<Nny; j++ ){ 
      Na[i+Nnx*j]=0; 
      Nb[i+Nnx*j]=0; 
      if(i>=int(Nnx/2-nx/2) && i<int(Nnx/2+nx/2) && j>=int(Nny/2-ny/2) && j<int(Nny/2+ny/2)){
	int ii = i-int(Nnx/2-nx/2); 
	int jj = j-int(Nny/2-ny/2); 

	if(ii>=nx || jj>=ny){ 
	  std::cout << " 1 error mapping " << ii << "  " << jj << std::endl; 
	  exit(1); 
	} 
	if(ii<0 || jj<0){ 
	  std::cout << " 2 error mapping " << ii << "  " << jj << std::endl; 
	  exit(1); 
	} 
	Na[i+Nnx*j]=a[ii+nx*jj];	     
	Nb[i+Nnx*j]=b[ii+nx*jj];	     
      }
    }

  // estimate the fourier transform
  // now we go in the Fourier Space
  double *dNa=new double[Nnx*Nny];
  double *dNb=new double[Nnx*Nny];
  double *input=new double[Nnx*Nny];

  fftw_complex *fNa=new fftw_complex[Nny*(Nnx/2+1)];
  fftw_complex *fNb=new fftw_complex[Nny*(Nnx/2+1)];
  fftw_complex *Nfcc=new fftw_complex[Nny*(Nnx/2+1)];
  size_t *ik=new size_t[Nny*(Nnx/2+1)];
  double *k=new double[Nny*(Nnx/2+1)];
  fftw_complex *output=new fftw_complex[Nny*(Nnx/2+1)];
  for(int i=0;i<Nnx;i++) for(int j=0;j<Nny; j++){ 
      dNa[i+Nnx*j] = double(Na[i+Nnx*j]);
      dNb[i+Nnx*j] = double(Nb[i+Nnx*j]);
    }
  fftw_plan p1;
  p1=fftw_plan_dft_r2c_2d(Nny,Nnx,dNa,fNa,FFTW_ESTIMATE);
  fftw_execute( p1 );
  fftw_destroy_plan(p1);
  fftw_plan p2;
  p2=fftw_plan_dft_r2c_2d(Nny,Nnx,dNb,fNb,FFTW_ESTIMATE);
  fftw_execute( p2 );
  fftw_destroy_plan(p2);

  // fb and fa are then used to compute the power spectrum
  // build modes
  for( int i=0; i<Nnx/2+1; i++ ){
    // kx = i if i<n/2 else i-n
    // double kx=(i<Nnx/2)?double(i):double(i-Nnx);
    double kx=double(i);
    for( int j=0; j<Nny; j++ ){
      // double ky=double(j);
      double ky=(j<Nny/2)?double(j):double(j-Nny);
      // rescale respect to the box size 
      k[i+(Nnx/2+1)*j] = sqrt(kx*kx/Nboxlx/Nboxlx + ky*ky/Nboxly/Nboxly)*2.*M_PI;      
      Nfcc[i+(Nnx/2+1)*j][0] =  (fNa[i+(Nnx/2+1)*j][0]*fNb[i+(Nnx/2+1)*j][0] + 
			   fNa[i+(Nnx/2+1)*j][1]*fNb[i+(Nnx/2+1)*j][1])*pow(1./2/M_PI/Nnx/Nny,2)*Nboxlx*Nboxly;
      Nfcc[i+(Nnx/2+1)*j][1] = -(fNa[i+(Nnx/2+1)*j][0]*fNb[i+(Nnx/2+1)*j][1] - 
			   fNa[i+(Nnx/2+1)*j][1]*fNb[i+(Nnx/2+1)*j][0])*pow(1./2/M_PI/Nnx/Nny,2)*Nboxlx*Nboxly;
    }
  }
  gsl_sort_index (ik,k,1,Nny*(Nnx/2+1));
  double *nk=new double[Nny*(Nnx/2+1)];
  double *nc=new double[Nny*(Nnx/2+1)];
  // sorted vectors
  for(int i=0; i<Nnx/2+1; i++) 
    for(int j=0; j<Nny;j++){
      nk[i+(Nnx/2+1)*j] = k[ik[i+(Nnx/2+1)*j]];
      nc[i+(Nnx/2+1)*j] = Nfcc[ik[i+(Nnx/2+1)*j]][0];
    }
  std:: vector<double> bink(nl);
  // build the binned power spectrum
  fill_linear(bink,nl,log10(nk[1]),log10(nk[Nny*(Nnx/2+1)-1]));
  double lk1,lk2;
  for(int i=0;i<nl;i++){
    if(i==0) lk1=bink[i];
    else lk1=bink[i]-0.5*(bink[i]-bink[i-1]);
    if(i==nl-1) lk2=bink[i];
    else lk2=bink[i]+0.5*(bink[i+1]-bink[i]);
    Pl[i]=0.;
    ll[i]=0.;
    int nin=0;
    // start from 1 because the first is 0
    for(int j=1;j<Nny*(Nnx/2+1)-1;j++){
      if(log10(nk[j])>=lk1 && log10(nk[j])<lk2){
	Pl[i]=Pl[i]+nc[j];
	ll[i]=ll[i]+nk[j];
	nin=nin+1;
      }
    }
    if(nin>0){
      Pl[i]=Pl[i]/double(nin)*zerosize*zerosize;
      ll[i]=ll[i]/double(nin);
    }
  }
  delete [] dNa;
  delete [] dNb;
  delete [] input;
  delete [] output;
  delete [] fNa;
  delete [] fNb;
  delete [] Nfcc;
  delete [] ik;
  delete [] k;
  delete [] nk;
  delete [] nc;
}
