#include <string>
#include <sstream>
#include "cMZhao.h"
#include "power2D.h"
#include "../Moka/utilities.h"
#include "../Moka/cosmology.h"
#include "../Moka/halo.h"
#include "../Moka/nfwHalo.h"
#include "../Moka/nfwLens.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <CCfits/CCfits>
#include <random>
#include <gsl/gsl_randist.h>

using namespace CCfits;

const gsl_rng_type * Th; 
gsl_rng * rh; // host halo concentration

// parameters in the input file
struct InputParams{
  double zs;                   // redshift of the sources
  std:: string pathcats;       // path fof and subs cat
  double omegam,omegal,h0,wq;  // cosmological parameters
  double cutR;                 // radius at which the profile is cut along the line of sight
  int nplanes;                 // number of planes i.e. files in which the catalogues are divided
  std:: string simtype;        // simulation type
  std:: string planelist;      // planelist produced by MapSim
  int nx,ny;                   // number of pixels in x and y directions
  double fx,fy;                // size of the field of view in x and y
};

// fof properties
struct FoF{
  int nfof;                    // total number of fof
  std:: vector<long> idfof;      
  std:: vector<double> mfof,m200,r200,ra,dec,z,Dc,c200,zsnap;
  std:: vector<int> plane;
};

// subs properties
struct Sub{
  int nsub;               
  std:: vector<long> idsub;
  std:: vector<long> idfof;
  std:: vector<double> ra,dec,z,Dc,msub,vdisp,vmax,rhalfmass,rmax,zsnap;
  std:: vector<int> plane;
};

// plane list
struct PlaneList{
  std:: vector<int> plane;               
  std:: vector<double> zl,zsnap;
  std:: vector<double> DlDOWN,DlUP;
  std:: vector<int> rep,snap;
};

// read fof files
void readFoF(struct FoF *f, struct InputParams p, struct PlaneList pl){
  for(int i=0;i<p.nplanes;i++){
    std:: ostringstream osid;
    std:: string sid;
    osid << i+1;
    if(i+1<10){
      sid = "0"+osid.str();      
    }else{
      sid = osid.str();
    }
    if(i+1>99){
      std:: cout << " I will stop here increase the number of digit in the file name " << std:: endl;
      exit(1);
    }	   
    std:: string filin = p.pathcats + "fofinfield_" + p.simtype + "_CoDECS.0" + sid + ".dat";
    std:: cout << " reading " << filin << std:: endl;
    std:: ifstream ifilin;
    ifilin.open(filin.c_str());
    if(ifilin.is_open()){
      // read the content
      long id;
      double mfof,ra,dec,z,Dc,m200,r200;
      while(ifilin >> id >> mfof >> ra >> dec >> z >> Dc >> m200 >> r200){
	f->idfof.push_back(id);
	f->mfof.push_back(mfof);
	f->ra.push_back(ra);
	f->dec.push_back(dec);
	f->z.push_back(z);
	f->Dc.push_back(Dc);
	f->m200.push_back(m200);
	f->r200.push_back(r200);
	f->plane.push_back(i+1);
	f->zsnap.push_back(pl.zsnap[i]);
      }
      ifilin.close();
    }else{
      std:: cout << " problem with FOF files " << std:: endl;
      std:: cout << " file = " << filin << std:: endl;
      std:: cout << " does not exist ... I will STOP here!!! " << std:: endl;
      exit(1);
    }
  }
  f->nfof = f->idfof.size();
  std:: cout << " total number of fof " << f->nfof << std:: endl;
}

// read subfind files
void readSubs(struct Sub *s,struct InputParams p, struct PlaneList pl){
  for(int i=0;i<p.nplanes;i++){
    std:: ostringstream osid;
    std:: string sid;
    osid << i+1;
    if(i+1<10){
      sid = "0"+osid.str();      
    }else{
      sid = osid.str();
    }
    if(i+1>99){
      std:: cout << " I will stop here increase the number of digit in the file name " << std:: endl;
      exit(1);
    }	   
    std:: string filin = p.pathcats + "subinfield_" + p.simtype + "_CoDECS.0" + sid + ".dat";
    std:: cout << " reading " << filin << std:: endl;
    std:: ifstream ifilin;
    ifilin.open(filin.c_str());
    if(ifilin.is_open()){
      // read the content
      long idsub,idfof;
      double ra,dec,z,Dc,msub,vdisp,vmax,rhalfmass,rmax;
      while(ifilin >> idsub >> idfof >> ra >> dec >> Dc >> z >> msub >> vdisp >> vmax >> rhalfmass >> rmax){
	s->idsub.push_back(idsub);
	s->idfof.push_back(idfof);
	s->ra.push_back(ra);
	s->dec.push_back(dec);
	s->Dc.push_back(Dc);
	s->z.push_back(z);
	s->msub.push_back(msub);
	s->vdisp.push_back(vdisp);
	s->vmax.push_back(vmax);
	s->rhalfmass.push_back(rhalfmass);
	s->rmax.push_back(rmax);
	s->plane.push_back(i+1);
	s->zsnap.push_back(pl.zsnap[i]);
      }
      ifilin.close();
    }else{
      std:: cout << " problem with Subfind files " << std:: endl;
      std:: cout << " file = " << filin << std:: endl;
      std:: cout << " does not exist ... I will STOP here!!! " << std:: endl;
      exit(1);
    }
  }
  s->nsub = s->idsub.size();
  std:: cout << " total number of subs " << s->nsub << std:: endl;
}

// read the input file "weaklensingMOKA.ini"
void readInput(struct InputParams *p){
  // read fileinput
  std:: ifstream fin ("weaklensingMOKA.ini");
  if(fin.is_open());
  else{
    std:: cout << " weaklensingMOKA.ini file does not exist where you are running the code " << std:: endl;
    std:: cout << " I will STOP here!!! " << std:: endl;
    exit(1);
  }
  std:: string str;
  fin >> str;
  fin >> p->omegam;          //  1. omega matter
  fin >> str; 
  fin >> p->omegal;          //  2. omega lambda
  fin >> str;
  fin >> p->h0;              //  3. hubble constant in unit of 100 
  fin >> str;
  fin >> p->wq;              //  4. dark energy equation of state parameter
  fin >> str;
  fin >> p->zs;              //  5. source redshift
  fin >> str;
  fin >> p->cutR;            //  6. radius at which cut the density profile in 1/Rvir
  fin >> str;
  fin >> p->pathcats;        //  7. path of files
  fin >> str;
  fin >> p->nplanes;         //  8. number of planes in which the catalogues are divided 
  fin >> str;
  fin >> p->simtype;         //  9. simulation type name
  fin >> str;
  fin >> p->planelist;       // 10. planelist file produced by MapSim
  fin >> str;
  fin >> p->nx;              // 11. number of pixels in x
  fin >> str;
  fin >> p->ny;              // 12. number of pixels in y
  fin >> str;
  fin >> p->fx;              // 13. size of the field of view in x
  fin >> str;
  fin >> p->fy;              // 14. size of the field of view in y

  /**
   *  check on cutR:  
   *  this value could be positive or negative ... not null
   *  if positive the integral in z of kappa is extended as the profile
   *  on the plane of the sky
   *  if negative on the plane of the sky is threated as it was positive
   *  while alone of the line-of-sight the integral of the density is 
   *  extended up to infinity
   **/
  if(fabs(p->cutR)<1e-4){
    std:: cout << " cutR parameter too small = " << p->cutR << std:: endl;
    std:: cout << " check this out ... for now I will STOP here!!! " << std:: endl;
    exit(1);
  }

  std:: cout << "  " << std:: endl;
  std:: cout << "  ___________________________________________"  << std:: endl;
  std:: cout << "  running with the following paramters:  " << std:: endl;
  std:: cout << "  Omegam           = " << p->omegam << std:: endl;
  std:: cout << "  Omegal           = " << p->omegal << std:: endl;
  std:: cout << "  Hubble parameter = " << p->h0 << std:: endl;
  std:: cout << "  w                = " << p->wq << std:: endl;
  std:: cout << "  source redshift  = " << p->zs << std:: endl;
  std:: cout << "  Rcut             = " << p->cutR << std:: endl;
  std:: cout << "  path catalogues  = " << p->pathcats << std:: endl;
  std:: cout << "  n planes of cats = " << p->nplanes << std:: endl;
  std:: cout << "  sim name         = " << p->simtype << std:: endl;
  std:: cout << "  plane list file  = " << p->planelist << std:: endl;
  std:: cout << "  nx = " << p->nx << "  ny = " << p->ny << std:: endl;
  std:: cout << "  field of view: x = " << p->fx << "  y = " << p->fy << std:: endl;
  std:: cout << "  ___________________________________________"  << std:: endl;  
  std:: cout << "  " << std:: endl;    
}

struct gslfparams_fc{
  double vmax,rmax,h;
};

double gslf_fc(double lc, void *ip){
  struct gslfparams_fc *fp  = (struct gslfparams_fc *) ip;
  double c = pow(10.,lc);
  return log(200./3.*gsl_pow_3(c)/(log(1.+c)-c/(1.+c))) - log(14.426*gsl_pow_2(fp->vmax/(fp->h*100.)/(fp->rmax/1000.)));
}

double getC(double vmax,double rmax,double h){
  int status;
  int iter = 0, max_iter = 512;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  struct  gslfparams_fc params = {vmax,rmax,h};
  F.function = &gslf_fc;
  F.params = &params;
  double lc;
  double esp = 1e-4;
  double x_lo = -3, x_hi = 4.; 
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      lc = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
				       0, esp);
      //if (status == GSL_SUCCESS)
      //	printf ("Converged:\n");	  
      //std:: cout << lc << std:: endl;
    }
  while (status == GSL_CONTINUE && iter < max_iter);      
  gsl_root_fsolver_free (s);
  return pow(10.,lc);
}

double getCNeto(double m, double z){
  double cm = 4.67*pow(m/1.e+14,-0.11)/(1+z); 
  return gsl_ran_lognormal (rh,log(cm),0.25);
}

// starting will show when the code has been run
void showStart(){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  std:: cout << "  " << std:: endl;
  std:: cout << " ********************************************************** " << std:: endl;
  std:: cout << " *                                                        * " << std:: endl;
  std:: cout << " *       Weak lensing Light-Cones from MOKA-Haloes        * " << std:: endl;
  std:: cout << " *                                                        * " << std:: endl;
  std:: cout << " *                                                        * " << std:: endl;
  std:: cout << "  " << std:: endl;
  std:: cout << "  running on: " << asctime(timeinfo) << std:: endl;
  std:: cout << " ********************************************************** " << std:: endl;
}

void readPlaneList(struct PlaneList *pl, struct InputParams p){
  std:: string filin = p.planelist;
  std:: cout << " reading planes_list file " << filin << std:: endl;
  std:: ifstream ifilin;
  ifilin.open(filin.c_str());
  if(ifilin.is_open()){
    int plane,rep,snap;
    double zl,dlDOWN,dlUP,zsnap;
    while(ifilin >> plane >> zl >> dlDOWN >> dlUP >> rep >> snap >> zsnap){
      pl->plane.push_back(plane);
      pl->zl.push_back(zl);
      pl->DlDOWN.push_back(dlDOWN);
      pl->DlUP.push_back(dlUP);
      pl->rep.push_back(rep);
      pl->snap.push_back(snap);
      pl->zsnap.push_back(zsnap);
    }
  }else{
    std:: cout << " planes_list file does not exist.. I will STOP here!!! " << std:: endl;
    exit(1);
  }
}

void writeFits(std:: string filename,std:: valarray<float> f, int npix, int npixy){
  long naxis=2;
  long naxes[2]={npix,npixy};
  std::auto_ptr<FITS> fout(new FITS(filename,FLOAT_IMG,naxis,naxes));
  std::vector<long> naxex(2);
  naxex[0]=npix;
  naxex[1]=npixy;
  PHDU *phout=&fout->pHDU();
  phout->write( 1, npix*npixy, f );
}

int main(){
  Th = gsl_rng_default;
  rh = gsl_rng_alloc (Th);
  gsl_rng_set(rh,1234); 
  showStart();
  struct InputParams p;
  readInput(&p);
  //************** C O S M O L O G Y **********************//
  //**************   using MOKA lib  **********************//
  cosmology co(p.omegam,p.omegal,p.h0,p.wq);
  //*******************************************************//  
  struct FoF fof;
  struct Sub sub;
  struct PlaneList pl;
  readPlaneList(&pl,p);
  readFoF(&fof,p,pl);
  readSubs(&sub,p,pl);
  // loop on haloes
  // allocate c200
  fof.c200.resize(fof.nfof);
  for(int i=0;i<fof.nfof;i++){
    if(fof.zsnap[i]<p.zs){
      // one plane per time
      int planefof = fof.plane[i];
      long idfof = fof.idfof[i];
      std:: vector<double> vmax,rmax,msub;
      // std:: cout << i << "  " << fof.mfof[i] << std:: endl;
      for(int j=0;j<sub.nsub;j++){
	if(sub.plane[j] == planefof && sub.idfof[j] == idfof){
	  msub.push_back(sub.msub[j]);
	  vmax.push_back(sub.vmax[j]);
	  rmax.push_back(sub.rmax[j]);
	  // std:: cout << sub.msub[j] << "   " << sub.vmax[j] << "  " << sub.rmax[j] << std:: endl;
	}
      }
      if(msub.size()>0){
	std:: vector<double> :: iterator itmostmassive = max_element(msub.begin(), msub.end());      
	double vmaxh = vmax[distance(msub.begin(),itmostmassive)];
	double rmaxh = rmax[distance(msub.begin(),itmostmassive)];
	// use the simulation snapshot 
	double hz = co.hubble(fof.zsnap[i]);
	// Net et al. 2008 cM relation evolving with z as in Bullock et al. 2010
	double c = getCNeto(fof.m200[i],fof.zsnap[i]);
	// use Hz
	// double c = getC(vmaxh,rmaxh,hz);
	// use H0 as in Springel et al. 2008 unrealistic redshift evolution!
	// double c = getC(vmaxh,rmaxh,p.h0);
	fof.c200[i] = c;
	if(fof.m200[i]>0){
	  std:: cout << i << "  " << fof.mfof[i] << "  " << fof.m200[i] << "  " << fof.zsnap[i] << "  " 
		     << fof.plane[i] << "  "  << c << std:: endl;      
	}
      }else{
	// not a sub associated!
	fof.c200[i] = 0;
      }
    }
  }
  std:: cout << " done with the concentrations " << std:: endl;
  double fxrad = p.fx*M_PI/180.;
  double fyrad = p.fy*M_PI/180.;
  double fxmin = -fxrad/2.;
  double fxmax =  fxrad/2.;
  double fymin = -fyrad/2.;
  double fymax =  fyrad/2.;
  // start to build the map
  std:: valarray<float> kappa(p.nx*p.ny);
  std:: vector<double> x,y;
  fill_linear(x,p.nx,-fxrad/2.,fxrad/2.);
  fill_linear(y,p.ny,-fyrad/2.,fyrad/2.);

  //for each lens plane
  //for(int j=0;j<p.nplanes;j++){
  //std:: valarray<float> kappai(p.nx*p.ny);
  for(int i=0;i<fof.nfof;i++){
    //if(fof.plane[i]>(j+1)) break;
    
    bool twohundredc=true;
    halo *ha;
    nfwHalo *nha;
    nfwLens *nLha;
    
    // if there is info on c200 (means a sub associated) and m200 (means a m200 halo associated)
    // we want also that the haloes has a lens redshift smaller than the source redshift
    // if(fof.m200[i]>0 && fof.c200[i]>0 && fof.z[i]<p.zs){
    // simulation redshift ... better comparison with the particle planes
    
    // if(fof.m200[i]>0 && fof.c200[i]>0 && fof.zsnap[i]<p.zs && fof.plane[i]==(j+1)){
    if(fof.m200[i]>0 && fof.c200[i]>0 && fof.zsnap[i]<p.zs){
      ha = new halo(&co,fof.m200[i],fof.z[i],twohundredc);
      std:: cout << fof.m200[i] << "  " << fof.z[i] << "  " << fof.c200[i] << std:: endl;
      nha = new nfwHalo(*ha,fof.c200[i]);
      nLha = new nfwLens(*nha,p.zs);
      double rs = nLha->getScaleRadius();
      double R200 = nLha->getRvir();
      // is in unit of the scale radius
      double Rz = R200*fabs(p.cutR)/rs;
      int xi0=locate(x,fof.ra[i]);
      xi0=GSL_MIN( GSL_MAX( xi0, 0 ), p.nx-1 );
      int yi0=locate(y,fof.dec[i]);
      yi0=GSL_MIN( GSL_MAX( yi0, 0 ), p.ny-1 );
      double Dl = (co.angularDist(0.,fof.z[i])*co.cn.lightspeed*1.e-7);
      double dr = R200*fabs(p.cutR)/Dl;
      double xmin = fof.ra[i] - dr;
      double xmax = fof.ra[i] + dr;
      double ymin = fof.dec[i] - dr;
      double ymax = fof.dec[i] + dr;
      int imin=0, imax=p.nx-1, jmin=0, jmax=p.ny-1;	  
      // check the buffer !!!
      // ----> ( CASE 1 )
      if(fof.ra[i]>=fxmin && fof.ra[i]<=fxmax && fof.dec[i]>=fymin && fof.dec[i]<=fymax){
	xmin = xmin>=fxmin ? xmin : fxmin;
	xmax = xmax<=fxmax ? xmax : fxmax;
	ymin = ymin>=fymin ? ymin : fymin;
	ymax = ymax<=fymax ? ymax : fymax;
	if(xmin == fxmin) imin = 0;
	else{
	  imin = locate(x,xmin);
	  imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
	}
	if(ymin == fymin) jmin = 0;
	else{
	  jmin = locate(y,ymin);
	  jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
	}
	if(xmax == fxmax) imax = p.nx-1; 
	else{
	  imax = locate(x,xmax);
	  imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
	}
	if(ymax == fymax) jmax = p.ny-1; 
	else{
	  jmax = locate(y,ymax);
	  jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
	}
      }
      // ----> ( CASE 2 ) 
      if(fof.ra[i] < fxmin && fof.dec[i] >= fymin && fof.dec[i] <= fymax){
	xmin = fxmin;
	imin = 0;
	xmax = fof.ra[i] + dr;
	double dx = fxmin - fof.ra[i];
	double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dx));
	ymin = (fof.dec[i] - dy)>=fymin ? (fof.dec[i] - dy) : fymin;
	ymax = (fof.dec[i] + dr)<=fymax ? (fof.dec[i] + dy) : fymax;
	if(ymin == fymin) jmin = 0;
	else{
	  jmin = locate(y,ymin);
	  jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
	}
	if(ymax == fymax) jmax = p.ny-1;
	else{
	  jmax = locate(y,ymax);
	  jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
	}
	if(xmax>=fxmin){
	  imax = locate(x,xmax);
	  imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
	}else{
	  imax=0;
	}
      }
      // ----> ( CASE 3 ) 	  
      if(fof.ra[i] > fxmax && fof.dec[i] >= fymin && fof.dec[i] <= fymax){
	xmin = fof.ra[i] - dr;
	xmax = fxmax;
	imax = p.nx - 1;
	double dx = fof.ra[i] - fxmax;
	double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dx));
	ymin = (fof.dec[i] - dy)>=fymin ? (fof.dec[i] - dy) : fymin ;
	ymax = (fof.dec[i] + dy)<=fymax ? (fof.dec[i] + dy) : fymax ;
	if(ymin == fymin) jmin = 0;
	else{
	  jmin = locate(y,ymin);
	  jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
	}
	if(ymax == fymax) jmax = p.ny-1;
	else{
	  jmax = locate(y,ymax);
	  jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
	}
	if(xmin<=fxmax){
	  imin = locate(x,xmin);
	  imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
	}else{
	  imin=p.nx-1;
	}
      } 
      // ----> ( CASE 4 ) 	        
      if(fof.dec[i] < fymin && fof.ra[i] >= fxmin && fof.ra[i] <=fxmax){
	ymin = fymin;
	jmin = 0;
	ymax = fof.dec[i] + dr;
	double dy = ymin - fof.dec[i];
	double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dy));
	xmin = (fof.ra[i] - dx)>=fxmin ? (fof.ra[i] - dx) : fxmin ;
	xmax = (fof.ra[i] + dx)<=fxmax ? (fof.ra[i] + dx) : fxmax ;
	if(xmin == fxmin) imin = 0;
	else{
	  imin = locate(x,xmin);
	  imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
	}
	if(xmax == fxmax) imax = p.nx-1;
	else{
	  imax = locate(x,xmax);
	  imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
	}
	if(ymax>=fymin){
	  jmax = locate(y,ymax);
	  jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
	}else{
	  jmax = 0;
	}
      }
      // ----> ( CASE 5 ) 	        
      if(fof.dec[i] > fymax && fof.ra[i] >= fxmin && fof.ra[i] <=fxmax){
	ymin = fof.dec[i] - dr;
	ymax = fymax;
	jmax = p.ny-1;
	double dy = fof.dec[i]-fymax;
	double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dy));
	xmin = (fof.ra[i] - dx)>=fxmin ? (fof.ra[i] - dx) : fxmin ;
	xmax = (fof.ra[i] + dx)<=fxmax ? (fof.ra[i] + dx) : fxmax ;
	if(xmin == fxmin) imin = 0;
	else{
	  imin = locate(x,xmin);
	  imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
	}
	if(xmax == fxmax) imax = p.nx-1;
	else{
	  imax = locate(x,xmax);
	  imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
	}
	if(ymin<=fymax){
	  jmin = locate(y,ymin);
	  jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
	}else{
	  jmin = p.ny-1;
	}
      }
      // ---> ( CASE 6 )
      if(fof.ra[i] > fxmax && fof.dec[i] > fymax){
	xmax = fxmax;
	ymax = fymax;
	imax = p.nx-1;
	jmax = p.ny-1;
	double dxcm = fof.ra[i] -fxmax;
	double dycm = fof.dec[i]-fymax;
	double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
	double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
	xmin = (fof.ra[i] - dx)>=fxmin ? (fof.ra[i] - dx) : fxmin ;
	ymin = (fof.dec[i] - dy)>=fymin ? (fof.dec[i] - dy) : fymin ;
	if(xmin == fxmin) imin = 0;
	else{
	  imin = locate(x,xmin);
	  imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
	}
	if(ymin == fymin) jmin = 0;
	else{
	  jmin = locate(y,ymin);
	  jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
	}
      }
      // --> ( CASE 7 )
      if(fof.ra[i] > fxmax && fof.dec[i] < fymin){
	xmax = fxmax;
	imax = p.nx-1;
	ymin = fymin;
	jmin = 0;
	double dycm = fymin-fof.dec[i];
	double dxcm = fof.ra[i]-fxmax;
	double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
	double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
	xmin = (fof.ra[i] - dx)>=fxmin ? (fof.ra[i] - dx) : fxmin ;
	ymax = (fof.dec[i] + dy)<=fymax ? (fof.dec[i] + dy) : fymax ;
	if(xmin == fxmin) imin = 0;
	else{
	  imin = locate(x,xmin);
	  imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
	}
	if(ymax == fxmax) jmax = p.ny-1;
	else{
	  jmax = locate(y,ymax);
	  jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
	}
      }
      // ---> ( CASE 8 )
      if(fof.ra[i] < fxmin && fof.dec[i] < fymin){
	xmin = fxmin;
	ymin = fymin;
	imin = 0;
	jmin = 0;
	double dxcm = fxmin-fof.ra[i];
	double dycm = fymin-fof.dec[i];
	double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
	double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
	xmax = (fof.ra[i] + dx)<=fxmax ? (fof.ra[i] + dx) : fxmax ;
	ymax = (fof.dec[i] + dy)<=fxmax ? (fof.dec[i] + dy) : fymax ;
	if(xmax == fxmax) imax = p.nx-1;
	else{
	  imax = locate(x,xmax);
	  imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
	}
	if(ymax == fymax) jmax = p.ny-1;
	else{
	  jmax = locate(y,ymax);
	  jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
	}
      }
      // ---> ( CASE 9 )
      if(fof.ra[i] < fxmin && fof.dec[i] > fymax){
	xmin = fxmin;
	ymax = fymax;
	imin = 0;
	jmax = p.ny-1;
	double dxcm = fxmin-fof.ra[i];
	double dycm = fof.dec[i]-fymax;
	double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
	double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
	xmax = (fof.ra[i] + dx)<=fxmax ? (fof.ra[i] + dx) : fxmax ;
	ymin = (fof.dec[i] - dy)>=fymin ? (fof.dec[i] - dy) : fymin ;
	if(xmax == fxmax) imax = p.nx-1;
	else{
	  imax = locate(x,xmax);
	  imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
	}
	if(ymax == fymin) jmin = 0;
	else{
	  jmin = locate(y,ymax);
	  jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
	}
      }
      if(imin<0) imin =0;
      if(jmin<0) jmin =0;
      for(int jj=jmin; jj<=jmax; jj++ ) for(int ii=imin; ii<=imax; ii++ ){
	  // this need to be in Mpc/h physical!!!
	  double dx=(x[ii]-fof.ra[i])*Dl;
	  double dy=(y[jj]-fof.dec[i])*Dl;
	  double r=sqrt( dx*dx+dy*dy ); 
	  if(r<=dr){
	    if(p.cutR>0){
	      kappa[ii+p.nx*jj]+=nLha->kappaz(r/rs,Rz); 
	    }else{
	      // if negative consider the integral up to infinity
	      kappa[ii+p.nx*jj]+=nLha->kappa(r/rs); 
	    }
	  }
	}       
      delete ha;
      delete nha;
      delete nLha;
    }
  }
  // km
  float km = kappa.sum()/(p.nx*p.ny);
  for(int jj=0; jj<p.ny; jj++ ) for(int ii=0; ii<p.nx; ii++ ){    
      kappa[ii+p.nx*jj] = (kappa[ii+p.nx*jj]-km);
    }
  std:: string filout = "kappa_test.fits";
  writeFits(filout,kappa,p.nx,p.ny);
  
  double *ll;
  double *Pl;
  int nb = 128;
  ll=new double[nb];
  Pl=new double[nb];
  powerl(kappa,kappa,p.nx,p.ny,fxrad,fyrad,ll,Pl,nb);  
  std:: ofstream outfile;
  outfile.open("mapPowerSpectrum.dat");
  for(int i=0;i<nb;i++){
    outfile << ll[i] << "  " << Pl[i] << std::endl;
  } 
  outfile.close();
  return 0;
}
