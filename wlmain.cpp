#include <string>
#include <sstream>
#include "cMZhao.h"
#include "power2D.h"
#include "../Moka/utilities.h"
#include "../Moka/cosmology.h"
#include "../Moka/halo.h"
#include "../Moka/nfwHalo.h"
#include "../Moka/nfwLens.h"
#include "../Moka/hernq_Halo.h"
#include "../Moka/hernq_Lens.h"
#include "../Moka/jaffe_Halo.h"
#include "../Moka/jaffe_Lens.h"
#include "../Moka/sisHalo.h"
#include "../Moka/sisLens.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <CCfits/CCfits>
#include <random>
#include "time.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>   

using namespace CCfits;

const gsl_rng_type * Th; 
gsl_rng * rh; // host halo concentration
const gsl_rng_type * Thfof; 
gsl_rng * rhfof; // host halo concentration fof

inline bool exists_file (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

// parameters in the input file
struct InputParams{
  int nzs,ncutR;               // number of source redshifts and number of cutRadius
  std:: vector<double> zs;     // redshifts of the sources
  std:: string pathcats;       // path fof and subs cat
  double omegam,omegal,h0,wq;  // cosmological parameters
  double cutR;                 // radius at which the profile is cut along the line of sight
  int nplanes;                 // number of planes i.e. files in which the catalogues are divided
  std:: string simtype;        // simulation type
  std:: string planelist;      // planelist produced by MapSim
  int nx,ny;                   // number of pixels in x and y directions
  double fx,fy;                // size of the field of view in x and y
  std:: string filsigma;       // file contaning lm s relation used for c-M Zhao
  double sigmalnC;             // scatter in the c-M relation if a c-M model is used
  std:: string vir;            // vir use FOF haloes
  std:: string simcase;        // MapSim files, Pinocchio PLC or Begogna CAT
  std:: string PinocchioFile;  // PLC Pinocchio file or Begogna CAT
  // other varialbles in the structure
  double xmin,xmax,ymin,ymax;
  std:: vector<double> x,y;
  double fxmin,fxmax,fymin,fymax;
  int Galaxies;
  std:: string GalModel;
  //  std:: string filecmrelation;
};

// fof properties
struct FoF{
  int nfof;                    // total number of fof
  std:: vector<long> idfof;      
  std:: vector<double> mfof,m200,r200,ra,dec,z,Dc,c200,zsnap,cfof,msmooth,csmooth;
  std:: vector<int> plane;
};

// fof properties
struct PinocchioPLC{
  int nhaloes;                 // pinocchio PLC file structure
  std:: vector<long> id;      
  std:: vector<double> mass,redshift,x_c,y_c,z_c,vx,vy,vz,theta,phi,pec_vel,obs_redshift;
  std:: vector<double> concentration,ra,dec;
};

// GalCat Begogna
struct GalaxyCats{
  int nhaloes;                 // Begonga file structure saved by Carlo in IDL
  std:: vector<long> idh,idg;
  std:: vector<int> iscentral;        
  std:: vector<double> masshalo,massgal,redshift,ra,dec,z;
  std:: vector<double> concentration;
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
      std:: cout << " I will STOP here increase the number of digit in the file name " << std:: endl;
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
	//f->ra.push_back(ra);
	//f->dec.push_back(dec);
	// to match the kappa maps from GLAMER-------
	f->ra.push_back(dec);
	f->dec.push_back(ra);
	// -------------------------------------------
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

// read PLC Pinocchio file
void readPinocchio(struct PinocchioPLC *pinPLC, struct InputParams p){
  std:: ifstream ifilinp;
    ifilinp.open(p.PinocchioFile.c_str());
    if(ifilinp.is_open()){
      for(int i=0;i<58;i++){
	string sbut;
	ifilinp >> sbut;
      }
      long id;
      double mass,redshift,x_c,y_c,z_c,vx,vy,vz,theta,phi,pec_vel,obs_redshift;
      int n=0;
      while(ifilinp >> id 
	    >> redshift
	    >> x_c >> y_c >> z_c 
	    >> vx >> vy >> vz 
	    >> mass 
	    >> theta 
	    >> phi 
	    >> pec_vel
	    >> obs_redshift){
	// std:: cout << id << "   " <<  mass << "   " <<  redshift << "   " <<  x_c << "   " 
	// <<  y_c << "   " <<  z_c << "   " <<  vx << "   " <<  vy << "   " 
	// <<  vz << "   " <<  theta << "   " <<  phi << "   " <<  pec_vel << "   "  << obs_redshift << std:: endl;
	pinPLC->id.push_back(id);
	pinPLC->mass.push_back(mass);
	pinPLC->redshift.push_back(redshift);
	pinPLC->x_c.push_back(x_c);
	pinPLC->y_c.push_back(y_c);
	pinPLC->z_c.push_back(z_c);
	pinPLC->vx.push_back(vx);
	pinPLC->vy.push_back(vy);
	pinPLC->vz.push_back(vz);
	pinPLC->theta.push_back(theta);
	pinPLC->phi.push_back(phi);
	pinPLC->pec_vel.push_back(pec_vel);
	pinPLC->obs_redshift.push_back(obs_redshift);
	double ri = sqrt(x_c*x_c + y_c*y_c + z_c*z_c);
	double xc = ri*sin((90-theta)*M_PI/180.)*cos(phi*M_PI/180.);
	double yc = ri*sin((90-theta)*M_PI/180.)*sin(phi*M_PI/180.);
	double zc = ri*cos((90-theta)*M_PI/180.);
	double dec = asin(xc/ri);
	double ra  = atan(yc/zc);
	pinPLC->ra.push_back(ra);
	pinPLC->dec.push_back(dec);
	n++;
      }
      pinPLC->nhaloes=n;
      std:: cout << " nhaloes in Pinocchio PLC file = " << n << std:: endl;
      std:: cout << "  " << std:: endl;
      //
      //std:: vector<long> id;      
      //std:: vector<double> mass,redshift,x_c,y_c,z_c,vx,vy,vz,theta,phi,pec_vel,obs_redshift;
    }else{
      std:: cout << " the PLC file " << p.PinocchioFile << " does not exist " << std:: endl;
      std:: cout << " I will STOP here!!! " << std:: endl;
      exit(1);
    }
}

// read Begogna Galaxy Cats
void readGalCat(struct GalaxyCats *GalCat, struct InputParams p){
  std:: ifstream ifilinp;
    ifilinp.open(p.PinocchioFile.c_str());
    if(ifilinp.is_open()){
      for(int i=0;i<17;i++){
	string sbut;
	ifilinp >> sbut;
      }
      long idh,idg;
      double masshalo,massgalaxy,redshift,ra,dec;
      int iscentral;
      double dbut;
      int n=0;
      while(ifilinp >> idh >> idg 
	    >> iscentral
	    >> masshalo >> massgalaxy
	    >> ra >> dec 
	    >> redshift
	    >> dbut
	    >> dbut	    
	    >> dbut	    
	    >> dbut
	    >> dbut
	    >> dbut
	    >> dbut
	    >> dbut){	    
	// select only the region interesting to me values are in radiants
	if(ra >= -8.8 && dec >= -8.8){
	  // std:: cout << idh << "   " <<  idg << "   " <<  iscentral << "   " <<  masshalo << "   " 
	  // <<  massgalaxy << "   " <<  ra << "   " <<  dec << "   " <<  redshift << "   " << std:: endl;
	  GalCat->idh.push_back(idh);
	  GalCat->idg.push_back(idg);
	  GalCat->iscentral.push_back(iscentral);	  	  
	  GalCat->masshalo.push_back(masshalo);
	  GalCat->massgal.push_back(massgalaxy);	  
	  GalCat->z.push_back(redshift);
	  GalCat->ra.push_back(ra);
	  GalCat->dec.push_back(dec);
	  n++;
	}
      }
      GalCat->nhaloes=n;
      double minra = *min_element(GalCat->ra.begin(),GalCat->ra.end());
      double maxra = *max_element(GalCat->ra.begin(),GalCat->ra.end());      
      double mindec = *min_element(GalCat->dec.begin(),GalCat->dec.end());
      double maxdec = *max_element(GalCat->dec.begin(),GalCat->dec.end());
      std:: cout << " " << std:: endl;
      std:: cout << " check those values they should be identical " << std:: endl;
      std:: cout << maxra-minra << "  " << maxdec-mindec << std:: endl;
      std:: cout << p.fx << "  " << p.fy << std:: endl;
      std:: cout << " " << std:: endl;      
      for(int i=0;i<n;i++){
	// convert ra and dec in radiants and set (0,0) in the center of the map	
	GalCat->ra[i] =  (GalCat->ra[i] - minra - 0.5*p.fx)*M_PI/180.;
	GalCat->dec[i] = (GalCat->dec[i] - mindec - 0.5*p.fy)*M_PI/180.;
      }
      std:: cout << " nhaloes and galaxies in the Galaxy Cat file = " << n << std:: endl;
      std:: cout << "  " << std:: endl;
    }else{
      std:: cout << " the Calaxy Cats file " << p.PinocchioFile << " does not exist " << std:: endl;
      std:: cout << " I will STOP here!!! " << std:: endl;
      exit(1);
    }
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
      std:: cout << " I will STOP here increase the number of digit in the file name " << std:: endl;
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
  fin >> p->nzs;             //  5. number of source redshifts to consider
  if(p->nzs<=0){
    std:: cout << " number of zs " << p->nzs << std:: endl;
    std:: cout << " check this out, I will STOP here!!! " << std:: endl;
  }
  fin >> str;
  p->zs.resize(p->nzs);
  for(int i=0;i<p->nzs;i++){
    fin >> p->zs[i];              //  5. source redshift (1:nzs)
  }
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
  fin >> str;
  fin >> p->filsigma;        // 15. file with lm s relation to be adopted
  fin >> str;
  fin >> p->sigmalnC;        // 16. scatter in the c-M relation if a c-M model is used
  fin >> str;
  fin >> p->vir;             // 17. vir use FOF haloes
  fin >> str;
  fin >> p->simcase;         // 18. MapSim files or Pinocchio PLC file
  fin >> str;
  fin >> p->PinocchioFile;   // 19. Pinocchio PLC file, needed for Pincchio and BegognaCAT
  fin >> str;
  fin >> p->GalModel;        // 20. model for the galaxies
  //fin >> str;
  //fin >> p->filecmrelation;  // 21. file with the cM relation ... if NO forced to not be considered 
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
  if(p->nzs==1){
    std:: cout << "  source redshift  = " << p->zs[0] << std:: endl;
  }else{
    std:: cout << "  source redshifts  = ";
    for(int i=0;i<p->nzs;i++){
      std:: cout << p->zs[i] << "  ";
    }
    std:: cout << "  " << std:: endl;
  }
  std:: cout << "  Rcut             = " << p->cutR << std:: endl;
  std:: cout << "  path catalogues  = " << p->pathcats << std:: endl;
  std:: cout << "  n planes of cats = " << p->nplanes << std:: endl;
  std:: cout << "  sim name         = " << p->simtype << std:: endl;
  std:: cout << "  plane list file  = " << p->planelist << std:: endl;
  std:: cout << "  nx = " << p->nx << "  ny = " << p->ny << std:: endl;
  std:: cout << "  field of view: x = " << p->fx << "  y = " << p->fy << std:: endl;
  std:: ifstream fillms;
  fillms.open(p->filsigma.c_str());
  if(fillms.is_open()){
    std:: cout << "  file lm s rel    = " << p->filsigma << std:: endl;    
  }else{
    if(p->filsigma=="sim"){
      std:: cout << " no file for the lm s relation ... I will use Vmax and Rmax from the simulation " << std:: endl;
      std:: cout << " since filsigma variable is set equal to sim " << std:: endl;
    }else{
      std:: cout << " no file for the lm s relation ... I will use Neto08+Bullock01 " << std:: endl;
    }
  }
  std:: cout << "  scatter lnC      = " << p->sigmalnC << std:: endl;
  if(p->simcase=="MapSim"){
    if(p->vir=="vir"){
      std:: cout << " I will build the map on the FOF haloes assuming them to be virial " << std:: endl;
    }else{
      std:: cout << " I will build the map on the 200rhoc haloes " << std:: endl;
    }    
    std:: cout <<  " Effective convergence map will be built using MapSim catalogues from cosmological simulations " << std:: endl;
  }else if(p->simcase=="Pinocchio"){
    if(p->vir=="vir"){
      std:: cout << " I will build the map assuming virial masses " << std:: endl;
    }else{
      std:: cout << " I will build the map assuming 200rhoc masses " << std:: endl;
    }    
    std:: cout <<  " Effective convergence map will be built using Pinocchio PLC output file " << std:: endl;
    std:: cout << p->PinocchioFile << std:: endl;
  }else if(p->simcase=="BegognaCAT"){
    if(p->vir=="vir"){
      std:: cout << " I will build the map assuming virial masses " << std:: endl;
    }else{
      std:: cout << " I will build the map assuming 200rhoc masses " << std:: endl;
    }    
    std:: cout <<  " Effective convergence map will be built using the galaxy catalogue file " << std:: endl;
    std:: cout << p->PinocchioFile << std:: endl;  
  }else{
    std:: cout << " wrong parameter in the input file for simcase " << p->simcase  << std:: endl;
    std:: cout << " I will STOP here!!! " << std:: endl;
    exit(1);
  }
  if(p->GalModel == "SIS" || p->GalModel == "HERNQUIST" || p->GalModel == "JAFFE"){
    std:: cout << " if galaxies will be treated I will assume a " << p->GalModel << " profile " << std:: endl;
  }
  /**
  if(p->filecmrelation != "NO"){
    // check first if the file exsist
    if(exists_file (p->filecmrelation)){
      std:: cout << " I will consider " << p->filecmrelation << " for the c-M relation of the haloes " << std:: endl;
    }else{
      std:: cout << p->filecmrelation << " does not exsist I will not consider it! " << std:: endl;
    }
  }
  */
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
  // Bullock+01 redshift evolution
  double cm = 4.67*pow(m/1.e+14,-0.11)/(1+z); 
  return cm;
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

void writeFits(std:: string filename,std:: valarray<float> f, int npix, int npixy,double zs){
  long naxis=2;
  long naxes[2]={npix,npixy};
  std::auto_ptr<FITS> fout(new FITS(filename,FLOAT_IMG,naxis,naxes));
  std::vector<long> naxex(2);
  naxex[0]=npix;
  naxex[1]=npixy;
  PHDU *phout=&fout->pHDU();
  phout->write( 1, npix*npixy, f );
  phout->addKey ("ZSOURCE",zs, "source redshift");  
}


double uNFW(double k, nfwHalo *ha, double r){
  double m200 = ha->getMass();
  double rs = ha->getScaleRadius();
  double m = ha->massin(r);
  double f1 = 4*M_PI*ha->getDensity()*gsl_pow_3(rs)/m;
  double cR = r/rs; // equal to c if r = Rvir
  double Si=gsl_sf_Si(k*rs);
  double Ci=gsl_sf_Ci(k*rs);
  double Sic=gsl_sf_Si((1+cR)*k*rs);
  double Cic=gsl_sf_Ci((1+cR)*k*rs);
  double uk=f1*(sin(k*rs)*(Sic-Si) - sin(cR*k*rs)/(1+cR)/k/rs + cos(k*rs)*(Cic-Ci));
  return uk;
}

void locateHalo(struct InputParams p, double ra, double dec, double dr, int &imin, int &imax, int &jmin, int &jmax){
  // ----> ( CASE 1 )
  if(ra>=p.fxmin && ra<=p.fxmax && dec>=p.fymin && dec<=p.fymax){
    p.xmin = p.xmin>=p.fxmin ? p.xmin : p.fxmin;
    p.xmax = p.xmax<=p.fxmax ? p.xmax : p.fxmax;
    p.ymin = p.ymin>=p.fymin ? p.ymin : p.fymin;
    p.ymax = p.ymax<=p.fymax ? p.ymax : p.fymax;
    if(p.xmin == p.fxmin) imin = 0;
    else{
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }
    if(p.ymin == p.fymin) jmin = 0;
    else{
      jmin = locate(p.y,p.ymin);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }
    if(p.xmax == p.fxmax) imax = p.nx-1; 
    else{
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }
    if(p.ymax == p.fymax) jmax = p.ny-1; 
    else{
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }
  }
  // ----> ( CASE 2 ) 
  if(ra < p.fxmin && dec >= p.fymin && dec <= p.fymax){
    p.xmin = p.fxmin;
    imin = 0;
    p.xmax = ra + dr;
    double dx = p.fxmin - ra;
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dx));
    p.ymin = (dec - dy)>=p.fymin ? (dec - dy) : p.fymin;
    p.ymax = (dec + dr)<=p.fymax ? (dec + dy) : p.fymax;
    if(p.ymin == p.fymin) jmin = 0;
    else{
      jmin = locate(p.y,p.ymin);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }
    if(p.ymax == p.fymax) jmax = p.ny-1;
    else{
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }
    if(p.xmax>=p.fxmin){
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }else{
      imax=0;
    }
  }
  // ----> ( CASE 3 ) 	  
  if(ra > p.fxmax && dec >= p.fymin && dec <= p.fymax){
    p.xmin = ra - dr;
    p.xmax = p.fxmax;
    imax = p.nx - 1;
    double dx = ra - p.fxmax;
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dx));
    p.ymin = (dec - dy)>=p.fymin ? (dec - dy) : p.fymin ;
    p.ymax = (dec + dy)<=p.fymax ? (dec + dy) : p.fymax ;
    if(p.ymin == p.fymin) jmin = 0;
    else{
      jmin = locate(p.y,p.ymin);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }
    if(p.ymax == p.fymax) jmax = p.ny-1;
    else{
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }
    if(p.xmin<=p.fxmax){
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }else{
      imin=p.nx-1;
    }
  } 
  // ----> ( CASE 4 ) 	        
  if(dec < p.fymin && ra >= p.fxmin && ra <=p.fxmax){
    p.ymin = p.fymin;
    jmin = 0;
    p.ymax = dec + dr;
    double dy = p.ymin - dec;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dy));
    p.xmin = (ra - dx)>=p.fxmin ? (ra - dx) : p.fxmin ;
    p.xmax = (ra + dx)<=p.fxmax ? (ra + dx) : p.fxmax ;
    if(p.xmin == p.fxmin) imin = 0;
    else{
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }
    if(p.xmax == p.fxmax) imax = p.nx-1;
    else{
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }
    if(p.ymax>=p.fymin){
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }else{
      jmax = 0;
    }
  }
  // ----> ( CASE 5 ) 	        
  if(dec > p.fymax && ra >= p.fxmin && ra <=p.fxmax){
    p.ymin = dec - dr;
    p.ymax = p.fymax;
    jmax = p.ny-1;
    double dy = dec-p.fymax;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dy));
    p.xmin = (ra - dx)>=p.fxmin ? (ra - dx) : p.fxmin ;
    p.xmax = (ra + dx)<=p.fxmax ? (ra + dx) : p.fxmax ;
    if(p.xmin == p.fxmin) imin = 0;
    else{
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }
    if(p.xmax == p.fxmax) imax = p.nx-1;
    else{
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }
    if(p.ymin<=p.fymax){
      jmin = locate(p.y,p.ymin);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }else{
      jmin = p.ny-1;
    }
  }
  // ---> ( CASE 6 )
  if(ra > p.fxmax && dec > p.fymax){
    p.xmax = p.fxmax;
    p.ymax = p.fymax;
    imax = p.nx-1;
    jmax = p.ny-1;
    double dxcm = ra -p.fxmax;
    double dycm = dec-p.fymax;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
    p.xmin = (ra - dx)>=p.fxmin ? (ra - dx) : p.fxmin ;
    p.ymin = (dec - dy)>=p.fymin ? (dec - dy) : p.fymin ;
    if(p.xmin == p.fxmin) imin = 0;
    else{
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }
    if(p.ymin == p.fymin) jmin = 0;
    else{
      jmin = locate(p.y,p.ymin);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }
  }
  // --> ( CASE 7 )
  if(ra > p.fxmax && dec < p.fymin){
    p.xmax = p.fxmax;
    imax = p.nx-1;
    p.ymin = p.fymin;
    jmin = 0;
    double dycm = p.fymin-dec;
    double dxcm = ra-p.fxmax;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
    p.xmin = (ra - dx)>=p.fxmin ? (ra - dx) : p.fxmin ;
    p.ymax = (dec + dy)<=p.fymax ? (dec + dy) : p.fymax ;
    if(p.xmin == p.fxmin) imin = 0;
    else{
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }
    if(p.ymax == p.fxmax) jmax = p.ny-1;
    else{
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }
  }
  // ---> ( CASE 8 )
  if(ra < p.fxmin && dec < p.fymin){
    p.xmin = p.fxmin;
    p.ymin = p.fymin;
    imin = 0;
    jmin = 0;
    double dxcm = p.fxmin-ra;
    double dycm = p.fymin-dec;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
    p.xmax = (ra + dx)<=p.fxmax ? (ra + dx) : p.fxmax ;
    p.ymax = (dec + dy)<=p.fxmax ? (dec + dy) : p.fymax ;
    if(p.xmax == p.fxmax) imax = p.nx-1;
    else{
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }
    if(p.ymax == p.fymax) jmax = p.ny-1;
    else{
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }
  }
  // ---> ( CASE 9 )
  if(ra < p.fxmin && dec > p.fymax){
    p.xmin = p.fxmin;
    p.ymax = p.fymax;
    imin = 0;
    jmax = p.ny-1;
    double dxcm = p.fxmin-ra;
    double dycm = dec-p.fymax;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
    p.xmax = (ra + dx)<=p.fxmax ? (ra + dx) : p.fxmax ;
    p.ymin = (dec - dy)>=p.fymin ? (dec - dy) : p.fymin ;
    if(p.xmax == p.fxmax) imax = p.nx-1;
    else{
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }
    if(p.ymax == p.fymin) jmin = 0;
    else{
      jmin = locate(p.y,p.ymax);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }
  }
}

int main(){
  time_t start;
  time (&start);  
  // check if there is a file with R,z
  // check if it exists
  std:: string fileRz = "Rz.txt";
  std:: ifstream infileRz;
  infileRz.open(fileRz.c_str());
  std:: vector<double> Ri,zi;
  if(infileRz.is_open()){
    double aa,bb;
    std:: cout << " Rz.txt exists I will read it! " << std:: endl;
    while(infileRz >> aa >> bb){
      Ri.push_back(aa);
      zi.push_back(bb);
    }
    infileRz.close();
    std:: cout << " file read and it has " << Ri.size() << " lines " << std:: endl;    
  }
  bool twohundredc=true;
  Th = gsl_rng_default;
  Thfof = gsl_rng_default;
  rh = gsl_rng_alloc (Th);
  rhfof = gsl_rng_alloc (Thfof);
  gsl_rng_set(rh,1234); 
  gsl_rng_set(rhfof,1234); 
  showStart();
  struct InputParams p;
  readInput(&p);
  time_t tread;
  time (&tread);
  std:: cout << " time to initialize and read the input file " << difftime(tread,start) << " sec" << std:: endl;
  //************** C O S M O L O G Y **********************//
  //**************   using MOKA lib  **********************//
  cosmology co(p.omegam,p.omegal,p.h0,p.wq);
  //*******************************************************//
  time (&tread);
  std:: cout << " time to initialize cosmological model " << difftime(tread,start) << " sec" << std:: endl;  
  struct FoF fof;
  struct Sub sub;
  struct PlaneList pl;
  struct PinocchioPLC pinPLC;
  struct GalaxyCats GalCat;  
  if(p.simcase=="MapSim"){
    readPlaneList(&pl,p);
    readFoF(&fof,p,pl);
    readSubs(&sub,p,pl);
    p.Galaxies=0;
  }
  if(p.simcase=="Pinocchio"){
    readPinocchio(&pinPLC,p);
    p.Galaxies=0;    
  }
  if(p.simcase=="BegognaCAT"){
    readGalCat(&GalCat,p);
    // here we have galaxies!!!
    p.Galaxies=1;    
  }  
  time (&tread);
  std:: cout << " time to read data " << difftime(tread,start) << " sec" << std:: endl;  
  std:: vector<double> lm,s,ls;
  std:: ifstream fillms;
  fillms.open(p.filsigma.c_str());
  if(fillms.is_open()){
    double a,b;
    while(fillms>> a >> b){
      lm.push_back(a);
      s.push_back(b);
      ls.push_back(log10(b));
    }
    iniTables(&co,lm,ls);
  }
  // compare the two cm_relations
  /*
  for(int i=0;i<pl.zl.size();i++){
    std:: vector<double> lx;
    int nx=128;
    fill_linear(lx,nx,12.,15.5);
    double zx = pl.zl[i];
    for(int j=0;j<nx;j++){
      double x = pow(10.,lx[j]);
      double c1 = getCNeto(x,zx);
      double c2 = getCZhao(&co,x,zx);
      std:: cout << i << "  " << zx << "  " << x << "  " << c1 << "  " << c2 << std:: endl;
    }
  }
  exit(1);
  */
  // loop on haloes with z<max zs
  double zsmax = *max_element(p.zs.begin(),p.zs.end());
  double dsmax = co.angularDist(0.0,zsmax);
  std:: vector<double> dlens(p.nzs);  
  //allocate c200
  std:: vector<double> SUBFINDmass, SUBFINDredshift, SUBFINDc;
  std:: vector<double> SUBFINDra, SUBFINDdec, SUBFINDzli;
  if(p.simcase=="MapSim"){
    fof.c200.resize(fof.nfof);
    fof.msmooth.resize(fof.nfof);
    fof.csmooth.resize(fof.nfof);
    fof.cfof.resize(fof.nfof);
    for(int i=0;i<fof.nfof;i++){
      if(fof.zsnap[i]<zsmax){
	// one plane per time
	int planefof = fof.plane[i];
	long idfof = fof.idfof[i];
	if(p.filsigma=="sim"){
	  std:: vector<double> vmax,rmax,msub,rasub,decsub;
	  for(int j=0;j<sub.nsub;j++){
	    if(sub.plane[j] == planefof && sub.idfof[j] == idfof){
	      msub.push_back(sub.msub[j]);
	      vmax.push_back(sub.vmax[j]);
	      rmax.push_back(sub.rmax[j]);
	      rasub.push_back(sub.ra[j]);
	      decsub.push_back(sub.dec[j]);
	      // std:: cout << sub.msub[j] << "   " << sub.vmax[j] << "  " << sub.rmax[j] << std:: endl;
	    }
	  }
	  if(msub.size()>0){
	    std:: vector<double> :: iterator itmostmassive = max_element(msub.begin(), msub.end());      
	    int ihost = distance(msub.begin(),itmostmassive);
	    double vmaxh = vmax[ihost];
	    double rmaxh = rmax[ihost];
	    // use the simulation snapshot 
	    double hz = co.hubble(fof.zsnap[i]);
	    // use Hz
	    double c = getC(vmaxh,rmaxh,hz);
	    // using H0 we have unrealistic redshift evolution ... c too high!
	    // double c = getC(vmaxh,rmaxh,p.h0);
	    fof.c200[i] = c;
	    // std:: cout << fof.m200[i] << "  " << fof.z[i] << "   " << fof.c200[i] << std:: endl;
	    // for the subs
	    for(int j=0;j<msub.size();j++){
	      if(j != ihost){
		SUBFINDmass.push_back(msub[j]);
		SUBFINDredshift.push_back(fof.z[i]);
		SUBFINDc.push_back(4); // not used now we use SIS model
		SUBFINDra.push_back(rasub[j]);
		SUBFINDdec.push_back(decsub[j]);
		SUBFINDzli.push_back(fof.zsnap[i]);
	      }
	    }
	    // smooth component ... first subhalo!
	    fof.msmooth[i] = msub[ihost];
	    // assign the fof mass 
	    fof.msmooth[i] = fof.mfof[i];
	    fof.csmooth[i] = c;
	  }else{
	    // not a sub associated!
	    fof.c200[i] = 0;	    
	  }
	}else{
	  // use the c-M relation models
	  double c,cfof;
	  // if(p.filecmrelation != "NO" && exists_file (p.filecmrelation)){
	  // define c and cfof
	  // means have read the lm s file set in the input paramater file
	  if(lm.size()>1){
	    // Zhat+09 with the parameters of Giocoli+13
	    // we assign the redshift of the snap for the c-M
	    c = getCZhao(&co,fof.m200[i],fof.zsnap[i]);
	    cfof = getCZhao(&co,fof.mfof[i],fof.zsnap[i],p.vir);
	  }else{
	    // Neto+08 cM relation evolving with z as in Bullock+10
	    // we assign the redshift of the snap for the c-M
	    c = getCNeto(fof.m200[i],fof.zsnap[i]);
	  }	  
	  if(p.sigmalnC>1e-2){
	    // add a log-normal scatter!!!	  
	    fof.c200[i] = gsl_ran_lognormal (rh,log(c),p.sigmalnC);	
	    fof.cfof[i] = gsl_ran_lognormal (rhfof,log(cfof),p.sigmalnC);	
	  }else{
	    fof.c200[i] = c;
	    fof.cfof[i] = cfof;
	  }
	  std:: vector<double> vmax,rmax,msub,rasub,decsub;
	  for(int j=0;j<sub.nsub;j++){
	    if(sub.plane[j] == planefof && sub.idfof[j] == idfof){
	      msub.push_back(sub.msub[j]);
	      vmax.push_back(sub.vmax[j]);
	      rmax.push_back(sub.rmax[j]);
	      rasub.push_back(sub.ra[j]);
	      decsub.push_back(sub.dec[j]);
	      // std:: cout << sub.msub[j] << "   " << sub.vmax[j] << "  " << sub.rmax[j] << std:: endl;
	    }
	  }
	  if(msub.size()>0){
	    std:: vector<double> :: iterator itmostmassive = max_element(msub.begin(), msub.end());      
	    int ihost = distance(msub.begin(),itmostmassive);
	    // for the subs
	    double massinsubs = 0;
	    for(int j=0;j<msub.size();j++){
	      if(j != ihost){
		SUBFINDmass.push_back(msub[j]);
		SUBFINDredshift.push_back(fof.z[i]);
		SUBFINDc.push_back(4); // not used now we use SIS model
		SUBFINDra.push_back(rasub[j]);
		SUBFINDdec.push_back(decsub[j]);
		SUBFINDzli.push_back(fof.zsnap[i]);
		massinsubs+=msub[j];
	      }
	    }
	    // smooth component ... first subhalo!...comment the second definition if you want to use this!!!
	    // fof.msmooth[i] = msub[ihost];
	    // smooth component ... mfof-massinsubs minus the smooth component 
	    // fof.msmooth[i] = fof.mfof[i]-massinsubs;
	    //  smooth component = massfof
	    fof.msmooth[i] = fof.mfof[i];
	    // assign the fof concentration 
	    // fof.csmooth[i] = getCZhao(&co,fof.msmooth[i],fof.zsnap[i],p.vir)*fof.cfof[i]/cfof;
	    fof.csmooth[i] = fof.cfof[i];
	  }else{
	    // use the fof if there is no associated subhalo
	    fof.msmooth[i] = fof.mfof[i];
            // assign the fof concentration  
            fof.csmooth[i] = fof.cfof[i];	    
	  }
	}
      }
    }
  }
  if(p.simcase=="Pinocchio"){
    for(int i=0;i<pinPLC.nhaloes;i++){
      // compute the halo concentration 
      // use the c-M relation models
      double c;
      // means have read the lm s file set in the input paramater file
      if(lm.size()>1){
	// Zhat+09 with the parameters of Giocoli+13
	// we assign the redshift of the snap for the c-M
	c = getCZhao(&co,pinPLC.mass[i],pinPLC.redshift[i],p.vir); // if set use virial definition
      }else{
	// Neto+08 cM relation evolving with z as in Bullock+10
	// we assign the redshift of the snap for the c-M
	c = getCNeto(pinPLC.mass[i],pinPLC.redshift[i]);
      }
      if(p.sigmalnC>1e-2){
	// add a log-normal scatter!!!	  
	pinPLC.concentration.push_back(gsl_ran_lognormal (rh,log(c),p.sigmalnC));	
      }else{
	pinPLC.concentration.push_back(c);	
      }
      // std:: cout << pinPLC.mass[i] << "  " << c << "  " << pinPLC.redshift[i] << std:: endl;
    }
  }
  if(p.simcase=="BegognaCAT"){
    for(int i=0;i<GalCat.nhaloes;i++){
      // compute the halo concentration 
      // use the c-M relation models
      double c;
      // means have read the lm s file set in the input paramater file
      if(lm.size()>1){
	// Zhat+09 with the parameters of Giocoli+13
	// we assign the redshift of the snap for the c-M
	c = getCZhao(&co,GalCat.masshalo[i],GalCat.z[i],p.vir); // if set use virial definition
      }else{
	// Neto+08 cM relation evolving with z as in Bullock+10
	// we assign the redshift of the snap for the c-M
	c = getCNeto(GalCat.masshalo[i],GalCat.z[i]);
      }
      if(p.sigmalnC>1e-2){
	// add a log-normal scatter!!!	  
	GalCat.concentration.push_back(gsl_ran_lognormal (rh,log(c),p.sigmalnC));	
      }else{
	GalCat.concentration.push_back(c);	
      }
      // std:: cout << GalCat.masshalo[i] << "  " << c << "  " << GalCat.z[i] << std:: endl;
    }
  }
  std:: cout << " done with the concentrations " << std:: endl;
  time (&tread);
  std:: cout << " until now I spent " << difftime(tread,start) << " sec" << std:: endl;

  double fxrad = p.fx*M_PI/180.;
  double fyrad = p.fy*M_PI/180.;
  p.fxmin = -fxrad/2.;
  p.fxmax =  fxrad/2.;
  p.fymin = -fyrad/2.;
  p.fymax =  fyrad/2.;
  // start to build the map
  std:: valarray<std:: valarray<float> > kappa(p.nzs),kappaGalaxies(p.nzs);
  for(int i=0;i<p.nzs;i++){
    kappa[i].resize(p.nx*p.ny);
    kappaGalaxies[i].resize(p.nx*p.ny);    
  }
  fill_linear(p.x,p.nx,-fxrad/2.,fxrad/2.);
  fill_linear(p.y,p.ny,-fyrad/2.,fyrad/2.);
  //for each lens plane
  //for(int j=0;j<p.nplanes;j++){
  //std:: valarray<float> kappai(p.nx*p.ny);
  int nh;
  if(p.simcase=="MapSim") nh = fof.nfof;
  if(p.simcase=="Pinocchio") nh = pinPLC.nhaloes;
  if(p.simcase=="BegognaCAT") nh = GalCat.nhaloes;  
  time_t loop;
  time (&loop);
  std:: vector<double> Radii(nh);
  for(int i=0;i<nh;i++){
    //if(fof.plane[i]>(j+1)) break;
    halo *ha;
    nfwHalo *nha;
    nfwLens *nLha;
    
    double mass, concentration, redshift;
    double zli;
    double ra,dec;

    double fact=1;
    int iscentral = 1;
    
    if(p.vir=="vir") twohundredc=false;
    if(p.vir=="subs") twohundredc=false;
    if(p.simcase=="MapSim"){
      if(p.vir=="vir"){
	mass = fof.mfof[i];
	concentration = fof.cfof[i];
      }else{
	// assign to the host halo mass m200
	mass = fof.m200[i];
	concentration = fof.c200[i];
	if(p.vir=="subs"){
	  // assign to the host halo mass the smooth component 
	  mass = fof.msmooth[i];
	  concentration = fof.csmooth[i];
	  // modify kappa
	  // mass = fof.mfof[i];
	  // fact = fof.msmooth[i]/fof.mfof[i];
	}
      }
      zli=fof.zsnap[i];
      redshift=fof.z[i];
      ra=fof.ra[i];
      dec=fof.dec[i];
    }
    if(p.simcase=="Pinocchio"){
      mass = pinPLC.mass[i];
      concentration = pinPLC.concentration[i];
      zli=pinPLC.redshift[i];
      redshift=zli;
      ra=pinPLC.ra[i];
      dec=pinPLC.dec[i];
    }
    if(p.simcase=="BegognaCAT"){
      mass = GalCat.masshalo[i];
      concentration = GalCat.concentration[i];
      zli=GalCat.z[i];
      redshift=zli;
      ra=GalCat.ra[i];
      dec=GalCat.dec[i];
      iscentral = GalCat.iscentral[i];
    }
    // if there is info on c200 (means a sub associated) and m200 (means a m200 halo associated)
    //if(mass>0 && concentration>0 && fof.zsnap[i]<p.zs){
    time (&tread);
    std:: cout << " start halo " << i << "/" <<  nh << "  " << difftime(tread,start) << " sec" << std:: endl;
    // std:: cout << " I should finish in " << difftime(tread,loop)/i*(nh-i)/60./60. << " hours" << std:: endl;    
    if(mass>0 && concentration>0 && zli<zsmax){
      ha = new halo(&co,mass,redshift,twohundredc);
      std:: cout << mass << "  " << redshift << "  " << concentration << "  " << ra << "  " << dec << std:: endl;
      nha = new nfwHalo(*ha,concentration);
      nLha = new nfwLens(*nha,zsmax);
      double rs = nLha->getScaleRadius();
      double Radius = nLha->getRvir();
      Radii[i] = Radius;
      // for the underling DM I use only central haloes ... no subhaloes for now
      if(iscentral == 1){
	// is in unit of the scale radius
	// TEST redshift evolution!!!
	// double Rz = Radius*fabs(p.cutR/pow(1+redshift,0.5))/rs;
	double Rz = Radius*fabs(p.cutR)/rs;      
	double fR=-1;
	if(Ri.size()>0){
	  fR = getY(zi,Ri,redshift);
	  Rz = Radius*fR/rs;      
	}
	int xi0=locate(p.x,ra);
	xi0=GSL_MIN( GSL_MAX( xi0, 0 ), p.nx-1 );
	int yi0=locate(p.y,dec);
	yi0=GSL_MIN( GSL_MAX( yi0, 0 ), p.ny-1 );
	double Dl = (co.angularDist(0.,redshift)*co.cn.lightspeed*1.e-7);
	// TEST redshift evolution!!!
	// double dr = Radius*fabs(p.cutR/pow(1+redshift,0.5))/Dl;
	double dr = Radius*fabs(p.cutR)/Dl;      
	if(fR>0){
	  dr = Radius*fR/Dl;      
	}
	p.xmin = ra - dr;
	p.xmax = ra + dr;
	p.ymin = dec - dr;
	p.ymax = dec + dr;
	int imin=0, imax=p.nx-1, jmin=0, jmax=p.ny-1;
	// std:: cout << "  Dl " << Dl << std:: endl;
	// check the buffer !!!      
	locateHalo(p,ra,dec,dr,imin,imax,jmin,jmax);
	// outputs should be jmin, jmax, imin, imax
	if(imin<0) imin =0;
	if(jmin<0) jmin =0;
	for(int ij=0;ij<p.nzs;ij++){
	  dlens[ij] = (co.angularDist(zli,p.zs[ij])/co.angularDist(0.0,p.zs[ij]))/(co.angularDist(zli,zsmax)/dsmax);
	}
	for(int jj=jmin; jj<=jmax; jj++ ) for(int ii=imin; ii<=imax; ii++ ){
	    double dx=(p.x[ii]-ra);
	    double dy=(p.y[jj]-dec);
	    double r=sqrt( dx*dx+dy*dy ); 
	    if(r<=dr){
	      if(p.cutR>0){
		for(int ij=0;ij<p.nzs;ij++){
		  // this need to be in Mpc/h physical!!!
		  if(zli<p.zs[ij]){
		    kappa[ij][ii+p.nx*jj]+=(nLha->kappaz(Dl*r/rs,Rz)*dlens[ij]);
		  }
		}
	      }else{
		for(int ij=0;ij<p.nzs;ij++){		
		  // this need to be in Mpc/h physical!!!
		  // if negative consider the integral up to infinity
		  if(zli<p.zs[ij]){		
		    kappa[ij][ii+p.nx*jj]+=(nLha->kappa(Dl*r/rs)*dlens[ij]);
		  }
		}
	      }
	    }
	  }       
	delete ha;
	delete nha;
	delete nLha;
      }
    }
  }

  if(p.GalModel == "SIS" || p.GalModel == "HERNQUIST" || p.GalModel == "JAFFE"){
    std:: cout << "  " << std:: endl;
    std:: cout << " I will now consider the galaxies " << std:: endl;
    std:: cout << "  " << std:: endl;    
  }else{
    p.Galaxies = 0;
  }
  
  if(p.simcase=="MapSim" && p.vir=="subs"){
    nh = SUBFINDmass.size();
    p.Galaxies=1;
    // only this implemented for MapSim now!!!
    p.GalModel = "SIS";
  }
  // I can create a separate map for the Galaxy components here!!!
  if(p.Galaxies==1){
    for(int i=0;i<nh;i++){
      
      hernq_Halo *hgalhalo;      
      hernq_Lens *hgallens;

      jaffe_Halo *jgalhalo;      
      jaffe_Lens *jgallens;
      
      halo *shalo;
      sisHalo *sgalhalo;      
      sisLens *sgallens;

      double mass, concentration, redshift, massgalaxy;
      double zli;
      double ra,dec;
      
      if(p.simcase=="MapSim"){
	mass = SUBFINDmass[i];
	concentration = SUBFINDc[i];
	redshift = SUBFINDredshift[i];
	ra=SUBFINDra[i];
	dec=SUBFINDdec[i];
	zli=SUBFINDzli[i];
	massgalaxy = mass;
      }
      if(p.simcase=="Pinocchio"){
	std:: cout << " galaxies in Pinocchio not implemented yet " << std:: endl;
	std:: cout << " I will STOP here!!! " << std:: endl;
	exit(1);
	mass = pinPLC.mass[i];
	concentration = pinPLC.concentration[i];
	zli=pinPLC.redshift[i];
	redshift=zli;
	ra=pinPLC.ra[i];
	dec=pinPLC.dec[i];
      }
      if(p.simcase=="BegognaCAT"){
	mass = GalCat.masshalo[i];
	concentration = GalCat.concentration[i];
	zli=GalCat.z[i];
	redshift=zli;
	ra=GalCat.ra[i];
	dec=GalCat.dec[i];
	massgalaxy = GalCat.massgal[i];	
      }
      
      time (&tread);
      std:: cout << " start Galaxy/Subs " << i << "/" <<  nh << "  " << difftime(tread,start) << " sec" << std:: endl;
      // std:: cout << " I should finish in " << difftime(tread,loop)/i*(nh-i)/60./60. << " hours" << std:: endl;
      std:: cout << mass << "  " << redshift << "  " << massgalaxy << "  " << ra << "  " << dec << std:: endl;
      
      if(mass>0 && concentration>0 && zli<zsmax){

	double Radius = Radii[i];	
	double rs;
	if(p.GalModel == "SIS"){
	  shalo = new halo(&co,massgalaxy,redshift,twohundredc);
	  bool indvir;	
	  if(twohundredc) indvir = false;
	  else indvir = true;
	  sgalhalo = new sisHalo(*shalo,indvir);
	  sgallens = new sisLens(*sgalhalo,zsmax);	  
	  rs = sgallens->getScale();
	  // Radius enclosing the mass ... important for SIS
	  Radius = sgallens->getRadius();  	  
	}
	if(p.GalModel == "HERNQUIST"){
	  hgalhalo = new hernq_Halo(&co,massgalaxy,mass,redshift);
	  hgallens = new hernq_Lens(*hgalhalo,zsmax);
	  rs = hgallens->getscaleRadius();
	}
	if(p.GalModel == "JAFFE"){
	  jgalhalo = new jaffe_Halo(&co,massgalaxy,mass,redshift);
	  jgallens = new jaffe_Lens(*jgalhalo,zsmax);
	  rs = jgallens->getscaleRadius();
	}
	
	// is in unit of the scale radius
	// TEST redshift evolution!!!
	// * * * Warning kappaz not implemented yet in the Galaxy convergence profile * * * //
	// double Rz = Radius*fabs(p.cutR/pow(1+redshift,0.5))/rs;
	double Rz = Radius*fabs(p.cutR)/rs;      
	double fR=-1;
	if(Ri.size()>0){
	  fR = getY(zi,Ri,redshift);
	  Rz = Radius*fR/rs;      
	}
	int xi0=locate(p.x,ra);
	xi0=GSL_MIN( GSL_MAX( xi0, 0 ), p.nx-1 );
	int yi0=locate(p.y,dec);
	yi0=GSL_MIN( GSL_MAX( yi0, 0 ), p.ny-1 );
	double Dl = (co.angularDist(0.,redshift)*co.cn.lightspeed*1.e-7);
	// TEST redshift evolution!!!
	// double dr = Radius*fabs(p.cutR/pow(1+redshift,0.5))/Dl;
	double dr = Radius*fabs(p.cutR)/Dl;      
	if(fR>0){
	  dr = Radius*fR/Dl;      
	}
	p.xmin = ra - dr;
	p.xmax = ra + dr;
	p.ymin = dec - dr;
	p.ymax = dec + dr;
	int imin=0, imax=p.nx-1, jmin=0, jmax=p.ny-1;
	// std:: cout << "  Dl " << Dl << std:: endl;
	// check the buffer !!!      
	locateHalo(p,ra,dec,dr,imin,imax,jmin,jmax);
	// outputs should be jmin, jmax, imin, imax
	if(imin<0) imin =0;
	if(jmin<0) jmin =0;
	for(int ij=0;ij<p.nzs;ij++){
	  dlens[ij] = (co.angularDist(zli,p.zs[ij])/co.angularDist(0.0,p.zs[ij]))/(co.angularDist(zli,zsmax)/dsmax);
	}
	for(int jj=jmin; jj<=jmax; jj++ ) for(int ii=imin; ii<=imax; ii++ ){
	    double dx=(p.x[ii]-ra);
	    double dy=(p.y[jj]-dec);
	    double r=sqrt( dx*dx+dy*dy ); 
	    if(r<=dr){
	      for(int ij=0;ij<p.nzs;ij++){
		// this need to be in Mpc/h physical!!!
		if(zli<p.zs[ij]){
		  if(p.GalModel == "SIS") kappaGalaxies[ij][ii+p.nx*jj]+=sgallens->kappa(Dl*r/rs)*dlens[ij];
		  if(p.GalModel == "HERNQUIST") kappaGalaxies[ij][ii+p.nx*jj]+=hgallens->kappa(Dl*r/rs)*dlens[ij];
		  if(p.GalModel == "JAFFE") kappaGalaxies[ij][ii+p.nx*jj]+=jgallens->kappa(Dl*r/rs)*dlens[ij];
		}
	      }
	    }
	  }
	if(p.GalModel == "SIS"){
	  delete sgalhalo;
	  delete sgallens;
	}
	if(p.GalModel == "HERNQUIST"){
	  delete hgalhalo;
	  delete hgallens;
	}
	if(p.GalModel == "JAFFE"){
	  delete jgalhalo;
	  delete jgallens;
	}
      }  
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int ij=0;ij<p.nzs;ij++){
    // kmGal
    float kmG = kappaGalaxies[ij].sum()/(p.nx*p.ny);
    for(int jj=0; jj<p.ny; jj++ ) for(int ii=0; ii<p.nx; ii++ ){    
	// test
	kappaGalaxies[ij][ii+p.nx*jj] = kappaGalaxies[ij][ii+p.nx*jj] -kmG;
      }    
    // km
    float km = (kappa[ij].sum() + kappaGalaxies[ij].sum())/(p.nx*p.ny);
    for(int jj=0; jj<p.ny; jj++ ) for(int ii=0; ii<p.nx; ii++ ){    
	// test
	kappa[ij][ii+p.nx*jj] = (kappa[ij][ii+p.nx*jj] + kappaGalaxies[ij][ii+p.nx*jj] -km);
      }
    std:: string szid;    
    std:: ostringstream oszid;
    oszid << ij;
    szid = oszid.str();
    std:: string filout = "kappa_test_" + szid + ".fits";
    writeFits(filout,kappa[ij],p.nx,p.ny,p.zs[ij]);
    int nb = 128;    
    if(p.Galaxies==1){
      filout = "kappaGalaxies_test_" + szid + ".fits";
      writeFits(filout,kappaGalaxies[ij],p.nx,p.ny,p.zs[ij]);
      double *ll;
      double *Pl;
      ll=new double[nb];
      Pl=new double[nb];
      powerl(kappaGalaxies[ij],kappaGalaxies[ij],p.nx,p.ny,fxrad,fyrad,ll,Pl,nb);  
      std:: ofstream outfile;
      outfile.open("mapPowerSpectrumGalaxies_" + szid+ ".dat");
      for(int i=0;i<nb;i++){
	outfile << ll[i] << "  " << Pl[i] << std::endl;
      } 
      outfile.close();      
    }
    double *ll;
    double *Pl;
    ll=new double[nb];
    Pl=new double[nb];
    powerl(kappa[ij],kappa[ij],p.nx,p.ny,fxrad,fyrad,ll,Pl,nb);  
    std:: ofstream outfile;
    outfile.open("mapPowerSpectrum_" + szid+ ".dat");
    for(int i=0;i<nb;i++){
      outfile << ll[i] << "  " << Pl[i] << std::endl;
    } 
    outfile.close();
  }
    return 0;
}
