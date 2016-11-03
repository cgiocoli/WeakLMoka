#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gsl/gsl_integration.h> 
#include "../Moka/utilities.h"
#include "../Moka/cosmology.h"

double getS(double mi);
void iniTables(cosmology *cosmology, std:: vector<double> lm, std:: vector<double> ls);
double getCZhao(cosmology *cosmology, double mass, double redshift, std:: string vir="200");
 
