#ifndef   	POWER_H_
# define   	POWER_H_
#include <valarray>
#include <fftw3.h>
#include <gsl/gsl_sort_double.h>
#include <vector>
#include <iostream>
#include "../Moka/utilities.h"

/** 
 * created by:  Carlo Giocoli, ZAH-ITA Heidelberg, 2010; INAF-OABO Bologna,  2011 - (carlo.giocoli@unibo.it)
 */
void powerl(std::valarray<float> a,std::valarray<float> b,
	int nx, int ny, double boxlx, double boxly, 
	double *ll, double *Pl, int nl);

#endif
