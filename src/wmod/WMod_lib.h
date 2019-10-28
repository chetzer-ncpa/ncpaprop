//utility functions
#ifndef __NCPA_WMOD_LIB_H__

#define __NCPA_WMOD_LIB_H__

#include "Atmosphere.h"
#include "ProcessOptionsNB.h"


double **dmatrix(long nr, long nc);
int      free_dmatrix(double**v, long nr, long nc);
int saveAtm_profile(NCPA::SampledProfile *p, std::string wind_units);


#endif