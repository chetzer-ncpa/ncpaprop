//utility functions
#include "Atmosphere.h"
#include "ProcessOptionsNB.h"
double **dmatrix(long nr, long nc);
int      free_dmatrix(double**v, long nr, long nc);
//int      plotwGNUplot_from_script();
//int      plotwGNUplot(double freq, bool write_2D_TLoss, bool Nby2Dprop, bool write_phase_speeds);
//int      plotwGNUplot(NCPA::ProcessOptionsNB *oNB);
//int saveAtm_profile(NCPA::SampledProfile *p);
int saveAtm_profile(NCPA::SampledProfile *p, std::string wind_units);


