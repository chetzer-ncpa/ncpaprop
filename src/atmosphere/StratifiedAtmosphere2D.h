#ifndef NCPAPROP_STRATIFIEDATMOSPHERE2D_H_INCLUDED
#define NCPAPROP_STRATIFIEDATMOSPHERE2D_H_INCLUDED

#include "Atmosphere2D.h"


namespace NCPA {
	class StratifiedAtmosphere2D : public Atmosphere2D {

	public:
		StratifiedAtmosphere2D( const Atmosphere1D *atm );
		StratifiedAtmosphere2D( const std::string &filename );
		//StratifiedAtmosphere2D( const StratifiedAtmosphere2D &atm );
		~StratifiedAtmosphere2D();

	};
}




#endif