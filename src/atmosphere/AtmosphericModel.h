#ifndef NCPAPROP_ATMOSPHERICMODEL_H_INCLUDED
#define NCPAPROP_ATMOSPHERICMODEL_H_INCLUDED

#include <string>

namespace NCPA {

	class AtmosphericModel {
	public:
		virtual void calculate_sound_speed_from_temperature( std::string new_key, std::string temperature_key ) = 0;
		virtual void calculate_sound_speed_from_pressure_and_density( std::string new_key, std::string pressure_key, std::string density_key ) = 0;
	};

}










#endif
