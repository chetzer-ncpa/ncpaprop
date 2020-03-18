#ifndef NCPAPROP_ATMOSPHERE_H_INCLUDED
#define NCPAPROP_ATMOSPHERE_H_INCLUDED



namespace NCPA {

	class Atmosphere {
	public:
		void calculate_sound_speed_from_temperature( std::string new_key, std::string temperature_key ) = 0;
		void calculate_sound_speed_from_pressure_and_density( std::string new_key, std::string pressure_key, std::string density_key ) = 0;
	};

}










#endif