CC= g++
CFLAGS=  -Wno-write-strings

MAINS_DIR= Code/
ATMO_DIR=  Code/Atmo/
GEOAC_DIR= Code/GeoAc/

INSTALL_DIR= /usr/local/bin

all: GeoAc2D GeoAc3D GeoAc3D.RngDep GeoAcGlobal GeoAcGlobal.RngDep

GeoAc2D:
	${CC}  ${ATMO_DIR}G2S_Spline1D.cpp ${ATMO_DIR}Atmo_State.Absorption.cpp ${GEOAC_DIR}GeoAc.Parameters.cpp ${GEOAC_DIR}GeoAc.EquationSets.2DStratified.cpp ${GEOAC_DIR}GeoAc.Solver.cpp ${GEOAC_DIR}GeoAc.Interface.cpp ${MAINS_DIR}GeoAc2D_main.cpp ${CFLAGS} -o GeoAc2D-ncpaprop

GeoAc3D:
	${CC}  ${ATMO_DIR}G2S_Spline1D.cpp ${ATMO_DIR}Atmo_State.Absorption.cpp ${GEOAC_DIR}GeoAc.Parameters.cpp ${GEOAC_DIR}GeoAc.Eigenray.cpp ${GEOAC_DIR}GeoAc.EquationSets.3DStratified.cpp ${GEOAC_DIR}GeoAc.Solver.cpp ${GEOAC_DIR}GeoAc.Interface.cpp ${MAINS_DIR}GeoAc3D_main.cpp ${CFLAGS} -o GeoAc3D-ncpaprop 

GeoAc3D.RngDep:
	${CC}  ${ATMO_DIR}G2S_MultiDimSpline3D.cpp ${ATMO_DIR}Atmo_State.Absorption.cpp ${GEOAC_DIR}GeoAc.Parameters.RngDep.cpp ${GEOAC_DIR}GeoAc.Eigenray.cpp ${GEOAC_DIR}GeoAc.EquationSets.3DRngDep.cpp ${GEOAC_DIR}GeoAc.Solver.cpp ${GEOAC_DIR}GeoAc.Interface.cpp ${MAINS_DIR}GeoAc3D.RngDep_main.cpp ${CFLAGS} -o GeoAc3D.RngDep-ncpaprop 

GeoAcGlobal:
	${CC}  ${ATMO_DIR}G2S_GlobalSpline1D.cpp ${ATMO_DIR}Atmo_State.Absorption.Global.cpp ${GEOAC_DIR}GeoAc.Parameters.Global.cpp ${GEOAC_DIR}GeoAc.Eigenray.Global.cpp ${GEOAC_DIR}GeoAc.EquationSets.Global.cpp ${GEOAC_DIR}GeoAc.Solver.cpp ${GEOAC_DIR}GeoAc.Interface.Global.cpp ${MAINS_DIR}GeoAcGlobal_main.cpp ${CFLAGS} -o GeoAcGlobal-ncpaprop

GeoAcGlobal.RngDep:
	${CC}  ${ATMO_DIR}G2S_GlobalMultiDimSpline3D.cpp ${ATMO_DIR}Atmo_State.Absorption.Global.cpp ${GEOAC_DIR}GeoAc.Parameters.Global.cpp ${GEOAC_DIR}GeoAc.Eigenray.Global.cpp ${GEOAC_DIR}GeoAc.EquationSets.GlobalRngDep.cpp ${GEOAC_DIR}GeoAc.Solver.cpp ${GEOAC_DIR}GeoAc.Interface.Global.cpp ${MAINS_DIR}GeoAcGlobal.RngDep_main.cpp ${CFLAGS} -o GeoAcGlobal.RngDep-ncpaprop

all: GeoAc2D GeoAc3D GeoAc3D.RngDep GeoAcGlobal GeoAcGlobal.RngDep

clean: 
	rm GeoAc2D-ncpaprop GeoAc3D-ncpaprop GeoAc3D.RngDep-ncpaprop GeoAcGlobal-ncpaprop GeoAcGlobal.RngDep-ncpaprop

install: GeoAc2D GeoAc3D GeoAc3D.RngDep GeoAcGlobal GeoAcGlobal.RngDep
	install -m 0755 GeoAc2D-ncpaprop ${INSTALL_DIR}
	install -m 0755 GeoAc3D-ncpaprop ${INSTALL_DIR}
	install -m 0755 GeoAc3D.RngDep-ncpaprop ${INSTALL_DIR}
	install -m 0755 GeoAcGlobal-ncpaprop ${INSTALL_DIR}
	install -m 0755 GeoAcGlobal.RngDep-ncpaprop ${INSTALL_DIR}

uninstall:
	rm ${INSTALL_DIR}/GeoAc2D-ncpaprop 
	rm ${INSTALL_DIR}/GeoAc3D-ncpaprop 
	rm ${INSTALL_DIR}/GeoAc3D.RngDep-ncpaprop 
	rm ${INSTALL_DIR}/GeoAcGlobal-ncpaprop 
	rm ${INSTALL_DIR}/GeoAcGlobal.RngDep-ncpaprop 
	rm GeoAc2D-ncpaprop GeoAc3D-ncpaprop GeoAc3D.RngDep-ncpaprop GeoAcGlobal-ncpaprop GeoAcGlobal.RngDep-ncpaprop
