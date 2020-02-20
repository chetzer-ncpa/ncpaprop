.PHONY: clean tidy install raytrace

#CXX=mpicxx
#CXX_LINKER=mpicxx
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

INCPATHS = -I. -I../common -I../atmosphere -I/usr/local/include -I/usr/include ${PETSC_CC_INCLUDES} ${SLEPC_CC_INCLUDES}


all: raytrace

raytrace: raytrace.2d raytrace.3d

raytrace.2d: Acoustic2DEquationSet.o AcousticEquationSet.o AtmosphericBreakConditions.o GeneralBreakConditions.o ODESystemBreakCondition.o ODESystem.o ReflectionCondition2D.o raytrace.2d.o
	${CXX_LINKER} -o $@ $^  -L/usr/lib ../../lib/libatmosphere.a ../../lib/libcommon.a -lgsl -lgslcblas -lm -lfftw3 
	cp $@ ../../bin

raytrace.3d: Acoustic3DEquationSet.o AcousticEquationSet.o AtmosphericBreakConditions.o GeneralBreakConditions.o ODESystemBreakCondition.o ODESystem.o ReflectionCondition3D.o raytrace.3d.o
	${CXX_LINKER} -o $@ $^  -L/usr/lib ../../lib/libatmosphere.a ../../lib/libcommon.a -lgsl -lgslcblas -lm -lfftw3 
	cp $@ ../../bin


%.o: %.cpp
	${CXX} ${INCPATHS} -fpic -c -ggdb -Wall -I. -I../common -I../atmosphere -I/usr/local/include -o $@ $<

clean::
	rm -rf raytrace.2d raytrace.3d Acoustic2DEquationSet.o AcousticEquationSet.o AtmosphericBreakConditions.o GeneralBreakConditions.o ODESystemBreakCondition.o ODESystem.o ReflectionCondition2D.o raytrace.2d.o Acoustic3DEquationSet.o ReflectionCondition3D.o raytrace.2d.o raytrace.3d.o

tidy:
	rm -rf Acoustic2DEquationSet.o AcousticEquationSet.o AtmosphericBreakConditions.o GeneralBreakConditions.o ODESystemBreakCondition.o ODESystem.o ReflectionCondition2D.o raytrace.2d.o Acoustic3DEquationSet.o ReflectionCondition3D.o raytrace.2d.o raytrace.3d.o
