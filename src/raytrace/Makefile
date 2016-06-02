.PHONY: clean tidy install raytrace

all: raytrace

raytrace: raytrace.2d raytrace.3d

raytrace.2d: Acoustic2DEquationSet.o AcousticEquationSet.o AtmosphericBreakConditions.o GeneralBreakConditions.o ODESystemBreakCondition.o ODESystem.o ReflectionCondition2D.o raytrace.2d.o
	g++ -o $@ $^  -L/usr/lib ../../lib/libatmosphere.a ../../lib/libcommon.a -lgsl -lgslcblas -lm -lfftw3 
	cp $@ ../../bin

raytrace.3d: Acoustic3DEquationSet.o AcousticEquationSet.o AtmosphericBreakConditions.o GeneralBreakConditions.o ODESystemBreakCondition.o ODESystem.o ReflectionCondition3D.o raytrace.3d.o
	g++ -o $@ $^  -L/usr/lib ../../lib/libatmosphere.a ../../lib/libcommon.a -lgsl -lgslcblas -lm -lfftw3 
	cp $@ ../../bin


%.o: %.cpp
	g++ -fpic -c -ggdb -Wall -I. -I../common -I../atmosphere -I/usr/local/include -o $@ $<

clean:
	rm -rf raytrace.2d raytrace.3d Acoustic2DEquationSet.o AcousticEquationSet.o AtmosphericBreakConditions.o GeneralBreakConditions.o ODESystemBreakCondition.o ODESystem.o ReflectionCondition2D.o raytrace.2d.o Acoustic3DEquationSet.o ReflectionCondition3D.o raytrace.2d.o raytrace.3d.o

tidy:
	rm -rf Acoustic2DEquationSet.o AcousticEquationSet.o AtmosphericBreakConditions.o GeneralBreakConditions.o ODESystemBreakCondition.o ODESystem.o ReflectionCondition2D.o raytrace.2d.o Acoustic3DEquationSet.o ReflectionCondition3D.o raytrace.2d.o raytrace.3d.o
