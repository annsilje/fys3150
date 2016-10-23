CC = g++
CFLAGS = -std=c++11 -Wall -Wextra -c 
LFLAGS = 
OBJS = main.o examples.o system.o vec3.o particle.o potential.o nopotential.o newtoniangravity.o modifiednewtoniangravity.o integrator.o velocityverlet.o initialcondition.o twobody.o threebody.o forwardeuler.o solarsystem.o
PROG = project3

all: ${PROG}
            
${PROG}: ${OBJS} 
	${CC} ${OBJS} -o ${PROG} ${LFLAGS} 

solarsystem.o: InitialConditions/solarsystem.cpp InitialConditions/solarsystem.h InitialConditions/initialcondition.h system.h vec3.h
	${CC} ${CFLAGS} -IInitialConditions InitialConditions/solarsystem.cpp

threebody.o: InitialConditions/threebody.cpp InitialConditions/threebody.h InitialConditions/initialcondition.h system.h vec3.h
	${CC} ${CFLAGS} -IInitialConditions InitialConditions/threebody.cpp

twobody.o: InitialConditions/twobody.cpp InitialConditions/twobody.h InitialConditions/initialcondition.h system.h vec3.h
	${CC} ${CFLAGS} -IInitialConditions InitialConditions/twobody.cpp

initialcondition.o: InitialConditions/initialcondition.cpp InitialConditions/initialcondition.h particle.h
	${CC} ${CFLAGS} -IInitialConditions InitialConditions/initialcondition.cpp

velocityverlet.o: Integrators/velocityverlet.cpp Integrators/velocityverlet.h Integrators/integrator.h particle.h system.h
	${CC} ${CFLAGS} -IIntegrators Integrators/velocityverlet.cpp

forwardeuler.o: Integrators/forwardeuler.cpp Integrators/forwardeuler.h Integrators/integrator.h particle.h system.h
	${CC} ${CFLAGS} -IIntegrators Integrators/forwardeuler.cpp
	
integrator.o: Integrators/integrator.cpp Integrators/integrator.h particle.h system.h
	${CC} ${CFLAGS} -IIntegrators Integrators/integrator.cpp	

modifiednewtoniangravity.o: Potentials/modifiednewtoniangravity.cpp Potentials/modifiednewtoniangravity.h Potentials/potential.h particle.h
	${CC} ${CFLAGS} -IPotentials Potentials/modifiednewtoniangravity.cpp

newtoniangravity.o: Potentials/newtoniangravity.cpp Potentials/newtoniangravity.h Potentials/potential.h particle.h
	${CC} ${CFLAGS} -IPotentials Potentials/newtoniangravity.cpp

nopotential.o: Potentials/nopotential.cpp Potentials/nopotential.h Potentials/potential.h
	${CC} ${CFLAGS} -IPotentials Potentials/nopotential.cpp

potential.o: Potentials/potential.cpp Potentials/potential.h particle.h
	${CC} ${CFLAGS} -IPotentials Potentials/potential.cpp

particle.o: particle.cpp particle.h vec3.h
	${CC} ${CFLAGS} particle.cpp

vec3.o: vec3.cpp vec3.h
	${CC} ${CFLAGS} vec3.cpp

system.o: system.cpp system.h particle.h Integrators/integrator.h Potentials/potential.h InitialConditions/initialcondition.h
	${CC} ${CFLAGS} system.cpp

examples.o: examples.cpp examples.h system.h particle.h Integrators/velocityverlet.h Potentials/newtoniangravity.h Potentials/nopotential.h InitialConditions/twobody.h InitialConditions/threebody.h InitialConditions/solarsystem.h
	${CC} ${CFLAGS} examples.cpp

main.o: main.cpp examples.h 
	${CC} ${CFLAGS} main.cpp

clean:
	rm -f *.o ${PROG} InitialConditions/*.o Potentials/*.o Integrators/*.o 
