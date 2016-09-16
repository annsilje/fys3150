CC = g++
CFLAGS = -Wall -c
LFLAGS = 
LIBDIR = lib
OBJS = project1.o
LIBOBJS = ${LIBDIR}/lib.o
PROG = project1

all: ${PROG}
            
${PROG}: ${OBJS} ${LIBOBJS}
	${CC} ${LFLAGS} ${OBJS} ${LIBOBJS} -o ${PROG}

lib.o: ${LIBDIR}/lib.cpp ${LIBDIR}/lib.h
	${CC} ${CFLAGS} -I${LIBDIR} ${LIBDIR}/lib.cpp 

project1.o: project1.cpp ${LIBDIR}/lib.h
	${CC} ${CFLAGS} -I${LIBDIR} project1.cpp

clean:
	rm -f *.o  ${LIBDIR}/*.o ${PROG} 
