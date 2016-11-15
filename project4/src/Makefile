CC = mpic++
CFLAGS = -Wall -O3 -c -std=c++11 
LFLAGS = 
LIBDIR = lib
OBJS = project4.o
LIBOBJS = ${LIBDIR}/lib.o
PROG = project4

all: ${PROG}
            
${PROG}: ${OBJS} ${LIBOBJS}
	${CC} ${LFLAGS} ${OBJS} ${LIBOBJS} -o ${PROG}

${LIBDIR}/lib.o: ${LIBDIR}/lib.cpp ${LIBDIR}/lib.h
	${CC} ${CFLAGS} -I${LIBDIR} ${LIBDIR}/lib.cpp -o ${LIBDIR}/lib.o

project4.o: project4.cpp project4.h ${LIBDIR}/lib.h
	${CC} ${CFLAGS} -I${LIBDIR} project4.cpp

clean:
	rm -f *.o  ${LIBDIR}/*.o ${PROG} 
