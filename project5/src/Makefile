CC = g++
CFLAGS = -Wall -O3 -c -std=c++11 
LFLAGS = 
LIBDIR = lib
OBJS = project5.o
LIBOBJS = ${LIBDIR}/lib.o
PROG = project5

all: ${PROG}
            
${PROG}: ${OBJS} ${LIBOBJS}
	${CC} ${LFLAGS} ${OBJS} ${LIBOBJS} -o ${PROG}

${LIBDIR}/lib.o: ${LIBDIR}/lib.cpp ${LIBDIR}/lib.h
	${CC} ${CFLAGS} -I${LIBDIR} ${LIBDIR}/lib.cpp -o ${LIBDIR}/lib.o

project5.o: project5.cpp project5.h ${LIBDIR}/lib.h
	${CC} ${CFLAGS} -I${LIBDIR} project5.cpp

clean:
	rm -f *.o  ${LIBDIR}/*.o ${PROG} 
