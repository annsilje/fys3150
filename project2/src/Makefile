CC = g++
CFLAGS = -Wall -c
LFLAGS = -O2 -larmadillo
OBJS = project2.o
PROG = project2

all: ${PROG}
            
${PROG}: ${OBJS} 
	${CC} ${OBJS} -o ${PROG} ${LFLAGS} 


project2.o: project2.cpp 
	${CC} ${CFLAGS} project2.cpp

clean:
	rm -f *.o ${PROG} 
