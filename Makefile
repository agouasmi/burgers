OBJECTS_MAIN = main.o utilities.o 
OBJECTS_CONVO = test_convo.o utilities.o
#SOURCES = main.cpp utilities.cpp utilities.h 
CFLAGS = -c -std=c++11 -I/usr/include
LFLAGS = -lcblas -lfftw3 -lm 
CC = g++


	
burgers: ${OBJECTS_MAIN}
	${CC} ${OBJECTS_MAIN} ${LFLAGS} -o burgers

test_convo: ${OBJECTS_CONVO}
	${CC} ${OBJECTS_CONVO}  ${LFLAGS} -o test_convo
	

main.o: main.cpp utilities.h
	${CC} ${CFLAGS} main.cpp utilities.h	

utilities.o: utilities.cpp utilities.h
	${CC} ${CFLAGS} utilities.cpp utilities.h	

test_convo.o: test_convo.cpp utilities.h
	${CC} ${CFLAGS} test_convo.cpp utilities.h
	
#config.o: test2.cpp utilities.h
#	${CC} ${CFLAGS} test3.cpp utilities.h

clean:
	rm -r ${OBJECTS}
	
run:
	./burgers

#test: test2.cpp utilities.cpp utilities.h 
#	g++ -std=c++11 -o test test2.cpp utilities.cpp utilities.h -lfftw3 -lm

#run: test
#	./test