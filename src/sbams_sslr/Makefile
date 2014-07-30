main: main.o controller.o parser.o SSLR.o
	g++  -O3 main.o parser.o controller.o SSLR.o -lm -L /usr/local/lib -lgsl -lgslcblas -o sbams_sslr
main.o: main.cc
	g++ -c main.cc
parser.o: parser.cc parser.h
	g++ -c parser.cc
controller.o: controller.cc controller.h
	g++ -c controller.cc
SSLR.o: SSLR.h SSLR.cc
	g++ -c SSLR.cc
clean:
	rm *.o sbams_sslr
