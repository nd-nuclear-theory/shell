#Simple makefile to keep track of all the files needed for the correct compilation of the code
CC = g++
CCFLAGS =  -I/usr/common/usg/gsl/1.16/gnu/include -L/usr/common/usg/gsl/1.16/gnu/lib -lgsl -lgslcblas 
#CCFLAGS = -I /usr/common/usg/gsl/1.16/pgi/include -L /usr/common/usg/gsl/1.16/pgi/lib -lgsl -lgslcblas

# Semptember 20. The code wasn't compiling so I did
# module swap PrgEnv-pgi PrgEnv-gnu
# module load gsl
# module show gsl
# This gives the correct path to the libraries..

# March 3rd. New gsl on Edison..
# -I/usr/common/usg/gsl/1.16/gnu/include -L/usr/common/usg/gsl/1.16/gnu/lib -lgsl -lgslcblas
# Same as old Edison..

all: emcalc emcalc_lj

emcalc_lj: halfint.o emcalc_lj.o halfint_bound.o angular_momentum.o mfdn_file_processing.o shell_radial.o shell_indexing_lj.o
	$(CC) -o emcalc_lj halfint.o emcalc_lj.o halfint_bound.o angular_momentum.o mfdn_file_processing.o shell_radial.o shell_indexing_lj.o $(CCFLAGS)

emcalc: halfint.o emcalc.o halfint_bound.o angular_momentum.o mfdn_file_processing.o shell_radial.o shell_indexing_lj.o
	$(CC) -o emcalc halfint.o emcalc.o halfint_bound.o angular_momentum.o mfdn_file_processing.o shell_radial.o shell_indexing_lj.o $(CCFLAGS) 

halfint.o: halfint.cpp halfint.h
	$(CC) -c halfint.cpp

halfint_bound.o: halfint_bound.cpp halfint_bound.h
	$(CC) -c halfint_bound.cpp

emcalc.o: emcalc.cpp shell_radial.h halfint.h wigner_gsl.h angular_momentum.h mfdn_file_processing.h
	$(CC) -c emcalc.cpp $(CCFLAGS) 

emcalc_lj.o: emcalc_lj.cpp shell_radial.h halfint.h wigner_gsl.h angular_momentum.h mfdn_file_processing.h
	$(CC) -c emcalc_lj.cpp $(CCFLAGS)

angular_momentum.o: angular_momentum.cpp angular_momentum.h
	$(CC) -c angular_momentum.cpp

mfdn_file_processing.o: mfdn_file_processing.cpp mfdn_file_processing.h halfint.h
	$(CC) -c mfdn_file_processing.cpp

shell_radial.o: shell_radial.cpp shell_radial.h shell_indexing_lj.h halfint.h angular_momentum.h
	$(CC) -c shell_radial.cpp 

shell_indexing_lj.o: shell_indexing_lj.cpp shell_indexing_lj.h 
	$(CC) -c shell_indexing_lj.cpp

clean:
	rm *.o emcalc emcalc_lj
