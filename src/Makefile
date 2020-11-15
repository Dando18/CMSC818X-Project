CC = gcc
CXX = g++
MPICC = mpicxx
WARNINGS = -w # -Wall -pedantic -Wno-unused-variable -g
STD = -std=c++11
OPT = -O0 -g
FLAGS = $(WARNINGS) $(STD) $(OPT)

METISLIB_INC = -I../metis-5.1.0/programs -I../metis-5.1.0/libmetis/ -I../metis-5.1.0/GKlib -I../metis-5.1.0/include
INC = $(METISLIB_INC)
LD = -L../metis/lib
METIS_SRC = ../metis-5.1.0/programs/io.c ../metis-5.1.0/programs/stat.c ../metis-5.1.0/programs/cmdline_gpmetis.c

# TARGET = ../bin/run
TARGET = partition

all: $(TARGET)

partition: partition.cpp 
	$(MPICC) $(FLAGS) -o $@ $< --std=c++11 $(METISLIB_INC) $(LD) -lmetis -lm

../bin/run: main.C
	if [ ! -d "../bin" ]; then mkdir -p ../bin; fi;
	$(MPICC) $(FLAGS) $(INC) -o $@ $< 


clean:
	rm $(TARGET)