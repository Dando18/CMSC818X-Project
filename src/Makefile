CC = gcc
MPICC = mpicxx
WARNINGS = -Wall -pedantic
STD = -std=c++11
OPT = -O3
FLAGS = $(WARNINGS) $(STD) $(OPT)

METIS_SRC = ../metis-5.1.0/programs/io.c ../metis-5.1.0/programs/stat.c ../metis-5.1.0/programs/cmdline_gpmetis.c

TARGET = ../bin/run

all: $(TARGET)

test: test.c 
	$(CC) -o $@ $< ${METIS_SRC} --std=c11 -I../metis-5.1.0/programs -I../metis/include -I../metis-5.1.0/GKlib -L../metis/lib -lmetis -lm

$(TARGET): main.C
	if [ ! -d "../bin" ]; then mkdir -p ../bin; fi;
	$(MPICC) $(FLAGS) -o $@ $< 


clean:
	rm $(TARGET)