METIS_SRC = ../metis-5.1.0/programs/io.c ../metis-5.1.0/programs/stat.c ../metis-5.1.0/programs/cmdline_gpmetis.c

all: test

test: test.c 
	$(CC) -o $@ $< ${METIS_SRC} --std=c11 -I../metis-5.1.0/programs -I../metis/include -I../metis-5.1.0/GKlib -L../metis/lib -lmetis -lm