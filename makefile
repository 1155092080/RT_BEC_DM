# BEC-GP-ROT-OMP programs are developed by:
#
# R. Kishor Kumar
# (Instituto de Fisica, Universidade de Sao Paulo, Sao Paulo, Brazil)
#
# Vladimir Loncar, Antun Balaz
# (Scientific Computing Laboratory, Center for the Study of Complex Systems, Institute of Physics Belgrade, Serbia)
#
# Paulsamy Muruganandam
# (Department of Physics, Bharathidasan University, Tiruchirappalli, Tamil Nadu, India)
#
# Sadhan K. Adhikari
# (Instituto de Fisica Teorica, UNESP - Sao Paulo State University, Brazil)
#
# Public use and modification of these codes are allowed provided that the following papers are cited:
# [1] R. Kishor Kumar et al., Comput. Phys. Commun. 240 (2019) 74.
# [2] P. Muruganandam and S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
# [3] L. E. Young-S. et al., Comput. Phys. Commun. 220 (2017) 503.
# [4] D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
# [5] R. Kishor Kumar et al., Comput. Phys. Commun. 195 (2015) 117.
# [6] B. Sataric et al., Comput. Phys. Commun. 200 (2016) 411.
# [7] V. Loncar et al., Comput. Phys. Commun. 200 (2016) 406.
# [8] L. E. Young-S. et al., Comput. Phys. Commun. 204 (2016) 209.
# [9] V. Loncar et al., Comput. Phys. Commun. 209 (2016) 190.
#
# The authors would be grateful for all information and/or comments
# regarding the use of the programs.

ifndef $(compiler)
    compiler = gnu
endif

ifeq ($(compiler), gnu)
   CC = gcc
   CCFLAGS = -std=c99 -O3
   OMPFLAGS = -fopenmp
endif

ifeq ($(compiler), intel)
   CC = icc
   CCFLAGS = -std=c99 -O3 -xHOST -no-prec-div -fno-alias
   OMPFLAGS = -qopenmp
#If older version of Intel C compiler is present, the option -qopenmp in the
#line above should be replaced by -openmp
endif

ifeq ($(compiler), pgi)
   CC = pgcc
   CCFLAGS = -c99 -fast -O3
   OMPFLAGS = -mp=allcores
endif

ifeq ($(compiler), clang)
   CC = clang
   CCFLAGS = -std=c99 -O3
   OMPFLAGS = -fopenmp
endif

ifeq ($(compiler), oracle)
   CC = suncc
   CCFLAGS = -std=c99 -fast
   OMPFLAGS = -fopenmp
#If suncc cannot find system header files in /usr/include (e.g., on Ubuntu),
#add -I/usr/include/x86_64-linux-gnu/ to the CCFLAGS line
endif

LIBS = -lm -lrt

UNAME := $(shell uname -s)

ifeq ($(UNAME), Darwin)
    LIBS = -lm
endif

UTIL = diffint mem bec_rand cfg out visit_writer timer
UTIL_OBJECTS=$(UTIL:=.o)

all: bec-gp-rot-2d-th bec-gp-rot-3d-th
	rm -rf *.o

help: README.md
	less $^

bec-gp-rot-2d-th: $(UTIL_OBJECTS)
	$(CC) $(CCFLAGS) $(OMPFLAGS) -c src/bec-gp-rot-2d-th.c -o bec-gp-rot-2d-th.o
	$(CC) $(CCFLAGS) $(OMPFLAGS) $(UTIL_OBJECTS) bec-gp-rot-2d-th.o -o bec-gp-rot-2d-th $(LIBS)
	rm bec-gp-rot-2d-th.o

bec-gp-rot-3d-th: $(UTIL_OBJECTS)
	$(CC) $(CCFLAGS) $(OMPFLAGS) -c src/bec-gp-rot-3d-th.c -o bec-gp-rot-3d-th.o
	$(CC) $(CCFLAGS) $(OMPFLAGS) $(UTIL_OBJECTS) bec-gp-rot-3d-th.o -o bec-gp-rot-3d-th $(LIBS)
	rm bec-gp-rot-3d-th.o

$(UTIL_OBJECTS): %.o:
	$(CC) $(CCFLAGS) -c src/util/$*.c -o $@

clean:
	rm -rf *.o bec-gp-rot-2d-th  bec-gp-rot-3d-th
