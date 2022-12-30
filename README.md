# BEC-GP-ROT-OMP

BEC-GP-ROT-OMP is a set of OpenMP-parallelized C (C99) and Fortran 90 programs, optimized for GNU, Intel, PGI, Clang and Oracle C and Fortran compilers, that solves the time-(in)dependent Gross-Pitaevskii nonlinear partial differential equation for rotating BECs with contact interaction in two and three spatial dimensions in a trap using imaginary-time and real-time propagations, chosen by selecting the parameter `OPTREIM = 1` (imaginary time propagation) and `OPTREIM = 2` (real time propagation). The Gross-Pitaevskii equation describes the properties of a dilute trapped Bose-Einstein condensate. The equation is solved using the split-step Crank-Nicolson method by discretizing space and time, as described in Ref. [1]. The discretized equation is then propagated in imaginary or real time over small time steps. This document describes the C programs of BEC-GP-ROT-OMP code distribution.

## Description of BEC-GP-ROT-OMP C code distribution

### I) Source codes

Programs are written in C99 revision of C programming language with OpenMP extensions and are located in the [src](src/) folder, which has the following structure:

 - [src/bec-gp-rot-2d-th.c](src/bec-gp-rot-2d-th.c) code solves the imaginary- and real-time Gross-Pitaevskii equations in two space dimensions for the dynamics in x-y plane in an anisotropic harmonic trap.

 - [src/bec-gp-rot-3d-th.c](src/bec-gp-rot-3d-th.c) code solves the imaginary- and real-time Gross-Pitaevskii equations in three space dimensions in an anisotropic harmonic trap.

### II) Input parameters

For each program, a specific parameter input file has to be given on the command line. Examples of input files with characteristic set of options and their detailed descriptions are provided in the [input](input/) folder.

### III) Examples of matching outputs

The [output](output/) folder contains examples of matching outputs for all programs and default inputs available in the BEC-GP-ROT-OMP distribution. Some large density files are omitted to save space.

### IV) Compilation

Programs are compiled via a provided `makefile`.

The use of the makefile:

    make <target> [compiler=<compiler>]

where possible targets are:

    all, clean, help

as well as program-specific targets, which compile only a specified program:

    bec-gp-rot-2d-th, bec-gp-rot-3d-th.

The provided makefile allows compilation of the BEC-GP-ROT-OMP programs, using either the GNU, Intel, PGI, Clang or Oracle C compiler.

Possible `<compiler>` values are:

    gnu, intel, pgi, clang, oracle

and should be chosen according to the available compilers on the compile host. If compiler is not specified on the command line, the default GNU C compiler is used. If an older version of the Intel compiler is used, compilation may fail due to the use of `-qopenmp` switch, which used to be `-openmp` in earlier versions. In that case, modify the makefile to use `-openmp` switch in the variable OMPFLAGS. If Oracle compiler is used on Ubuntu Linux, `-I/usr/include/x86_64-linux-gnu/` should be added to CCFLAGS variable.

The default optimization option set in the makefile is -O3, and can be adjusted if needed in the makefile variable OMPFLAGS.

**Examples compiling the programs:**

1. Compile all BEC-GP-ROT-OMP programs with the GNU C compiler:

        make all
   or

        make all compiler=gnu

2. Compile `bec-gp-rot-2d-th` program with the Intel C compiler (if available):

        make bec-gp-rot-2d-th compiler=intel

### V) Running the compiled programs

To run any of the programs compiled with the make command, you need to use the syntax:

    ./<programname> -i <parameterfile>

where `<programname>` is a name of the compiled executable, and `<parameterfile>` is a parameter input file prepared by the user. Examples of parameter input files are described in section II above, and are provided in the folder [input](input/). Matching output of the principal output files are given in folder [output](output/); very large density output files are omitted.

Example of running a program:

    ./bec-gp-rot-3d-th -i input/imag3d-input

Note that you can specify the number of OpenMP threads by setting the `OMP_NUM_THREADS` environment variable:

    export OMP_NUM_THREADS=<number>

where `<number>` is the number of threads wanted, before the execution of the program. If the `OMP_NUM_THREAD` is not specified, the OpenMP runtime will try to use the same number of threads as there are available CPU cores on the running host.

### VI) Authors

**BEC-GP-ROT-OMP** programs are developed by:

R. Kishor Kumar *(Instituto de Fisica,  Universidade de Sao Paulo, Sao Paulo, Brazil)*  
Vladimir Loncar and Antun Balaz *(Scientific Computing Laboratory, Center for the Study of Complex Systems, Institute of Physics Belgrade, Serbia)*  
Paulsamy Muruganandam *(Bharathidasan University, Tamil Nadu, India)*  
Sadhan K. Adhikari *(UNESP - Sao Paulo State University, Brazil)*  

Public use and modification of these codes are allowed provided that the following papers are cited:  
[1] [R. Kishor Kumar et al., Comput. Phys. Commun. 240 (2019) 74.](https://doi.org/10.1016/j.cpc.2019.03.004)  
[2] [P. Muruganandam and S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.](https://doi.org/10.1016/j.cpc.2009.04.015)  
[3] [L. E. Young-S. et al., Comput. Phys. Commun. 220 (2017) 503.](https://doi.org/10.1016/j.cpc.2017.07.013)  
[4] [D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.](https://doi.org/10.1016/j.cpc.2012.03.022)  
[5] [R. Kishor Kumar et al., Comput. Phys. Commun. 195 (2015) 117.](https://doi.org/10.1016/j.cpc.2015.03.024)  
[6] [B. Sataric et al., Comput. Phys. Commun. 200 (2016) 411.](https://doi.org/10.1016/j.cpc.2015.12.006)  
[7] [V. Loncar et al., Comput. Phys. Commun. 200 (2016) 406.](https://doi.org/10.1016/j.cpc.2015.11.014)  
[8] [L. E. Young-S. et al., Comput. Phys. Commun. 204 (2016) 209.](https://doi.org/10.1016/j.cpc.2016.03.015)  
[9] [V. Loncar et al., Comput. Phys. Commun. 209 (2016) 190.](https://doi.org/10.1016/j.cpc.2016.07.029)  

The authors would be grateful for all information and/or comments regarding the use of the programs.

### VII) Licence

**BEC-GP-ROT-OMP** code distribution is licensed under Apache License, Version 2.0. See [LICENCE](LICENCE) for details. Portions of code related to writing the output in VTK file format are copyrighted by Lawrence Livermore National Security, LLC, see [visit_writer.h](src/util/visit_writer.h)/[.c](src/util/visit_writer.c) for details.
