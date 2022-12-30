/**
 * BEC-GP-ROT-OMP programs are developed by:
 *
 * R. Kishor Kumar
 * (Instituto de Fisica, Universidade de Sao Paulo, Sao Paulo, Brazil)
 *
 * Vladimir Loncar, Antun Balaz
 * (Scientific Computing Laboratory, Center for the Study of Complex Systems, Institute of Physics Belgrade, Serbia)
 *
 * Paulsamy Muruganandam
 * (Department of Physics, Bharathidasan University, Tiruchirappalli, Tamil Nadu, India)
 *
 * Sadhan K. Adhikari
 * (Instituto de Fisica Teorica, UNESP - Sao Paulo State University, Brazil)
 *
 * Public use and modification of these codes are allowed provided that the following papers are cited:
 * [1] R. Kishor Kumar et al., Comput. Phys. Commun. 240 (2019) 74.
 * [2] P. Muruganandam and S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
 * [3] L. E. Young-S. et al., Comput. Phys. Commun. 220 (2017) 503.
 * [4] D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
 * [5] R. Kishor Kumar et al., Comput. Phys. Commun. 195 (2015) 117.
 * [6] B. Sataric et al., Comput. Phys. Commun. 200 (2016) 411.
 * [7] V. Loncar et al., Comput. Phys. Commun. 200 (2016) 406.
 * [8] L. E. Young-S. et al., Comput. Phys. Commun. 204 (2016) 209.
 * [9] V. Loncar et al., Comput. Phys. Commun. 209 (2016) 190.
 *
 * The authors would be grateful for all information and/or comments
 * regarding the use of the programs.
 */
#ifndef OUT_H
#define OUT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "mem.h"
#include "visit_writer.h"

#define VTK_BINARY         1

void out_curve_txt(char *filename, int N, int step, double *f, double *curve);
void out_curve_bin(char *filename, int N, int step, double *f, double *curve);
void out_curve_vis(char *filename, int N, int step, double *f, double *curve, char *name);

void out_matrix_txt(char *filename, int Ny, int Nx, int stepy, int stepx, double *y, double *x, double **matrix);
void out_matrix_bin(char *filename, int Ny, int Nx, int stepy, int stepx, double *y, double *x, double **matrix);
void out_matrix_vis(char *filename, int Ny, int Nx, int stepy, int stepx, double *y, double *x, double **matrix, char *name);

void out_tensor_txt(char *filename, int Nz, int Ny, int Nx, int stepz, int stepy, int stepx, double *z, double *y, double *x, double ***tensor);
void out_tensor_bin(char *filename, int Nz, int Ny, int Nx, int stepz, int stepy, int stepx, double *z, double *y, double *x, double ***tensor);
void out_tensor_vis(char *filename, int Nz, int Ny, int Nx, int stepz, int stepy, int stepx, double *z, double *y, double *x, double ***tensor, char *name);

void out_complex_matrix_txt(char *filename, int Ny, int Nx, double complex **matrix);
void out_complex_matrix_bin(char *filename, int Ny, int Nx, double complex **matrix);

void out_complex_tensor_txt(char *filename, int Nz, int Ny, int Nx, double complex ***tensor);
void out_complex_tensor_bin(char *filename, int Nz, int Ny, int Nx, double complex ***tensor);

#endif /* OUT_H */
