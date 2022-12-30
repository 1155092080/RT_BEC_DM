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
#include "mem.h"

/**
 *    Double vector allocation
 */
double *alloc_double_vector(int Nx) {
   double *vector;

   if ((vector = (double *) malloc((size_t) (Nx * sizeof(double)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the vector.\n");
      exit(EXIT_FAILURE);
   }

   return vector;
}

/**
 *    Complex vector allocation
 */
double complex *alloc_complex_vector(int Nx) {
   double complex *vector;

   if ((vector = (double complex *) malloc((size_t) (Nx * sizeof(double complex)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the vector.\n");
      exit(EXIT_FAILURE);
   }

   return vector;
}

/**
 *    Double matrix allocation
 */
double **alloc_double_matrix(int Nx, int Ny) {
   double **matrix;

   if ((matrix = (double **) malloc((size_t) (Nx * sizeof(double *)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the matrix.\n");
      exit(EXIT_FAILURE);
   }
   if ((matrix[0] = (double *) malloc((size_t) (Nx * Ny * sizeof(double)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the matrix.\n");
      exit(EXIT_FAILURE);
   }
   for (int i = 1; i < Nx; i ++)
      matrix[i] = matrix[i - 1] + Ny;

   return matrix;
}

/**
 *    Complex matrix allocation
 */
double complex **alloc_complex_matrix(int Nx, int Ny) {
   double complex **matrix;

   if ((matrix = (double complex **) malloc((size_t) (Nx * sizeof(double complex *)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the matrix.\n");
      exit(EXIT_FAILURE);
   }
   if ((matrix[0] = (double complex *) malloc((size_t) (Nx * Ny * sizeof(double complex)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the matrix.\n");
      exit(EXIT_FAILURE);
   }
   for (int i = 1; i < Nx; i ++)
      matrix[i] = matrix[i - 1] + Ny;

   return matrix;
}

/**
 *    Double tensor allocation
 */
double ***alloc_double_tensor(int Nx, int Ny, int Nz) {
   double ***tensor;

   if ((tensor = (double ***) malloc((size_t) (Nx * sizeof(double **)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the tensor.\n");
      exit(EXIT_FAILURE);
   }
   if ((tensor[0] = (double **) malloc((size_t) (Nx * Ny * sizeof(double *)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the tensor.\n");
      exit(EXIT_FAILURE);
   }
   if ((tensor[0][0] = (double *) malloc((size_t) (Nx * Ny * Nz * sizeof(double)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the tensor.\n");
      exit(EXIT_FAILURE);
   }
   for (int j = 1; j < Ny; j ++)
      tensor[0][j] = tensor[0][j - 1] + Nz;
   for (int i = 1; i < Nx; i ++) {
      tensor[i] = tensor[i - 1] + Ny;
      tensor[i][0] = tensor[i - 1][0] + Ny * Nz;
      for (int j = 1; j < Ny; j ++)
         tensor[i][j] = tensor[i][j - 1] + Nz;
   }

   return tensor;
}

/**
 *    Complex tensor allocation
 */
double complex ***alloc_complex_tensor(int Nx, int Ny, int Nz) {
   double complex ***tensor;

   if ((tensor = (double complex ***) malloc((size_t) (Nx * sizeof(double complex **)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the tensor.\n");
      exit(EXIT_FAILURE);
   }
   if ((tensor[0] = (double complex **) malloc((size_t) (Nx * Ny * sizeof(double complex *)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the tensor.\n");
      exit(EXIT_FAILURE);
   }
   if ((tensor[0][0] = (double complex *) malloc((size_t) (Nx * Ny * Nz * sizeof(double complex)))) == NULL) {
      fprintf(stderr, "Failed to allocate memory for the tensor.\n");
      exit(EXIT_FAILURE);
   }
   for (int j = 1; j < Ny; j ++)
      tensor[0][j] = tensor[0][j - 1] + Nz;
   for (int i = 1; i < Nx; i ++) {
      tensor[i] = tensor[i - 1] + Ny;
      tensor[i][0] = tensor[i - 1][0] + Ny * Nz;
      for (int j = 1; j < Ny; j ++)
         tensor[i][j] = tensor[i][j - 1] + Nz;
   }

   return tensor;
}

/**
 *    Free double vector
 */
void free_double_vector(double *vector) {
   free((char *) vector);
}

/**
 *    Free complex vector
 */
void free_complex_vector(double complex *vector) {
   free((char *) vector);
}

/**
 *    Free double matrix
 */
void free_double_matrix(double **matrix) {
   free((char *) matrix[0]);
   free((char *) matrix);
}

/**
 *    Free complex matrix
 */
void free_complex_matrix(double complex **matrix) {
   free((char *) matrix[0]);
   free((char *) matrix);
}

/**
 *    Free double tensor
 */
void free_double_tensor(double ***tensor) {
   free((char *) tensor[0][0]);
   free((char *) tensor[0]);
   free((char *) tensor);
}

/**
 *    Free complex tensor
 */
void free_complex_tensor(double complex ***tensor) {
   free((char *) tensor[0][0]);
   free((char *) tensor[0]);
   free((char *) tensor);
}
