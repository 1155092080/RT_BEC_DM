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
#include "out.h"

static char full_filename[FILENAME_MAX];

static void add_ext(char *filename, char *ext) {
    if (strstr(filename, ext) != NULL) strcpy(full_filename, filename);
    else sprintf(full_filename, "%s.%s", filename, ext);
}

void out_curve_txt(char *filename, int N, int step, double *f, double *curve) {
    FILE *file;

    add_ext(filename, "txt");
    file = fopen(full_filename, "w");

    for(int i = 0; i < N; i += step) {
        fprintf(file, "%8le %8le\n", f[i], curve[i]);
    }

    fflush(file);
    fclose(file);
}

void out_curve_bin(char *filename, int N, int step, double *f, double *curve) {
    FILE *file;
    double tmp[2];

    add_ext(filename, "bin");
    file = fopen(full_filename, "wb");

    for(int i = 0; i < N; i += step) {
        tmp[0] = f[i];
        tmp[1] = curve[i];
        fwrite(&tmp, sizeof(double), 2, file);
    }

    fflush(file);
    fclose(file);
}

void out_curve_vis(char *filename, int N, int step, double *f, double *curve, char *name) {
    FILE *file;

    add_ext(filename, "ult");
    file = fopen(full_filename, "w");

    fprintf(file, "# %s\n", name);

    for(int i = 0; i < N; i += step) {
        fprintf(file, "%8le %8le\n", f[i], curve[i]);
    }

    fflush(file);
    fclose(file);
}

void out_matrix_txt(char *filename, int Ny, int Nx, int stepy, int stepx, double *y, double *x, double **matrix) {
    FILE *file;

    add_ext(filename, "txt");
    file = fopen(full_filename, "w");

    for(int j = 0; j < Ny; j += stepy) {
        for(int i = 0; i < Nx; i += stepx) {
            fprintf(file, "%8le %8le %8le\n", x[j], y[i], matrix[j][i]);
        }
        fprintf(file, "\n");
        fflush(file);
    }

    fclose(file);
}

void out_matrix_bin(char *filename, int Ny, int Nx, int stepy, int stepx, double *y, double *x, double **matrix) {
    FILE *file;
    double tmp[3];

    add_ext(filename, "bin");
    file = fopen(full_filename, "wb");

    for(int j = 0; j < Ny; j += stepy) {
        tmp[0] = y[j];
        for(int i = 0; i < Nx; i += stepx) {
            tmp[1] = x[i];
            tmp[2] = matrix[j][i];
            fwrite(&tmp, sizeof(double), 3, file);
        }
        fflush(file);
    }

    fclose(file);
}

void out_matrix_vis(char *filename, int Ny, int Nx, int stepy, int stepx, double *y, double *x, double **matrix, char *name) {
    int vardims[] = { 1 };
    int centering[] = { 1 };
    const char *varnames[] = { name };
    double axis0[] = { 0 };

    add_ext(filename, "vtk");

    if (stepx > 1 || stepy > 1) {
        double *sx, *sy, **smatrix;
        int modx, mody;

        modx = Nx % stepx;
        mody = Ny % stepy;

        sx = alloc_double_vector(Nx / stepx);
        sy = alloc_double_vector(Ny / stepy);
        smatrix = alloc_double_matrix(Ny / stepy, Nx / stepx);

        for (int i = 0; i < Nx - modx; i += stepx) {
            sx[i / stepx] = x[i];
        }

        for(int j = 0; j < Ny - mody; j += stepy) {
            sy[j / stepy] = y[j];
        }

        for(int j = 0; j < Ny - mody; j += stepy) {
            for (int i = 0; i < Nx - modx; i += stepx) {
                smatrix[j / stepy][i / stepx] = matrix[j][i];
            }
        }

        int dimxy[] = { Nx / stepx, Ny / stepy, 1 };
        double *varxy[] = { *smatrix };
        write_rectilinear_mesh(full_filename, VTK_BINARY, dimxy, sx, sy, axis0, 1, vardims, centering, varnames, varxy);

        free_double_vector(sx);
        free_double_vector(sy);
        free_double_matrix(smatrix);
    } else {
        int dimxy[] = { Nx, Ny, 1 };
        double *varxy[] = { *matrix };
        write_rectilinear_mesh(full_filename, VTK_BINARY, dimxy, x, y, axis0, 1, vardims, centering, varnames, varxy);
    }
}

void out_tensor_txt(char *filename, int Nz, int Ny, int Nx, int stepz, int stepy, int stepx, double *z, double *y, double *x, double ***tensor) {
    FILE *file;

    add_ext(filename, "txt");
    file = fopen(full_filename, "w");

    for(int k = 0; k < Nz; k += stepz) {
        for(int j = 0; j < Ny; j += stepy) {
            for(int i = 0; i < Nx; i += stepx) {
                fprintf(file, "%8le %8le %8le %8le\n", z[k], y[j], x[i], tensor[k][j][i]);
                //fprintf(file, "%8le\n", tensor[k][j][i]);
            }
            fprintf(file, "\n");
            fflush(file);
        }
        fprintf(file, "\n");
        fflush(file);
    }

    fclose(file);
}

void out_tensor_bin(char *filename, int Nz, int Ny, int Nx, int stepz, int stepy, int stepx, double *z, double *y, double *x, double ***tensor) {
    FILE *file;
    double tmp[4];

    add_ext(filename, "bin");
    file = fopen(full_filename, "wb");

    for(int k = 0; k < Nz; k += stepz) {
        tmp[0] = z[k];
        for(int j = 0; j < Ny; j += stepy) {
            tmp[1] = y[j];
            for(int i = 0; i < Nx; i += stepx) {
                tmp[2] = x[i];
                tmp[3] = tensor[k][j][i];
                fwrite(&tmp, sizeof(double), 4, file);
                //fwrite(&tmp[3], sizeof(double), 1, file);
            }
            fflush(file);
        }
    }

    fclose(file);
}

void out_tensor_vis(char *filename, int Nz, int Ny, int Nx, int stepz, int stepy, int stepx, double *z, double *y, double *x, double ***tensor, char *name) {
    int vardims[] = { 1 };
    int centering[] = { 1 };
    const char *varnames[] = { name };

    add_ext(filename, "vtk");

    if (stepx > 1 || stepy > 1 || stepz > 1) {
        double *sx, *sy, *sz, ***stensor;
        int modx, mody, modz;

        modx = Nx % stepx;
        mody = Ny % stepy;
        modz = Nz % stepz;

        sx = alloc_double_vector(Nx / stepx);
        sy = alloc_double_vector(Ny / stepy);
        sz = alloc_double_vector(Nz / stepz);
        stensor = alloc_double_tensor(Nz / stepz, Ny / stepy, Nx / stepx);

        for (int i = 0; i < Nx - modx; i += stepx) {
            sx[i / stepx] = x[i];
        }

        for(int j = 0; j < Ny - mody; j += stepy) {
            sy[j / stepy] = y[j];
        }

        for(int k = 0; k < Nz - modz; k += stepz) {
            sz[k / stepz] = z[k];
        }

        for(int k = 0; k < Nz - modz; k += stepz) {
            for(int j = 0; j < Ny - mody; j += stepy) {
                for (int i = 0; i < Nx - modx; i += stepx) {
                    stensor[k / stepz][j / stepy][i / stepx] = tensor[k][j][i];
                }
            }
        }

        int dimxyz[] = { Nx / stepx, Ny / stepy, Nz / stepz };
        double *varxyz[] = { **stensor };
        write_rectilinear_mesh(full_filename, VTK_BINARY, dimxyz, sx, sy, sz, 1, vardims, centering, varnames, varxyz);

        free_double_vector(sx);
        free_double_vector(sy);
        free_double_vector(sz);
        free_double_tensor(stensor);
    } else {
        int dimxyz[] = { Nx, Ny, Nz };
        double *varxyz[] = { **tensor };
        write_rectilinear_mesh(full_filename, VTK_BINARY, dimxyz, x, y, z, 1, vardims, centering, varnames, varxyz);
    }
}

void out_complex_matrix_txt(char *filename, int Ny, int Nx, double complex **matrix) {
    FILE *file;

    add_ext(filename, "txt");
    file = fopen(full_filename, "w");

    for(int j = 0; j < Ny; j ++) {
        for(int i = 0; i < Nx; i ++) {
            fprintf(file, "%1.9e %1.9e\n", creal(matrix[j][i]), cimag(matrix[j][i]));
        }
        fprintf(file, "\n");
        fflush(file);
    }

    fclose(file);
}

void out_complex_matrix_bin(char *filename, int Ny, int Nx, double complex **matrix) {
    FILE *file;
    double tmp[2];

    add_ext(filename, "bin");
    file = fopen(full_filename, "wb");

    for(int j = 0; j < Ny; j ++) {
        for(int i = 0; i < Nx; i ++) {
            tmp[0] = creal(matrix[j][i]);
            tmp[1] = cimag(matrix[j][i]);
            fwrite(&tmp, sizeof(double), 2, file);
        }
        fflush(file);
    }

    fclose(file);
}

void out_complex_tensor_txt(char *filename, int Nz, int Ny, int Nx, double complex ***tensor) {
    FILE *file;

    add_ext(filename, "txt");
    file = fopen(full_filename, "w");

    for(int k = 0; k < Nz; k ++) {
        for(int j = 0; j < Ny; j ++) {
            for(int i = 0; i < Nx; i ++) {
                fprintf(file, "%1.9e %1.9e\n", creal(tensor[k][j][i]), cimag(tensor[k][j][i]));
            }
            fprintf(file, "\n");
            fflush(file);
        }
        fprintf(file, "\n");
        fflush(file);
    }

    fclose(file);
}

void out_complex_tensor_bin(char *filename, int Nz, int Ny, int Nx, double complex ***tensor) {
    FILE *file;
    double tmp[2];

    add_ext(filename, "bin");
    file = fopen(full_filename, "wb");

    for(int k = 0; k < Nz; k ++) {
        for(int j = 0; j < Ny; j ++) {
            for(int i = 0; i < Nx; i ++) {
                tmp[0] = creal(tensor[k][j][i]);
                tmp[1] = cimag(tensor[k][j][i]);
                fwrite(&tmp, sizeof(double), 2, file);
            }
            fflush(file);
        }
    }

    fclose(file);
}
