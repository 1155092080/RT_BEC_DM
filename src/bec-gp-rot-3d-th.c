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
#define _XOPEN_SOURCE 600

// Unfortunately, PGI complains when using collapse(2), so we ommit this clause when using PGI compilers
#ifdef __PGI
#define omp_collapse
#else
#define omp_collapse collapse(2)
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <omp.h>
#include <assert.h>
#include "util/cfg.h"
#include "util/mem.h"
#include "util/diffint.h"
#include "util/exp.h"
#include "util/bec_rand.h"
#include "util/out.h"
#include "util/timer.h"

#define MAX(a, b, c)       (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c)
#define BOHR_RADIUS  5.2917720859e-11
#define PI           3.14159265358979

static int Nx, Ny, Nz;
static int ADD_RANDOM_PHASE, ADD_ONE_VORTEX;

static double dx, dy, dz;
static double dt;

static double vnu, vgamma, vlambda;
static double g, g0;
static double aho, as;
static int Natoms;

static int Nstp, Npas, Nrun, iter;

static double *x, *y, *z, *x2, *y2, *z2;
static double complex ***psi;
static double ***psi2;
static double ***pot;
static double ***dpsir, ***dpsii;
static double **dpsixr, **dpsixi, **dpsiyr, **dpsiyi, **dpsizr, **dpsizi;
static double *rms;

static double complex Ax0r, Ay0r, Az, Az0r;
static double complex *Axm, *Axp, *Aym, *Ayp;
static double complex **alphax, **alphay, *alphaz;
static double complex **beta;
static double complex **gammax, **gammay, *gammaz;

static double complex **ctmpxy;
static double **tmpx, **tmpy, **tmpz;
static double **tmpxy, **tmpxz, **tmpyz;
static double ***tmpxyz;

static double omega;

static int optscale;
static int optreim;
static int seed;
static double complex C;

static char *input, *output, *rmsout, *initout, *Nstpout, *Npasout, *Nrunout, *finalout;
static enum OutputType { TEXTUAL = 1, BINARY = 2, VISUAL = 3 } outtype;
static enum OutputFlags { DEN_X = 1 << 0, DEN_Y = 1 << 1, DEN_Z = 1 << 2, DEN_XY = 1 << 3, DEN_XZ = 1 << 4, DEN_YZ = 1 << 5, DEN_XY0 = 1 << 6, DEN_X0Z = 1 << 7, DEN_0YZ = 1 << 8, DEN_XYZ = 1 << 9 } outflags;
static int outstpx, outstpy, outstpz, outstpt;

static int nthreads;

void alloc();
void freemem();
void readparam(char *filename);

void initspace(double *x, double *x2, double *y, double *y2, double *z, double *z2);
void initpsi(double complex ***psi);
void initpot(double ***pot);
void initcoef(double complex *Axm, double complex *Ax0r, double complex *Axp, double complex **alphax, double complex **gammax,
    double complex *Aym, double complex *Ay0r, double complex *Ayp, double complex **alphay, double complex **gammay,
    double complex *Az, double complex *Az0r, double complex *alphaz, double complex *gammaz);

double complex (*e)(double x);

void calcnu(double complex ***psi, double ***pot);
void calclux(double complex ***psi, double complex **beta);
void calcluy(double complex ***psi, double complex **beta);
void calcluz(double complex ***psi, double complex **beta);

void calcpsi2(double complex ***psi, double ***psi2);
void calcnorm(double *norm, double complex ***psi, double **tmpx, double **tmpy, double **tmpz);
void calcmuen(double *mu, double *en, double complex ***psi, double ***dpsir, double ***dpsii,
    double **dpsixr, double **dpsixi, double **dpsiyr, double **dpsiyi, double **dpsizr, double **dpsizi);
void calcrms(double *rms, double ***psi2, double **tmpx, double **tmpy, double **tmpz);

void writeout(char *prefix, int step);

void outdenx(double ***psi2, double **tmpy, double **tmpz, double *den, char *filename);
void outdeny(double ***psi2, double **tmpx, double **tmpz, double *den, char *filename);
void outdenz(double ***psi2, double **tmpx, double **tmpy, double *den, char *filename);

void outdenxy(double ***psi2, double **tmpz, double **den, char *filename);
void outdenxz(double ***psi2, double **tmpy, double **den, char *filename);
void outdenyz(double ***psi2, double **tmpx, double **den, char *filename);

void outpsi2xy(double ***psi2, double **outxy, char *filename);
void outpsi2xz(double ***psi2, double **outxz, char *filename);
void outpsi2yz(double ***psi2, double **outyz, char *filename);

void outpsi2xyz(double ***psi2, double ***outxyz, char *filename);

void outpsi(double complex ***psi, char *filename);

int main(int argc, char **argv) {
    FILE *fileout;
    FILE *filerms;
    char filename[FILENAME_MAX];
    double norm, mu, en;
    struct timespec progstart, progend, cpustart, cpuend;
    double walltime, cputime;

    if ((argc != 3) || (strcmp(*(argv + 1), "-i") != 0)) {
        fprintf(stderr, "Usage: %s -i <input parameter file> \n", *argv);
        exit(EXIT_FAILURE);
    }

    readparam(argv[2]);
    bec_srand(seed);

    get_wall_time(&progstart);
    get_cpu_time(&cpustart);

    #pragma omp parallel
    #pragma omp master
    nthreads = omp_get_num_threads();

    // Round up Nx and Ny to the nearest larger even number and allocate memory
    Nx = (Nx + 2) & ~1;
    Ny = (Ny + 2) & ~1;

    alloc();

    // After allocation, use the nearest odd number as the grid size
    Nx--;
    Ny--;

    if (optreim == 1) {
        C = 1. + I * 0.;
        e = &d_exp;
    } else {
        C = 0. + I * 1.;
        e = &c_exp;
    }

    if (Nstp > 0) {
        g = 0.;
    } else {
        g = optscale * g0;
    }

    if (output != NULL) {
        sprintf(filename, "%s.txt", output);
        fileout = fopen(filename, "w");
    } else fileout = stdout;

    if (rmsout != NULL) {
        sprintf(filename, "%s.txt", rmsout);
        filerms = fopen(filename, "w");
    } else filerms = NULL;

    if (optreim == 1) {
        fprintf(fileout, "Imaginary time propagation 3d, OPTION = %d, NUM_THREADS = %d\n\n", optscale, nthreads);
    } else {
        fprintf(fileout, "Real time propagation 3d, OPTION = %d, NUM_THREADS = %d\n\n", optscale, nthreads);
    }

    if (cfg_defined("G0")) {
        fprintf(fileout, "Nonlinearity G_3D = %.6f, G_2D = %.6f\n", g0, g0 * sqrt(vlambda / (2 * PI)));
    } else {
        fprintf(fileout, "Number of atoms N = %d, Unit of length AHO = %.8f m\n", Natoms, aho);
        fprintf(fileout, "Scattering length a = %.6f*a0\n", as);
        fprintf(fileout, "Calculated nonlinearity G_3D = %.6f, G_2D = %.6f\n", g0, g0 * sqrt(vlambda / (2 * PI)));
    }

    fprintf(fileout, "Parameters of trap: GAMMA = %.2f, NU = %.2f, LAMBDA = %.2f\n", vgamma, vnu, vlambda);
    fprintf(fileout, "Parameters of rotation: ANG VEL = %.4f\n\n", omega);
    fprintf(fileout, "Space steps: NX = %d, NY = %d, NZ = %d\n", Nx, Ny, Nz);
    fprintf(fileout, "             DX = %.6f, DY = %.6f, DZ = %.6f\n", dx, dy, dz);
    fprintf(fileout, "Time steps:  NSTP = %d, NPAS = %d, NRUN = %d\n", Nstp, Npas, Nrun);
    fprintf(fileout, "             DT = %.6f\n",  dt);
    if (Nstp != 0) fprintf(fileout, "RNG seed:    SEED = %d\n",  seed);
    fprintf(fileout, "\n             -----------------------------------------------------\n");
    fprintf(fileout, "               Iter     Norm       Chem       Ener/N      <rho>\n");
    fprintf(fileout, "             -----------------------------------------------------\n");
    fflush(fileout);

    initspace(x, x2, y, y2, z, z2);
    initpsi(psi);
    initpot(pot);
    initcoef(Axm, &Ax0r, Axp, alphax, gammax, Aym, &Ay0r, Ayp, alphay, gammay, &Az, &Az0r, alphaz, gammaz);

    calcnorm(&norm, psi, tmpx, tmpy, tmpz);
    calcmuen(&mu, &en, psi, dpsir, dpsii, dpsixr, dpsixi, dpsiyr, dpsiyi, dpsizr, dpsizi);
    calcpsi2(psi, psi2);
    calcrms(rms, psi2, tmpx, tmpy, tmpz);

    fprintf(fileout, "Initial:  %19.4f %11.5f %11.5f %10.5f\n", norm, mu / optscale, en / optscale, *rms);
    fflush(fileout);
    writeout(initout, 0);

    if (rmsout != NULL) {
        if (optreim == 1) {
            fprintf(filerms, "Imaginary time propagation 3d, OPTION = %d\n\n", optscale);
        } else {
            fprintf(filerms, "Real time propagation 3d, OPTION = %d\n\n", optscale);
        }
        fprintf(filerms, "           ------------------------------------------------------------\n");
        fprintf(filerms, "RMS size:    Iter      <r>          <x>          <y>          <z>\n");
        fprintf(filerms, "           ------------------------------------------------------------\n");
        fprintf(filerms, "Initial: %21.5f %12.5f %12.5f %12.5f\n", rms[0], rms[1], rms[2], rms[3]);
        fflush(filerms);
    }

    if (Nstp > 0) {
        double g_stp = optscale * g0 / (double) Nstp;
        g = 0.;
        for (int i = 1; i <= Nstp; i ++) {
            g += g_stp;
            calcnu(psi, pot);
            calclux(psi, beta);
            calcluy(psi, beta);
            calcluz(psi, beta);
            if (optreim == 1) calcnorm(&norm, psi, tmpx, tmpy, tmpz);
        }
        if (optreim == 2) calcnorm(&norm, psi, tmpx, tmpy, tmpz);
        calcmuen(&mu, &en, psi, dpsir, dpsii, dpsixr, dpsixi, dpsiyr, dpsiyi, dpsizr, dpsizi);
        calcpsi2(psi, psi2);
        calcrms(rms, psi2, tmpx, tmpy, tmpz);

        fprintf(fileout, "NSTP iter.:  %16.4f %11.5f %11.5f %10.5f\n", norm, mu / optscale, en / optscale, *rms);
        fflush(fileout);

        writeout(Nstpout, 0);
        if (rmsout != NULL) {
            fprintf(filerms, "NSTP iter.: %18.5f %12.5f %12.5f %12.5f\n", rms[0], rms[1], rms[2], rms[3]);
            fflush(filerms);
        }
    }

    int step = 1;
    for (int i = 1; i <= Npas; i ++) {
        calcnu(psi, pot);
        calclux(psi, beta);
        calcluy(psi, beta);
        calcluz(psi, beta);
        if (optreim == 1) calcnorm(&norm, psi, tmpx, tmpy, tmpz);
        if (i % iter == 0) {
            if (optreim == 2) calcnorm(&norm, psi, tmpx, tmpy, tmpz);
            calcmuen(&mu, &en, psi, dpsir, dpsii, dpsixr, dpsixi, dpsiyr, dpsiyi, dpsizr, dpsizi);
            calcpsi2(psi, psi2);
            calcrms(rms, psi2, tmpx, tmpy, tmpz);

            fprintf(fileout, "NPAS iter.:  %5d %10.4f %11.5f %11.5f %10.5f\n", step, norm, mu / optscale, en / optscale, *rms);
            fflush(fileout);
            writeout(Npasout, step);

            if (rmsout != NULL) {
                fprintf(filerms, "NPAS iter.: %5d %12.5f %12.5f %12.5f %12.5f\n", step, rms[0], rms[1], rms[2], rms[3]);
                fflush(filerms);
            }

            step++;
        }
    }

    for (int i = 1; i <= Nrun; i ++) {
        calcnu(psi, pot);
        calclux(psi, beta);
        calcluy(psi, beta);
        calcluz(psi, beta);
        if (optreim == 1) calcnorm(&norm, psi, tmpx, tmpy, tmpz);
    }
    if (optreim == 2) calcnorm(&norm, psi, tmpx, tmpy, tmpz);
    calcmuen(&mu, &en, psi, dpsir, dpsii, dpsixr, dpsixi, dpsiyr, dpsiyi, dpsizr, dpsizi);
    calcpsi2(psi, psi2);
    calcrms(rms, psi2, tmpx, tmpy, tmpz);

    fprintf(fileout, "NRUN iter.:  %16.4f %11.5f %11.5f %10.5f\n", norm, mu / optscale, en / optscale, *rms);
    fflush(fileout);
    writeout(Nrunout, 0);
    if (rmsout != NULL) {
        fprintf(filerms, "NRUN iter.: %18.5f %12.5f %12.5f %12.5f\n", rms[0], rms[1], rms[2], rms[3]);
        fprintf(filerms, "           ------------------------------------------------------------\n");
    }

    fprintf(fileout, "            -----------------------------------------------------\n\n");

    if (finalout != NULL) {
        outpsi(psi, finalout);
    }

    freemem();

    get_wall_time(&progend);
    get_cpu_time(&cpuend);
    walltime = elapsed_time(&progstart, &progend);
    cputime = elapsed_time(&cpustart, &cpuend);

    fprintf(fileout, " Clock Time: %.f seconds\n", walltime);
    fprintf(fileout, " CPU Time: %.f seconds\n", cputime);

    if (output != NULL) fclose(fileout);
    if (rmsout != NULL) fclose(filerms);

    return EXIT_SUCCESS;
}

/**
 *    Allocation of memory.
 */
void alloc() {
    x = alloc_double_vector(Nx);
    y = alloc_double_vector(Ny);
    z = alloc_double_vector(Nz);
    x2 = alloc_double_vector(Nx);
    y2 = alloc_double_vector(Ny);
    z2 = alloc_double_vector(Nz);

    pot = alloc_double_tensor(Nz, Ny, Nx);

    Axm = alloc_complex_vector(Ny);
    Axp = alloc_complex_vector(Ny);
    alphax = alloc_complex_matrix(Ny, Nx);
    gammax = alloc_complex_matrix(Ny, Nx);

    Aym = alloc_complex_vector(Nx);
    Ayp = alloc_complex_vector(Nx);
    alphay = alloc_complex_matrix(Ny, Nx);
    gammay = alloc_complex_matrix(Ny, Nx);

    alphaz = alloc_complex_vector(Nz);
    gammaz = alloc_complex_vector(Nz);

    beta = alloc_complex_matrix(nthreads, MAX(Nx, Ny, Nz) - 1);

    psi = alloc_complex_tensor(Nz, Ny, Nx);
    psi2 = alloc_double_tensor(Nz, Ny, Nx);

    dpsir = alloc_double_tensor(Nz, Ny, Nx);
    dpsii = alloc_double_tensor(Nz, Ny, Nx);
    dpsixr = alloc_double_matrix(nthreads, Nx);
    dpsixi = alloc_double_matrix(nthreads, Nx);
    dpsiyr = alloc_double_matrix(nthreads, Ny);
    dpsiyi = alloc_double_matrix(nthreads, Ny);
    dpsizr = alloc_double_matrix(nthreads, Nz);
    dpsizi = alloc_double_matrix(nthreads, Nz);

    rms = alloc_double_vector(4);

    tmpx = dpsixr;
    tmpy = dpsiyr;
    tmpz = dpsizr;

    ctmpxy = alloc_complex_matrix(Ny, Nx);
    tmpxy = alloc_double_matrix(Ny, Nx);
    tmpxz = alloc_double_matrix(Nz, Nx);
    tmpyz = alloc_double_matrix(Nz, Ny);

    tmpxyz = dpsir;

    memset(**psi, 0, sizeof(double complex) * Nz * Ny * Nx);
}

/**
 *    Free all dynamically allocated memory.
 */
void freemem() {
    free_double_vector(x);
    free_double_vector(y);
    free_double_vector(z);
    free_double_vector(x2);
    free_double_vector(y2);
    free_double_vector(z2);

    free_double_tensor(pot);

    free_complex_vector(Axm);
    free_complex_vector(Axp);
    free_complex_matrix(alphax);
    free_complex_matrix(gammax);

    free_complex_vector(Aym);
    free_complex_vector(Ayp);
    free_complex_matrix(alphay);
    free_complex_matrix(gammay);

    free_complex_vector(alphaz);
    free_complex_vector(gammaz);

    free_complex_matrix(beta);

    free_complex_tensor(psi);
    free_double_tensor(psi2);

    free_double_tensor(dpsir);
    free_double_tensor(dpsii);
    free_double_matrix(dpsixr);
    free_double_matrix(dpsixi);
    free_double_matrix(dpsiyr);
    free_double_matrix(dpsiyi);
    free_double_matrix(dpsizr);
    free_double_matrix(dpsizi);

    free_double_vector(rms);

    free_complex_matrix(ctmpxy);
    free_double_matrix(tmpxy);
    free_double_matrix(tmpxz);
    free_double_matrix(tmpyz);
}

/**
 *    Reading input parameters from the configuration file.
 *    filename - name of file with input parameters
 */
void readparam(char *filename) {
    if (! cfg_init(filename)) {
        fprintf(stderr, "Wrong input parameter file.\n");
        exit(EXIT_FAILURE);
    }

    if (cfg_defined("G0")) {
        g0 = cfg_read_double("G0");
    } else {
        Natoms = cfg_read_int("NATOMS");
        aho = cfg_read_double("AHO");
        as = cfg_read_double("AS");

        g0 = 4. * PI * as * Natoms * BOHR_RADIUS / aho;
    }

    omega = cfg_read_double("OMEGA");
    seed = (cfg_defined("SEED")) ? cfg_read_int("SEED") : 1;

    if (cfg_defined("ADD_RANDOM_PHASE")) {
        ADD_RANDOM_PHASE = cfg_read_int("ADD_RANDOM_PHASE");
    } else {
        ADD_RANDOM_PHASE = 1;
    }

    if (cfg_defined("ADD_ONE_VORTEX")) {
        ADD_ONE_VORTEX = cfg_read_int("ADD_ONE_VORTEX");
    } else {
        ADD_ONE_VORTEX = 1;
    }

    Nx = cfg_read_long("NX");
    Ny = cfg_read_long("NY");
    Nz = cfg_read_long("NZ");

    Nstp = cfg_read_int("NSTP");
    Npas = cfg_read_int("NPAS");
    Nrun = cfg_read_int("NRUN");
    iter = (cfg_defined("ITER")) ? Npas / cfg_read_int("ITER") : Npas;

    dx = cfg_read_double("DX");
    dy = cfg_read_double("DY");
    dz = cfg_read_double("DZ");
    dt = cfg_read_double("DT");

    vgamma = cfg_read_double("GAMMA");
    vnu = cfg_read_double("NU");
    vlambda = cfg_read_double("LAMBDA");

    optscale = cfg_read_int("OPTSCALE");
    optreim = cfg_read_int("OPTREIM");

    if (Nstp == 0) input = cfg_read_string("INPUT");
    if (cfg_defined("OUTPUT")) output = cfg_read_string("OUTPUT");
    if (cfg_defined("RMSOUT")) rmsout = cfg_read_string("RMSOUT");
    if (cfg_defined("INITOUT")) initout = cfg_read_string("INITOUT");
    if (cfg_defined("NSTPOUT")) Nstpout = cfg_read_string("NSTPOUT");
    if (cfg_defined("NPASOUT")) Npasout = cfg_read_string("NPASOUT");
    if (cfg_defined("NRUNOUT")) Nrunout = cfg_read_string("NRUNOUT");
    if (cfg_defined("FINALOUT")) finalout = cfg_read_string("FINALOUT");

    if (initout != NULL || Nstpout != NULL || Npasout != NULL || Nrunout != NULL) {
        outtype = cfg_read_int("OUTTYPE");
        outflags = cfg_read_int("OUTFLAGS");

        outstpx = cfg_read_int("OUTSTPX");
        outstpy = cfg_read_int("OUTSTPY");
        outstpz = cfg_read_int("OUTSTPZ");
    }
}

/**
 *    Initialization of the space mesh.
 */
void initspace(double *x, double *x2, double *y, double *y2, double *z, double *z2) {
    for (int i = 0; i < Nx; i ++) {
        x[i] = (i - Nx / 2) * dx;
        x2[i] = x[i] * x[i];
    }

    for (int j = 0; j < Ny; j ++) {
        y[j] = (j - Ny / 2) * dy;
        y2[j] = y[j] * y[j];
    }

    for (int k = 0; k < Nz; k ++) {
        z[k] = (k - Nz / 2) * dz;
        z2[k] = z[k] * z[k];
    }
}

/**
 *    Initialization of the wave function.
 *    psi - array with the wave function values
 */
void initpsi(double complex ***psi) {
    double cpsi;
    double tmpr, tmpi, ctmp[2];
    FILE *file;

    cpsi = sqrt(PI * sqrt(PI / (vgamma * vnu * vlambda)));

    if (Nstp == 0) {
        char *dot = strrchr(input, '.');
        bool textual;

        if (dot && !strcmp(dot, ".txt")) {
            file = fopen(input, "r");
            textual = true;
        } else {
            file = fopen(input, "rb");
            textual = false;
        }

        if (file == NULL) {
            fprintf(stderr, "Could not open the requested file.\n");
            exit(EXIT_FAILURE);
        }

        for (int k = 0; k < Nz; k ++) {
            for (int j = 0; j < Ny; j ++) {
                for (int i = 0; i < Nx; i ++) {
                    if (textual) {
                        if(fscanf(file,"%le %le\n", &tmpr, &tmpi) == 2) {
                            psi[k][j][i] = tmpr + I * tmpi;
                        }
                    } else {
                        if (fread(ctmp, sizeof(double), 2, file) == 2) {
                            psi[k][j][i] = ctmp[0] + I * ctmp[1];
                        }
                    }
                }
            }
        }
        fclose(file);
    } else {
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                ctmpxy[j][i] = 1.;
                if(ADD_RANDOM_PHASE == 1) ctmpxy[j][i] *= cexp(2 * PI * I * bec_rand());
                if(ADD_ONE_VORTEX == 1) ctmpxy[j][i] *= (x[i] + I * y[j]);
            }
        }

        for (int k = 0; k < Nz; k ++) {
            for (int j = 0; j < Ny; j ++) {
                for (int i = 0; i < Nx; i ++) {
                    psi[k][j][i] = ctmpxy[j][i] * exp(- 0.5 * (vgamma * x2[i] + vnu * y2[j] + vlambda * z2[k])) / cpsi;
                }
            }
        }
    }
}

/**
 *    Initialization of the potential.
 *    pot - array with the trap potential
 */
void initpot(double ***pot) {
    double vgamma2, vnu2, vlambda2;

    vgamma2 = vgamma * vgamma;
    vnu2 = vnu * vnu;
    vlambda2 = vlambda * vlambda;

    #pragma omp parallel for omp_collapse
    for (int k = 0; k < Nz; k ++) {
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                pot[k][j][i] = optscale * 0.5 * (vgamma2 * x2[i] + vnu2 * y2[j] + vlambda2 * z2[k]);
            }
        }
    }
}

/**
 *    Initialization of Crank-Nicolson scheme coefficients.
 */
void initcoef(double complex *Axm, double complex *Ax0r, double complex *Axp, double complex **alphax, double complex **gammax,
    double complex *Aym, double complex *Ay0r, double complex *Ayp, double complex **alphay, double complex **gammay,
    double complex *Az, double complex *Az0r, double complex *alphaz, double complex *gammaz) {

    double dx2, dy2, dz2;
    double complex cdt;
    double complex A0, Ax, Ay, Ax0, Ay0, Az0, Atmpx, Atmpy;

    dx2 = dx * dx;
    dy2 = dy * dy;
    dz2 = dz * dz;

    cdt = C * dt;

    Ax0 = 1. + cdt / dx2;
    *Ax0r = 1. - cdt / dx2;
    Ax = 0.5 * cdt / dx2;
    Ay0 = 1. + cdt / dy2;
    *Ay0r = 1. - cdt / dy2;
    Ay = 0.5 * cdt / dy2;
    Az0 = 1. + cdt / dz2;
    *Az0r = 1. - cdt / dz2;
    *Az = - 0.5 * cdt / dz2;

    A0 = 0.5 * I * cdt * optscale * omega;

    Atmpx = A0 / (2. * dx);
    Atmpy = A0 / (2. * dy);

    for (int j = 0; j < Ny; j ++) {
        Axm[j] = - (Ax - Atmpx * y[j]);
        Axp[j] = - (Ax + Atmpx * y[j]);
    }

    for (int j = 0; j < Ny; j ++) {
        alphax[j][Nx - 2] = 0.;
        gammax[j][Nx - 2] = - 1. / Ax0;
        for (int i = Nx - 2; i > 0; i --) {
            alphax[j][i - 1] = Axm[j] * gammax[j][i];
            gammax[j][i - 1] = - 1. / (Ax0 + Axp[j] * alphax[j][i - 1]);
        }
    }

    for (int i = 0; i < Nx; i ++) {
        Aym[i] = - (Ay + Atmpy * x[i]);
        Ayp[i] = - (Ay - Atmpy * x[i]);
    }

    for (int i = 0; i < Nx; i ++) {
        alphay[Ny - 2][i] = 0.;
        gammay[Ny - 2][i] = - 1. / Ay0;
        for (int j = Ny - 2; j > 0; j --) {
            alphay[j - 1][i] = Aym[i] * gammay[j][i];
            gammay[j - 1][i] = - 1. / (Ay0 + Ayp[i] * alphay[j - 1][i]);
        }
    }

    alphaz[Nz - 2] = 0.;
    gammaz[Nz - 2] = - 1. / Az0;
    for (int k = Nz - 2; k > 0; k --) {
        alphaz[k - 1] = *Az * gammaz[k];
        gammaz[k - 1] = - 1. / (Az0 + *Az * alphaz[k - 1]);
    }
}

/**
 *    Time propagation with respect to H1 (part of the Hamiltonian without spatial derivatives).
 *    psi - array with the wave function values
 *    pot - array with the trap potential
 */
void calcnu(double complex ***psi, double ***pot) {
    double psi2, tmp;

    #pragma omp parallel for private(psi2, tmp) omp_collapse
    for (int k = 0; k < Nz; k ++) {
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                psi2 = psi[k][j][i] * conj(psi[k][j][i]);
                tmp = dt * (pot[k][j][i] + g * psi2);
                psi[k][j][i] *= (*e)(- tmp);
            }
        }
    }
}

/**
 *    Time propagation with respect to H2 (x-part of the Laplacian).
 *    psi  - array with the wave function values
 *    beta - Crank-Nicolson scheme coefficients
 */
void calclux(double complex ***psi, double complex **beta) {
    double complex c;

    #pragma omp parallel private(c)
    {
        int threadid = omp_get_thread_num();

        #pragma omp for omp_collapse
        for (int k = 0; k < Nz; k ++) {
            for (int j = 0; j < Ny; j ++) {
                beta[threadid][Nx - 2] = (optreim == 1) ? 0. : psi[k][j][Nx - 1];
                for (int i = Nx - 2; i > 0; i --) {
                    c = - Axp[j] * psi[k][j][i + 1] + Ax0r * psi[k][j][i] - Axm[j] * psi[k][j][i - 1];
                    beta[threadid][i - 1] = gammax[j][i] * (Axp[j] * beta[threadid][i] - c);
                }
                psi[k][j][0] = 0.;
                for (int i = 0; i < Nx - 2; i ++) {
                    psi[k][j][i + 1] = alphax[j][i] * psi[k][j][i] + beta[threadid][i];
                }
                psi[k][j][Nx - 1] = 0.;
            }
        }
    }
}

/**
 *    Time propagation with respect to H3 (y-part of the Laplacian).
 *    psi  - array with the wave function values
 *    beta - Crank-Nicolson scheme coefficients
 */
void calcluy(double complex ***psi, double complex **beta) {
    double complex c;

    #pragma omp parallel private(c)
    {
        int threadid = omp_get_thread_num();

        #pragma omp for omp_collapse
        for (int k = 0; k < Nz; k ++) {
            for (int i = 0; i < Nx; i ++) {
                beta[threadid][Ny - 2] = (optreim == 1) ? 0. : psi[k][Ny - 1][i];
                for (int j = Ny - 2; j > 0; j --) {
                    c = - Ayp[i] * psi[k][j + 1][i] + Ay0r * psi[k][j][i] - Aym[i] * psi[k][j - 1][i];
                    beta[threadid][j - 1] = gammay[j][i] * (Ayp[i] * beta[threadid][j] - c);
                }
                psi[k][0][i] = 0.;
                for (int j = 0; j < Ny - 2; j ++) {
                    psi[k][j + 1][i] = alphay[j][i] * psi[k][j][i] + beta[threadid][j];
                }
                psi[k][Ny - 1][i] = 0.;
            }
        }
    }
}

/**
 *    Time propagation with respect to H4 (z-part of the Laplacian).
 *    psi  - array with the wave function values
 *    beta - Crank-Nicolson scheme coefficients
 */
void calcluz(double complex ***psi, double complex **beta) {
    double complex c;

    #pragma omp parallel private(c)
    {
        int threadid = omp_get_thread_num();

        #pragma omp for omp_collapse
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                beta[threadid][Nz - 2] = (optreim == 1) ? 0. : psi[Nz - 1][j][i];
                for (int k = Nz - 2; k > 0; k --) {
                    c = - Az * psi[k + 1][j][i] + Az0r * psi[k][j][i] - Az * psi[k - 1][j][i];
                    beta[threadid][k - 1] = gammaz[k] * (Az * beta[threadid][k] - c);
                }
                psi[0][j][i] = 0.;
                for (int k = 0; k < Nz - 2; k ++) {
                    psi[k + 1][j][i] = alphaz[k] * psi[k][j][i] + beta[threadid][k];
                }
                psi[Nz - 1][j][i] = 0.;
            }
        }
    }
}

/**
 *    Calculation of squared wave function values.
 *    psi  - array with the wave function values
 *    psi2 - array with the squared wave function values
 */
void calcpsi2(double complex ***psi, double ***psi2) {
    #pragma omp parallel for omp_collapse
    for (int k = 0; k < Nz; k ++) {
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                psi2[k][j][i] = psi[k][j][i] * conj(psi[k][j][i]);
            }
        }
    }
}

/**
 *    Calculation of the wave function norm and normalization.
 *    norm - wave function norm
 *    psi  - array with the wave function values
 *    tmpx - temporary array
 *    tmpy - temporary array
 *    tmpz - temporary array
 */
void calcnorm(double *norm, double complex ***psi, double **tmpx, double **tmpy, double **tmpz) {
    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int k = 0; k < Nz; k ++) {
            for (int j = 0; j < Ny; j ++) {
                for (int i = 0; i < Nx; i ++) {
                    tmpx[threadid][i] = psi[k][j][i] * conj(psi[k][j][i]);
                }
                tmpy[threadid][j] = simpint(dx, tmpx[threadid], Nx);
            }
            (*tmpz)[k] = simpint(dy, tmpy[threadid], Ny);
        }
    }

    *norm = sqrt(simpint(dz, *tmpz, Nz));

    if (optreim == 1) {
        double tmp = 1. / *norm;

        #pragma omp parallel for omp_collapse
        for (int k = 0; k < Nz; k ++) {
            for (int j = 0; j < Ny; j ++) {
                for (int i = 0; i < Nx; i ++) {
                    psi[k][j][i] *= tmp;
                }
            }
        }
    }
}

/**
 *    Calculation of the chemical potential and energy.
 *    mu     - chemical potential
 *    en     - energy
 *    psi    - array with the wave function values
 *    dpsir  - temporary array
 *    dpsii  - temporary array
 *    dpsixr - temporary array
 *    dpsixi - temporary array
 *    dpsiyr - temporary array
 *    dpsiyi - temporary array
 *    dpsizr - temporary array
 *    dpsizi - temporary array
 */
void calcmuen(double *mu, double *en, double complex ***psi, double ***dpsir, double ***dpsii, double **dpsixr, double **dpsixi, double **dpsiyr, double **dpsiyi, double **dpsizr, double **dpsizi) {
    int threadid;
    double psi2, psi2r, dpsiz;
    double norm;

    #pragma omp parallel private(threadid, psi2, psi2r, dpsiz)
    {
        threadid = omp_get_thread_num();

        #pragma omp for omp_collapse
        for (int k = 0; k < Nz; k ++) {
            for (int i = 0; i < Nx; i ++) {
                cdiff(dy, psi[k][0] + i, dpsiyr[threadid], dpsiyi[threadid], Ny, Nx + 1); // Nx + 1 is due to padded allocation
                for (int j = 0; j < Ny; j ++) {
                    dpsir[k][j][i] = dpsiyr[threadid][j] * dpsiyr[threadid][j];
                    dpsii[k][j][i] = x[i] * dpsiyi[threadid][j];
                }
            }
        }

        #pragma omp for omp_collapse
        for (int k = 0; k < Nz; k ++) {
            for (int j = 0; j < Ny; j ++) {
                cdiff(dx, psi[k][j], dpsixr[threadid], dpsixi[threadid], Nx, 1);
                for (int i = 0; i < Nx; i ++) {
                    dpsir[k][j][i] += dpsixr[threadid][i] * dpsixr[threadid][i];
                    dpsii[k][j][i] -= y[j] * dpsixi[threadid][i];
                }
            }
        }

        #pragma omp for omp_collapse
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                for (int k = 0; k < Nz; k ++) {
                    dpsizi[threadid][k] = creal(psi[k][j][i]);
                }
                diff(dz, dpsizi[threadid], dpsizr[threadid], Nz);
                for (int k = 0; k < Nz; k ++) {
                    dpsir[k][j][i] += dpsizr[threadid][k] * dpsizr[threadid][k];
                }
            }
        }

        #pragma omp for
        for (int k = 0; k < Nz; k ++) {
            for (int j = 0; j < Ny; j ++) {
                for (int i = 0; i < Nx; i ++) {
                    psi2 = psi[k][j][i] * conj(psi[k][j][i]);
                    psi2r = creal(psi[k][j][i]);
                    psi2r *= psi2r;
                    dpsiz = creal(psi[k][j][i]) * dpsii[k][j][i];
                    dpsixr[threadid][i] = (pot[k][j][i] + psi2 * g) * psi2r + dpsir[k][j][i] - optscale * (omega * dpsiz);
                    dpsixi[threadid][i] = (pot[k][j][i] + 0.5 * psi2 * g) * psi2r + dpsir[k][j][i] - optscale * (omega * dpsiz);
                }
                dpsiyr[threadid][j] = simpint(dx, dpsixr[threadid], Nx);
                dpsiyi[threadid][j] = simpint(dx, dpsixi[threadid], Nx);
            }
            (*dpsizr)[k] = simpint(dy, dpsiyr[threadid], Ny);
            (*dpsizi)[k] = simpint(dy, dpsiyi[threadid], Ny);
        }
    }

    *mu = simpint(dz, *dpsizr, Nz);
    *en = simpint(dz, *dpsizi, Nz);

    #pragma omp parallel private (threadid)
    {
        threadid = omp_get_thread_num();

        #pragma omp for
        for (int k = 0; k < Nz; k ++) {
            for (int j = 0; j < Ny; j ++) {
                for (int i = 0; i < Nx; i ++) {
                    dpsixr[threadid][i] = creal(psi[k][j][i]) * creal(psi[k][j][i]);
                }
                dpsiyr[threadid][j] = simpint(dx, dpsixr[threadid], Nx);
            }
            (*dpsizr)[k] = simpint(dy, dpsiyr[threadid], Ny);
        }
    }

    norm = simpint(dz, *dpsizr, Nz);

    *mu /= norm;
    *en /= norm;
}

/**
 *    Calculation of the root mean square radius.
 *    rms  - root mean square radius
 *    psi2 - array with the squared wave function values
 *    tmpx - temporary array
 *    tmpy - temporary array
 *    tmpz - temporary array
 */
void calcrms(double *rms, double ***psi2, double **tmpx, double **tmpy, double **tmpz) {
    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int k = 0; k < Nz; k ++) {
            for (int j = 0; j < Ny; j ++) {
                for (int i = 0; i < Nx; i ++) {
                    tmpx[threadid][i] = x2[i] * psi2[k][j][i];
                }
                tmpy[threadid][j] = simpint(dx, tmpx[threadid], Nx);
            }
            (*tmpz)[k] = simpint(dy, tmpy[threadid], Ny);
        }

        #pragma omp single
        rms[1] = sqrt(simpint(dz, *tmpz, Nz));

        #pragma omp for
        for (int k = 0; k < Nz; k ++) {
            for (int i = 0; i < Nx; i ++) {
                for (int j = 0; j < Ny; j ++) {
                    tmpy[threadid][j] = y2[j] * psi2[k][j][i];
                }
                tmpx[threadid][i] = simpint(dy, tmpy[threadid], Ny);
            }
            (*tmpz)[k] = simpint(dx, tmpx[threadid], Nx);
        }

        #pragma omp single
        rms[2] = sqrt(simpint(dz, *tmpz, Nz));

        #pragma omp for
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                for (int k = 0; k < Nz; k ++) {
                    tmpz[threadid][k] = z2[k] * psi2[k][j][i];
                }
                tmpx[threadid][i] = simpint(dz, tmpz[threadid], Nz);
            }
            (*tmpy)[j] = simpint(dx, tmpx[threadid], Nx);
        }

        #pragma omp single
        rms[3] = sqrt(simpint(dy, *tmpy, Ny));
    }

    rms[0] = sqrt(rms[1] * rms[1] + rms[2] * rms[2] + rms[3] * rms[3]);
}

/**
 *    Calculate and write density profiles.
 *    prefix - Prefix for all file names used
 *    step   - Optional step used to construct file name
 */
void writeout(char *prefix, int step) {
    char filename[FILENAME_MAX];
    char iter[10];

    if (prefix != NULL) {
        calcpsi2(psi, psi2);

        if (step > 0) {
            sprintf(iter, "-%d", step);
        } else {
            iter[0] = '\0';
        }

        if (outflags & DEN_X) {
            sprintf(filename, "%s%s-1d_x", prefix, iter);
            outdenx(psi2, tmpy, tmpz, *tmpxy, filename);
        }

        if (outflags & DEN_Y) {
            sprintf(filename, "%s%s-1d_y", prefix, iter);
            outdeny(psi2, tmpx, tmpz, *tmpxy, filename);
        }

        if (outflags & DEN_Z) {
            sprintf(filename, "%s%s-1d_z", prefix, iter);
            outdenz(psi2, tmpx, tmpy, *tmpxy, filename);
        }

        if (outflags & DEN_XY) {
            sprintf(filename, "%s%s-2d_xy", prefix, iter);
            outdenxy(psi2, tmpz, tmpxy, filename);
        }

        if (outflags & DEN_XZ) {
            sprintf(filename, "%s%s-2d_xz", prefix, iter);
            outdenxz(psi2, tmpy, tmpxz, filename);
        }

        if (outflags & DEN_YZ) {
            sprintf(filename, "%s%s-2d_yz", prefix, iter);
            outdenyz(psi2, tmpx, tmpyz, filename);
        }

        if (outflags & DEN_XY0) {
            sprintf(filename, "%s%s-2d_xy0", prefix, iter);
            outpsi2xy(psi2, tmpxy, filename);
        }

        if (outflags & DEN_X0Z) {
            sprintf(filename, "%s%s-2d_x0z", prefix, iter);
            outpsi2xz(psi2, tmpxz, filename);
        }

        if (outflags & DEN_0YZ) {
            sprintf(filename, "%s%s-2d_0yz", prefix, iter);
            outpsi2yz(psi2, tmpyz, filename);
        }

        if (outflags & DEN_XYZ) {
            sprintf(filename, "%s%s-3d", prefix, iter);
            outpsi2xyz(psi2, tmpxyz, filename);
        }
    }
}

void outdenx(double ***psi2, double **tmpy, double **tmpz, double *den, char *filename) {
    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int i = 0; i < Nx; i += outstpx) {
            for (int k = 0; k < Nz; k ++) {
                for (int j = 0; j < Ny; j ++) {
                    tmpy[threadid][j] = psi2[k][j][i];
                }
                tmpz[threadid][k] = simpint(dy, tmpy[threadid], Ny);
            }
            den[i] = simpint(dz, tmpz[threadid], Nz);
        }
    }

    switch(outtype) {
        case TEXTUAL: out_curve_txt(filename, Nx, outstpx, x, den); break;
        case BINARY:  out_curve_bin(filename, Nx, outstpx, x, den); break;
        case VISUAL:  out_curve_vis(filename, Nx, outstpx, x, den, "Integrated density X"); break;
    }
}

void outdeny(double ***psi2, double **tmpx, double **tmpz, double *den, char *filename) {
    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int j = 0; j < Ny; j += outstpy) {
            for (int k = 0; k < Nz; k ++) {
                for (int i = 0; i < Nx; i ++) {
                    tmpx[threadid][i] = psi2[k][j][i];
                }
                tmpz[threadid][k] = simpint(dx, tmpx[threadid], Nx);
            }
            den[j] = simpint(dz, tmpz[threadid], Nz);
        }
    }

    switch(outtype) {
        case TEXTUAL: out_curve_txt(filename, Ny, outstpy, y, den); break;
        case BINARY:  out_curve_bin(filename, Ny, outstpy, y, den); break;
        case VISUAL:  out_curve_vis(filename, Ny, outstpy, y, den, "Integrated density Y"); break;
    }
}

void outdenz(double ***psi2, double **tmpx, double **tmpy, double *den, char *filename) {
    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int k = 0; k < Nz; k += outstpz) {
            for (int j = 0; j < Ny; j ++) {
                for (int i = 0; i < Nx; i ++) {
                    tmpx[threadid][i] = psi2[k][j][i];
                }
                tmpy[threadid][j] = simpint(dx, tmpx[threadid], Nx);
            }
            den[k] = simpint(dy, tmpy[threadid], Ny);
        }
    }

    switch(outtype) {
        case TEXTUAL: out_curve_txt(filename, Nz, outstpz, z, den); break;
        case BINARY:  out_curve_bin(filename, Nz, outstpz, z, den); break;
        case VISUAL:  out_curve_vis(filename, Nz, outstpz, z, den, "Integrated density Z"); break;
    }
}

void outdenxy(double ***psi2, double **tmpz, double **den, char *filename) {
    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int j = 0; j < Ny; j += outstpy) {
            for (int i = 0; i < Nx; i += outstpx) {
                for (int k = 0; k < Nz; k ++) {
                    tmpz[threadid][k] = psi2[k][j][i];
                }
                den[j][i] = simpint(dz, tmpz[threadid], Nz);
            }
        }
    }

    switch(outtype) {
        case TEXTUAL: out_matrix_txt(filename, Ny, Nx, outstpy, outstpx, y, x, den); break;
        case BINARY:  out_matrix_bin(filename, Ny, Nx, outstpy, outstpx, y, x, den); break;
        case VISUAL:  out_matrix_vis(filename, Ny, Nx, outstpy, outstpx, y, x, den, "Integrated_density_XY"); break;
    }
}

void outdenxz(double ***psi2, double **tmpy, double **den, char *filename) {
    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int k = 0; k < Nz; k += outstpz) {
            for (int i = 0; i < Nx; i += outstpx) {
                for (int j = 0; j < Ny; j ++) {
                    tmpy[threadid][j] = psi2[k][j][i];
                }
                den[k][i] = simpint(dy, tmpy[threadid], Ny);
            }
        }
    }

    switch(outtype) {
        case TEXTUAL: out_matrix_txt(filename, Nz, Nx, outstpz, outstpx, z, x, den); break;
        case BINARY:  out_matrix_bin(filename, Nz, Nx, outstpz, outstpx, z, x, den); break;
        case VISUAL:  out_matrix_vis(filename, Nz, Nx, outstpz, outstpx, z, x, den, "Integrated_density_XZ"); break;
    }
}

void outdenyz(double ***psi2, double **tmpx, double **den, char *filename) {
    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int k = 0; k < Nz; k += outstpz) {
            for (int j = 0; j < Ny; j += outstpy) {
                for (int i = 0; i < Nx; i ++) {
                    tmpx[threadid][i] = psi2[k][j][i];
                }
                den[k][j] = simpint(dx, tmpx[threadid], Nx);
            }
        }
    }

    switch(outtype) {
        case TEXTUAL: out_matrix_txt(filename, Nz, Ny, outstpz, outstpy, z, y, den); break;
        case BINARY:  out_matrix_bin(filename, Nz, Ny, outstpz, outstpy, z, y, den); break;
        case VISUAL:  out_matrix_vis(filename, Nz, Ny, outstpz, outstpy, z, y, den, "Integrated_density_YZ"); break;
    }
}

void outpsi2xy(double ***psi2, double **outxy, char *filename) {
    #pragma omp parallel for
    for (int j = 0; j < Ny; j += outstpy) {
        for (int i = 0; i < Nx; i += outstpx) {
            outxy[j][i] = psi2[Nz / 2][j][i];
        }
    }

    switch(outtype) {
        case TEXTUAL: out_matrix_txt(filename, Ny, Nx, outstpy, outstpx, y, x, outxy); break;
        case BINARY:  out_matrix_bin(filename, Ny, Nx, outstpy, outstpx, y, x, outxy); break;
        case VISUAL:  out_matrix_vis(filename, Ny, Nx, outstpy, outstpx, y, x, outxy, "Density_XY0"); break;
    }
}

void outpsi2xz(double ***psi2, double **outxz, char *filename) {
    #pragma omp parallel for
    for (int k = 0; k < Nz; k += outstpz) {
        for (int i = 0; i < Nx; i += outstpx) {
            outxz[k][i] = psi2[k][Ny / 2][i];
        }
    }

    switch(outtype) {
        case TEXTUAL: out_matrix_txt(filename, Nz, Nx, outstpz, outstpx, z, x, outxz); break;
        case BINARY:  out_matrix_bin(filename, Nz, Nx, outstpz, outstpx, z, x, outxz); break;
        case VISUAL:  out_matrix_vis(filename, Nz, Nx, outstpz, outstpx, z, x, outxz, "Density_X0Z"); break;
    }
}

void outpsi2yz(double ***psi2, double **outyz, char *filename) {
    #pragma omp parallel for
    for (int k = 0; k < Nz; k += outstpz) {
        for (int j = 0; j < Ny; j += outstpy) {
            outyz[k][j] = psi2[k][j][Nx / 2];
        }
    }

    switch(outtype) {
        case TEXTUAL: out_matrix_txt(filename, Nz, Ny, outstpz, outstpy, z, y, outyz); break;
        case BINARY:  out_matrix_bin(filename, Nz, Ny, outstpz, outstpy, z, y, outyz); break;
        case VISUAL:  out_matrix_vis(filename, Nz, Ny, outstpz, outstpy, z, y, outyz, "Density_0YZ"); break;
    }
}

void outpsi2xyz(double ***psi2, double ***outxyz, char *filename) {
    #pragma omp parallel for
    for (int k = 0; k < Nz; k += outstpz) {
        for (int j = 0; j < Ny; j += outstpy) {
            for (int i = 0; i < Nx; i += outstpx) {
                outxyz[k][j][i] = psi2[k][j][i];
            }
        }
    }

    switch(outtype) {
        case TEXTUAL: out_tensor_txt(filename, Nz, Ny, Nx, outstpz, outstpy, outstpx, z, y, x, outxyz); break;
        case BINARY:  out_tensor_bin(filename, Nz, Ny, Nx, outstpz, outstpy, outstpx, z, y, x, outxyz); break;
        case VISUAL:  out_tensor_vis(filename, Nz, Ny, Nx, outstpz, outstpy, outstpx, z, y, x, outxyz, "Density"); break;
    }
}

void outpsi(double complex ***psi, char *filename) {
    switch(outtype) {
        case TEXTUAL: out_complex_tensor_txt(filename, Nz, Ny, Nx, psi); break;
        case BINARY:  out_complex_tensor_bin(filename, Nz, Ny, Nx, psi); break;
        case VISUAL:  out_complex_tensor_bin(filename, Nz, Ny, Nx, psi); break;
    }
}
