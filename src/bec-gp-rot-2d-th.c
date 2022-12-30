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

#define MAX(a, b)    (a > b) ? a : b
#define BOHR_RADIUS  5.2917720859e-11
#define PI           3.14159265358979

static int Nx, Ny;
static int ADD_RANDOM_PHASE, ADD_ONE_VORTEX;

static double dx, dy;
static double dt;

static double vnu, vgamma;
static double g, g0, g_2d, d_z;
static double aho, as;
static int Natoms;

static int Nstp, Npas, Nrun, iter;

static double *x, *y, *x2, *y2;
static double complex **psi;
static double **psi2;
static double **pot;
static double **dpsir, **dpsii;
static double **dpsixr, **dpsixi, **dpsiyr, **dpsiyi;
static double *rms;

static double complex Ax0r, Ay0r;
static double complex *Axm, *Axp, *Aym, *Ayp;
static double complex **alphax, **alphay;
static double complex **beta;
static double complex **gammax, **gammay;

static double **tmpx, **tmpy;
static double **tmpxy;

static double omega;

static int optscale;
static int optreim;
static int seed;
static double complex C;

static char *input, *output, *rmsout, *initout, *Nstpout, *Npasout, *Nrunout, *finalout;
static enum OutputType { TEXTUAL = 1, BINARY = 2, VISUAL = 3 } outtype;
static enum OutputFlags { DEN_X = 1 << 0, DEN_Y = 1 << 1, DEN_XY = 1 << 2 } outflags;
static int outstpx, outstpy, outstpt;

static int nthreads;

void alloc();
void freemem();
void readparam(char *filename);

void initspace(double *x, double *x2, double *y, double *y2);
void initpsi(double complex **psi);
void initpot(double **pot);
void initcoef(double complex *Axm, double complex *Ax0r, double complex *Axp, double complex **alphax, double complex **gammax,
    double complex *Aym, double complex *Ay0r, double complex *Ayp, double complex **alphay, double complex **gammay);

double complex (*e)(double x);

void calcnu(double complex **psi, double **pot);
void calclux(double complex **psi, double complex **beta);
void calcluy(double complex **psi, double complex **beta);

void calcpsi2(double complex **psi, double **psi2);
void calcnorm(double *norm, double complex **psi, double **tmpx, double **tmpy);
void calcmuen(double *mu, double *en, double complex **psi, double **dpsir, double **dpsii,
    double **dpsixr, double **dpsixi, double **dpsiyr, double **dpsiyi);
void calcrms(double *rms, double **psi2, double **tmpx, double **tmpy);

void writeout(char *prefix, int step);

void outdenx(double **psi2, double **tmpy, double *den, char *filename);
void outdeny(double **psi2, double **tmpx, double *den, char *filename);

void outpsi2xy(double **psi2, double **den, char *filename);

void outpsi(double complex **psi, char *filename);

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

    g_2d = g0 / (sqrt(2. * PI) * d_z);

    if (Nstp > 0) {
        g = 0.;
    } else {
        g = optscale * g_2d;
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
        fprintf(fileout, "Imaginary time propagation 2d, OPTION = %d, NUM_THREADS = %d\n\n", optscale, nthreads);
    } else {
        fprintf(fileout, "Real time propagation 2d, OPTION = %d, NUM_THREADS = %d\n\n", optscale, nthreads);
    }

    if (cfg_defined("G0")) {
        fprintf(fileout, "Nonlinearity G_3D = %.6f, G_2D = %.6f\n", g0, g_2d);
    } else {
        fprintf(fileout, "Number of atoms N = %d, Unit of length AHO = %.8f m\n", Natoms, aho);
        fprintf(fileout, "Scattering length a = %.6f*a0\n", as);
        fprintf(fileout, "Calculated nonlinearity G_3D = %.6f, G_2D = %.6f\n", g0, g_2d);
    }

    fprintf(fileout, "Parameters of trap: GAMMA = %.2f, NU = %.2f\n", vgamma, vnu);
    fprintf(fileout, "Parameters of rotation: ANG VEL = %.4f\n", omega);
    fprintf(fileout, "Axial trap parameter = %.2f\n\n", d_z);
    fprintf(fileout, "Space steps: NX = %d, NY = %d\n", Nx, Ny);
    fprintf(fileout, "             DX = %.6f, DY = %.6f\n", dx, dy);
    fprintf(fileout, "Time steps:  NSTP = %d, NPAS = %d, NRUN = %d\n", Nstp, Npas, Nrun);
    fprintf(fileout, "             DT = %.6f\n",  dt);
    if (Nstp != 0) fprintf(fileout, "RNG seed:    SEED = %d\n",  seed);
    fprintf(fileout, "\n             -----------------------------------------------------\n");
    fprintf(fileout, "               Iter     Norm       Chem       Ener/N      <rho>\n");
    fprintf(fileout, "             -----------------------------------------------------\n");
    fflush(fileout);

    initspace(x, x2, y, y2);
    initpsi(psi);
    initpot(pot);
    initcoef(Axm, &Ax0r, Axp, alphax, gammax, Aym, &Ay0r, Ayp, alphay, gammay);

    calcnorm(&norm, psi, tmpx, tmpy);
    calcmuen(&mu, &en, psi, dpsir, dpsii, dpsixr, dpsixi, dpsiyr, dpsiyi);
    calcpsi2(psi, psi2);
    calcrms(rms, psi2, tmpx, tmpy);

    fprintf(fileout, "Initial:  %19.4f %11.5f %11.5f %10.5f\n", norm, mu / optscale, en / optscale, *rms);
    fflush(fileout);
    writeout(initout, 0);

    if (rmsout != NULL) {
        if (optreim == 1) {
            fprintf(filerms, "Imaginary time propagation 2d, OPTION = %d\n\n", optscale);
        } else {
            fprintf(filerms, "Real time propagation 2d, OPTION = %d\n\n", optscale);
        }
        fprintf(filerms, "           -----------------------------------------------\n");
        fprintf(filerms, "RMS size:    Iter      <r>          <x>          <y>\n");
        fprintf(filerms, "           -----------------------------------------------\n");
        fprintf(filerms, "Initial: %21.5f %12.5f %12.5f\n", rms[0], rms[1], rms[2]);
        fflush(filerms);
    }

    if (Nstp > 0) {
        double g_stp = optscale * g_2d / (double) Nstp;
        g = 0.;
        for (int i = 1; i <= Nstp; i ++) {
            g += g_stp;
            calcnu(psi, pot);
            calclux(psi, beta);
            calcluy(psi, beta);
            if (optreim == 1) calcnorm(&norm, psi, tmpx, tmpy);
        }
        if (optreim == 2) calcnorm(&norm, psi, tmpx, tmpy);
        calcmuen(&mu, &en, psi, dpsir, dpsii, dpsixr, dpsixi, dpsiyr, dpsiyi);
        calcpsi2(psi, psi2);
        calcrms(rms, psi2, tmpx, tmpy);

        fprintf(fileout, "NSTP iter.:  %16.4f %11.5f %11.5f %10.5f\n", norm, mu / optscale, en / optscale, *rms);
        fflush(fileout);

        writeout(Nstpout, 0);
        if (rmsout != NULL) {
            fprintf(filerms, "NSTP iter.: %18.5f %12.5f %12.5f\n", rms[0], rms[1], rms[2]);
            fflush(filerms);
        }
    }

    int step = 1;
    for (int i = 1; i <= Npas; i ++) {
        calcnu(psi, pot);
        calclux(psi, beta);
        calcluy(psi, beta);
        if (optreim == 1) calcnorm(&norm, psi, tmpx, tmpy);
        if (i % iter == 0) {
            if (optreim == 2) calcnorm(&norm, psi, tmpx, tmpy);
            calcmuen(&mu, &en, psi, dpsir, dpsii, dpsixr, dpsixi, dpsiyr, dpsiyi);
            calcpsi2(psi, psi2);
            calcrms(rms, psi2, tmpx, tmpy);

            fprintf(fileout, "NPAS iter.:  %5d %10.4f %11.5f %11.5f %10.5f\n", step, norm, mu / optscale, en / optscale, *rms);
            fflush(fileout);
            writeout(Npasout, step);

            if (rmsout != NULL) {
                fprintf(filerms, "NPAS iter.: %5d %12.5f %12.5f %12.5f\n", step, rms[0], rms[1], rms[2]);
                fflush(filerms);
            }

            step++;
        }
    }

    for (int i = 1; i <= Nrun; i ++) {
        calcnu(psi, pot);
        calclux(psi, beta);
        calcluy(psi, beta);
        if (optreim == 1) calcnorm(&norm, psi, tmpx, tmpy);
    }
    if (optreim == 2) calcnorm(&norm, psi, tmpx, tmpy);
    calcmuen(&mu, &en, psi, dpsir, dpsii, dpsixr, dpsixi, dpsiyr, dpsiyi);
    calcpsi2(psi, psi2);
    calcrms(rms, psi2, tmpx, tmpy);

    fprintf(fileout, "NRUN iter.:  %16.4f %11.5f %11.5f %10.5f\n", norm, mu / optscale, en / optscale, *rms);
    fflush(fileout);
    writeout(Nrunout, 0);
    if (rmsout != NULL) {
        fprintf(filerms, "NRUN iter.: %18.5f %12.5f %12.5f\n", rms[0], rms[1], rms[2]);
        fprintf(filerms, "           -----------------------------------------------\n");
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
    x2 = alloc_double_vector(Nx);
    y2 = alloc_double_vector(Ny);

    pot = alloc_double_matrix(Ny, Nx);

    Axm = alloc_complex_vector(Ny);
    Axp = alloc_complex_vector(Ny);
    alphax = alloc_complex_matrix(Ny, Nx);
    gammax = alloc_complex_matrix(Ny, Nx);

    Aym = alloc_complex_vector(Nx);
    Ayp = alloc_complex_vector(Nx);
    alphay = alloc_complex_matrix(Ny, Nx);
    gammay = alloc_complex_matrix(Ny, Nx);

    beta = alloc_complex_matrix(nthreads, MAX(Nx, Ny) - 1);

    psi = alloc_complex_matrix(Ny, Nx);
    psi2 = alloc_double_matrix(Ny, Nx);

    dpsir = alloc_double_matrix(Ny, Nx);
    dpsii = alloc_double_matrix(Ny, Nx);
    dpsixr = alloc_double_matrix(nthreads, Nx);
    dpsixi = alloc_double_matrix(nthreads, Nx);
    dpsiyr = alloc_double_matrix(nthreads, Ny);
    dpsiyi = alloc_double_matrix(nthreads, Ny);

    rms = alloc_double_vector(3);

    tmpx = dpsixr;
    tmpy = dpsiyr;
    tmpxy = dpsir;

    memset(*psi, 0, sizeof(double complex) * Ny * Nx);
}

/**
 *    Free all dynamically allocated memory.
 */
void freemem() {
    free_double_vector(x);
    free_double_vector(y);
    free_double_vector(x2);
    free_double_vector(y2);

    free_double_matrix(pot);

    free_complex_vector(Axm);
    free_complex_vector(Axp);
    free_complex_matrix(alphax);
    free_complex_matrix(gammax);

    free_complex_vector(Aym);
    free_complex_vector(Ayp);
    free_complex_matrix(alphay);
    free_complex_matrix(gammay);

    free_complex_matrix(beta);

    free_complex_matrix(psi);
    free_double_matrix(psi2);

    free_double_matrix(dpsir);
    free_double_matrix(dpsii);
    free_double_matrix(dpsixr);
    free_double_matrix(dpsixi);
    free_double_matrix(dpsiyr);
    free_double_matrix(dpsiyi);

    free_double_vector(rms);
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

    Nstp = cfg_read_int("NSTP");
    Npas = cfg_read_int("NPAS");
    Nrun = cfg_read_int("NRUN");
    iter = (cfg_defined("ITER")) ? Npas / cfg_read_int("ITER") : Npas;

    dx = cfg_read_double("DX");
    dy = cfg_read_double("DY");
    dt = cfg_read_double("DT");

    vgamma = cfg_read_double("GAMMA");
    vnu = cfg_read_double("NU");
    d_z = cfg_read_double("D_Z");

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
    }
}

/**
 *    Initialization of the space mesh.
 */
void initspace(double *x, double *x2, double *y, double *y2) {
    for (int i = 0; i < Nx; i ++) {
        x[i] = (i - Nx / 2) * dx;
        x2[i] = x[i] * x[i];
    }

    for (int j = 0; j < Ny; j ++) {
        y[j] = (j - Ny / 2) * dy;
        y2[j] = y[j] * y[j];
    }
}

/**
 *    Initialization of the wave function.
 *    psi - array with the wave function values
 */
void initpsi(double complex **psi) {
    double cpsi;
    double tmpr, tmpi, ctmp[2];
    double complex tmp;
    FILE *file;

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

        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                if (textual) {
                    if(fscanf(file,"%le %le\n", &tmpr, &tmpi) == 2) {
                        psi[j][i] = tmpr + I * tmpi;
                    }
                } else {
                    if (fread(ctmp, sizeof(double), 2, file) == 2) {
                        psi[j][i] = ctmp[0] + I * ctmp[1];
                    }
                }
            }
        }
        fclose(file);
    } else {
        cpsi = sqrt(sqrt(vgamma * vnu) / PI);

        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                psi[j][i] = cpsi * cexp(- 0.5 * (vgamma * x2[i] + vnu * y2[j]));
                if(ADD_RANDOM_PHASE == 1) psi[j][i] *= cexp(2 * PI * I * bec_rand());
                if(ADD_ONE_VORTEX == 1) psi[j][i] *= (x[i] + I * y[j]);
            }
        }
    }
}

/**
 *    Initialization of the potential.
 *    pot - array with the trap potential
 */
void initpot(double **pot) {
    double vgamma2, vnu2;

    vgamma2 = vgamma * vgamma;
    vnu2 = vnu * vnu;

    #pragma omp parallel for
    for (int j = 0; j < Ny; j ++) {
        for (int i = 0; i < Nx; i ++) {
            pot[j][i] = optscale * 0.5 * (vgamma2 * x2[i] + vnu2 * y2[j]);
        }
    }
}

/**
 *    Initialization of Crank-Nicolson scheme coefficients.
 */
void initcoef(double complex *Axm, double complex *Ax0r, double complex *Axp, double complex **alphax, double complex **gammax,
    double complex *Aym, double complex *Ay0r, double complex *Ayp, double complex **alphay, double complex **gammay) {

    double dx2, dy2;
    double complex cdt;
    double complex A0, Ax, Ay, Ax0, Ay0, Atmpx, Atmpy;

    dx2 = dx * dx;
    dy2 = dy * dy;

    cdt = C * dt;

    Ax0 = 1. + cdt / dx2;
    *Ax0r = 1. - cdt / dx2;
    Ax = 0.5 * cdt / dx2;
    Ay0 = 1. + cdt / dy2;
    *Ay0r = 1. - cdt / dy2;
    Ay = 0.5 * cdt / dy2;

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
}

/**
 *    Time propagation with respect to H1 (part of the Hamiltonian without spatial derivatives).
 *    psi - array with the wave function values
 *    pot - array with the trap potential
 */
void calcnu(double complex **psi, double **pot) {
    double psi2, tmp;

    #pragma omp parallel for private(psi2, tmp)
    for (int j = 0; j < Ny; j ++) {
        for (int i = 0; i < Nx; i ++) {
            psi2 = psi[j][i] * conj(psi[j][i]);
            tmp = dt * (pot[j][i] + g * psi2);
            psi[j][i] *= (*e)(- tmp);
        }
    }
}

/**
 *    Time propagation with respect to H2 (x-part of the Laplacian).
 *    psi  - array with the wave function values
 *    beta - Crank-Nicolson scheme coefficients
 */
void calclux(double complex **psi, double complex **beta) {
    double complex c;

    #pragma omp parallel private(c)
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int j = 0; j < Ny; j ++) {
            beta[threadid][Nx - 2] = (optreim == 1) ? 0. : psi[j][Nx - 1];
            for (int i = Nx - 2; i > 0; i --) {
                c = - Axp[j] * psi[j][i + 1] + Ax0r * psi[j][i] - Axm[j] * psi[j][i - 1];
                beta[threadid][i - 1] = gammax[j][i] * (Axp[j] * beta[threadid][i] - c);
            }
            psi[j][0] = 0.;
            for (int i = 0; i < Nx - 2; i ++) {
                psi[j][i + 1] = alphax[j][i] * psi[j][i] + beta[threadid][i];
            }
            psi[j][Nx - 1] = 0.;
        }
    }
}

/**
 *    Time propagation with respect to H3 (y-part of the Laplacian).
 *    psi  - array with the wave function values
 *    beta - Crank-Nicolson scheme coefficients
 */
void calcluy(double complex **psi, double complex **beta) {
    double complex c;

    #pragma omp parallel private(c)
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int i = 0; i < Nx; i ++) {
            beta[threadid][Ny - 2] = (optreim == 1) ? 0. : psi[Ny - 1][i];
            for (int j = Ny - 2; j > 0; j --) {
                c = - Ayp[i] * psi[j + 1][i] + Ay0r * psi[j][i] - Aym[i] * psi[j - 1][i];
                beta[threadid][j - 1] = gammay[j][i] * (Ayp[i] * beta[threadid][j] - c);
            }
            psi[0][i] = 0.;
            for (int j = 0; j < Ny - 2; j ++) {
                psi[j + 1][i] = alphay[j][i] * psi[j][i] + beta[threadid][j];
            }
            psi[Ny - 1][i] = 0.;
        }
    }
}

/**
 *    Calculation of squared wave function values.
 *    psi  - array with the wave function values
 *    psi2 - array with the squared wave function values
 */
void calcpsi2(double complex **psi, double **psi2) {
    #pragma omp parallel for
    for (int j = 0; j < Ny; j ++) {
        for (int i = 0; i < Nx; i ++) {
            psi2[j][i] = psi[j][i] * conj(psi[j][i]);
        }
    }
}

/**
 *    Calculation of the wave function norm and normalization.
 *    norm - wave function norm
 *    psi  - array with the wave function values
 *    tmpx - temporary array
 *    tmpy - temporary array
 */
void calcnorm(double *norm, double complex **psi, double **tmpx, double **tmpy) {
    double tmp;

    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                tmpx[threadid][i] = psi[j][i] * conj(psi[j][i]);
            }
            (*tmpy)[j] = simpint(dx, tmpx[threadid], Nx);
        }
    }

    *norm = sqrt(simpint(dy, *tmpy, Ny));

    if (optreim == 1) {
        tmp = 1. / *norm;

        #pragma omp parallel for
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                psi[j][i] *= tmp;
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
 */
void calcmuen(double *mu, double *en, double complex **psi, double **dpsir, double **dpsii,
    double **dpsixr, double **dpsixi, double **dpsiyr, double **dpsiyi) {

    int threadid;
    double psi2, psi2r, dpsiz;
    double norm;

    #pragma omp parallel private(threadid, psi2, psi2r, dpsiz)
    {
        threadid = omp_get_thread_num();

        #pragma omp for
        for (int i = 0; i < Nx; i ++) {
            cdiff(dy, psi[0] + i, dpsiyr[threadid], dpsiyi[threadid], Ny, Nx + 1); // Nx + 1 is due to padded allocation
            for (int j = 0; j < Ny; j ++) {
                dpsir[j][i] = dpsiyr[threadid][j] * dpsiyr[threadid][j];
                dpsii[j][i] = x[i] * dpsiyi[threadid][j];
            }
        }

        #pragma omp for
        for (int j = 0; j < Ny; j ++) {
            cdiff(dx, psi[j], dpsixr[threadid], dpsixi[threadid], Nx, 1);
            for (int i = 0; i < Nx; i ++) {
                dpsir[j][i] += dpsixr[threadid][i] * dpsixr[threadid][i];
                dpsii[j][i] -= y[j] * dpsixi[threadid][i];
            }
        }

        #pragma omp for
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                psi2 = psi[j][i] * conj(psi[j][i]);
                psi2r = creal(psi[j][i]);
                psi2r *= psi2r;
                dpsiz = creal(psi[j][i]) * dpsii[j][i];
                dpsixr[threadid][i] = (pot[j][i] + psi2 * g) * psi2r + dpsir[j][i] - optscale * (omega * dpsiz);
                dpsixi[threadid][i] = (pot[j][i] + 0.5 * psi2 * g) * psi2r + dpsir[j][i] - optscale * (omega * dpsiz);
            }
            (*dpsiyr)[j] = simpint(dx, dpsixr[threadid], Nx);
            (*dpsiyi)[j] = simpint(dx, dpsixi[threadid], Nx);
        }
    }

    *mu = simpint(dy, *dpsiyr, Ny);
    *en = simpint(dy, *dpsiyi, Ny);

    #pragma omp parallel private (threadid)
    {
        threadid = omp_get_thread_num();

        #pragma omp for
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                dpsixr[threadid][i] = creal(psi[j][i]) * creal(psi[j][i]);
            }
            (*dpsiyr)[j] = simpint(dx, dpsixr[threadid], Nx);
        }
    }

    norm = simpint(dy, *dpsiyr, Ny);

    *mu /= norm;
    *en /= norm;
}

/**
 *    Calculation of the root mean square radius.
 *    rms  - root mean square radius
 *    psi2 - array with the squared wave function values
 *    tmpx - temporary array
 *    tmpy - temporary array
 */
void calcrms(double *rms, double **psi2, double **tmpx, double **tmpy) {
    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                tmpx[threadid][i] = x2[i] * psi2[j][i];
            }
            (*tmpy)[j] = simpint(dx, tmpx[threadid], Nx);
        }

        #pragma omp single
        rms[1] = sqrt(simpint(dy, *tmpy, Ny));

        #pragma omp for
        for (int i = 0; i < Nx; i ++) {
            for (int j = 0; j < Ny; j ++) {
                tmpy[threadid][j] = y2[j] * psi2[j][i];
            }
            (*tmpx)[i] = simpint(dy, tmpy[threadid], Ny);
        }

        #pragma omp single
        rms[2] = sqrt(simpint(dx, *tmpx, Nx));
    }

    rms[0] = sqrt(rms[1] * rms[1] + rms[2] * rms[2]);
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

        if (outflags & DEN_XY) {
            sprintf(filename, "%s%s-2d", prefix, iter);
            outpsi2xy(psi2, tmpxy, filename);
        }

        if (outflags & DEN_X) {
            sprintf(filename, "%s%s-1d_x", prefix, iter);
            outdenx(psi2, tmpy, *tmpxy, filename);
        }

        if (outflags & DEN_Y) {
            sprintf(filename, "%s%s-1d_y", prefix, iter);
            outdeny(psi2, tmpx, *tmpxy, filename);
        }
    }
}

void outdenx(double **psi2, double **tmpy, double *den, char *filename) {
    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int i = 0; i < Nx; i += outstpx) {
            for (int j = 0; j < Ny; j ++) {
                tmpy[threadid][j] = psi2[j][i];
            }
            den[i] = simpint(dy, tmpy[threadid], Ny);
        }
    }

    switch(outtype) {
        case TEXTUAL: out_curve_txt(filename, Nx, outstpx, x, den); break;
        case BINARY:  out_curve_bin(filename, Nx, outstpx, x, den); break;
        case VISUAL:  out_curve_vis(filename, Nx, outstpx, x, den, "Integrated density X"); break;
    }
}

void outdeny(double **psi2, double **tmpx, double *den, char *filename) {
    #pragma omp parallel
    {
        int threadid = omp_get_thread_num();

        #pragma omp for
        for (int j = 0; j < Ny; j += outstpy) {
            for (int i = 0; i < Nx; i ++) {
                tmpx[threadid][i] = psi2[j][i];
            }
            den[j] = simpint(dx, tmpx[threadid], Nx);
        }
    }

    switch(outtype) {
        case TEXTUAL: out_curve_txt(filename, Ny, outstpy, y, den); break;
        case BINARY:  out_curve_bin(filename, Ny, outstpy, y, den); break;
        case VISUAL:  out_curve_vis(filename, Ny, outstpy, y, den, "Integrated density Y"); break;
    }
}

void outpsi2xy(double **psi2, double **den, char *filename) {
    #pragma omp parallel for
    for (int j = 0; j < Ny; j += outstpy) {
        for (int i = 0; i < Nx; i += outstpx) {
            den[j][i] = psi2[j][i];
        }
    }

    switch(outtype) {
        case TEXTUAL: out_matrix_txt(filename, Ny, Nx, outstpy, outstpx, y, x, den); break;
        case BINARY:  out_matrix_bin(filename, Ny, Nx, outstpy, outstpx, y, x, den); break;
        case VISUAL:  out_matrix_vis(filename, Ny, Nx, outstpy, outstpx, y, x, den, "Density"); break;
    }
}

void outpsi(double complex **psi, char *filename) {
    switch(outtype) {
        case TEXTUAL: out_complex_matrix_txt(filename, Ny, Nx, psi); break;
        case BINARY:  out_complex_matrix_bin(filename, Ny, Nx, psi); break;
        case VISUAL:  out_complex_matrix_bin(filename, Ny, Nx, psi); break;
    }
}
