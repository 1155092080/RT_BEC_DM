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
#include "diffint.h"

/**
 *    Spatial 1D integration with Simpson's rule.
 *    h - space step
 *    f - array with the function values
 *    N - number of integration points
 */
double simpint(double h, double *f, int N) {
    /*int c;
    double sum;

    sum = f[0];
    c = 4;
    for (int i = 1; i < N - 1; i ++) {
        //c = 2 + 2 * (i % 2);
        sum += c * f[i];
        c = 6 - c;
    }
    sum += f[N - 1];

    return sum * h / 3.;*/

    double sum;
    double f1, f2;

    f1 = f[1] + f[N - 2];
    f2 = f[2];

    for (int i = 3; i < N - 2; i += 2) {
        f1 += f[i];
        f2 += f[i + 1];
    }

    return h * (f[0] + 4. * f1 + 2. * f2 + f[N - 1]) / 3.;
}

/**
 *    Richardson extrapolation formula for calculation of space
 *    derivatives.
 *    h  - space step
 *    f  - array with the function values
 *    df - array with the first derivatives of the function
 *    N  - number of space mesh points
 */
void diff(double h, double *f, double *df, int N) {
    df[0] = 0.;
    df[1] = (f[2] - f[0]) / (2. * h);

    for (int i = 2; i < N - 2; i ++) {
        df[i] = (f[i - 2] - 8. * f[i - 1] + 8. * f[i + 1] - f[i + 2]) / (12. * h);
    }

    df[N - 2] = (f[N - 1] - f[N - 3]) / (2. * h);
    df[N - 1] = 0.;
}

/**
 *    Richardson extrapolation formula for calculation of space
 *    derivatives of a complex-valued function.
 *    h   - space step
 *    f   - complex array with the function values
 *    dfr - array with the first derivatives of the real part of function
 *    dfi - array with the first derivatives of the imaginary part of function
 *    N   - number of space mesh points
 */
void cdiff(double h, double complex *f, double *dfr, double *dfi, int N, int stride) {
    dfr[0] = 0.;
    dfi[0] = 0.;

    dfr[1] = (creal(f[stride * 2]) - creal(f[stride * 0])) / (2. * h);
    dfi[1] = (cimag(f[stride * 2]) - cimag(f[stride * 0])) / (2. * h);

    for (int i = 2; i < N - 2; i ++) {
        dfr[i] = (creal(f[stride * (i - 2)]) - 8. * creal(f[stride * (i - 1)]) + 8. * creal(f[stride * (i + 1)]) - creal(f[stride * (i + 2)])) / (12. * h);
        dfi[i] = (cimag(f[stride * (i - 2)]) - 8. * cimag(f[stride * (i - 1)]) + 8. * cimag(f[stride * (i + 1)]) - cimag(f[stride * (i + 2)])) / (12. * h);
    }

    dfr[N - 2] = (creal(f[stride * (N - 1)]) - creal(f[stride * (N - 3)])) / (2. * h);
    dfi[N - 2] = (cimag(f[stride * (N - 1)]) - cimag(f[stride * (N - 3)])) / (2. * h);
    
    dfr[N - 1] = 0.;
    dfi[N - 1] = 0.;
}
