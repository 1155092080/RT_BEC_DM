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
#ifndef BEC_RAND_H
#define BEC_RAND_H

static unsigned long int rand_next = 1;

double bec_rand(void);
void bec_srand(unsigned int);

#endif /* BEC_RAND_H */
