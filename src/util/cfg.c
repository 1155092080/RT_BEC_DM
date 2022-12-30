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
#include "cfg.h"

static int cfg_size;
static char cfg_key[256][256], cfg_val[256][256];

/**
 *    Configuration property value.
 *    key      - property
 *    required - is property required
 */
static char *cfg_read(char *key, bool required) {
    for (int i = 0; i < cfg_size; i ++) {
        if (! strcmp(key, cfg_key[i])) {
            return cfg_val[i];
        }
    }
    
    if (required) {
        fprintf(stderr, "%s is not defined in the configuration file.\n", key);
        exit(EXIT_FAILURE);
    }

    return NULL;
}

/**
 *    Configuration file parsing.
 *    cfg_file - stands for a configuration file, which is supplied on a command line
 */
int cfg_init(char *cfg_file) {
    FILE *file;
    char buf[256];
    
    file = fopen(cfg_file, "r");
    if (! file) return 0;
    
    cfg_size = 0;
    while (fgets(buf, 256, file) != NULL) {
        if (sscanf(buf, "%s = %s", cfg_key[cfg_size], cfg_val[cfg_size]) == 2) {
            cfg_size ++;
        }
    }
    
    fclose(file);
    
    return cfg_size;
}

/**
 *    Read integer value from configuration file.
 *    key - property
 */
int cfg_read_int(char *key) {
    return atoi(cfg_read(key, true));
}

/**
 *    Read long value from configuration file.
 *    key - property
 */
long cfg_read_long(char *key) {
    return atol(cfg_read(key, true));
}

/**
 *    Read double value from configuration file.
 *    key - property
 */
double cfg_read_double(char *key) {
    return atof(cfg_read(key, true));
}

/**
 *    Read string from configuration file.
 *    key - property
 */
char * cfg_read_string(char *key) {
    return cfg_read(key, true);
}

/**
 *    Check if property exists in configuration file.
 *    key - property
 */
bool cfg_defined(char *key) {
    return cfg_read(key, false) != NULL;
}
