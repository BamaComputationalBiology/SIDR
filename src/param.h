

#ifndef param_h
#define param_h

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdint.h>

typedef struct {

    uint32_t n_threads;

    // filepaths
    char *assembly;
    char *alignment;
    char *blast;
    char *delnodes;
    char *merged;
    char *nodes;
    char *names;

    // taxonomy info
    char *rank;
    char *classification;

    // kmer distribution info
    uint32_t max_kmer_freq;
    uint32_t n_kmers;
    uint32_t *kmer_list;

    PyObject *module;

} PARAM;

extern void initPARAM(PARAM *, int, char **);
extern void freePARAM(PARAM *);

#endif
