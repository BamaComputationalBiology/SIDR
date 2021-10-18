
#ifndef sequence_h
#define sequence_h

#include <stdint.h>
#include <stdbool.h>


typedef struct sequence_encoding SEQCODE;

extern SEQCODE *encode(const char *, uint32_t, uint64_t *);
extern void insertSEQCODE(SEQCODE *, uint64_t);
extern uint64_t seed0(SEQCODE *);
extern uint64_t seed1(SEQCODE *);
extern void freeSEQCODE(SEQCODE *);
extern uint64_t get_kmer(SEQCODE *, uint32_t, uint32_t);

extern void printBits(uint64_t);

#endif
