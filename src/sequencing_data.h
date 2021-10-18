

#ifndef sequencing_data_h
#define sequencing_data_h

#include <stdbool.h>
#include <stdint.h>


typedef struct sequencing_data SEQDATA;

extern SEQDATA *newSEQDATA(uint32_t);
extern void freeSEQDATA(SEQDATA *, uint32_t);
extern void displaySEQDATA(SEQDATA *, uint32_t, uint32_t);

extern char *get_name(SEQDATA *);
extern uint32_t get_length(SEQDATA *);
extern float get_gc(SEQDATA *);
extern float get_cov(SEQDATA *);
extern bool get_blast(SEQDATA *);
extern bool get_tax(SEQDATA *);
extern uint32_t get_kpoint(SEQDATA *, uint32_t, uint32_t);

extern void update_name(SEQDATA *, const char *);
extern void update_length(SEQDATA *, uint32_t);
extern void update_cov(SEQDATA *, uint64_t);
extern void update_gc(SEQDATA *, uint64_t);
extern void update_blast(SEQDATA *, bool);
extern void update_tax(SEQDATA *, bool);
extern void update_kdist(SEQDATA *, uint32_t, uint32_t *);



#endif
