



#ifndef taxonomy_h
#define taxonomy_h

#include "sequencing_data.h"

typedef struct taxonomy_container TAX;

extern TAX *newTAX(uint32_t);
extern bool checkTAX(TAX *, uint32_t);
extern void insertTAX(TAX *, uint32_t);
extern void freeTAX(TAX *);
extern void displayTAX(TAX *);
extern uint32_t numHits(TAX *);
extern uint32_t getTaxID(TAX *, uint32_t);

extern TAX **parse_blast(char *, char **, uint32_t);

extern bool check_delnodes(const char *, const uint32_t);
extern void check_merged(const char *, uint32_t *);
extern void check_nodes(const char *, uint32_t *, const uint32_t, const char *);
extern bool check_names(const char *, const uint32_t, const char *);

extern char *strlower(char *);


#endif
