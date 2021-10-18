

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "sequencing_data.h"



struct sequencing_data {

    char *name;                           // region name
    uint32_t length;                      // region length
    float gc;                             // gc content %
    float cov;                            // read coverage
    bool blast;                           // does region have blast hit?
    bool target;                          // does blast hit match target classification?
    uint32_t **k_dists;                   // kmer coverage distributions

};



SEQDATA *newSEQDATA(uint32_t n_kmers) {
//printf("inside newseqdata seqdata.c\n");

    SEQDATA *seqdata = malloc(sizeof(SEQDATA));
    assert(seqdata != 0);
    seqdata->k_dists = malloc(n_kmers * sizeof(uint32_t *));
    assert(seqdata->k_dists != 0);
    return seqdata;

}



void freeSEQDATA(SEQDATA *seqdata, uint32_t n_kmers) {
//printf("inside freeseqdata seqdata.c\n");

    if (seqdata == 0) return;
    if (seqdata->name) {
        for (uint32_t k = 0; k < n_kmers; k++) free(seqdata->k_dists[k]);
        free(seqdata->k_dists);
        free(seqdata->name);
    }

    free(seqdata);

}



void displaySEQDATA(SEQDATA *seqdata, uint32_t n_kmers, uint32_t max_freq) {
//printf("inside displayseqdata sequencingdata.c\n");
    printf("Region: %s\n", seqdata->name);
    printf("Length: %d\n", seqdata->length);
    printf("GC Content: %f\n", seqdata->gc);
    printf("Read Coverage: %f\n", seqdata->cov);
    printf("BLAST Hit: %s\n", seqdata->blast ? "true" : "false");
    printf("Target taxon: %s\n", seqdata->target ? "true" : "false");
    printf("Kmer dists:\n");
    for(uint32_t i = 0; i < n_kmers; i++){
        for(uint32_t j = 0; j < max_freq; j++) printf("%lu ", seqdata->k_dists[i][j]);
        printf("\n");
    }
    printf("\n");

}



char *get_name(SEQDATA *seqdata) {
//printf("inside get_name sequencingdata.c\n");
    
    return seqdata->name;

}



uint32_t get_length(SEQDATA *seqdata) {
//printf("inside get_length sequencingdata.c\n");
    return seqdata->length;

}



float get_cov(SEQDATA *seqdata) {

//printf("inside get_cov sequencingdata.c\n");
    return seqdata->cov;

}



float get_gc(SEQDATA *seqdata) {
//printf("inside get_gc sequencingdata.c\n");
    return seqdata->gc;

}



bool get_blast(SEQDATA *seqdata) {
//printf("inside get_blast sequencingdata.c\n");
    return seqdata->blast;
}



bool get_tax(SEQDATA *seqdata) {
//printf("inside get_Tax sequencingdata.c\n");
    return seqdata->target;

}



extern uint32_t get_kpoint(SEQDATA *seqdata, uint32_t k, uint32_t i) {
//printf("inside get_kpoint sequencingdata.c\n");

    return seqdata->k_dists[k][i];

}


// called by process_taxonomy in pipeline.c
void update_name(SEQDATA *seqdata, const char *name) {
//printf("inside update_name sequencingdata.c\n");

    seqdata->name = malloc(strlen(name+1));
    strcpy(seqdata->name, name);

}



void update_length(SEQDATA *seqdata, uint32_t length) {
//printf("inside update_length sequencingdata.c\n");

    seqdata->length = length;

}



void update_gc(SEQDATA *seqdata, uint64_t gc_count) {
//printf("inside update_gc sequencingdata.c\n");

    assert(seqdata->length != 0);
    float gc = ((double)gc_count / (float)seqdata->length) * 100;
    seqdata->gc = gc;

}



void update_cov(SEQDATA *seqdata, uint64_t bp_count) {
//printf("inside update_cov sequencingdata.c\n");

    assert(seqdata->length != 0);
    float cov = (double)bp_count / (float)seqdata->length;
    seqdata->cov = cov;

}



void update_blast(SEQDATA *seqdata, bool has_blast) {
//printf("inside update_blast sequencingdata.c\n");

    seqdata->blast = has_blast;

}



void update_tax(SEQDATA *seqdata, bool is_target) {
//printf("inside update_tax sequencingdata.c\n");

    seqdata->target = is_target;

}



void update_kdist(SEQDATA *seqdata, uint32_t k, uint32_t *dist) {
//printf("inside updatekdist sequencingdata.c\n");
    seqdata->k_dists[k] = dist;
    
}
