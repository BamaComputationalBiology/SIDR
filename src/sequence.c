


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>

#include "sequence.h"




struct sequence_encoding {

    uint32_t capacity;
    uint64_t *code;

};



static SEQCODE *newSEQCODE(uint32_t seq_len) {
//  printf("inside newSEQCODE in sequence.c\n");

  SEQCODE *seqcode = malloc(sizeof(SEQCODE));
  assert(seqcode != 0);
  seqcode->capacity = seq_len / 32;
  if(seq_len % 32 != 0) seqcode->capacity++;
  seqcode->code = malloc(seqcode->capacity * sizeof(uint64_t));
  assert(seqcode->code != 0);
  return seqcode;

}



uint64_t seed0(SEQCODE *seqcode) {
//    printf("inside seed0 sequence.c\n");

    return seqcode->code[0];

}


uint64_t seed1(SEQCODE *seqcode) {
//    printf("inside seed1 sequence.c\n");


    return ~seqcode->code[0];

}



void freeSEQCODE(SEQCODE *seqcode) {
//    printf("inside freeSEQCODE sequence.c\n");
    free(seqcode->code);
    free(seqcode);

}



SEQCODE *encode(const char *seq, uint32_t seq_len, uint64_t *gc_count) {
//    printf("inside seqcode encode sequence.c\n");

    SEQCODE *seqcode = newSEQCODE(seq_len);

    for(uint32_t i = 0; i < seqcode->capacity; i++) {

        uint64_t code = 0;
        char seq_buff[32];
        uint32_t buff_len;

        if(i < seqcode->capacity - 1) buff_len = 32;
        else buff_len = seq_len % 32;

        memcpy(seq_buff, &seq[i * 32], buff_len);

        for(uint32_t j = 0; j < buff_len; j++) {
            switch(seq_buff[j]) {
            //    case 'A':
            //    case 'a':
            //        code |= ((uint64_t)0x0 << (j * 2));
            //        break;
                  case 'C':
                  case 'c':
                      code |= ((uint64_t)0x1 << (j * 2));
                      (*gc_count)++;
                      break;
                  case 'G':
                  case 'g':
                      code |= ((uint64_t)0x2 << (j * 2));
                      (*gc_count)++;
                      break;
                  case 'T':
                  case 't':
                      code |= ((uint64_t)0x3 << (j * 2));
                      break;
            }
        }

        seqcode->code[i] = code;

    }

    return seqcode;

}


// zero-based indexing
uint64_t get_kmer(SEQCODE *seqcode, uint32_t k, uint32_t pos) {
//    printf("inside get_kmer sequence.c\n");

    uint64_t kmer = 0;
    uint32_t div = pos / 32;
    uint32_t mod = pos % 32;
    uint64_t buffer = seqcode->code[div];

    if(32 - mod < k) {
        kmer |= (buffer >> (mod * 2));
        uint64_t overflow = seqcode->code[div + 1];
        kmer |= (overflow << ((64 - mod - k) * 2) >> ((32 - k) * 2));
    }
    else if(32 - mod == k) {
        kmer |= (buffer >> (mod * 2));
    }
    else {
        kmer |= (buffer << ((32 - mod - k) * 2) >> ((32 - k) * 2));
    }

    return kmer;

}


// https://stackoverflow.com/questions/52845040/printing-a-long-in-binary-64-bit-representation
void printBits(uint64_t n){
//    printf("inside printBits sequence.c\n");
    uint64_t i = 1UL<<(sizeof(n)*CHAR_BIT-1);
    while(i>0){
         if(n&i)
              printf("1");
         else
              printf("0");
         i >>= 1;
    }
    printf("\n\n");
}
