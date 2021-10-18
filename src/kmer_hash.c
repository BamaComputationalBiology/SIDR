// based on https://github.com/tidwall/hashmap.c

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "kmer_hash.h"



typedef struct data_bucket {    // 20 bytes

    uint64_t hash:48;           // 8 bytes
    uint16_t dib;               // 2 bytes
    uint64_t kmer;              // 8 bytes
    uint16_t count;             // 2 bytes

} BUCKET;



struct kmer_hash {

    uint64_t seed0;
    uint64_t seed1;
    size_t bucketsz;
    size_t nbuckets;
    size_t count;
    size_t mask;
    BUCKET *buckets;
    BUCKET *temp;
    BUCKET *entry;

};



HASH *newHASH(uint32_t min, uint64_t seed0, uint64_t seed1) {
//    printf("inside HASH *newHASH kmer_hash.c\n");
    int cap = 256;
    if(cap < min) {
        while (cap < min) cap *= 2;
    }

    size_t bucketsz = sizeof(BUCKET);
    while (bucketsz & (sizeof(uintptr_t) - 1)) bucketsz++;

    // hashmap + (spare + edata)
    size_t size = sizeof(HASH) + (bucketsz * 2);
    HASH *map = malloc(size);
    assert(map != 0);
    memset(map, 0, sizeof(HASH));

    map->bucketsz = bucketsz;
    map->seed0 = seed0;
    map->seed1 = seed1;
    map->temp = (BUCKET *)((char *)map + sizeof(HASH));
    map->entry = (BUCKET *)((char *)map->temp + bucketsz);
    map->nbuckets = cap;
    map->mask = map->nbuckets-1;
    map->buckets = malloc(map->bucketsz * map->nbuckets);
    assert(map->buckets != 0);
    memset(map->buckets, 0, map->bucketsz * map->nbuckets);

    return map;

}



void clearHASH(HASH *map) {
//    printf("inside clearHASHi kmerhas.c \n");

    map->count = 0;
    memset(map->buckets, 0, map->bucketsz * map->nbuckets);

}



void insertHASH(HASH *map, uint64_t kmer) {
//    printf("inside insertHASH kmerhash.c\n");    

    BUCKET *entry = map->entry;
    entry->hash = hash_func(&kmer, sizeof(uint64_t), map->seed0, map->seed1);
    entry->dib = 1;
    entry->kmer = kmer;
    entry->count = 1;

    size_t i = entry->hash & map->mask;

	  for (;;) {

        BUCKET *bucket = &map->buckets[i];

        if (bucket->dib == 0) {
            memcpy(bucket, entry, map->bucketsz);
            map->count++;
            break;
        }
        if (entry->hash == bucket->hash && entry->kmer == bucket->kmer){
            bucket->count++;
            break;
        }
        if (bucket->dib < entry->dib) {
            memcpy(map->temp, bucket, map->bucketsz);
            memcpy(bucket, entry, map->bucketsz);
            memcpy(entry, map->temp, map->bucketsz);
  		  }

        entry->dib++;
		    i = (i + 1) & map->mask;

	  }

}



void displayHASH(HASH *map) {
//    printf("inside displayHASH\n");

    for(size_t i = 0; i < map->nbuckets; i++) {
        BUCKET *bucket = &map->buckets[i];
        printf("[%lu, %lu]\n", bucket->kmer, bucket->count);
    }

}



void freeHASH(HASH *map) {
//    printf("inside freehash kmerhash.c\n");

    if (!map) return;
    free(map->buckets);
    free(map);

}



uint32_t *distHASH(HASH *map, uint32_t max_freq){
//    printf("inside distHASH kmer_hash.c\n");

    uint32_t *dist = malloc(max_freq * sizeof(uint32_t));
    memset(dist, 0, max_freq * sizeof(uint32_t));
    if(map->count == 0) return dist;
    for(size_t i = 0; i < map->nbuckets; i++) {
      BUCKET *bucket = &map->buckets[i];
      if (bucket->dib) {
          if (bucket->count < max_freq) dist[bucket->count]++;
          else dist[max_freq-1]++;
      }
    }

    return dist;

}



//-----------------------------------------------------------------------------
// SipHash reference C implementation
//
// Copyright (c) 2012-2016 Jean-Philippe Aumasson
// <jeanphilippe.aumasson@gmail.com>
// Copyright (c) 2012-2014 Daniel J. Bernstein <djb@cr.yp.to>
//
// To the extent possible under law, the author(s) have dedicated all copyright
// and related and neighboring rights to this software to the public domain
// worldwide. This software is distributed without any warranty.
//
// You should have received a copy of the CC0 Public Domain Dedication along
// with this software. If not, see
// <http://creativecommons.org/publicdomain/zero/1.0/>.
//
// default: SipHash-2-4
//-----------------------------------------------------------------------------
static uint64_t SIP64(const uint8_t *in, const size_t inlen, uint64_t seed0, uint64_t seed1) {
 //    printf("inside uint64_t kmer_hash.c\n");
#define U8TO64_LE(p) \
    {  (((uint64_t)((p)[0])) | ((uint64_t)((p)[1]) << 8) | \
        ((uint64_t)((p)[2]) << 16) | ((uint64_t)((p)[3]) << 24) | \
        ((uint64_t)((p)[4]) << 32) | ((uint64_t)((p)[5]) << 40) | \
        ((uint64_t)((p)[6]) << 48) | ((uint64_t)((p)[7]) << 56)) }
#define U64TO8_LE(p, v) \
    { U32TO8_LE((p), (uint32_t)((v))); \
      U32TO8_LE((p) + 4, (uint32_t)((v) >> 32)); }
#define U32TO8_LE(p, v) \
    { (p)[0] = (uint8_t)((v)); \
      (p)[1] = (uint8_t)((v) >> 8); \
      (p)[2] = (uint8_t)((v) >> 16); \
      (p)[3] = (uint8_t)((v) >> 24); }
#define ROTL(x, b) (uint64_t)(((x) << (b)) | ((x) >> (64 - (b))))
#define SIPROUND \
    { v0 += v1; v1 = ROTL(v1, 13); \
      v1 ^= v0; v0 = ROTL(v0, 32); \
      v2 += v3; v3 = ROTL(v3, 16); \
      v3 ^= v2; \
      v0 += v3; v3 = ROTL(v3, 21); \
      v3 ^= v0; \
      v2 += v1; v1 = ROTL(v1, 17); \
      v1 ^= v2; v2 = ROTL(v2, 32); }
    uint64_t k0 = U8TO64_LE((uint8_t*)&seed0);
    uint64_t k1 = U8TO64_LE((uint8_t*)&seed1);
    uint64_t v3 = UINT64_C(0x7465646279746573) ^ k1;
    uint64_t v2 = UINT64_C(0x6c7967656e657261) ^ k0;
    uint64_t v1 = UINT64_C(0x646f72616e646f6d) ^ k1;
    uint64_t v0 = UINT64_C(0x736f6d6570736575) ^ k0;
    const uint8_t *end = in + inlen - (inlen % sizeof(uint64_t));
    for (; in != end; in += 8) {
        uint64_t m = U8TO64_LE(in);
        v3 ^= m;
        SIPROUND; SIPROUND;
        v0 ^= m;
    }
    const int left = inlen & 7;
    uint64_t b = ((uint64_t)inlen) << 56;
    switch (left) {
    case 7: b |= ((uint64_t)in[6]) << 48;
    case 6: b |= ((uint64_t)in[5]) << 40;
    case 5: b |= ((uint64_t)in[4]) << 32;
    case 4: b |= ((uint64_t)in[3]) << 24;
    case 3: b |= ((uint64_t)in[2]) << 16;
    case 2: b |= ((uint64_t)in[1]) << 8;
    case 1: b |= ((uint64_t)in[0]); break;
    case 0: break;
    }
    v3 ^= b;
    SIPROUND; SIPROUND;
    v0 ^= b;
    v2 ^= 0xff;
    SIPROUND; SIPROUND; SIPROUND; SIPROUND;
    b = v0 ^ v1 ^ v2 ^ v3;
    uint64_t out = 0;
    U64TO8_LE((uint8_t*)&out, b);
    return out;
}



uint64_t hash_func(const uint64_t *item, size_t len, uint64_t seed0, uint64_t seed1){
//    printf("inside hash_func kmerhash.c\n");

    return (SIP64((uint8_t *)item, len, seed0, seed1) << 16 >> 16);

}
