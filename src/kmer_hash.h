// Copyright 2020 Joshua J Baker. All rights reserved.
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file.

#ifndef HASHMAP_H
#define HASHMAP_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>


typedef struct kmer_hash HASH;

extern HASH *newHASH(uint32_t, uint64_t, uint64_t);
extern void clearHASH(HASH *);
extern void insertHASH(HASH *, uint64_t);
extern void displayHASH(HASH *);
extern void freeHASH(HASH *);
extern uint32_t *distHASH(HASH *, uint32_t);

extern uint64_t hash_func(const uint64_t *, size_t, uint64_t, uint64_t);


#endif
