
#ifndef KSSDHASH_H
#define KSSDHASH_H
#include <stdint.h>

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
  key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
  key = key ^ key >> 28;
  key = (key + (key << 31)) & mask;
  return key;
};
// hash 4: Cheap but high-quality 64-bit mix (SplitMix64 finalizer).
static inline uint64_t mix64(uint64_t x){
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
};

#endif