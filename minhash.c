#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "minhash.h"
#include "xxHash/xxhash.c"

Min minhash (char *s, int k, int seed) {
  int i = 0;
  Min m;
  m.hash = UINT32_MAX;
  m.pos = 0;
  int slen = strlen(s);
  char kmer[k];
  for(; i < slen-k; i++) {
    memcpy(kmer, &s[i], k); // get k-mer as a substring of s
    uint32_t hsh = XXH32(kmer, k, seed);
    //printf("kmer: %s, hash value: %u\n", kmer, hsh);
    if(hsh < m.hash) {
      m.hash = hsh;
      m.pos = i;
    }
  }
  return m;
}

Min minhash_bitpack (char *s, int k, unsigned char reverse) {

  int i;
  int slen = strlen(s);

  Min m;
  m.hash = UINT32_MAX;
  m.pos = 0;

  // set up kmer uint32
  uint32_t kmer = 0;
  int rev_shift = 2*(k-1);
  for(i = 0; i < k-1; i++) {
    if(reverse == 0) {
      kmer = (kmer << 2) + ((s[i] >> 1) & 3);
    } else {
      kmer = (kmer >> 2) + ((((s[i] >> 1) & 3) ^ 2) << rev_shift);
    }
  }
  for(i = 0; i <= slen-k; i++) {

    // just use uint32 kmer representation, with random xor
    if(reverse == 0) {
      kmer = (kmer << 2) + ((s[i+k-1] >> 1) & 3);
      if(k < 16) {
        kmer = kmer % (1 << (2*k));
      }
    } else {
      kmer = (kmer >> 2) + ((((s[i+k-1] >> 1) & 3) ^ 2) << rev_shift);
    }
    printf("pos %i, val: %u\n", i, kmer);

    if(kmer < m.hash) {
      m.hash = kmer;
      m.pos = i;
    }
  }

  return m;
}

int main(int argc, char *argv[]) {
  if(argc < 4) {
    printf("Usage: minhash <string> <k> <seed>\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *s = argv[1];
  int k = atoi(argv[2]);
  /*
  int seed = atoi(argv[3]);
  Min m = minhash(s, k, seed);
  */
  Min m = minhash_bitpack(s, k, 0);
  printf("Min hash value: %u\n", m.hash);
  printf("Min hash pos: %u\n", m.pos);

  m = minhash_bitpack(s, k, 1);
  printf("RC Min hash value: %u\n", m.hash);
  printf("RC Min hash pos: %u\n", m.pos);
  return 0;
}
