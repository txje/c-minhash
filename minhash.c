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

int main(int argc, char *argv[]) {
  if(argc < 4) {
    printf("Usage: minhash <string> <k> <seed>\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *s = argv[1];
  int k = atoi(argv[2]);
  int seed = atoi(argv[3]);
  Min m = minhash(s, k, seed);
  printf("Min hash value: %i\n", m.hash);
  printf("Min hash pos: %i\n", m.pos);
  return 0;
}
