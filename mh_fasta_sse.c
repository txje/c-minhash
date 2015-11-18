#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <emmintrin.h>
#include "../klib/kseq.h" // FASTA/Q parser
#include "../klib/kvec.h" // C dynamic vector
#include "xxHash/xxhash.c" // fast hashing

typedef struct {
  uint32_t hash;
  uint32_t pos;
} Min;

typedef kvec_t(Min*) signature_list;

// init kseq struct
KSEQ_INIT(FILE*, fileread);

Min* minhash (char *s, int k, int h, uint32_t hash_seeds[], unsigned char reverse) {

  int i, j;
  int slen = strlen(s);

  // initialize minimums, one per hash seed
  Min* m = (Min*)malloc(sizeof(Min)*h);

  for(j = 0; j < h; j++) {
    m[j].hash = UINT32_MAX;
  }

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

    //kmer = XXH32(&s[i], k, 24242424); // this seed should not be random because hash_seeds is filled with random

    // just use uint32 kmer representation, with random xor
    if(reverse == 0) {
      kmer = (kmer << 2) + ((s[i+k-1] >> 1) & 3);
      if(k < 16) {
        kmer = kmer % (1 << (2*k));
      }
    } else {
      kmer = (kmer >> 2) + ((((s[i+k-1] >> 1) & 3) ^ 2) << rev_shift);
    }

    for(j = 0; j < h; j++) {
      kmer = kmer ^ hash_seeds[j];
      if(kmer < m[j].hash) {
        m[j].hash = kmer;
        m[j].pos = i;
      }
    }
  }

  //printf("first min: %d at %d\n", m[0].hash, m[0].pos);

  return m;
}

// if reverse is true (1), forward and reverse signatures will be adjacent such that
// signature of read X forward is at index (2*X), and reverse is (2*X + 1)
signature_list get_fasta_signatures(char* fa, int k, int h, uint32_t hash_seeds[], unsigned char reverse) {

  FILE* fp;
  kseq_t* seq;
  int l;
  signature_list signatures;
  kv_init(signatures);
  //kv_resize(Min*, signatures, <num_reads>); // resize all at once, assuming we know the total # reads

  fp = fopen(fa, "r");
  seq = kseq_init(fp);
  //printf("Reading fasta file: %s\n", fa);

  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l

    Min *m = minhash(seq->seq.s, k, h, hash_seeds, 0);
    kv_push(Min*, signatures, m);

    if(reverse == 1) {
      m = minhash(seq->seq.s, k, h, hash_seeds, 1);
      kv_push(Min*, signatures, m);
    }
  }

  kseq_destroy(seq);
  fclose(fp);

  return signatures;
}

// have to reorder params to make this work with kseq
// so stupid
int fileread(FILE* f, char* buffer, int size) {
  return fread(buffer, 1, size, f);
}

int mh_fasta(char *query_fa, char *target_fa, int k, int h, int seed, int threshold) {

  // construct array of hash seeds
  srand(seed); // seed random number generator
  // RAND_MAX must be at least 2**15, mine is 2**31
  // as long as we don't duplicate seeds, it should be fine
  uint32_t hash_seeds[h];
  int i;
  for(i = 0; i < h; i++) {
    hash_seeds[i] = rand(); // 0 to RAND_MAX
  }

  time_t t0 = time(NULL);
  // read FASTA file
  signature_list query_signatures = get_fasta_signatures(query_fa, k, h, hash_seeds, 1);
  time_t t1 = time(NULL);
  printf("# Loaded and processed %d reads in %d seconds\n", query_signatures.n, (t1-t0));
  t0 = t1;

  __m128i qhashbin, thashbin, qposbin, tposbin, matchMask, matches, offset;
  matchMask = _mm_set_epi32(1,1,1,1);
  uint32_t matchArr[4]; // to collect parts of SSE results at the end
  uint32_t offsetArr[4];
  uint32_t matchFinal;
  uint32_t offsetFinal;
  int q, qidx, qrev, t;

  // if target fasta is different
  if(strcmp(query_fa, target_fa) != 0) {
    signature_list target_signatures = get_fasta_signatures(target_fa, k, h, hash_seeds, 0);
    t1 = time(NULL);
    printf("# Loaded and processed %d (target) reads in %d seconds\n", target_signatures.n, (t1-t0));
    t0 = t1;
    for(q = 0; q < query_signatures.n/2; q++) {
      for(qrev = 0; qrev <= 1; qrev++) {
        if(qrev == 1) {
          qidx = q*2 + 1;
        } else {
          qidx = q*2;
        }
        Min* qmin = (Min*)(query_signatures.a[qidx]);
        for(t = 0; t < target_signatures.n; t++) {
          matches = _mm_set_epi32(0,0,0,0);
          offset = _mm_set_epi32(0,0,0,0);
          Min* tmin = (Min*)(target_signatures.a[t]);
          for(i = 0; i < h; i=i+4) { // do in sets of 4 uint32's w/SSE2 (128bits)
            qhashbin = _mm_loadu_si128(qhash+i);
            thashbin = _mm_loadu_si128(thash+i);
            qposbin = _mm_loadu_si128(qpos+i);
            tposbin = _mm_loadu_si128(tpos+i);
            hashmatch = _mm_cmpeq_epi32(qhashbin, thashbin);
            offset = _mm_add_epi32(offset, _mm_and_si128(hashmatch, _mm_sub_epi32(tposbin, qposbin)));
            matches = _mm_add_epi32(matches, _mm_and_si128(hashmatch, matchMask));
            // increment matches
          }
          // collapse matches and offset
          _mm_stream_si128(offsetArr, offset); // copy offset contents into array
          _mm_stream_si128(matchArr, matches);
          matchFinal = 0;
          offsetFinal = 0;
          for(i = 0; i < 4; i++) {
            matchFinal = matchFinal + matchArr[i];
            offsetFinal = offsetFinal + offsetArr[i];
          }

          if(matchFinal >= threshold) {
            offsetFinal = offsetFinal / matchesFinal;
            printf("%d,%d,%d,%d,%d\n", q, qrev, t, matchFinal, offsetFinal);
          }
        }
      }
    }
  } else {
    int q, qidx, qrev, q2, q2idx, matches, offset; // in this case, tpos refers to q2
    for(q = 0; q < query_signatures.n/2; q++) {
      for(qrev = 0; qrev <= 1; qrev++) {
        if(qrev == 1) {
          qidx = q*2 + 1;
        } else {
          qidx = q*2;
        }
        Min* qmin = (Min*)(query_signatures.a[qidx]);
        for(q2 = q+1; q2 < query_signatures.n/2; q2++) {
          q2idx = q2 * 2;
          matches = _mm_set_epi32(0,0,0,0);
          offset = _mm_set_epi32(0,0,0,0);
          Min* tmin = (Min*)(query_signatures.a[q2idx]);
          for(i = 0; i < h; i=i+4) { // do in sets of 4 uint32's w/SSE2 (128bits)
            qhashbin = _mm_loadu_si128(qhash+i);
            thashbin = _mm_loadu_si128(thash+i);
            qposbin = _mm_loadu_si128(qpos+i);
            tposbin = _mm_loadu_si128(tpos+i);
            hashmatch = _mm_cmpeq_epi32(qhashbin, thashbin);
            offset = _mm_add_epi32(offset, _mm_and_si128(hashmatch, _mm_sub_epi32(tposbin, qposbin)));
            matches = _mm_add_epi32(matches, _mm_and_si128(hashmatch, matchMask));
            // increment matches
          }
          // collapse matches and offset
          _mm_stream_si128(offsetArr, offset); // copy offset contents into array
          _mm_stream_si128(matchArr, matches);
          matchFinal = 0;
          offsetFinal = 0;
          for(i = 0; i < 4; i++) {
            matchFinal = matchFinal + matchArr[i];
            offsetFinal = offsetFinal + offsetArr[i];
          }

          if(matchFinal >= threshold) {
            offsetFinal = offsetFinal / matchesFinal;
            printf("%d,%d,%d,%d,%d\n", q, qrev, q2, matchFinal, offsetFinal);
          }
        }
      }
    }
  }

  t1 = time(NULL);
  printf("# Compared and output in %d seconds\n", (t1-t0));

  return 0;
}

int main(int argc, char *argv[]) {
  if(argc < 6) {
    printf("Usage: mh_fasta <query_fasta> <target_fasta> <k> <h> <seed> <threshold>\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *query_fasta = argv[1];
  char *target_fasta = argv[2];
  int k = atoi(argv[3]);
  if(k > 16) {
    printf("Current optimizations do not allow for k > 16");
    return -1;
  }
  int h = atoi(argv[4]);
  if(h % 4 != 0) {
    printf("[h] must be a multiple of 4 for SSE");
    return -1;
  }
  int seed = atoi(argv[5]);
  int threshold = atoi(argv[6]);
  return mh_fasta(query_fasta, target_fasta, k, h, seed, threshold);
}
