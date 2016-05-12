/* ***************************************************************************** */
/* A C-program for MEMT607-64                                                    */
/* Copyright:      Shin Harase, Ritsumeikan University                           */
/*                 Takamitsu Kimoto, Recruit Holdings Co., Ltd.                  */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact S. Harase at: harase @ fc.ritsumei.ac.jp       */
/*                                                                               */
/* Remark:         We recommend using the most significant bits (not taking the  */
/*                 least significant bits) because our generators are optimized  */
/*                 preferentially from the most significant bits,                */
/*                 see Remark 4.1 for details.                                   */
/* ***************************************************************************** */

#include <stdio.h>

#define NN 9
#define MM 5
#define MATRIX_A 0x81f1fd68012348bcULL
#define P 31
#define W 64
#define MASKU (0xffffffffffffffffULL << (W-P))
#define MASKL (~MASKU)
#define MAT3NEG(t, v) (v ^ (v << ((t))))
#define MAT3POS(t, v) (v ^ (v >> ((t))))
#define LAG1 3
#define SHIFT1 30
#define MASK1 0x66edc62a6bf8c826ULL
#define LAG1over 6

static unsigned long long mt[NN]; 
static int mti;
static unsigned long long lung;
static unsigned long long mag01[2]={0ULL, MATRIX_A};
static unsigned long long x;

static unsigned long long case_1(void);
static unsigned long long case_2(void);
static unsigned long long case_3(void);
static unsigned long long case_4(void);
unsigned long long (*genrand64_int64)(void);

/* initializes mt[NN] and lung with a seed */
void init_genrand64(unsigned long long seed)
{
  mt[0] = seed;
  for (mti=1; mti<NN; mti++) {
    mt[mti] = (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
  }
  lung = (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
  mti = 0;
  genrand64_int64 = case_1;
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array64(unsigned long long init_key[],
		     unsigned long long key_length)
{
	unsigned long long i, j, k;
    init_genrand64(19650218ULL);
    i=1; j=0;
    k = (NN>key_length ? NN : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 3935559000370003845ULL))
          + init_key[j] + j; /* non linear */
        i++; j++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=NN-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 2862933555777941757ULL))
          - i; /* non linear */
        i++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
    }

    mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */
	lung = (6364136223846793005ULL * (mt[NN-2] ^ (mt[NN-2] >> 62)) + NN-1);
	mti = 0;
}

static unsigned long long case_1(void) {
    x = (mt[mti] & MASKU) | (mt[mti+1] & MASKL);
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ mt[mti+MM] ^ MAT3NEG(13, lung);
    mt[mti] = x ^ MAT3POS(35, lung);
    x = mt[mti] ^ (mt[mti] << SHIFT1);
    x = x ^ (mt[mti + LAG1] & MASK1);
    ++mti;
    if (mti == NN - MM) genrand64_int64 = case_2;
    return x;
}

static unsigned long long case_2(void) {
    x = (mt[mti] & MASKU) | (mt[mti+1] & MASKL);
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ mt[mti+(MM-NN)] ^ MAT3NEG(13, lung);
    mt[mti] = x ^ MAT3POS(35, lung);
    x = mt[mti] ^ (mt[mti] << SHIFT1);
    x = x ^ (mt[mti + LAG1] & MASK1);
    ++mti;
    if (mti == LAG1over) genrand64_int64 = case_3;
    return x;
}

static unsigned long long case_3(void) {
    x = (mt[mti] & MASKU) | (mt[mti+1] & MASKL);
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ mt[mti+(MM-NN)] ^ MAT3NEG(13, lung);
    mt[mti] = x ^ MAT3POS(35, lung);
    x = mt[mti] ^ (mt[mti] << SHIFT1);
    x = x ^ (mt[mti - LAG1over] & MASK1);
    ++mti;
    if (mti == NN-1) genrand64_int64 = case_4;
    return x;
}

static unsigned long long case_4(void) {
    x = (mt[NN-1] & MASKU) | (mt[0] & MASKL);
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ mt[MM-1] ^ MAT3NEG(13, lung);
    mt[NN-1] = x ^ MAT3POS(35, lung);
    x = mt[mti] ^ (mt[mti] << SHIFT1);
    x = x ^ (mt[mti - LAG1over] & MASK1);
    mti = 0;
    genrand64_int64 = case_1;
    return x;
}

/* generates a random number on [0, 2^63-1]-interval */
long long genrand64_int63(void)
{
    return (long long)(genrand64_int64() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double genrand64_real1(void)
{
    return (genrand64_int64() >> 11) * (1.0/9007199254740991.0);
}

/* generates a random number on [0,1)-real-interval */
double genrand64_real2(void)
{
    return (genrand64_int64() >> 11) * (1.0/9007199254740992.0);
}

/* generates a random number on (0,1)-real-interval */
double genrand64_real3(void)
{
    return ((genrand64_int64() >> 12) + 0.5) * (1.0/4503599627370496.0);
}

/* generates a random number on [0,1)-real-interval using a union trick */
double genrand64_res53(void)
{
    union {
	unsigned long long u;
	double d;
    } conv;
	
	conv.u = (genrand64_int64() >> 12) | 0x3FF0000000000000ULL;
	
	return (conv.d - 1.0);
}

/* generates a random number on (0,1)-real-interval using a union trick */
double genrand64_res53_open(void)
{
    union {
	unsigned long long u;
	double d;
    } conv;
	
	conv.u = (genrand64_int64() >> 12) | 0x3FF0000000000001ULL;
	
	return (conv.d - 1.0);
}

int main(void)
{
    int i;
    unsigned long long init[4]={0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL}, length=4;
    init_by_array64(init, length);
    printf("1000 outputs of genrand64_int64()\n");
    for (i=0; i<1000; i++) {
      printf("%20llu ", genrand64_int64());
      if (i%5==4) printf("\n");
    }
    printf("\n1000 outputs of genrand64_real2()\n");
    for (i=0; i<1000; i++) {
      printf("%10.15f ", genrand64_real2());
      if (i%5==4) printf("\n");
    }
    printf("\n1000 outputs of genrand64_res53()\n");
    for (i=0; i<1000; i++) {
      printf("%10.15f ", genrand64_res53());
      if (i%5==4) printf("\n");
    }
    return 0;
}
