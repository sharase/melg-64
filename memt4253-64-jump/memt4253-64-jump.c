/* ***************************************************************************** */
/* A C-program for MEMT4253-64-JUMP                                              */
/* Copyright:      Shin Harase, Ritsumeikan University                           */
/*                 Takamitsu Kimoto, Recruit Holdings Co., Ltd.                  */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact S. Harase at: harase @ fc.ritsumei.ac.jp       */
/*                                                                               */
/* Remark:         We recommend using the most significant bits (not taking the  */
/*                 least significant bits) because our generators are optimized  */
/*				   preferentially from the most significant bits,                */
/*                 see Remark 4.1 for details.                                   */
/* ***************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NN 66 //N-1
#define MM 29
#define MATRIX_A 0xfac1e8c56471d722ULL
#define P 29 //W-r
#define W 64
#define MASKU (0xffffffffffffffffULL << (W-P))
#define MASKL (~MASKU)
#define MAT3NEG(t, v) (v ^ (v << ((t))))
#define MAT3POS(t, v) (v ^ (v >> ((t))))
#define LAG1 9
#define SHIFT1 5
#define MASK1 0xcb67b0c18fe14f4dULL
#define LAG1over 57 //NN-LAG1

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

struct memt_state{
	unsigned long long int lung;
	unsigned long long int mt[NN];
	int mti;
	unsigned long long int (*function_p)(void);
};

static void add(struct memt_state *state);

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
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ mt[mti+MM] ^ MAT3NEG(30, lung);
    mt[mti] = x ^ MAT3POS(20, lung);
    x = mt[mti] ^ (mt[mti] << SHIFT1);
    x = x ^ (mt[mti + LAG1] & MASK1);
    ++mti;
    if (mti == NN - MM) genrand64_int64 = case_2;
    return x;
}

static unsigned long long case_2(void) {
    x = (mt[mti] & MASKU) | (mt[mti+1] & MASKL);
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ mt[mti+(MM-NN)] ^ MAT3NEG(30, lung);
    mt[mti] = x ^ MAT3POS(20, lung);
    x = mt[mti] ^ (mt[mti] << SHIFT1);
    x = x ^ (mt[mti + LAG1] & MASK1);
    ++mti;
    if (mti == LAG1over) genrand64_int64 = case_3;
    return x;
}

static unsigned long long case_3(void) {
    x = (mt[mti] & MASKU) | (mt[mti+1] & MASKL);
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ mt[mti+(MM-NN)] ^ MAT3NEG(30, lung);
    mt[mti] = x ^ MAT3POS(20, lung);
	x = mt[mti] ^ (mt[mti] << SHIFT1);
    x = x ^ (mt[mti - LAG1over] & MASK1);
    ++mti;
    if (mti == NN-1) genrand64_int64 = case_4;
    return x;
}

static unsigned long long case_4(void) {
    x = (mt[NN-1] & MASKU) | (mt[0] & MASKL);
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ mt[MM-1] ^ MAT3NEG(30, lung);
    mt[NN-1] = x ^ MAT3POS(20, lung);
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
	unsigned long long int u;
	double d;
    } conv;
	
	conv.u = (genrand64_int64() >> 12) | 0x3FF0000000000000ULL;
	
	return (conv.d - 1.0);
}

/* generates a random number on (0,1)-real-interval using a union trick */
double genrand64_res53_open(void)
{
    union {
	unsigned long long int u;
	double d;
    } conv;
	
	conv.u = (genrand64_int64() >> 12) | 0x3FF0000000000001ULL;
	
	return (conv.d - 1.0);
}

/* This is a jump function for the generator. It is equivalent
   to 2^256 calls to genrand64_int64(). */
void memt_jump(void)
{
	struct memt_state *memt_state_init;
	int i, j;
	int bits, mask;
	
	//jump size 2^256
	char jump_string[] = 
        "514609396aa32e1815afd614eabdd3ebacb4868f08cfed4ca1"
        "e27c40ed5a24db338fe372795db756f0f632ce67327a5e61e5"
        "53e9248920f860bf759719e5db8ace1d5334763fa5df0e92dc"
        "9e78719aae25aa0e8125b0a63fc035b9605c185a4fe35cc18f"
        "98210fb398dc6cd68932a6c4dead9efd6f410086ccbe8d2518"
        "9be700ed70bd07af780e7cbda0172647d929221aec90cd5bc7"
        "1d52b673c34edf12ab6fa5b72cb466b514dec1695e3aafbc15"
        "6e1a4c7d289d7644359e108ccad0247e120f3ca0d7d5007776"
        "f2df463f383eaa3abe97e4248764e79e8219ac22b00c622376"
        "d0d17dcb3d280de5e87c3b0b826a65c36c84704026ef8351df"
        "3e7f428d113311ae397d20fe518709867ae8f076ae58cf2498"
        "945fa9fc5dfa12d6db79078d3ad42c07655feb5a7846af5d6d"
        "1422db5ae9dcc999418ac24a1f4e4c2145fa7a7c74631de210"
        "6b284f0f26377cca29a1740104cdb723b8c907d50204da74d2"
        "d3ebe9fa9eed13e21ed507151567b864798ac67aa55dec472f"
        "2907042795f242448c0b9772d51b18fd7ce9b2eed2effbd069"
        "417d2d1ea1b14d2b5753d3a11295aa1de6b0d9e7172646c86f"
        "610133672dd7a1014e2916c5dbdff1040034e07d52e90edb36"
        "59e7b612a1c45d49280072c8da8a0f5e5497d8eebde67d6a9a"
        "1b29375720036d84245ff9b96c670c555610f35ed15889e97a"
        "cc8a7812be1fd09440714e353d1c197a9addf30acc7a0bc90f"
        "b1e354125eb388";
	
	/*allocates memt_state_init*/
	memt_state_init = (struct memt_state *)malloc(sizeof(struct memt_state));
	
	/*initializes memt_state_init*/
	memt_state_init->lung = 0ULL;
	for(i = 0; i < NN; i++) memt_state_init->mt[i] = 0ULL;
	memt_state_init->mti = mti;
	memt_state_init->function_p = genrand64_int64;
	
	for (i = 0; i < ceil((double)(NN*W+P)/4); i++) {
	bits = jump_string[i];
	if (bits >= 'a' && bits <= 'f') {
	    bits = bits - 'a' + 10;
	} else {
	    bits = bits - '0';
	}
	bits = bits & 0x0f;
	mask = 0x08;
	for (j = 0; j < 4; j++) {
	    if ((bits & mask) != 0) {
			add(memt_state_init);
			}
			genrand64_int64();
			mask = mask >> 1;
		}
	}
	
	/*updates the new initial state*/
	lung = memt_state_init->lung;
	for(i = 0; i < NN; i++) mt[i] = memt_state_init->mt[i];
	mti = memt_state_init->mti;
	genrand64_int64 = memt_state_init->function_p;
	
	free(memt_state_init);
}

static void add(struct memt_state *state)
{
	int i;
	int n1, n2;
	int diff1, diff2;
	
	/*adds the lung*/
	state->lung ^= lung;
	
	n1 = state->mti;
	n2 = mti;

	/*adds the states*/
	if(n1 <= n2)
	{
		diff1 = NN - n2 + n1;
		diff2 = n2 - n1;
		
		for(i = n1; i < diff1; i++)
			state->mt[i] ^= mt[i + diff2];
		
		for(; i < NN; i++)
			state->mt[i] ^= mt[i - diff1];

		for(i = 0; i < n1; i++)
			state->mt[i] ^= mt[i + diff2];
	} else {
		diff1 = NN - n1 + n2;
		diff2 = n1 - n2;
		
		for(i = n1; i < NN; i++)
			state->mt[i] ^= mt[i - diff2];
		
		for(i = 0; i < diff2; i++)
			state->mt[i] ^= mt[i + diff1];
	
		for(; i < n1; i++)
			state->mt[i] ^= mt[i - diff2];
	}
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
    printf("\n1000 outputs of genrand64_res53()\n");
    for (i=0; i<1000; i++) {
      printf("%10.15f ", genrand64_res53());
      if (i%5==4) printf("\n");
    }
    printf("\njump ahead by 2^256 steps");
    memt_jump(); // It is equivalent to 2^256 calls to genrand64_int64()
    printf("\n1000 outputs of genrand64_int64()\n");
    for (i=0; i<1000; i++) {
      printf("%20llu ", genrand64_int64());
      if (i%5==4) printf("\n");
    }
	
    return 0;
}