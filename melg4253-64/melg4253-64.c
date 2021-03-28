/* ***************************************************************************** */
/* A C-program for MELG4253-64                                                   */
/* Copyright:   Shin Harase, Ritsumeikan University                              */
/*              Takamitsu Kimoto                                                 */
/* Notice:      This code can be used freely for personal, academic,             */
/*              or non-commercial purposes. For commercial purposes,             */
/*              please contact S. Harase at: harase @ fc.ritsumei.ac.jp          */
/* Reference:   S. Harase and T. Kimoto, "Implementing 64-bit maximally          */ 
/*              equidistributed F2-linear generators with Mersenne prime period",*/ 
/*              ACM Transactions on Mathematical Software, Volume 44, Issue 3,   */ 
/*              April 2018, Article No. 30, 11 Pages.                            */
/* Remark:      We recommend using the most significant bits (not taking the     */
/*              least significant bits) because our generators are optimized     */
/*              preferentially from the most significant bits,                   */
/*              see Remark 4.1 in the above paper for details.                   */
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

static unsigned long long melg[NN]; 
static int melgi;
static unsigned long long lung; //extra state variable
static unsigned long long mag01[2]={0ULL, MATRIX_A};
static unsigned long long x;

static unsigned long long case_1(void);
static unsigned long long case_2(void);
static unsigned long long case_3(void);
static unsigned long long case_4(void);
unsigned long long (*genrand64_int64)(void);

struct melg_state{
	unsigned long long lung;
	unsigned long long melg[NN];
	int melgi;
	unsigned long long (*function_p)(void);
};

void melg_jump(void); //jump ahead by 2^256 steps
static void add(struct melg_state *state);

/* initializes melg[NN] and lung with a seed */
void init_genrand64(unsigned long long seed)
{
    melg[0] = seed;
    for (melgi=1; melgi<NN; melgi++) {
        melg[melgi] = (6364136223846793005ULL * (melg[melgi-1] ^ (melg[melgi-1] >> 62)) + melgi);
    }
    lung = (6364136223846793005ULL * (melg[melgi-1] ^ (melg[melgi-1] >> 62)) + melgi);
    melgi = 0;
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
        melg[i] = (melg[i] ^ ((melg[i-1] ^ (melg[i-1] >> 62)) * 3935559000370003845ULL))
          + init_key[j] + j; /* non linear */
        i++; j++;
        if (i>=NN) { melg[0] = melg[NN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=NN-1; k; k--) {
        melg[i] = (melg[i] ^ ((melg[i-1] ^ (melg[i-1] >> 62)) * 2862933555777941757ULL))
          - i; /* non linear */
        i++;
        if (i>=NN) { melg[0] = melg[NN-1]; i=1; }
    }
    lung = (lung ^ ((melg[NN-1] ^ (melg[NN-1] >> 62)) * 2862933555777941757ULL))
	  - NN; /* non linear */
    melg[0] = (melg[0] | (1ULL << 63)); /* MSB is 1; assuring non-zero initial array. Corrected.  */
    melgi = 0;
}

static unsigned long long case_1(void) {
    x = (melg[melgi] & MASKU) | (melg[melgi+1] & MASKL);
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ melg[melgi+MM] ^ MAT3NEG(30, lung);
    melg[melgi] = x ^ MAT3POS(20, lung);
    x = melg[melgi] ^ (melg[melgi] << SHIFT1);
    x = x ^ (melg[melgi + LAG1] & MASK1);
    ++melgi;
    if (melgi == NN - MM) genrand64_int64 = case_2;
    return x;
}

static unsigned long long case_2(void) {
    x = (melg[melgi] & MASKU) | (melg[melgi+1] & MASKL);
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ melg[melgi+(MM-NN)] ^ MAT3NEG(30, lung);
    melg[melgi] = x ^ MAT3POS(20, lung);
    x = melg[melgi] ^ (melg[melgi] << SHIFT1);
    x = x ^ (melg[melgi + LAG1] & MASK1);
    ++melgi;
    if (melgi == LAG1over) genrand64_int64 = case_3;
    return x;
}

static unsigned long long case_3(void) {
    x = (melg[melgi] & MASKU) | (melg[melgi+1] & MASKL);
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ melg[melgi+(MM-NN)] ^ MAT3NEG(30, lung);
    melg[melgi] = x ^ MAT3POS(20, lung);
	x = melg[melgi] ^ (melg[melgi] << SHIFT1);
    x = x ^ (melg[melgi - LAG1over] & MASK1);
    ++melgi;
    if (melgi == NN-1) genrand64_int64 = case_4;
    return x;
}

static unsigned long long case_4(void) {
    x = (melg[NN-1] & MASKU) | (melg[0] & MASKL);
    lung = (x >> 1) ^ mag01[(int)(x & 1ULL)] ^ melg[MM-1] ^ MAT3NEG(30, lung);
    melg[NN-1] = x ^ MAT3POS(20, lung);
	x = melg[melgi] ^ (melg[melgi] << SHIFT1);
    x = x ^ (melg[melgi - LAG1over] & MASK1);
    melgi = 0;
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
double genrand64_fast_res52(void)
{
    union {
	unsigned long long u;
	double d;
    } conv;
	
	conv.u = (genrand64_int64() >> 12) | 0x3FF0000000000000ULL;
	return (conv.d - 1.0);
}

/* generates a random number on (0,1)-real-interval using a union trick */
double genrand64_fast_res52_open(void)
{
    union {
	unsigned long long u;
	double d;
    } conv;
	
	conv.u = (genrand64_int64() >> 12) | 0x3FF0000000000001ULL;
	return (conv.d - 1.0);
}

/* generates a random number on [0,1)-real-interval with 53-bit significant bits */
double genrand64_res53(void)
{
	return (genrand64_int64() >> 11) * 0x1.0p-53;
}

/* This is a jump function for the generator. It is equivalent
   to 2^256 calls to genrand64_int64(). */
void melg_jump(void)
{
	struct melg_state *melg_state_init;
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
	
	/*allocates melg_state_init*/
	melg_state_init = (struct melg_state *)malloc(sizeof(struct melg_state));
	
	/*initializes melg_state_init*/
	melg_state_init->lung = 0ULL;
	for(i = 0; i < NN; i++) melg_state_init->melg[i] = 0ULL;
	melg_state_init->melgi = melgi;
	melg_state_init->function_p = genrand64_int64;
	
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
			add(melg_state_init);
			}
			genrand64_int64();
			mask = mask >> 1;
		}
	}
	
	/*updates the new initial state*/
	lung = melg_state_init->lung;
	for(i = 0; i < NN; i++) melg[i] = melg_state_init->melg[i];
	melgi = melg_state_init->melgi;
	genrand64_int64 = melg_state_init->function_p;
	
	free(melg_state_init);
}

static void add(struct melg_state *state)
{
	int i;
	int n1, n2;
	int diff1, diff2;
	
	/*adds the lung*/
	state->lung ^= lung;
	
	n1 = state->melgi;
	n2 = melgi;

	/*adds the states*/
	if(n1 <= n2)
	{
		diff1 = NN - n2 + n1;
		diff2 = n2 - n1;
		
		for(i = n1; i < diff1; i++)
			state->melg[i] ^= melg[i + diff2];
		
		for(; i < NN; i++)
			state->melg[i] ^= melg[i - diff1];

		for(i = 0; i < n1; i++)
			state->melg[i] ^= melg[i + diff2];
	} else {
		diff1 = NN - n1 + n2;
		diff2 = n1 - n2;
		
		for(i = n1; i < NN; i++)
			state->melg[i] ^= melg[i - diff2];
		
		for(i = 0; i < diff2; i++)
			state->melg[i] ^= melg[i + diff1];
	
		for(; i < n1; i++)
			state->melg[i] ^= melg[i - diff2];
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
    melg_jump(); // It is equivalent to 2^256 calls to genrand64_int64()
    printf("\n1000 outputs of genrand64_int64()\n");
    for (i=0; i<1000; i++) {
      printf("%20llu ", genrand64_int64());
      if (i%5==4) printf("\n");
    }
	
    return 0;
}
