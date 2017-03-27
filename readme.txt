Implementing 64-bit Maximally Equidistributed F2-Linear Generators with Mersenne Prime Period

S. Harase and T. Kimoto, "Implementing 64-bit Maximally Equidistributed F2-Linear Generators with Mersenne Prime Period", submitted. http://arxiv.org/abs/1505.06582

Abstract:
CPUs and operating systems are moving from 32 to 64 bits, and hence it is important to have good pseudorandom number generators designed to fully exploit these word lengths. However, existing 64-bit very long period generators based on linear recurrences modulo 2 are not completely optimized in terms of the equidistribution properties. Here we develop 64-bit maximally equidistributed pseudorandom number generators that are optimal in this respect and have speeds equivalent to 64-bit Mersenne Twisters. We provide a table of specific parameters with period lengths from $2^{607}-1$ to $2^{44497}-1$.

(27 March 2017)
The paper was revised.
We renamed the proposed generators.
We attached the acronym 64-bit MELGs, 
for 64-bit ``Maximally Equidistributed $\mathbb{F}_2$-Linear Generators", 
to the proposed generators.

We modified a initialization schime in init_by_array64.
Thus, the output values have been chenaged.

(12 May 2016)
Minor changes: unsigned long long int -> unsigned long long
               void memt_jump(void); //jump ahead by 2^256 steps

(10 May 2016)
We implemented the following functions:
  genrand64_res53()
  genrand64_res53_open()
to generate double-precision floating point numbers on [0,1)- and (0,1)-real-interval by using union tricks.

We implemented a jump ahead algorithm for parallel computing. 
The default skip size is 2^256.

We also changed some examples in the main functions.

(24 May 2015)
The code in C was released.
