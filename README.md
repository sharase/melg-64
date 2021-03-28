# Implementing 64-bit Maximally Equidistributed  **F**<sub>2</sub>-Linear Generators with Mersenne Prime Period

## What is MELG-64?

The 64-bit **M**aximally **E**quidistributed **F**<sub>2</sub>-**L**inear **G**enerators with Mersenne Prime Period (**MELG-64**) are 64-bit Mersenne-Twister-type pseudorandom number generators developed between 2014 and 2017, and the corresponding paper was published on ACM TOMS in 2018.

> S. Harase and T. Kimoto, "Implementing 64-bit maximally equidistributed **F**<sub>2</sub>-linear generators with Mersenne prime period", ACM Transactions on Mathematical Software, Volume 44, Issue 3, April 2018, Article No. 30, 11 pp. <a href="http://doi.acm.org/10.1145/3159444">Artcle</a>

## Background
CPUs and operating systems are moving from 32 to 64 bits, and hence it is important to have good 64-bit pseudorandom number generators (PRNGs) designed to fully exploit these word lengths. 

The 32-bit Mersenne Twister (MT) MT19937 (Matsumoto--Nishimura, 1998) is one of the most widely used PRNGs, but it is _not completely optimized_ in terms of _high-dimensional uniformity_, which is a theoretical criterion of PRNGs. The 32-bit WELL generators (Panneton et. al., 2006) was developed in order to overcome this weakness. 

However, for 64-bit PRNGs, MT19937-64 (Nishimura, 2000) and SFMT19937 using SIMD (Saito and Matsumoto, 2008), etc., have been proposed, but there exists no 64-bit MT-type long-period linear PRNG completely optimized for high-dimensional uniformity, such as a variant of WELL generators.

In this page, we introduce 64-bit maximally equidistributed **F**<sub>2</sub>-linear PRNGs (**MELG-64**) that are optimal in this respect and have speeds equivalent to 64-bit Mersenne Twisters. 

## Feature

MELG19937-64 has the following properties:
-  Very long period 2<sup>19937</sup>-1 <span></span> <span>&asymp;</span> 10<sup>6000</sup>;
- High-dimensional uniformity completely optimized;
- Fast generation competitive with MT19937-64;
- Memory size requiring only 312 words (similarly to MT19937-64).

We provide the codes for **MELG-64** with various period lengths from 2<sup>607</sup>-1 to 2<sup>44497</sup>-1. Also, we have implemented the jump-ahead algorithm in order to obtain disjoint streams in parallel computing.

## Usage

Please click the "Code" button in the upper right of the content pane, and clone this page: 

> clone https://github.com/sharase/melg-64.git

or download the zip file. If you want to use **MELG19937-64**, 
> cd melg-64 </br>
> cd melg19937-64 </br>
> gcc melg19937-64.c -o melg19937-64 -O3 -Wall </br>
> ./melg19937-64

Before using, please initialize the state by the function init_by_array64(init, length) in a similar way as Mersenne Twisters.

## High-dimensional uniformity
The high-dimensional uniformity is a theoretical criterion for PRNGs, which is assessed via the _dimension of equidistribution with v-bit accuracy_ as follows:

Let

> **x**<sub>0</sub>, **x**<sub>1</sub>, ..., **x**<sub>P-1</sub>, **x**<sub>P</sub> = **x**<sub>0</sub>, ...

be a _unsigned w-bit binary integer sequence_ with period P, 
where w is the word size of the intended machine. Let trunc<sub>v</sub>(**x**<sub>i</sub>)  
denote the number formed by the v most significant bits of **x**<sub>i</sub>. Consider the kv-bit vectors for the entire period:

> (trunc<sub>v</sub>(**x**<sub>i</sub>), trunc<sub>v</sub>(**x**<sub>i+1</sub>), ..., trunc<sub>v</sub>(**x**<sub>i+k-1</sub>)), &emsp; i = 0, ..., P-1 .

A pseudorandom sequence **x**<sub>i</sub> of w-bit integers of period P is said to be _k-dimensionally equidistributed with v-bit accuracy_ if each of the 2<sup>kv</sup> possible combinations of bits occurs the same number of times over the whole period P, except for the all-zero combination that occurs once less often. 

The largest value of k with this property is called the _dimension of equidistribution 
with v-bit accuracy_, denoted by k(v).

This definition is based on the assumption that the higher digits are large numbers. In particular, the dimension of equidistribution assures that the output values with the v most significant bits are uniformly distributed up to dimension k(v). Thus, as a criterion of uniformity, larger values of k(v) for each 1
<span>&le;</span> v <span>&le;</span>  w is desirable.

Now we have a trivial upper bound 

> k(v) <span>&le;</span> <span>&lfloor;</span> log<sub>2</sub> (P+1) /v <span>&rfloor;</span>

for each v = 1,2, ...,w. Define the sum of the gaps 

> <span>&Delta;</span> := <span>&sum;</span>&ensp;<span>&lfloor;</span> \log <sub>2</sub> (P+1)/v <span>&rfloor;</span> - k(v)), 

where the sum is over all 1 <span>&le;</span> v <span>&le;</span> w.

If <span>&Delta;</span> = 0, the generator is said to be _maximally equidistributed (ME)_.  

The aim of our study is to design maximally equidistributed **F**<sub>2</sub>-linear PRNGs with similar speed as 64-bit Mersenne Twisters.

## Performance

We compare the following MT-type PRNGs corresponding to 64-bit output sequences:
- MELG19937-64: the 64-bit integer output of our proposed generator;
- MT19937-64: the 64-bit integer output of the 64-bit Mersenne Twister (downloaded from
http://www.math.sci.hiroshima-u.ac.jp/∼m-mat/MT/emt64.html);
- MT19937-64 (ID3): the 64-bit integer output of a 64-bit Mersenne Twister based on a five-term
recursion (ID3) (Nishimura 2000);
- SFMT19937-64 (without SIMD): the 64-bit integer output of the SIMD-oriented Fast
Mersenne Twister SFMT19937 without SIMD (Saito and Matsumoto 2008);
- SFMT19937-64 (with SIMD): the 64-bit integer output of the foregoing with SIMD (Saito
and Matsumoto 2008).

We measure the CPU time (in seconds) taken to generate 10<sup>9</sup> 64-bit unsigned integers. 
The following table summarizes the timing and the figures of merit <span>&Delta;</span> and _N_<sub>1</sub>. 
See the remark below for the definition of _N_<sub>1</sub>.

| Generators | CPU time (Intel) | CPU time (AMD) | <span>&Delta;</span> | _N_<sub>1</sub> |
| :--------- | ---------------: | -------------: | -------: | ----: |
| **MELG19937-64** | **4.2123** | **6.2920** | **0** | **9603** |
| MT19937-64 | 5.1002 | 6.6490 | 7820 | 285 |
| MT19937-64 (ID3) | 4.8993 | 6.7930 | 7940 | 5795 |
| SFMT19937-64 (without SIMD) | 4.2654 | 5.6123 | **14095** | 6711 |
| SFMT19937-64 (with SIMD) | 1.8457 | 2.8806 | **14095** | 6711 |

</br>
Platforms (64-bit CPUs and OSs):

- CPU time (Intel): Intel Core i7-3770 (3.40GHz) Linux gcc compiler with -O3
- CPU time (AMD): AMD Phenom II X6 1045T (2.70 GHz) Linux gcc compiler with -O3

**SFMT19937 is very fast but <span>&Delta;</span> for SFMT19937 is large.** 
(In fact, the SFMT
generators are optimized under the assumption that one will mainly be using 32-bit output sequences. For double-precision floating-point numbers, **dSFMT** is faster than SFMT and is also improved
from the viewpoint of the dimensions of equidistribution with v-bit accuracy.)

Remark:
>For MT-type PRNGs, each bit of the sequence obeys a linear feedback shift register generator with the characteristic polynomial P(x) over the two-element field **F**<sub>2</sub>.
Let _N_<sub>1</sub> be the number of nonzero coefficients of P(x). 
In addition to the excellent equidistribution, as a secondary criterion, _N_<sub>1</sub> should be large enough. 
This criterion implies that the PRNG avoids a long-lasting impact for poor initialization, such as 0-excess states (Panneton et al, 2006).

## Further readings

- S. Harase and T. Kimoto, "Implementing 64-bit maximally equidistributed **F**<sub>2</sub>-linear generators with Mersenne prime period", ACM Transactions on Mathematical Software, Volume 44, Issue 3, April 2018, Article No. 30, 11 Pages. <a href="http://doi.acm.org/10.1145/3159444">Artcle</a>
- S. Harase, "Conversion of Mersenne Twister to double-precision floating-point numbers", Mathematics and Computers in Simulation, Volume 161, July 2019, Pages 76-83. <a href="https://doi.org/10.1016/j.matcom.2018.08.006"> Article</a>

## Acknowledgments
This work was partially supported by his work was supported by JSPS KAKENHI Grant Numbers JP18K18016, JP26730015, JP26310211, JP15K13460, JP12J07985.
