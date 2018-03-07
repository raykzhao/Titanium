/* ****************************** *
 * Titanium_CPA_super             *
 * Implemented by Raymond K. ZHAO *
 *                                *
 * Samplers                       *
 * ****************************** */
 
#include "sampler.h"
#include "param.h"
#include "fastrandombytes.h"
#include "littleendian.h"
#include "fastmodulo.h"
#include <stdint.h>
#include <x86intrin.h>

/* rejection sampling in the range of ZQ_T * Q < 2^24, 
 * then mod the sample by Q 
 * (to reduce the rejection rate) */
#define ZQ_T BARRETT_FACTOR_ZQ
#define ZQ_Q (Q * ZQ_T)

#define B1_T (NUM_B1 / (K1 + 1)) /* the number of r_i exclusively sampling from B_1 */
#define B1_REMAINING (NUM_B1 - B1_T * (K1 + 1)) /* the number of B_1 samples in r_t which has both B_1 and B_2 */

static const __m256i V_1_1_1_1 = {1, 1, 1, 1};
static const __m256i V_0_0_0_0 = {0, 0, 0, 0};
static const __m256i V_fe_fe_fe_fe = {0xfffffffffffffffeLL, 0xfffffffffffffffeLL, 0xfffffffffffffffeLL, 0xfffffffffffffffeLL};
static const __m256i V_Q_Q_Q_Q = {Q, Q, Q, Q};
static const __m256i V_0_1_2_3 = {0, 1, 2, 3};

/* Sampler of r_i
 * generate Ndec1 samples by B_1 and T*(K+1)-Ndec1 samples by B_2
 * then sample and apply the sign bit to each sample */
void sampler_zb(uint64_t sample[T][N + K1 + 1])
{
	unsigned char r[B_BYTE * T * (K1 + 1) + T * (K1 + 1) / 8];
	uint32_t i;
	uint32_t t;
	unsigned char *rr;
	__m256i y1, x1, z1;
	uint32_t y;

	fastrandombytes(r, B_BYTE * T * (K1 + 1) + T * (K1 + 1) / 8);
	
	/* sample B1_T r_i full of B_1 */
	for (t = 0; t < B1_T; t++)
	{
		rr = r + B_BYTE * t * (K1 + 1);
		for (i = 0; i < K1 + 1; i++)
		{
			y = LOAD_B(rr + i * B_BYTE);
			
			sample[t][i] = (y & B_1_BITMASK) + 1;
		}
	}
	
	/* for r_t, sample the B_1 samples first */
	rr = r + B_BYTE * B1_T * (K1 + 1);
	for (i = 0; i < B1_REMAINING; i++)
	{
		y = LOAD_B(rr + i * B_BYTE);
		
		sample[B1_T][i] = (y & B_1_BITMASK) + 1;
	}
	
	/* sample the remaining part of B_2 for r_t */
	for (i = B1_REMAINING; i < K1 + 1; i++)
	{
		y = LOAD_B(rr + i * B_BYTE);
		
		sample[B1_T][i] = (y & B_2_BITMASK) + 1;
	}
	
	/* sample those r_i full of B_2 */
	for (t = B1_T + 1; t < T; t++)
	{
		rr = r + B_BYTE * t * (K1 + 1);
		for (i = 0; i < K1 + 1; i++)
		{
			y = LOAD_B(rr + i * B_BYTE);
			
			sample[t][i] = (y & B_2_BITMASK) + 1;
		}
	}

	/* apply the sign bits in constant time
	 * each random byte will generate 8 sign bits
	 * each sample will become [-B,-1] U [1,B] after applying the sign bit */
	for (t = 0; t < T; t++)
	{
		rr = r + (B_BYTE * T * (K1 + 1)) + (t * (K1 + 1) / 8);
		for (i = 0; i < K1 + 1; i += 8)
		{
			x1 = _mm256_set1_epi64x(rr[i / 8]);
			x1 = _mm256_srlv_epi64(x1, V_0_1_2_3);
			z1 = _mm256_and_si256(x1, V_1_1_1_1);
			z1 = _mm256_sub_epi64(V_0_0_0_0, z1);
			z1 = _mm256_and_si256(z1, V_fe_fe_fe_fe);
			z1 = _mm256_xor_si256(z1, V_1_1_1_1);
			y1 = _mm256_loadu_si256((__m256i *)(sample[t] + i));
			y1 = _mm256_mul_epi32(y1, z1);
			y1 = _mm256_add_epi64(V_Q_Q_Q_Q, y1);
			_mm256_storeu_si256((__m256i *)(sample[t] + i), y1);
			
			x1 = _mm256_srli_epi64(x1, 4);
			z1 = _mm256_and_si256(x1, V_1_1_1_1);
			z1 = _mm256_sub_epi64(V_0_0_0_0, z1);
			z1 = _mm256_and_si256(z1, V_fe_fe_fe_fe);
			z1 = _mm256_xor_si256(z1, V_1_1_1_1);
			y1 = _mm256_loadu_si256((__m256i *)(sample[t] + i + 4));
			y1 = _mm256_mul_epi32(y1, z1);
			y1 = _mm256_add_epi64(V_Q_Q_Q_Q, y1);
			_mm256_storeu_si256((__m256i *)(sample[t] + i + 4), y1);
		}	
	}
}

/* return (x >= y) */
static inline int ct_ge_u32(uint32_t x, uint32_t y)
{
    return 1 ^ ((x - y) >> 31);
}

/* Sample slen uniforms in Z_q 
 * 
 * Generate a little bit more random bytes, since re-run the PRG if running out of the randomness would be costly (also for the constant time implementation)
 * 
 * also, to reduce the rejection rate, do the rejection in the range of ZQ_T * Q, which is "just" smaller than 2^24, 
 * and mod each sample by Q */
void sampler_zq(uint64_t *sample, uint32_t slen, uint32_t bytpc)
{
	unsigned char r[ZQ_BYTES * ZQ_BYTPCS];
	uint32_t i = 0, j = 0;
	uint32_t x;
	
	fastrandombytes(r, ZQ_BYTES * bytpc);
	
	while (j < slen)
	{
		do
		{
			if (i == bytpc)
			{
				fastrandombytes(r, ZQ_BYTES * bytpc);
				i = 0;
			}
			
			x = LOAD_ZQ(r + ZQ_BYTES * (i++));
		} while (ct_ge_u32(x, ZQ_Q)); /* rejection in the range of ZQ_T * Q */
		
		sample[j++] = barrett_short(x); /* barrett reduction here */
	}
}

/* Binomial Sampler */
void sampler_binomial(uint64_t sample[T][D + K1 + 1])
{
	unsigned char r[T * (D + K1 + 1) * BINOMIAL_BYTE];
	uint32_t i, j;
	uint32_t x, s;
	uint32_t t;
	unsigned char *rr;
	
	fastrandombytes(r, T * (D + K1 + 1) * BINOMIAL_BYTE);
	
	for (t = 0; t < T; t++)
	{
		rr = r + t * (D + K1 + 1) * BINOMIAL_BYTE;
		
		/* x = (r1,r2)
		 * s1 = bit_sum(r1 & bitmask), s2 = bit_sum(r2 & bitmask), where bit_sum is the number of "1" bits, and bitmask is the mask to select k bits
		 * sample <-- s1 - s2 */
		for (i = 0; i < D + K1; i++)
		{
			x = LOAD_BINOMIAL(rr + i * BINOMIAL_BYTE);
		
			s = 0;
			for (j = 0; j < BINOMIAL_K; j++)
			{
				s += (x >> j) & BINOMIAL_ADDMASK;
			}
		
			sample[t][i] = Q + (s & BINOMIAL_MASK) - (s >> BINOMIAL_SHIFT);
		}
	}
}
