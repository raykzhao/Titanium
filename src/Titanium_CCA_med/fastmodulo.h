/* ****************************** *
 * Titanium_CCA_med               *
 * Implemented by Raymond K. ZHAO *
 *                                *
 * Modulo reduction               *
 * ****************************** */
 
#ifndef FASTMODULO_H
#define FASTMODULO_H

#include "param.h"
#include <stdint.h>

/* Montgomery reduction
 * Input: x < Q*R, where R=2^k and Q<R
 * Output: m = x*R^{-1} % Q
 * 
 * b = -Q^{-1} % R
 * t = ((x % R)*b) % R
 * m = (x + t * Q) / R */

#define MONTGOMERY_FACTOR 430079
#define MONTGOMERY_SHIFT 23
#define MONTGOMERY_MASK ((1 << MONTGOMERY_SHIFT) - 1)

inline uint32_t montgomery(uint64_t t)
{
	return (t + ((((t & MONTGOMERY_MASK) * MONTGOMERY_FACTOR) & MONTGOMERY_MASK) * Q)) >> MONTGOMERY_SHIFT;
}

/* Input: x < 2^k
 * Output m = x % Q in [0, 2Q)
 * 
 * b = floor(2^k/Q)
 * t = floor((x * b) / 2^k), where t is an estimation of x / Q
 * m = x - t * Q */

#define BARRETT_BITSHIFT_4Q 21 
#define BARRETT_BITSHIFT_8Q 22 
#define BARRETT_BITSHIFT_16Q 23 

#define BARRETT_BITSHIFT_2Q2 39
#define BARRETT_BITSHIFT_4Q2 40 
#define BARRETT_BITSHIFT_8Q2 41 

#define BARRETT_BITSHIFT_ZQ (ZQ_BYTES * 8) 

#define BARRETT_FACTOR_4Q 4
#define BARRETT_FACTOR_8Q 9
#define BARRETT_FACTOR_16Q 19

#define BARRETT_FACTOR_2Q2 1278261
#define BARRETT_FACTOR_4Q2 2556522
#define BARRETT_FACTOR_8Q2 5113044

#define BARRETT_FACTOR_ZQ 39

inline uint32_t barrett_4q(uint32_t t)
{
	return t - (((t * BARRETT_FACTOR_4Q) >> BARRETT_BITSHIFT_4Q) * Q);
}

inline uint32_t barrett_8q(uint32_t t)
{
	return t - (((t * BARRETT_FACTOR_8Q) >> BARRETT_BITSHIFT_8Q) * Q);
}

inline uint32_t barrett_16q(uint32_t t)
{
	return t - (((t * BARRETT_FACTOR_16Q) >> BARRETT_BITSHIFT_16Q) * Q);
}

inline uint32_t barrett_2q2(uint64_t t)
{
	return t - (((t * BARRETT_FACTOR_2Q2) >> BARRETT_BITSHIFT_2Q2) * Q);
}

inline uint32_t barrett_4q2(uint64_t t)
{
	return t - (((t * BARRETT_FACTOR_4Q2) >> BARRETT_BITSHIFT_4Q2) * Q);
}

inline uint32_t barrett_8q2(uint64_t t)
{
	return t - (((t * BARRETT_FACTOR_8Q2) >> BARRETT_BITSHIFT_8Q2) * Q);
}

inline uint32_t barrett_zq(uint32_t t)
{
	return t - (((t * BARRETT_FACTOR_ZQ) >> BARRETT_BITSHIFT_ZQ) * Q);
}

#endif
