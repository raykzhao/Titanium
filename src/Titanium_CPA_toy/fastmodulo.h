/* ****************************** *
 * Titanium_CPA_toy               *
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

#define MONTGOMERY_FACTOR 1148431359
#define MONTGOMERY_SHIFT 32

/* Input: x < 2^k
 * Output m = x % Q in [0, 2Q)
 * 
 * b = floor(2^k/Q)
 * t = floor((x * b) / 2^k), where t is an estimation of x / Q
 * m = x - t * Q */

#define BARRETT_BITSHIFT_SHORT 32
#define BARRETT_FACTOR_SHORT 17848

#define BARRETT_BITSHIFT_4Q2 38 
#define BARRETT_FACTOR_4Q2 1142273

#define BARRETT_BITSHIFT_32Q2 41 
#define BARRETT_FACTOR_32Q2 9138190

#define BARRETT_BITSHIFT_ZQ (ZQ_BYTES * 8) 
#define BARRETT_FACTOR_ZQ 69

inline uint32_t barrett_short(uint64_t t)
{
	return t - (((t * BARRETT_FACTOR_SHORT) >> BARRETT_BITSHIFT_SHORT) * Q);
}

inline uint32_t barrett_4q2(uint64_t t)
{
	return t - (((t * BARRETT_FACTOR_4Q2) >> BARRETT_BITSHIFT_4Q2) * Q);
}

inline uint32_t barrett_32q2(uint64_t t)
{
	return t - (((t * BARRETT_FACTOR_32Q2) >> BARRETT_BITSHIFT_32Q2) * Q);
}

#endif
