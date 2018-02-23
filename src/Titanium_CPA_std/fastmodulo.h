/* ****************************** *
 * Titanium_CPA_std               *
 * Implemented by Raymond K. ZHAO *
 *                                *
 * Modulo reduction               *
 * ****************************** */
 
#ifndef FASTMODULO_H
#define FASTMODULO_H

#include "param.h"
#include <stdint.h>

/* Input: x < 2^k
 * Output m = x % Q in [0, 2Q)
 * 
 * b = floor(2^k/Q)
 * t = floor((x * b) / 2^k), where t is an estimation of x / Q
 * m = x - t * Q */

#define BARRETT_BITSHIFT_4Q 19 
#define BARRETT_BITSHIFT_8Q 20 
#define BARRETT_BITSHIFT_16Q 21 

#define BARRETT_BITSHIFT_2Q2 35 
#define BARRETT_BITSHIFT_4Q2 36 
#define BARRETT_BITSHIFT_8Q2 37 
#define BARRETT_BITSHIFT_16Q2 38 

#define BARRETT_BITSHIFT_ZQ (ZQ_BYTES * 8) 

#define BARRETT_FACTOR_4Q 6
#define BARRETT_FACTOR_8Q 12
#define BARRETT_FACTOR_16Q 24

#define BARRETT_FACTOR_2Q2 399452
#define BARRETT_FACTOR_4Q2 798905
#define BARRETT_FACTOR_8Q2 1597811
#define BARRETT_FACTOR_16Q2 3195623

#define BARRETT_FACTOR_ZQ 195

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

inline uint32_t barrett_16q2(uint64_t t)
{
	return t - (((t * BARRETT_FACTOR_16Q2) >> BARRETT_BITSHIFT_16Q2) * Q);
}

inline uint32_t barrett_zq(uint32_t t)
{
	return t - (((t * BARRETT_FACTOR_ZQ) >> BARRETT_BITSHIFT_ZQ) * Q);
}

#endif
