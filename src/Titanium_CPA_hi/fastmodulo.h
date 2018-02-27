/* ****************************** *
 * Titanium_CPA_hi                *
 * Implemented by Raymond K. ZHAO *
 *                                *
 * Modulo reduction               *
 * ****************************** */
 
#ifndef FASTMODULO_H
#define FASTMODULO_H

#include "param.h"
#include <stdint.h>

#define MONTGOMERY_FACTOR 737279
#define MONTGOMERY_SHIFT 22
#define MONTGOMERY_MASK ((1 << MONTGOMERY_SHIFT) - 1)

/* Input: x < 2^k
 * Output m = x % Q in [0, 2Q)
 * 
 * b = floor(2^k/Q)
 * t = floor((x * b) / 2^k), where t is an estimation of x / Q
 * m = x - t * Q */

#define BARRETT_BITSHIFT_4Q 22 
#define BARRETT_BITSHIFT_8Q 23 
#define BARRETT_BITSHIFT_16Q 24 
#define BARRETT_BITSHIFT_32Q 25 

#define BARRETT_BITSHIFT_ZQ (ZQ_BYTES * 8) 

#define BARRETT_FACTOR_4Q 5
#define BARRETT_FACTOR_8Q 11
#define BARRETT_FACTOR_16Q 22
#define BARRETT_FACTOR_32Q 45

#define BARRETT_FACTOR_ZQ 22

inline uint32_t barrett_zq(uint32_t t)
{
	return t - (((t * BARRETT_FACTOR_ZQ) >> BARRETT_BITSHIFT_ZQ) * Q);
}

#endif
