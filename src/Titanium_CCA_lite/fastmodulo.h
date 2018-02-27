/* ****************************** *
 * Titanium_CCA_lite              *
 * Implemented by Raymond K. ZHAO *
 *                                *
 * Modulo reduction               *
 * ****************************** */
 
#ifndef FASTMODULO_H
#define FASTMODULO_H

#include "param.h"
#include <stdint.h>

#define MONTGOMERY_FACTOR 377343
#define MONTGOMERY_SHIFT 19
#define MONTGOMERY_MASK ((1 << MONTGOMERY_SHIFT) - 1)

/* Input: x < 2^k
 * Output m = x % Q in [0, 2Q)
 * 
 * b = floor(2^k/Q)
 * t = floor((x * b) / 2^k), where t is an estimation of x / Q
 * m = x - t * Q */

#define BARRETT_BITSHIFT_4Q 19 
#define BARRETT_BITSHIFT_8Q 20 
#define BARRETT_BITSHIFT_16Q 21 

#define BARRETT_BITSHIFT_ZQ (ZQ_BYTES * 8) 

#define BARRETT_FACTOR_4Q 4
#define BARRETT_FACTOR_8Q 9
#define BARRETT_FACTOR_16Q 18

#define BARRETT_FACTOR_ZQ 145

inline uint32_t barrett_zq(uint32_t t)
{
	return t - (((t * BARRETT_FACTOR_ZQ) >> BARRETT_BITSHIFT_ZQ) * Q);
}

#endif
