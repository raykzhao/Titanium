/* ****************************** *
 * Titanium                       *
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

inline uint32_t montgomery(uint64_t t)
{
	uint32_t x, y;
	
	x = t;
	y = ((uint64_t)x) * MONTGOMERY_FACTOR;
	
	return (t + ((uint64_t)y) * Q) >> MONTGOMERY_SHIFT;
}

/* Barrett reduction
 * Input: x < 2^k
 * Output m = x % Q in [0, 2Q)
 * 
 * b = floor(2^k/Q)
 * t = floor((x * b) / 2^k), where t is an estimation of x / Q
 * m = x - t * Q */

inline uint32_t barrett(uint64_t t, const uint64_t barrett_factor, const uint64_t barrett_bitshift)
{
	return t - (((t * barrett_factor) >> barrett_bitshift) * Q);
}

#endif
