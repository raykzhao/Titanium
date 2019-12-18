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
