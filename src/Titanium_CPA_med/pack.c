/* ****************************** *
 * Titanium_CPA_med               *
 * Implemented by Raymond K. ZHAO *
 *                                *
 * Packing/Unpacking              *
 * ****************************** */

#include "pack.h"
#include "param.h"
#include "littleendian.h"
#include <stdint.h>
#include <x86intrin.h>

static const __m256i V_Q_Q_Q_Q = {Q, Q, Q, Q};
static const __m256i V_1_1_1_1 = {1, 1, 1, 1};

/* convert a polynomial to a binary string */
void poly_encode(unsigned char *b, const uint64_t *p, uint32_t len)
{
	uint32_t i;
	unsigned char *bb;
	uint64_t pp[8] __attribute__ ((aligned (32)));
	__m256i u, t;
	
	/* pack 8 19-bit coordinates to 19 bytes */
	for (i = 0; i < len; i += 8)
	{
		/* make sure each coordinate is smaller than Q */
		u = _mm256_load_si256((__m256i *)(p + i));
		t = _mm256_sub_epi64(u, V_Q_Q_Q_Q);
		t = _mm256_srli_epi64(t, 63);
		t = _mm256_xor_si256(t, V_1_1_1_1);
		t = _mm256_mul_epu32(t, V_Q_Q_Q_Q);
		t = _mm256_sub_epi64(u, t);
		_mm256_store_si256((__m256i *)pp, t);

		u = _mm256_load_si256((__m256i *)(p + i + 4));
		t = _mm256_sub_epi64(u, V_Q_Q_Q_Q);
		t = _mm256_srli_epi64(t, 63);
		t = _mm256_xor_si256(t, V_1_1_1_1);
		t = _mm256_mul_epu32(t, V_Q_Q_Q_Q);
		t = _mm256_sub_epi64(u, t);
		_mm256_store_si256((__m256i *)(pp + 4), t);
		
		bb = b + (i / 8) * Q_BITS;
		bb[0] = pp[0];
		bb[1] = pp[0] >> 8;
		bb[2] = (pp[0] >> 16) | ((pp[1] & 0x1f) << 3);
		bb[3] = pp[1] >> 5;
		bb[4] = (pp[1] >> 13) | ((pp[2] & 0x03) << 6);
		bb[5] = pp[2] >> 2;
		bb[6] = pp[2] >> 10;
		bb[7] = (pp[2] >> 18) | ((pp[3] & 0x7f) << 1);
		bb[8] = pp[3] >> 7;
		bb[9] = (pp[3] >> 15) | ((pp[4] & 0x0f) << 4);
		bb[10] = pp[4] >> 4;
		bb[11] = (pp[4] >> 12) | ((pp[5] & 0x01) << 7);
		bb[12] = pp[5] >> 1;
		bb[13] = pp[5] >> 9;
		bb[14] = (pp[5] >> 17) | ((pp[6] & 0x3f) << 2);
		bb[15] = pp[6] >> 6;
		bb[16] = (pp[6] >> 14) | ((pp[7] & 0x07) << 5);
		bb[17] = pp[7] >> 3;
		bb[18] = pp[7] >> 11;
	}
}

/* convert a binary string to a polynomial */
void poly_decode(uint64_t *p, const unsigned char *b, uint32_t len)
{
	uint32_t i;
	unsigned char *bb;
	
	/* unpack 19 bytes to 8 19-bit coordinates */
	for (i = 0; i < len; i += 8)
	{
		bb = b + (i / 8) * Q_BITS;

		p[i] = ((uint64_t)bb[0]) | (((uint64_t)bb[1]) << 8) | ((((uint64_t)bb[2]) & 0x07) << 16);
		p[i + 1] = (((uint64_t)bb[2]) >> 3) | (((uint64_t)bb[3]) << 5) | ((((uint64_t)bb[4]) & 0x3f) << 13);
		p[i + 2] = (((uint64_t)bb[4]) >> 6) | (((uint64_t)bb[5]) << 2) | (((uint64_t)bb[6]) << 10) | ((((uint64_t)bb[7]) & 0x01) << 18);
		p[i + 3] = (((uint64_t)bb[7]) >> 1) | (((uint64_t)bb[8]) << 7) | ((((uint64_t)bb[9]) & 0x0f) << 15);
		p[i + 4] = (((uint64_t)bb[9]) >> 4) | (((uint64_t)bb[10]) << 4) | ((((uint64_t)bb[11]) & 0x7f) << 12);
		p[i + 5] = (((uint64_t)bb[11]) >> 7) | (((uint64_t)bb[12]) << 1) | (((uint64_t)bb[13]) << 9) | ((((uint64_t)bb[14]) & 0x03) << 17);
		p[i + 6] = (((uint64_t)bb[14]) >> 2) | (((uint64_t)bb[15]) << 6) | ((((uint64_t)bb[16]) & 0x1f) << 14);
		p[i + 7] = (((uint64_t)bb[16]) >> 5) | (((uint64_t)bb[17]) << 3) | (((uint64_t)bb[18]) << 11);
	}
}

/* convert a polynomial to a binary string (with compression) */
void poly_encode_c2(unsigned char *b, const uint64_t *p, uint32_t len)
{
	uint32_t i;
	uint64_t pp[4] __attribute__ ((aligned (32)));
	__m256i u, t;
	
	for (i = 0; i < len; i += 4)
	{
		u = _mm256_load_si256((__m256i *)(p + i));
		t = _mm256_sub_epi64(u, V_Q_Q_Q_Q);
		t = _mm256_srli_epi64(t, 63);
		t = _mm256_xor_si256(t, V_1_1_1_1);
		t = _mm256_mul_epu32(t, V_Q_Q_Q_Q);
		t = _mm256_sub_epi64(u, t);
		t = _mm256_srli_epi64(t, C2_COMPRESSION_BITS);
		_mm256_store_si256((__m256i *)pp, t);
		
		STORE_C2(b + i * C2_COMPRESSION_BYTE, pp[0]);
		STORE_C2(b + (i + 1) * C2_COMPRESSION_BYTE, pp[1]);
		STORE_C2(b + (i + 2) * C2_COMPRESSION_BYTE, pp[2]);
		STORE_C2(b + (i + 3) * C2_COMPRESSION_BYTE, pp[3]);
	}
}

/* convert a binary string to a polynomial (with compression) */
void poly_decode_c2(uint64_t *p, const unsigned char *b, uint32_t len)
{
	uint32_t i;
	__m256i t;
	
	for (i = 0; i < len; i += 4)
	{
		t = _mm256_set_epi64x(LOAD_C2(b + (i + 3) * C2_COMPRESSION_BYTE), LOAD_C2(b + (i + 2) * C2_COMPRESSION_BYTE), LOAD_C2(b + (i + 1) * C2_COMPRESSION_BYTE), LOAD_C2(b + i * C2_COMPRESSION_BYTE));
		t = _mm256_slli_epi64(t, C2_COMPRESSION_BITS);
		_mm256_store_si256((__m256i *)(p + i), t);
	}
}

