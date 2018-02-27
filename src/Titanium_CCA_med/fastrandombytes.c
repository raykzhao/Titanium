/* AES-256 CTR PRNG by using the OpenSSL hardware accelerated subroutine */

#include "fastrandombytes.h"
#include "littleendian.h"
#include <string.h>
#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>

static unsigned char iv[16];
static EVP_CIPHER_CTX *ctx;

/* round_key <-- aes256_key_expansion(randomness), iv <-- 0 */
void fastrandombytes_setseed(const unsigned char *randomness)
{
	ctx = EVP_CIPHER_CTX_new();
	EVP_EncryptInit_ex(ctx, EVP_aes_256_ecb(), NULL, randomness, NULL);

	memset(iv, 0, 16);
}

/* r <-- aes256_ctr(round_key, iv, rlen) */
void fastrandombytes(unsigned char *r, unsigned long long rlen)
{
	unsigned char ct[16];
	unsigned long long num_of_blocks = rlen >> 4;
	unsigned long long i;
	int len;
	
	for (i = 0; i < num_of_blocks; i++)
	{
		EVP_EncryptUpdate(ctx, r + (i << 4), &len, iv, 16);
		
		store_32(iv, load_32(iv) + 1);
	}
	
	if (rlen & 0x0f)
	{
		EVP_EncryptUpdate(ctx, ct, &len, iv, 16);
		store_32(iv, load_32(iv) + 1);

		memcpy(r + (i << 4), ct, rlen & 0x0f);
	}	
}
