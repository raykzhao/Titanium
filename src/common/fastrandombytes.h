/* ****************************** *
 * Titanium                       *
 * Implemented by Raymond K. ZHAO *
 *                                *
 * Wrapper of AES256 PRG          *
 * ****************************** */
 
#ifndef FASTRANDOMBYTES_H
#define FASTRANDOMBYTES_H

void fastrandombytes(unsigned char *r, unsigned long long rlen);
void fastrandombytes_setseed(const unsigned char *randomness);

#endif
