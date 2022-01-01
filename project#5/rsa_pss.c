/*
 * Copyright 2020,2021. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었습니다.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "rsa_pss.h"

#if defined(SHA224)
void (*sha)(const unsigned char *, unsigned int, unsigned char *) = sha224;
#elif defined(SHA256)
void (*sha)(const unsigned char *, unsigned int, unsigned char *) = sha256;
#elif defined(SHA384)
void (*sha)(const unsigned char *, unsigned int, unsigned char *) = sha384;
#else
void (*sha)(const unsigned char *, unsigned int, unsigned char *) = sha512;
#endif

/*
 * Copyright 2020, 2021. Heekuck Oh, all rights reserved
 * rsa_generate_key() - generates RSA keys e, d and n in octet strings.
 * If mode = 0, then e = 65537 is used. Otherwise e will be randomly selected.
 * Carmichael's totient function Lambda(n) is used.
 */
void rsa_generate_key(void *_e, void *_d, void *_n, int mode)
{
    mpz_t p, q, lambda, e, d, n, gcd;
    gmp_randstate_t state;

    /*
     * Initialize mpz variables
     */
    mpz_inits(p, q, lambda, e, d, n, gcd, NULL);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, arc4random());
    /*
     * Generate prime p and q such that 2^(RSAKEYSIZE-1) <= p*q < 2^RSAKEYSIZE
     */
    do {
        do {
            mpz_urandomb(p, state, RSAKEYSIZE/2);
            mpz_setbit(p, 0);
            mpz_setbit(p, RSAKEYSIZE/2-1);
       } while (mpz_probab_prime_p(p, 50) == 0);
        do {
            mpz_urandomb(q, state, RSAKEYSIZE/2);
            mpz_setbit(q, 0);
            mpz_setbit(q, RSAKEYSIZE/2-1);
        } while (mpz_probab_prime_p(q, 50) == 0);
        mpz_mul(n, p, q);
    } while (!mpz_tstbit(n, RSAKEYSIZE-1));
    /*
     * Generate e and d using Lambda(n)
     */
    mpz_sub_ui(p, p, 1);
    mpz_sub_ui(q, q, 1);
    mpz_lcm(lambda, p, q);
    if (mode == 0)
        mpz_set_ui(e, 65537);
    else do {
        mpz_urandomb(e, state, RSAKEYSIZE);
        mpz_gcd(gcd, e, lambda);
    } while (mpz_cmp(e, lambda) >= 0 || mpz_cmp_ui(gcd, 1) != 0);
    mpz_invert(d, e, lambda);
    /*
     * Convert mpz_t values into octet strings
     */
    mpz_export(_e, NULL, 1, RSAKEYSIZE/8, 1, 0, e);
    mpz_export(_d, NULL, 1, RSAKEYSIZE/8, 1, 0, d);
    mpz_export(_n, NULL, 1, RSAKEYSIZE/8, 1, 0, n);
    /*
     * Free the space occupied by mpz variables
     */
    mpz_clears(p, q, lambda, e, d, n, gcd, NULL);
}

/*
 * Copyright 2020. Heekuck Oh, all rights reserved
 * rsa_cipher() - compute m^k mod n
 * If m >= n then returns EM_MSG_OUT_OF_RANGE, otherwise returns 0 for success.
 */
static int rsa_cipher(void *_m, const void *_k, const void *_n)
{
    mpz_t m, k, n;

    /*
     * Initialize mpz variables
     */
    mpz_inits(m, k, n, NULL);
    /*
     * Convert big-endian octets into mpz_t values
     */
    mpz_import(m, RSAKEYSIZE/8, 1, 1, 1, 0, _m);
    mpz_import(k, RSAKEYSIZE/8, 1, 1, 1, 0, _k);
    mpz_import(n, RSAKEYSIZE/8, 1, 1, 1, 0, _n);
    /*
     * Compute m^k mod n
     */
    if (mpz_cmp(m, n) >= 0) {
        mpz_clears(m, k, n, NULL);
        return EM_MSG_OUT_OF_RANGE;
    }
    mpz_powm(m, m, k, n);
    /*
     * Convert mpz_t m into the octet string _m
     */
    mpz_export(_m, NULL, 1, RSAKEYSIZE/8, 1, 0, m);
    /*
     * Free the space occupied by mpz variables
     */
    mpz_clears(m, k, n, NULL);
    return 0;
}

/*
 * Copyright 2020. Heekuck Oh, all rights reserved
 * A mask generation function based on a hash function
 */
static unsigned char *mgf(const unsigned char *mgfSeed, size_t seedLen, unsigned char *mask, size_t maskLen)
{
    uint32_t i, count, c;
    size_t hLen;
    unsigned char *mgfIn, *m;

    /*
     * Check if maskLen > 2^32*hLen
     */
    hLen = SHASIZE/8;
    if (maskLen > 0x0100000000*hLen)
        return NULL;
    /*
     * Generate octet string mask
     */
    if ((mgfIn = (unsigned char *)malloc(seedLen+4)) == NULL)
        return NULL;;
    memcpy(mgfIn, mgfSeed, seedLen);
    count = maskLen/hLen + (maskLen%hLen ? 1 : 0);
    if ((m = (unsigned char *)malloc(count*hLen)) == NULL)
        return NULL;
    /*
     * Convert i to an octet string C of length 4 octets
     * Concatenate the hash of the seed mgfSeed and C to the octet string T:
     *       T = T || Hash(mgfSeed || C)
     */
    for (i = 0; i < count; i++) {
        c = i;
        mgfIn[seedLen+3] = c & 0x000000ff; c >>= 8;
        mgfIn[seedLen+2] = c & 0x000000ff; c >>= 8;
        mgfIn[seedLen+1] = c & 0x000000ff; c >>= 8;
        mgfIn[seedLen] = c & 0x000000ff;
        (*sha)(mgfIn, seedLen+4, m+i*hLen);
    }
    /*
     * Copy the mask and free memory
     */
    memcpy(mask, m, maskLen);
    free(mgfIn); free(m);
    return mask;
}

/*
 * rsassa_pss_sign - RSA Signature Scheme with Appendix
 */
int rsassa_pss_sign(const void *m, size_t mLen, const void *d, const void *n, void *s)
{
    uint8_t M[2 * SHA_DIGEST_SIZE + 8] = {};
    uint8_t mHash[SHA_DIGEST_SIZE] = {};
    uint8_t salt[SHA_DIGEST_SIZE] = {};
    uint8_t DB[DBSIZE] = {};
    uint8_t H[SHA_DIGEST_SIZE] = {};
    uint8_t EM[RSAKEYSIZE / 8] = {};
    uint8_t MGFH[DBSIZE] = {};
    uint8_t TF = 0xBC;
    uint8_t MSB = 0x01;

    if(mLen > MAX_MSG_SIZE) return EM_MSG_TOO_LONG; // 해시함수의 입력 데이터가 너무 길어 한도를 초과함
    if((SHA_DIGEST_SIZE / 8) * 2 + 2 > RSAKEYSIZE / 8) return EM_HASH_TOO_LONG; // 해시의 길이가 너무 커서 EM에 수용불가

    arc4random_buf(salt, SHA_DIGEST_SIZE); // make salt
    sha(m, mLen, mHash); // make mHash
    memcpy(M + 8, mHash, SHA_DIGEST_SIZE); // append mHash
    memcpy(M + 8 + SHA_DIGEST_SIZE, salt, SHA_DIGEST_SIZE); // append salt

    memcpy(DB + DBSIZE - SHA_DIGEST_SIZE - 1, &MSB, 1); // append 0x01
    memcpy(DB + DBSIZE - SHA_DIGEST_SIZE, salt, SHA_DIGEST_SIZE); // append salt

    sha(M, 2 * SHA_DIGEST_SIZE + 8, H); // make H

    mgf(H, SHA_DIGEST_SIZE, MGFH, DBSIZE); // make mgf(H)

    for(int i = 0 ; i < DBSIZE ; i++) {
        EM[i] = DB[i] ^ MGFH[i];  // append maskedDB
    }
    memcpy(EM + DBSIZE, H, SHA_DIGEST_SIZE); // append H
    memcpy(EM + DBSIZE + SHA_DIGEST_SIZE, &TF, 1); // TF는 1 바이트이며 0xBC로 채운다.
    if(EM[0] & 0x80) EM[0] &= 0x78; // EM의 가장 왼쪽 비트 (MSB)가 1이면 강제로 0으로 바꾼다.
    memcpy(s, EM, RSAKEYSIZE / 8);

    if(rsa_cipher(s, d, n)) return EM_MSG_OUT_OF_RANGE;

    return 0;
}

/*
 * rsassa_pss_verify - RSA Signature Scheme with Appendix
 */
int rsassa_pss_verify(const void *m, size_t mLen, const void *e, const void *n, const void *s)
{
    uint8_t M[2 * SHA_DIGEST_SIZE + 8] = {};
    uint8_t mHash[SHA_DIGEST_SIZE] = {};
    uint8_t salt[SHA_DIGEST_SIZE] = {};
    uint8_t DB[DBSIZE] = {};
    uint8_t maskedDB[DBSIZE] = {};
    uint8_t H[SHA_DIGEST_SIZE] = {};
    uint8_t EM[RSAKEYSIZE / 8] = {};
    uint8_t MGFH[DBSIZE] = {};
    uint8_t TF = 0xBC;

    memcpy(EM, s, RSAKEYSIZE / 8);
    if(rsa_cipher(EM, e, n)) return EM_MSG_OUT_OF_RANGE;
    if(EM[RSAKEYSIZE / 8 - 1] ^ TF) return EM_INVALID_LAST; // EM의 LSB 바이트가 0xbc가 아님
    if(EM[0] & 0x80) return EM_INVALID_INIT; // EM의 MSB 비트가 0이 아님

    memcpy(maskedDB, EM, DBSIZE); // make maskedDB
    memcpy(H, EM + DBSIZE, SHA_DIGEST_SIZE); // make H
    mgf(H, SHA_DIGEST_SIZE, MGFH, DBSIZE); // make mgf(H)

    sha(m, mLen, mHash); // make mHash

    for(int i = 1 ; i < DBSIZE ; i++) {
        DB[i] = maskedDB[i] ^ MGFH[i]; // make DB

        /*
            DB의 앞부분이 0x000000..|0x01 인지 확인.
         */
        if(i < DBSIZE - SHA_DIGEST_SIZE - 1) {
            if(DB[i] ^ 0x00) return EM_INVALID_PD2; // 00이 아니면 if문 발생
        }
        if(i == DBSIZE - SHA_DIGEST_SIZE - 1) {
            if(DB[i] ^ 0x01) return EM_INVALID_PD2; // 01이 아니면 if문 발생
        }
    }

    memcpy(salt, DB + DBSIZE - SHA_DIGEST_SIZE, SHA_DIGEST_SIZE); // make salt

    memcpy(M + 8, mHash, SHA_DIGEST_SIZE); // append mHash to M'
    memcpy(M + 8 + SHA_DIGEST_SIZE, salt, SHA_DIGEST_SIZE); // append salt to M'

    if(mLen > MAX_MSG_SIZE) return EM_MSG_TOO_LONG; // 해시함수의 입력 데이터가 너무 길어 한도를 초과함
    if((SHA_DIGEST_SIZE / 8) * 2 + 2 > RSAKEYSIZE / 8) return EM_HASH_TOO_LONG; // 해시의 길이가 너무 커서 EM에 수용불가

    sha(M, 8 + 2 * SHA_DIGEST_SIZE, mHash);

    for(int i = 0 ; i < SHA_DIGEST_SIZE ; i++) {
        if(H[i] ^ mHash[i]) return EM_HASH_MISMATCH; // 해시값이 일치하지 않음.
    }
    return 0;

}
