/*
 * Copyright 2020, 2021. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었습니다.
 */
#ifndef RSA_PSS_H
#define RSA_PSS_H

#include "sha2.h"

#define RSAKEYSIZE 2048
#define SHA256

#if defined(SHA224)
#define SHASIZE 224
#define SHA_DIGEST_SIZE (224 / 8)
#define SHA_BLOCK_SIZE (512 / 8)
#define MAX_MSG_SIZE 0x1FFFFFFFFFFFFFFF
#define DBSIZE RSAKEYSIZE / 8 - SHA_DIGEST_SIZE - 1
#elif defined(SHA256)
#define SHASIZE 256
#define SHA_DIGEST_SIZE (256 / 8)
#define SHA_BLOCK_SIZE (512 / 8)
#define MAX_MSG_SIZE 0x1FFFFFFFFFFFFFFF
#define DBSIZE RSAKEYSIZE / 8 - SHA_DIGEST_SIZE - 1
#elif defined(SHA384)
#define SHASIZE 384
#define SHA_DIGEST_SIZE (384 / 8)
#define SHA_BLOCK_SIZE (1024 / 8)
#define MAX_MSG_SIZE 0xFFFFFFFFFFFFFFFF
#define DBSIZE RSAKEYSIZE / 8 - SHA_DIGEST_SIZE - 1
#else
#define SHASIZE 512
#define SHA_DIGEST_SIZE (512 / 8)
#define SHA_BLOCK_SIZE (1024 / 8)
#define MAX_MSG_SIZE 0xFFFFFFFFFFFFFFFF
#define DBSIZE RSAKEYSIZE / 8 - SHA_DIGEST_SIZE - 1
#endif

#define EM_MSG_OUT_OF_RANGE 1
#define EM_MSG_TOO_LONG 2
#define EM_HASH_TOO_LONG 3
#define EM_INVALID_LAST 4
#define EM_INVALID_INIT 5
#define EM_INVALID_PD2 6
#define EM_HASH_MISMATCH 7

void rsa_generate_key(void *e, void *d, void *n, int mode);
int rsassa_pss_sign(const void *m, size_t mLen, const void *d, const void *n, void *s);
int rsassa_pss_verify(const void *m, size_t mLen, const void *e, const void *n, const void *s);

#endif
