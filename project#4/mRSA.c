/*
 * Copyright 2020, 2021. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었습니다.
 */
#include <stdio.h>
#include <stdlib.h>
#include "mRSA.h"

const uint64_t a[ALEN] = {2,3,5,7,11,13,17,19,23,29,31,37};

static uint64_t gcd(uint64_t a, uint64_t b)
{
    if(!a || !b) return a ? a : b;
    while(b != 0)
    {
      uint64_t tmp;
      tmp = a;
      a = b;
      b = tmp % b;
    }
    return a;
}

static uint64_t mul_inv(uint64_t a, uint64_t m)
{
    uint64_t d0 = a, d1 = m;
    uint64_t x0 = 1, x1 = 0, q;
    uint64_t tmp;
    int64_t sign = -1;;
    while(d1 > 1)
    {
      q = d0 / d1;
      tmp = d0 - q * d1; d0 = d1; d1 = tmp;
      tmp = x0 + q * x1; x0 = x1; x1 = tmp;
      sign = ~sign;
    }
    if(d1 == 1) return(sign ? m-x1 : x1);
    else return 0;
}

static uint64_t mod_add(uint64_t a, uint64_t b, uint64_t m)
{
    if(a >= m) a %= m;
    if(b >= m) b %= m;
    return (a >= m - b) ? a - (m - b) : a + b;
}


static uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t r = 0;
    while(b > 0) {
        if (b & 1)
            r = mod_add(r, a, m);
        b = b >> 1;
        a = mod_add(a, a, m);
    }
    return r;
}

static uint64_t mod_pow(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t r = 1;
    while (b > 0) {
        if (b & 1)
            r = mod_mul(r, a, m);
        b = b >> 1;
        a = mod_mul(a, a, m);
    }
    return r;
}

static int miller_rabin(uint64_t n)
{
    if(n == 2) return PRIME; // 2이면 소수
    if(n % 2 == 0 && n != 2) return COMPOSITE; // 그 외 짝수의 경우 합성수
    uint64_t tmp = n - 1;
    uint64_t k = 0;
    uint64_t q;
    while(tmp % 2 == 0) {
        tmp /= 2;
        k++;
    }
    q = tmp; // k와 q를 찾는 과정
    for(int i = 0 ; i < ALEN && a[i] < n - 1 ; i++) { // 모든 a배열 내에 존재하는 수에 대해 테스트
        uint64_t flag = 0;
        uint64_t tmp = mod_pow(a[i], q, n);
        if(tmp == 1) continue; // a[i]^q mod n 이 1일 경우 소수일 수도 있음. continue
        for(int j = 0 ; j < k ; j++) { // 0 to k - 1
            if(tmp % n == n - 1) {
                flag = 1;
                break; // (a[i]^q)^(2^j) mod n 이 n - 1일 경우 소수일 수도 있음. flag = 1 설정.
            }
            tmp = mod_mul(tmp, tmp, n);
        }
        if(flag) continue; // flag가 1이므로 소수일 수도 있음. continue;
        return COMPOSITE; // 아닐 경우 소수가 아님.
    }
    return PRIME;
}



void mRSA_generate_key(uint64_t *e, uint64_t *d, uint64_t *n)
{
    uint64_t p, q, lambda;
    p = 0;
    q = 0;
    lambda = 0;
    while (1)
    {
        while (1)
        {
            p = arc4random_uniform(4294967296 - 2147483648) + 2147483648;
            if (miller_rabin(p)) break;
        }

        while (1)
        {
            q = arc4random_uniform(4294967296 - 2147483648) + 2147483648;
            if (miller_rabin(q)) break;
        }
        
        if(p * q >= MINIMUM_N && p != q) break;
    }
    *n = p*q;
    lambda = ((p - 1) * (q - 1)) / gcd(p - 1, q - 1);
    if(lambda > 65537) *e = 65537;
    else
    {
        while(1)
        {
            *e = arc4random_uniform(lambda - 2) + 2;
            if(gcd(*e, lambda) == 1) break;
        }
    }
    *d = mul_inv(*e,lambda);
}


int mRSA_cipher(uint64_t *m, uint64_t k, uint64_t n)
{
    if(*m >= n)
    {
        return 1;
    }
    else
    {
        *m = mod_pow(*m, k, n);
        return 0;
    }
}
