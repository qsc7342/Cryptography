/*
 * Copyright 2020, 2021. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었습니다.
 */
#include "miller_rabin.h"
/*
 * Miller-Rabin Primality Testing against small sets of bases
 *
 * if n < 2^64,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37.
 *
 * if n < 3,317,044,064,679,887,385,961,981,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, and 41.
 */
const uint64_t a[ALEN] = {2,3,5,7,11,13,17,19,23,29,31,37};

/*
 * miller_rabin() - Miller-Rabin Primality Test (deterministic version)
 *
 * n > 3, an odd integer to be tested for primality
 * It returns 1 if n is prime, 0 otherwise.
 */ // a is cur n is n
int miller_rabin(uint64_t n)
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

