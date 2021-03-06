/*
 * Copyright 2020. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었습니다.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "aes.h"
#define LRotWord(word, type) type tmp = word[0]; word[0] = word[1]; word[1] = word[2];  word[2] = word[3]; word[3] = tmp;

static const uint8_t sbox[256] = {
  0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };

static const uint8_t isbox[256] = {
  0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
  0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
  0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
  0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
  0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
  0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
  0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
  0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
  0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
  0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
  0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
  0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
  0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
  0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
  0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
  0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d };

static const uint8_t Rcon[11] = {0x8d, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36};

static const uint8_t M[16] = {2, 3, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 3, 1, 1, 2};

static const uint8_t IM[16] = {0x0e, 0x0b, 0x0d, 0x09, 0x09, 0x0e, 0x0b, 0x0d, 0x0d, 0x09, 0x0e, 0x0b, 0x0b, 0x0d, 0x09, 0x0e};


/*
 MixColumns 과정을 위한 gf_mul 함수
 */
uint8_t gf8_mul(uint8_t a, uint8_t b)
{
    uint8_t r = 0;
    while(b > 0)
    {
      if(b & 1) r = r ^ a;
      b = b >> 1;
      a = XTIME(a);
    }
    return r;
}

/*
 라운드 키와 state를 xor 연산
 uint32_t로 자료형 통일
 */
void AddRoundKey(uint8_t *state, const uint32_t *roundKey)
{
    uint32_t *tmp = (uint32_t *) state;
    for(int i = 0 ; i < Nb ; i++)
    {
        tmp[i] ^= roundKey[i];
    }
}

/*
 Substitution 과정을 위한 함수.
 ENCRYPT면 sbox와, DECRYPT면 isbox값으로 치환
*/
 
void SubBytes(uint8_t *state, int mode)
{
    for(int i = 0 ; i < BLOCKLEN ; i++)
    {
        state[i] = mode ? sbox[state[i]] : isbox[state[i]];
    }
}

/*
 ShiftRows 연산을 위한 함수
 ENCRYPT면 23,4열을 왼쪽으로 1,2,3번
 DECRYPT면 2,3,4열을 오른쪽으로 1,2,3번 Shift
 */
void ShiftRows(uint8_t *state, int mode)
{
    uint8_t tmp[16] = {};
    if(mode)
    {
        for(int i = 0 ; i < 16 ; i++)
        {
            if(i%4 == 1)
            {
                if(i - 4 > 0) tmp[i - 4] = state[i];
                else tmp[i - 4 + 16] = state[i];
            }
            else if(i%4 == 2)
            {
                if(i - 8 > 0) tmp[i - 8] = state[i];
                else tmp[i - 8 + 16] = state[i];
            }
            else if(i%4 == 3)
            {
                if(i - 12 > 0) tmp[i - 12] = state[i];
                else tmp[i - 12 + 16] = state[i];
            }
            else
            {
                tmp[i] = state[i];
            }
        }
    }
    
    else
    {
        for(int i = 0 ; i < 16 ; i++)
        {
            if(i%4 == 1)
            {
                if(i + 4 <= 15) tmp[i + 4] = state[i];
                else tmp[i + 4 - 16] = state[i];
            }
            else if(i%4 == 2)
            {
                if(i + 8 <= 15) tmp[i + 8] = state[i];
                else tmp[i + 8 - 16] = state[i];
            }
            else if(i%4 == 3)
            {
                if(i + 12 <= 15) tmp[i + 12] = state[i];
                else tmp[i + 12 - 16] = state[i];
            }
            else
            {
                tmp[i] = state[i];
            }
        }
    }
    for(int i = 0; i < 16 ; i++) state[i] = tmp[i];
}

/*
 MixColumns 과정을 위한 함수
 ENCRYPT면 M 행렬과, DECRYPT면 IM 행렬과 연산 진행
 */
void MixColumns(uint8_t *state, int mode)
{
    uint8_t tmp[Nb*Nb] = {0,};
    
    for(int i = 0 ; i < Nb ; i++)
    {
        for(int j = 0 ; j < Nb ; j++)
        {
            for(int k = 0 ; k < Nb ; k++)
            {
                if(mode) tmp[j * Nb + i] ^= gf8_mul(M[i * Nb + k], state[j * Nb + k]);
                else tmp[j * Nb + i] ^= gf8_mul(IM[i * Nb + k], state[j * Nb + k]);
            }
        }
    }
    for(int i = 0 ; i < Nb * Nb ; i++) state[i] = tmp[i];
}
    
/*
    각 라운드에서 사용할 라운드 키 생성
    설명이 필요한 세부 과정은 함수 내부에 따로 주석 표기
*/
void KeyExpansion(const uint8_t *key, uint32_t *roundKey)
{
    uint8_t p[Nb* RNDKEYSIZE];
    for(int i = 0 ; i < KEYLEN ; i++)
    {
        p[i] = key[i];
    }
    
    for(int i = 1 ; i <= Nr ; i++)
    {
        uint8_t w3[KEYLEN/4] = {p[i * KEYLEN - 4], p[i * KEYLEN - 3], p[i * KEYLEN - 2], p[i * KEYLEN - 1]};
        LRotWord(w3, uint8_t);
        for(int j = 0 ; j < KEYLEN/4 ; j++)
        {
            w3[j] = sbox[w3[j]];
            if(!j) w3[j] ^= Rcon[i]; // g(w3) 생성
        }
        
        for(int j = KEYLEN * i ; j < KEYLEN * (i + 1) ; j++) // 현재 라운드부터 다음 라운드 전 까지
        {
            if(j - KEYLEN * i < KEYLEN / 4) p[j] = p[j - KEYLEN] ^ w3[j - KEYLEN * i]; // (w4의 경우 g(w3)과 XOR연산)
            else p[j] = p[j - KEYLEN] ^ p[j - 4]; // (w5, w6, w7의 경우 바로 이전 것과 이전 라운드와 XOR
        }
    }
    memcpy((uint8_t *)roundKey, p, sizeof(uint8_t) * Nb * RNDKEYSIZE);
}

/*
    mode에 맞게 AES 진행
    DECRYPT일 경우 라운드키의 주소값을 이동시킨 상태에서 시작
*/
void Cipher(uint8_t *state, const uint32_t *roundKey, int mode)
{
    if(mode)
    {
        AddRoundKey(state, roundKey);
        roundKey += Nk;
        for(int i = 1; i <= Nr - 1 ; i++){
            SubBytes(state, mode);
            ShiftRows(state, mode);
            MixColumns(state, mode);
            AddRoundKey(state, roundKey);
            roundKey += Nk;
            
        }
        SubBytes(state, mode);
        ShiftRows(state, mode);
        AddRoundKey(state, roundKey);
    }
    else
    {
        roundKey += Nk * Nr;
        AddRoundKey(state, roundKey);
        roundKey -= Nk;
        for(int i = Nr - 1 ; i >= 1 ; i--)
        {
            ShiftRows(state, mode);
            SubBytes(state, mode);
            AddRoundKey(state, roundKey);
            MixColumns(state, mode);
            roundKey -= Nk;
        }
        ShiftRows(state, mode);
        SubBytes(state, mode);
        AddRoundKey(state, roundKey);
    }
}


