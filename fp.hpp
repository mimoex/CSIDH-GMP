#pragma once
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <chrono>
#include <random>

struct Fp {
    int v;
    Fp(int _v=0) : v(_v) {}
    static const int N = 8; //64bit * 8 =512bit
    uint64_t buf[N]{};

    struct FpDbl {
        uint64_t buf[N * 2]{};
    };

    //montgomery剰余乗算
    static Fp inv4;
    static mpz_class sqrt4;    //4*(sqrt(mod))
    static size_t nbit;
    static Fp R;
    static Fp R2;
    static Fp nr;

    static Fp mrR2;
    static Fp MR2;  //2* para.R2
    static Fp MR4;  //4* para.R2
    static mpz_class modmpz;


    //足し算 Mod p OK
    static void add(Fp& z, const Fp& x, const Fp& y)
    {
        mpn_add_n(z.buf, x.buf, y.buf, N);
        if (mpn_cmp(z.buf, p.buf, N) >= 0)
            mpn_sub_n(z.buf, z.buf, p.buf, N);
    }

    //引き算 Mod p OK
    static void sub(Fp& z, const Fp& x, const Fp& y)
    {
        if (mpn_cmp(x.buf, y.buf, N) >= 0) {
            mpn_sub_n(z.buf, x.buf, y.buf, N);
        } else {
            mpn_add_n(z.buf, x.buf, p.buf, N);
            mpn_sub_n(z.buf, z.buf, y.buf, N);
        }
    }

    //MR    In:1024bit, Out:512bit
    static void MR(Fp& z, const FpDbl& x)
    {
        uint64_t temp1[N * 2]{}; //初期化しろ！
        uint64_t temp2[N * 2]{};
        uint64_t* tempH = temp2 + N;
        uint64_t temp512[N]{};

        mpn_mul_n(temp1, x.buf, p.nr.buf, N);
        mpn_mul_n(temp2, temp1, p.p.buf, N);
        mpn_add_n(temp2, temp2, x.buf, N * 2);

        if (mpn_cmp(tempH, p.p.buf, N) >= 0)
            mpn_sub_n(z.buf, tempH, p.p.buf, N);
        else {
            for (int i = 0; i < N; i++) {
                z.buf[i] = tempH[i];
            }
        }
    }

    //MR    In:512bit, Out:512bit
    static void MR512(Fp& z, const Fp& x)
    {
        FpDbl temp; //初期化しろ！

        for (int i = 0; i < 8; i++) {
            temp.buf[i] = x.buf[i];
        }

        MR(z, temp);
    }

    //掛け算 Mod p(1回で完結)
    static void mul_test(Fp& z, const Fp& x, const Fp& y)
    {
        mpz_class test;
        Fp temp512{};
        FpDbl temp{};
        mpn_mul_n(temp.buf, x.buf, y.buf,N);
        MR(temp512,temp);
        mpn_mul_n(temp.buf, temp512.buf, p.R2.buf, N);
        MR(z, temp);
    }

    //掛け算 Mod p(Montgomeryで返す)
    static void mul(Fp& z, const Fp& x, const Fp& y)
    {
        FpDbl temp;
        mpn_mul_n(temp.buf, x.buf, y.buf, N);
        MR(z, temp);
    }

    static void sqr(Fp& z, const Fp& x)
    {
        Fp temp512;
        FpDbl temp;
        mpn_sqr(temp.buf, x.buf, N);
        MR(z, temp);
    }

    static void pow(Fp& z, const Fp& x, const mpz_class& y)
    {
        Fp result=Fp::mrR2;
        
        for (int i = nbit-1; i >= 0; i--) {
            sqr(result, result);
            if (mpz_tstbit(y.get_mpz_t(), i) == 1) mul(result, result, x);
        }
        z = result;
    }

    static void inv(Fp& z, const Fp& y)
    {
        if (mpn_zero_p(y.buf, N))
            throw std::range_error("Divided by zero.");
        else
            pow(z, y, modmpz-2);
    }

    static void div(Fp& z, const Fp& x, const Fp& y)
    {
        FpDbl temp1024;
        Fp temp,xmon,ymon,invy;

        mul(xmon, x, p.R2);
        mul(ymon, y, p.R2);
        inv(invy, ymon);
        mul(temp, xmon, invy);

        MR512(z, temp);
    }

    //有限体pから乱数生成
    static void random_fp(Fp& z)
    {
        Fp temp;

        while (1) {
            mpn_random(temp.buf, N);

            if (mpn_cmp(temp.buf, p.buf, N) < 0) {
                z = temp;
                break;
            }
        }
    }

    static Fp p;
};

struct Point {
    Fp X, Z;
};
