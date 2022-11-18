#pragma once
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <chrono>
#include <random>
#include "mcl/bint.hpp"

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
        mcl::bint::addT<N>(z.buf, x.buf, y.buf);
        if (mpn_cmp(z.buf, p.buf, N) >= 0)
            mcl::bint::subT<N>(z.buf, z.buf, p.buf);
    }

    //引き算 Mod p OK
    static void sub(Fp& z, const Fp& x, const Fp& y)
    {
        if (mcl::bint::subT<N>(z.buf, x.buf, y.buf)) {
            mcl::bint::addT<N>(z.buf, z.buf, p.buf);
        }
    }

    //MR    In:1024bit, Out:512bit
    static void MR(Fp& z, const FpDbl& x)
    {
        uint64_t temp[N * 2];
        uint64_t* tempH = temp + N;
        uint64_t temp512[N];

        mcl::bint::mulLowT<N>(temp512, x.buf, p.nr.buf);
        mcl::bint::mulT<N>(temp, temp512, p.p.buf);
        mcl::bint::addT<N*2>(temp, temp, x.buf);

        if (mpn_cmp(tempH, p.p.buf, N) >= 0)
            mcl::bint::subT<N>(z.buf, tempH, p.p.buf);
        else {
            for (int i = 0; i < N; i++) {
                z.buf[i] = tempH[i];
            }
        }
    }

    //MR    In:512bit, Out:512bit
    static void MR512(Fp& z, const Fp& x)
    {
        FpDbl temp;

        for (int i = 0; i < 8; i++) {
            temp.buf[i] = x.buf[i];
        }

        MR(z, temp);
    }

    //掛け算 Mod p(1回で完結)
    static void mul_test(Fp& z, const Fp& x, const Fp& y)
    {
        Fp temp512;
        FpDbl temp;
        mcl::bint::mulT<N>(temp.buf, x.buf, y.buf);
        MR(temp512,temp);
        mcl::bint::mulT<N>(temp.buf, temp512.buf, p.R2.buf);
        MR(z, temp);
    }

    //掛け算 Mod p(Montgomeryで返す)
    static void mul(Fp& z, const Fp& x, const Fp& y)
    {
        FpDbl temp;
        mcl::bint::mulT<N>(temp.buf, x.buf, y.buf);
        MR(z, temp);
    }

    static void sqr(Fp& z, const Fp& x)
    {
        FpDbl temp;
        mcl::bint::sqrT<N>(temp.buf, x.buf);
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
