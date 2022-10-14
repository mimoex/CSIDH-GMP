#pragma once
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>

using namespace std;

struct Fp {
    int v;
    Fp(int _v) : v(_v) {};
    static const int N = 8;
    uint64_t buf[N];


    //足し算 Mod p
    static void add(Fp& z, const Fp& x, const Fp& y)
    {
        mpn_add_n(z.buf, x.buf, y.buf, N);
        if (mpn_cmp(z.buf, p.buf, N) >= 0) {
            mpn_sub_n(z.buf, z.buf, p.buf, N);
        }
    }

    //引き算 Mod p
    static void sub(Fp& z, const Fp& x, const Fp& y)
    {
        
        if (mpn_cmp(x.buf, y.buf, N) >= 0) {
            mpn_sub_n(z.buf, x.buf, y.buf, N);
        } else {
            mpn_add_n(z.buf, x.buf, p.buf, N);
            mpn_sub_n(z.buf, z.buf, y.buf, N);
        }
    }

    //MR
    static void MR(Fp& z, const Fp& x)
    {
        Fp temp;
        mpn_mul_n(temp.buf, x.buf, p.nr, N);
        mpn_and_n(temp.buf, temp.buf, p.MR, N);
        mpn_mul_n(temp.buf, temp.buf, p.p, N);
        mpn_add_n(temp.buf, temp.buf, x.buf, N);
        mpn_rshift(temp.buf, temp.buf, p.nbit, N);
        if (mpn_cmp(temp.buf, p.p, N) >= 0)
            mpn_sub_n(temp.buf, temp.buf, p.p, N);
        z = temp;
    }

    //掛け算 Mod p
    static void mul(Fp& z, const Fp& x, const Fp& y)
    {
        Fp temp;
        mpn_mul_n(temp.buf, x.buf, y.buf, N);
        MR(z, temp);
    }

    static Fp p;
};