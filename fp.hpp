#pragma once
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>

using namespace std;

struct Fp {
    int v;
    Fp(int _v=0) : v(_v) {}
    static const int N = 8;
    uint64_t buf[N]={};

    struct FpDbl {
        uint64_t buf[N * 2];
    };

    //montgomery剰余乗算
    static Fp inv4;
    //4*(sqrt(mod))
    static Fp sqrt4;

    static size_t nbit;
    static Fp R;
    static Fp R2;
    static Fp nr;
    static Fp mrR2;


    //足し算 Mod p OK
    static void add(Fp& z, const Fp& x, const Fp& y)
    {
        mpz_class result;
        mpn_add_n(z.buf, x.buf, y.buf, N);
        mpz_import(result.get_mpz_t(), Fp::N, -1, 8, 0, 0, z.buf);
        cout << result << endl;

        if (mpn_cmp(z.buf, p.buf, N) >= 0) {
            cout << "test" << endl;
            mpn_sub_n(z.buf, z.buf, p.buf, N);


            mpz_import(result.get_mpz_t(), Fp::N, -1, 8, 0, 0, z.buf);
            cout << result << endl;
        }
    }

    //引き算 Mod p OK
    static void sub(Fp& z, const Fp& x, const Fp& y)
    {

        if (mpn_cmp(x.buf, y.buf, N) >= 0) {
            mpn_sub_n(z.buf, x.buf, y.buf, N);
        }
        else {
            mpn_add_n(z.buf, x.buf, p.buf, N);
            mpn_sub_n(z.buf, z.buf, y.buf, N);
        }
    }

    //MR    In:1024bit, Out:512bit
    static void MR(Fp& z, const FpDbl& x)
    {
        mpz_class test;
        uint64_t temp[N * 2];
        uint64_t temp512[N]{};
        mpn_and_n(temp512, x.buf, p.R.buf, N * 2);
        mpn_mul_n(temp, temp512, p.nr.buf, N);
        mpn_and_n(temp512, temp, p.R.buf, N * 2);
        mpn_mul_n(temp, temp512, p.p.buf, N);
        mpn_add_n(temp, temp, x.buf, N * 2);
        mpn_rshift(temp, temp, N*2, p.nbit);

        for (int i = 0; i < 8; i++) {
            temp512[i] = temp[i+7];
        }

        if (mpn_cmp(temp512, p.p.buf, N) >= 0)
            mpn_sub_n(temp512, temp512, p.p.buf, N);

        *z.buf = *temp512;
    }

    //掛け算 Mod p
    static void mul(Fp& z, const Fp& x, const Fp& y)
    {
        Fp temp512{};
        FpDbl temp{};
        mpn_mul_n(temp.buf, x.buf, y.buf, N);
        MR(temp512,temp);
        mpn_mul_n(temp.buf, z.buf, p.R2.buf, N);
        MR(z,temp);
    }

    static Fp p;
};