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

    //MR
    static void MR(Fp& z, const FpDbl& x)
    {
        uint64_t temp[N * 2];
        mpn_mul_n(temp, x.buf, para.nr, N);
        mpn_and_n(temp, temp, p.MR, N);
        mpn_mul_n(temp, temp, p.p, N);
        mpn_add_n(temp, temp, x.buf, N);
        mpn_rshift(temp, temp, p.nbit, N);
        if (mpn_cmp(temp, p.p, N) >= 0)
            mpn_sub_n(temp, temp, p.p, N);
        z.buf = temp;
    }

    //掛け算 Mod p
    static void mul(Fp& z, const Fp& x, const Fp& y)
    {
        uint64_t temp[N * 2];;
        mpn_mul_n(temp, x.buf, y.buf, N);
        MR(z, temp);
    }

    static Fp p;
};