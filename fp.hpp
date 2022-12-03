#pragma once
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <chrono>
#include <random>
#include "mcl/bint.hpp"

void const_set();

struct Fp {
    //int v;
    //Fp(int _v = 0) : v(_v) {}
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
    static Fp fpone; //fptwo.buf[0]=1;
    static Fp fptwo; //fptwo.buf[0]=2;
    static mpz_class p_plus_1_quarte;

    static void set_mpz(Fp& z, const mpz_class& x)
    {
        mpz_export(&z, NULL, -1, 8, 0, 0, x.get_mpz_t());
    }

    static void get_mpz(mpz_class& z, const Fp& x)
    {
        mpz_import(z.get_mpz_t(), Fp::N, -1, 8, 0, 0, x.buf);
    }


    //足し算 Mod p OK
    static void add(Fp& z, const Fp& x, const Fp& y)
    {
        mcl::bint::addT<N>(z.buf, x.buf, y.buf);
        if (mcl::bint::cmpGeT<N>(z.buf, p.buf))
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
        mcl::bint::addT<N * 2>(temp, temp, x.buf);

        if(mcl::bint::cmpGeT<N>(tempH, p.p.buf))
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
        MR(temp512, temp);
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
        Fp result = Fp::mrR2;

        for (int i = nbit - 1; i >= 0; i--) {
            sqr(result, result);
            if (mpz_tstbit(y.get_mpz_t(), i) == 1) mul(result, result, x);
        }
        z = result;
    }

    //x^-1 = x^(p-2)
    static void inv(Fp& z, const Fp& x)
    {
        if (mcl::bint::isZeroT<N>(x.buf))
            throw std::range_error("Divided by zero.");
        else
            pow(z, x, modmpz - 2);
    }

    static int binary_inv(Fp& z, const Fp& y)
    {
        Fp u, v, r, s;
        u = p; v = y; s = fpone;

        if (mcl::bint::isZeroT<N>(y.buf)) {
            throw std::range_error("Divided by zero.");
            return -1;
        } else {
            while (!(mcl::bint::cmpEqT<N>(u.buf, fpone.buf)) && !(mcl::bint::cmpEqT<N>(v.buf, fpone.buf))) {

                while ((u.buf[0] & 1) == 0) {
                    mcl::bint::shrT<N>(u.buf, u.buf, 1);        //u >>= 1;
                    if ((r.buf[0] & 1) == 1) {
                        mcl::bint::addT<N>(r.buf, r.buf, p.buf);    //r += p;
                    }
                    mcl::bint::shrT<N>(r.buf, r.buf, 1);        //r >>= 1;
                }
                while ((v.buf[0] & 1) == 0) {
                    mcl::bint::shrT<N>(v.buf, v.buf, 1);        //v >>= 1;
                    if ((s.buf[0] & 1) == 1) {
                        mcl::bint::addT<N>(s.buf, s.buf, p.buf);    //s += p;
                    }
                    mcl::bint::shrT<N>(s.buf, s.buf, 1);        //s >>= 1;
                }
                if (mcl::bint::cmpGeT<N>(u.buf, v.buf)) {       //if (u >= v)
                    mcl::bint::subT<N>(u.buf, u.buf, v.buf);    //u -= v;
                    if (mcl::bint::subT<N>(r.buf, r.buf, s.buf))//r -= s;
                        mcl::bint::addT<N>(r.buf, r.buf, p.buf);
                }
                else {
                    mcl::bint::subT<N>(v.buf, v.buf, u.buf);    //v -= u;
                    if (mcl::bint::subT<N>(s.buf, s.buf, r.buf))    //s -= r;
                        mcl::bint::addT<N>(s.buf, s.buf, p.buf);
                }

            }
            if (mcl::bint::cmpEqT<N>(u.buf, fpone.buf)) {
                sub(z, r, p);
            }
            else {
                sub(z, s, p);
            }
        }
        return 0;
    }

    static void div_minus2(Fp& z, const Fp& x, const Fp& y)
    {
        Fp temp, xmon, ymon, invy;

        mul(xmon, x, p.R2);
        mul(ymon, y, p.R2);
        inv(invy, ymon);
        mul(temp, xmon, invy);

        MR512(z, temp);
    }

    static void div(Fp& z, const Fp& x, const Fp& y)
    {
        Fp temp, xmon, ymon, yinv;

        mul(xmon, x, p.R2);
        binary_inv(yinv, y);
        mul(ymon, yinv, p.R2);  //yinv*R

        mul(temp, xmon, ymon);  // x * y^-1

        MR512(z, temp);
    }

    // x^3+ Ax^2 + x (p mod 4 ≡ 3) is square => return +1
    // else => return -1
    static int isSquare(const Fp& x)
    {
        Fp temp,check;
        pow(temp, x, p_plus_1_quarte);     //x^((p + 1) / 4)
        pow(check, temp, 2);
        if (mcl::bint::cmpEqT<N>(x.buf, check.buf))
            return 1;
        return -1;
    }

    //有限体pから乱数生成
    static void RandomGen(Fp& z)
    {
        while (1) {
            mpn_random(z.buf, N);

            if (mcl::bint::cmpGtT<N>(p.buf, z.buf))   //if(p > z)
                break;
        }
    }

    static Fp p;
};

struct Point {
    Fp X, Z;
};
