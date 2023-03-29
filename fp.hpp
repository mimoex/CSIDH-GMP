#pragma once
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <chrono>
#include <random>
#include "mcl/bint.hpp"

#define USE_NEW_MCL
#ifdef USE_NEW_MCL
#include "mcl.h"
#endif

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
    static Fp MRinv4;

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
#ifdef USE_NEW_MCL
        mcl_add(z.buf, x.buf, y.buf);
#else
        mcl::bint::addT<N>(z.buf, x.buf, y.buf);
        if (mcl::bint::cmpGeT<N>(z.buf, p.buf))
            mcl::bint::subT<N>(z.buf, z.buf, p.buf);
#endif
    }

    //引き算 Mod p OK
    static void sub(Fp& z, const Fp& x, const Fp& y)
    {
#ifdef USE_NEW_MCL
        mcl_sub(z.buf, x.buf, y.buf);
#else
        if (mcl::bint::subT<N>(z.buf, x.buf, y.buf))
            mcl::bint::addT<N>(z.buf, z.buf, p.buf);
#endif
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
            for (int i = 0; i < N; ++i) {
                z.buf[i] = tempH[i];
            }
        }
    }

    //MR    In:512bit, Out:512bit
    static void MR512(Fp& z, const Fp& x)
    {
        FpDbl temp;

        for (int i = 0; i < 8; ++i) {
            temp.buf[i] = x.buf[i];
        }
        //MR(z, temp);
        mcl_mod(z.buf, temp.buf);
    }

    //掛け算 Mod p(1回で完結)
    static void mul_test(Fp& z, const Fp& x, const Fp& y)
    {
        Fp temp512;
        FpDbl temp;
#ifdef USE_NEW_MCL
        mcl_mont(temp512.buf, x.buf, y.buf);
#else
        mcl::bint::mulT<N>(temp.buf, x.buf, y.buf);
        MR(temp512, temp);
#endif

#ifdef USE_NEW_MCL
        mcl_mont(z.buf, temp512.buf, R2.buf);
#else
        mcl::bint::mulT<N>(temp.buf, temp512.buf, p.R2.buf);
        MR(z, temp);
#endif
       
    }

    //掛け算 Mod p(Montgomeryで返す)
    static void mul(Fp& z, const Fp& x, const Fp& y)
    {
#ifdef USE_NEW_MCL
        mcl_mont(z.buf, x.buf, y.buf);
#else
        FpDbl temp;
        mcl::bint::mulT<N>(temp.buf, x.buf, y.buf);
        MR(z, temp);
#endif
    }

    static void sqr(Fp& z, const Fp& x)
    {
#ifdef USE_NEW_MCL
        mcl_mont(z.buf, x.buf, x.buf);
#else
        FpDbl temp;
        mcl::bint::sqrT<N>(temp.buf, x.buf);
        MR(z, temp);
#endif
    }

    static void pow(Fp& z, const Fp& x, const mpz_class& y)
    {
        Fp result = Fp::mrR2;

        for (int i = nbit - 1; i >= 0; --i) {
            sqr(result, result);
            if (mpz_tstbit(y.get_mpz_t(), i) == 1) mul(result, result, x);
        }
        z = result;
    }

        static void pow_window(Fp& z, const Fp& x, const mpz_class& y)
    {
        Fp exp, t;
        set_mpz(exp, y);
        Fp pre[16];

        pre[0] = fpone;
        pre[1] = x;
        for (int i = 2; i < 16; i=i+2) {
            sqr(pre[i], pre[i / 2]);
            mul(pre[i + 1], pre[i], x);
        }
        z = fpone;
        for (int i = nbit / 4 - 1; i >= 0; --i) {
            for (int j = 0; j < 4; ++j) {
                sqr(z, z);
            }
            int idx = (exp.buf[i / 16] >> uint64_t((i % 16) * 4)) & 15;
            mul(z, z, pre[idx]);
        }

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


    //constant-time binaryGCD
    static const size_t inv_loop_N = 54;

    typedef struct {
        int64_t u, v, q, r;
    } trans2x2;

    static void divsteps(int64_t& zeta, int64_t& f0, int64_t& g0, trans2x2& t)
    {
        int64_t u = 1, v = 0, q = 0, r = 1;
        int64_t c1, c2, x, y, z;
        int i;

        for (i = 0; i < inv_loop_N; ++i) {
            c1 = zeta >> 63;
            c2 = -(g0 & 1);

            x = (f0 ^ c1) - c1;
            y = (u ^ c1) - c1;
            z = (v ^ c1) - c1;

            g0 += x & c2;
            q += y & c2;
            r += z & c2;
            c1 &= c2;
            zeta ^= c1;
            zeta -= 1;
            f0 += g0 & c1;
            u += q & c1;
            v += r & c1;
            /* Shifts */
            g0 >>= 1;
            u <<= 1;
            v <<= 1;
        }
        t.u = u;
        t.v = v;
        t.q = q;
        t.r = r;
    }

    static void update_fg(mpz_class& f, mpz_class& g, const trans2x2& t)
    {
        mpz_class cf, cg;
        mpz_class temp, temp1;

        temp = t.u * f;
        temp1 = t.v * g;
        cf = temp + temp1;

        temp = t.q * f;
        temp1 = t.r * g;
        cg = temp + temp1;

        f = cf >> inv_loop_N;
        g = cg >> inv_loop_N;
    }

    static const uint64_t mask = (1ull << inv_loop_N) - 1;

    static void update_de(mpz_class& d, mpz_class& e, const trans2x2& t, const uint64_t& Mi)
    {
        mpz_class cd, ce, uv, qr, temp;
        int64_t d_sign, e_sign, md, me;
        uv = t.u * d + t.v * e;
        qr = t.q * d + t.r * e;

        //d_sign = d >> 512;
        d_sign = sgn(d);
        //test = e >> 512;
        e_sign = sgn(e);

        md = (t.u & d_sign) + (t.v & e_sign);

        me = (t.q & d_sign) + (t.r & e_sign);

        cd = uv & mask;

        ce = qr & mask;

        //md -= (Mi * cd + md) & mask;
        temp = Mi * cd;
        temp += md;
        md -= temp.get_si() & mask;

        //me -= (Mi * ce + me) % 2 **N; Pythonの場合
        //me -= (Mi * ce + me) & mask;
        temp = Mi * ce;
        temp += me;
        me -= temp.get_si() & mask;

        cd = uv + modmpz * md;

        ce = qr + modmpz * me;

        d = cd >> inv_loop_N;
        e = ce >> inv_loop_N;
    }

    static void normalize(const mpz_class& sign, mpz_class& v)
    {
        mpz_class v_sign,c;

        v_sign = v >> 512;
        v += modmpz & v_sign;

        c = (sign - 1) >> 1;

        v ^= c;
        v -= c;
        v_sign = v >> 512;

        v += modmpz & v_sign;
    }


    static void my_const_gcd(Fp& z, const Fp& x)
    {
        uint64_t Mi = 17680012166682291;
        mpz_class d = 0, e = 1, f = Fp::modmpz, g;
        get_mpz(g, x);
        int i;
        int64_t zeta = -1;

        uint64_t shift=1;
        shift <<=  inv_loop_N;
        --shift;

        trans2x2 t;
        int64_t fs, gs;
        mpz_class temp;


        for (i = 0; i < 20; ++i) {

            temp = f & shift;
            fs = temp.get_si();
            temp = g & shift;
            gs = temp.get_si();

            divsteps(zeta, fs, gs, t);

            update_fg(f, g, t);

            update_de(d, e, t, Mi);

        }
        normalize(f, d);
        set_mpz(z, d);
    }

    // //fermat inv
    // static void div(Fp& z, const Fp& x, const Fp& y)
    // {
    //     Fp temp, xmon, ymon, invy;

    //     mul(xmon, x, p.R2);
    //     mul(ymon, y, p.R2);
    //     inv(invy, ymon);
    //     mul(temp, xmon, invy);

    //     MR512(z, temp);
    // }


    // //binary gcd inv
    // static void div(Fp& z, const Fp& xmon, const Fp& ymon)
    // {
    //     Fp yinv, y, ytemp;

    //     MR512(y, ymon);
    //     binary_inv(yinv, y);
    //     mul(ytemp, yinv, p.R2);  //yinv*R

    //     mul(z, xmon, ytemp);  // x * y^-1
    // }

    //const-time binary gcd inv
    static void div(Fp& z, const Fp& xmon, const Fp& ymon)
    {
        Fp yinv, y, ytemp;

        MR512(y, ymon);
        my_const_gcd(yinv, y);
        mul(ytemp, yinv, p.R2);  //yinv*R

        mul(z, xmon, ytemp);  // x * y^-1
    }


    // x^3+ Ax^2 + x (p mod 4 ≡ 3) is square => return +1
    // else => return -1
    static int isSquare(const Fp& x)
    {
        Fp temp,check;
        pow(temp, x, p_plus_1_quarte);     //x^((p + 1) / 4)
        pow(check, temp, 2);
        if (mcl::bint::cmpEqT<N>(x.buf, check.buf))
            return 0;
        return 1;
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
