#pragma once
#include "fp.hpp"

    //constant-time binaryGCD
    static const size_t inv_loop_N = 10;

    typedef struct {
        mpz_class u, v, q, r;
    } trans2x2;

    static mpz_class divsteps(mpz_class zeta, const mpz_class& f0, const mpz_class& g0, trans2x2* t)
    {
        mpz_class u = 1, v = 0, q = 0, r = 1;
        mpz_class c1, c2, f = f0, g = g0, x, y, z;
        int i;

        for (i = 0; i < inv_loop_N; ++i) {
            c1 = zeta >> 63;
            c2 = -(g & 1);

            x = (f ^ c1) - c1;
            y = (u ^ c1) - c1;
            z = (v ^ c1) - c1;

            g += x & c2;
            q += y & c2;
            r += z & c2;
            c1 &= c2;
            zeta = (zeta ^ c1) - 1;
            f += g & c1;
            u += q & c1;
            v += r & c1;
            /* Shifts */
            g >>= 1;
            u <<= 1;
            v <<= 1;
        }
        t->u = u;
        t->v = v;
        t->q = q;
        t->r = r;

        return zeta;
    }

    static void update_fg(mpz_class& f, mpz_class& g, const trans2x2* t)
    {
        mpz_class cf, cg;
        mpz_class u, v, q, r;
        mpz_class temp, temp1;
        u = t->u; v = t->v; q = t->q; r = t->r;

        temp = u * f;
        temp1 = v * g;
        cf = temp + temp1;

        temp = q * f;
        temp1 = r * g;
        cg = temp + temp1;

        f = cf >> inv_loop_N;
        g = cg >> inv_loop_N;
    }

    static void update_de(mpz_class& d, mpz_class& e, trans2x2 t, mpz_class& Mi)
    {
        mpz_class d_sign, e_sign, md, me, cd, ce;
        mpz_class u, v, q, r;
        mpz_class temp, temp1, temp2;
        u = t.u; v = t.v; q = t.q; r = t.r;

        d_sign = d >> 512;
        e_sign = e >> 512;

        md = (t.u & d_sign) + (t.v & e_sign);

        me = (t.q & d_sign) + (t.r & e_sign);

        cd = (t.u * d + t.v * e) & (1ull << inv_loop_N) - 1;

        ce = (t.q * d + t.r * e) & (1ull << inv_loop_N) - 1;

        md -= (Mi * cd + md) & ((1ULL << inv_loop_N) - 1);

        //me -= (Mi * ce + me) % 2 * *N;
        me -= (Mi * ce + me) & ((1ULL << inv_loop_N) - 1);

        cd = t.u * d + t.v * e + modmpz * md;

        ce = t.q * d + t.r * e + modmpz * me;

        d = cd >> inv_loop_N;
        e = ce >> inv_loop_N;
    }

    static mpz_class normalize(mpz_class& sign, mpz_class& v)
    {
        mpz_class v_sign, temp;
        mpz_class c;

        v_sign = v >> 512;
        v += modmpz & v_sign;

        c = (sign - 1) >> 1;

        v = (v ^ c) - c;
        v_sign = v >> 512;

        v += modmpz & v_sign;

        return v;
    }


    static void my_const_gcd(Fp& z, const Fp& x)
    {
        mpz_class Mi = 691;
        mpz_class d = 0, e = 1, f = Fp::modmpz, g;
        get_mpz(g, x);
        int i;
        mpz_class zeta = -1;
        mpz_class result;

        mpz_class shift;
        shift = 1 << inv_loop_N;
        shift -= 1;

        trans2x2 t;
        mpz_class fs, gs;

        for (i = 0; i < 108; ++i) {

            fs = f&shift;
            gs = g&shift;

            zeta = divsteps(zeta, fs, gs, &t);

            update_fg(f, g, &t);

            update_de(d, e, t, Mi);

        }
        result = normalize(f, d);
        set_mpz(z, result);
    }
