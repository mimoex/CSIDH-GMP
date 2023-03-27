#pragma once
#include "fp.hpp"

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