#pragma once
#include "fp.hpp"

size_t inv_loop_N = 59;

typedef struct {
    int64_t u, v, q ,r;
} trans2x2;

static int64_t divsteps(int64_t zeta, uint64_t f0, uint64_t g0, trans2x2* t)
{
    uint64_t u = 8, v = 0, q = 0, r = 8;
    uint64_t c1, c2, f = f0, g = g0, x, y, z;
    int i;

    for(i=3;i<62;++i) {
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
    t->u = (int64_t)u;
    t->v = (int64_t)v;
    t->q = (int64_t)q;
    t->r = (int64_t)r;

    return zeta;
}

void update_fg(Fp& f, Fp& g, const trans2x2 t)
{
    Fp cf, cg;
    Fp u, v, q, r;
    Fp temp, temp1;
    mpz_class mpz_u, mpz_v, mpz_q, mpz_r;
    mpz_u = t.u; mpz_v = t.v; mpz_q = t.q; mpz_r = t.r;
    Fp::set_mpz(u, mpz_u);
    Fp::set_mpz(v, mpz_v);
    Fp::set_mpz(q, mpz_q);
    Fp::set_mpz(r, mpz_r);

    mcl::bint::mulT<Fp::N>(temp.buf, u.buf, f.buf);
    mcl::bint::mulT<Fp::N>(temp1.buf, v.buf, g.buf);
    mcl::bint::addT<Fp::N>(cf.buf, temp.buf, temp1.buf);

    mcl::bint::mulT<Fp::N>(temp.buf, q.buf, f.buf);
    mcl::bint::mulT<Fp::N>(temp1.buf, r.buf, g.buf);
    mcl::bint::addT<Fp::N>(cg.buf, temp.buf, temp1.buf);

    mcl::bint::shiftRight(f.buf, cf.buf, inv_loop_N, Fp::N);
    mcl::bint::shiftRight(g.buf, cg.buf, inv_loop_N, Fp::N);
}

void update_de(Fp& d, Fp& e, trans2x2 t, Fp& M, Fp& Mi )
{
    Fp d_sign, e_sign, md, me, cd, ce;
    Fp u, v, q, r;
    Fp temp, temp1, temp2;
    mpz_class mpz_u, mpz_v, mpz_q, mpz_r;
    mpz_u = t.u; mpz_v = t.v; mpz_q = t.q; mpz_r = t.r;
    Fp::set_mpz(u, mpz_u);
    Fp::set_mpz(v, mpz_v);
    Fp::set_mpz(q, mpz_q);
    Fp::set_mpz(r, mpz_r);

    mcl::bint::shiftRight(d_sign.buf, d.buf, 512, Fp::N);
    mcl::bint::shiftRight(e_sign.buf, e.buf, 512, Fp::N);



    //md = (t.u & d_sign) + (t.v & e_sign);
    for (int i = 0; i < Fp::N; ++i) {
        temp.buf[i] = u.buf[i] & d_sign.buf[i];
        temp1.buf[i] = v.buf[i] & e_sign.buf[i];
    }
    mcl::bint::addT<Fp::N>(md.buf, temp.buf, temp1.buf);

    //me = (t.q & d_sign) + (t.r & e_sign);
    for (int i = 0; i < Fp::N; ++i) {
        temp.buf[i] = q.buf[i] & d_sign.buf[i];
        temp1.buf[i] = r.buf[i] & e_sign.buf[i];
    }
    mcl::bint::addT<Fp::N>(me.buf, temp.buf, temp1.buf);

    //cd = (t.u * d + t.v * e) & (1ull << inv_loop_N) - 1;
    Fp shift;
    mcl::bint::shiftLeft(shift.buf, Fp::fpone.buf, inv_loop_N,Fp::N);
    shift.buf[0] -= 1;
    mcl::bint::mulT<Fp::N>(temp.buf, u.buf, d.buf);
    mcl::bint::mulT<Fp::N>(temp1.buf, v.buf, e.buf);
    mcl::bint::addT<Fp::N>(temp2.buf, temp.buf, temp1.buf);
    for (int i = 0; i < Fp::N; ++i) {
        cd.buf[i] = temp2.buf[i] & shift.buf[i];
    }

    //ce = (t.q * d + t.r * e) & (1ull << inv_loop_N) - 1;
    mcl::bint::mulT<Fp::N>(temp.buf, q.buf, d.buf);
    mcl::bint::mulT<Fp::N>(temp1.buf, r.buf, e.buf);
    mcl::bint::addT<Fp::N>(temp2.buf, temp.buf, temp1.buf);
    for (int i = 0; i < Fp::N; ++i) {
        ce.buf[i] = temp2.buf[i] & shift.buf[i];
    }

    //md -= (Mi * cd + md) % 2 * *N;
    //md -= ((uint64_t)Mi * cd + md) & ((1ULL << inv_loop_N) - 1);
    mcl::bint::mulT<Fp::N>(temp.buf, Mi.buf, cd.buf);
    mcl::bint::addT<Fp::N>(temp2.buf, temp.buf, md.buf);
    for (int i = 0; i < Fp::N; ++i) {
        md.buf[i] = temp2.buf[i] & shift.buf[i];
    }

    //me -= (Mi * ce + me) % 2 * *N;
    //me -= ((uint64_t)Mi * ce + me) & ((1ULL << inv_loop_N) - 1);
    mcl::bint::mulT<Fp::N>(temp.buf, Mi.buf, ce.buf);
    mcl::bint::addT<Fp::N>(temp2.buf, temp.buf, me.buf);
    for (int i = 0; i < Fp::N; ++i) {
        me.buf[i] = temp2.buf[i] & shift.buf[i];
    }

    //cd = t.u * d + t.v * e + M * md;
    mcl::bint::mulT<Fp::N>(temp.buf, u.buf, d.buf);
    mcl::bint::mulT<Fp::N>(temp1.buf, v.buf, e.buf);
    mcl::bint::mulT<Fp::N>(temp2.buf, M.buf, md.buf);
    mcl::bint::addT<Fp::N>(temp1.buf, temp.buf, temp1.buf);
    mcl::bint::addT<Fp::N>(cd.buf, temp1.buf, temp2.buf);

    //ce = t.q * d + t.r * e + M * me;
    mcl::bint::mulT<Fp::N>(temp.buf, q.buf, d.buf);
    mcl::bint::mulT<Fp::N>(temp1.buf, r.buf, e.buf);
    mcl::bint::mulT<Fp::N>(temp2.buf, M.buf, me.buf);
    mcl::bint::addT<Fp::N>(temp1.buf, temp.buf, temp1.buf);
    mcl::bint::addT<Fp::N>(ce.buf, temp1.buf, temp2.buf);

    //d = cd >> inv_loop_N;
    //e = ce >> inv_loop_N;
    mcl::bint::shiftRight(d.buf, cd.buf, inv_loop_N, Fp::N);
    mcl::bint::shiftRight(e.buf, ce.buf, inv_loop_N, Fp::N);
}

Fp normalize(Fp& sign, Fp& v, Fp& M)
{
    Fp v_sign, temp;
    Fp c;
    mcl::bint::shiftRight(v_sign.buf, v.buf, 512, Fp::N);
    //v += M & v_sign
    for(int i=0;i<Fp::N;++i) 
        temp.buf[i]=M.buf[i] & v_sign.buf[i];
    mcl::bint::addT<Fp::N>(v.buf, v.buf, temp.buf);

    //c.buf[0] = (sign - 1) >> 1;
    mcl::bint::subT<Fp::N>(temp.buf, sign.buf, Fp::fpone.buf);
    mcl::bint::shiftRight(c.buf, temp.buf, 1, Fp::N);

    for (int i = 0; i < Fp::N; ++i)
        temp.buf[i] = v.buf[i] ^ c.buf[i];
    mcl::bint::subT<Fp::N>(v.buf, temp.buf, c.buf);

    mcl::bint::shiftRight(v_sign.buf, v.buf, 512, Fp::N);

    for (int i = 0; i < Fp::N; ++i)
        temp.buf[i] = M.buf[i] & v_sign.buf[i];
    mcl::bint::addT<Fp::N>(v.buf, v.buf, temp.buf);

    return v;
}


static Fp my_const_gcd(Fp& z, const Fp& x)
{
    Fp Mi;
    Mi.buf[0] = 691;
    Fp d, e, f, g;
    e = Fp::fpone; f = Fp::p; g = x;
    int i;
    int64_t zeta = -1;

    for (i = 0; i < 10; ++i) {
        trans2x2 t;
        zeta = divsteps(zeta, f.buf[0], g.buf[0], &t);
        update_de(d, e, t,Fp::p, Mi);
        update_fg(f, g, t);
    }
    return(normalize(f, d, Fp::p) );
