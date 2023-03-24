N=54

def divsteps_n_matrix(zeta, f, g):
    """Compute zeta and transition matrix t after N divsteps (multiplied by 2^N)."""
    u, v, q, r = 1, 0, 0, 1 # start with identity matrix
    for _ in range(N):
        c1 = zeta >> 63
        # Compute x, y, z as conditionally-negated versions of f, u, v.
        x, y, z = (f ^ c1) - c1, (u ^ c1) - c1, (v ^ c1) - c1


        c2 = -(g & 1)
        # Conditionally add x, y, z to g, q, r.
        g, q, r = g + (x & c2), q + (y & c2), r + (z & c2)
        c1 &= c2                     # reusing c1 here for the earlier c3 variable
        zeta = (zeta ^ c1) - 1       # inlining the unconditional zeta decrement here
        # Conditionally add g, q, r to f, u, v.
        f, u, v = f + (g & c1), u + (q & c1), v + (r & c1)
        # When shifting g down, don't shift q, r, as we construct a transition matrix multiplied
        # by 2^N. Instead, shift f's coefficients u and v up.
        g, u, v = g >> 1, u << 1, v << 1

    return zeta, (u, v, q, r)

def update_fg(f, g, t):
    """Multiply matrix t/2^N with [f, g]."""
    u, v, q, r = t
    cf, cg = u*f + v*g, q*f + r*g
    return cf >> N, cg >> N

def update_de(d, e, t, M, Mi):
    """Multiply matrix t/2^N with [d, e], modulo M."""
    u, v, q, r = t
    d_sign, e_sign = d >> 512, e >> 512
    md, me = (u & d_sign) + (v & e_sign), (q & d_sign) + (r & e_sign)
    cd, ce = (u*d + v*e) % 2**N, (q*d + r*e) % 2**N
    md -= (Mi*cd + md) % 2**N
    me -= (Mi*ce + me) % 2**N
    cd, ce = u*d + v*e + M*md, q*d + r*e + M*me
    return cd >> N, ce >> N

def normalize(sign, v, M):
    """Compute sign*v mod M, where v in (-2*M,M); output in [0,M)."""
    v_sign = v >> 512
    # Conditionally add M to v.
    v += M & v_sign
    c = (sign - 1) >> 1
    # Conditionally negate v.
    v = (v ^ c) - c
    v_sign = v >> 512
    # Conditionally add M to v again.
    v += M & v_sign
    return v

def modinv(M, Mi, x):
    """Compute the modular inverse of x mod M, given Mi=1/M mod 2^N."""
    zeta, f, g, d, e = -1, M, x, 0, 1
    i=0
    for _ in range((1041  + N - 1) // N):
        print("i {}".format(i))
        i+=1
        zeta, t = divsteps_n_matrix(zeta, f & ((1 << N) - 1), g & ((1 << N) - 1))
        f, g = update_fg(f, g, t)
        print(f)
        print(g)
        print("\n")
        d, e = update_de(d, e, t, M, Mi)
    return normalize(f, d, M)

def inversion(M, x):
    return pow(x,-1,M)

import random

mod=5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659
test=1331684699081905773686966904488651388517342873708180584403111660513502390006644134406723028256595313406156735410987361198165720310405343322235720072016415
for i in range(1):
    #test = random.randint(0, mod - 1)
    ans1=modinv(mod,17680012166682291,test)
    ans2=inversion(mod,test)

    if ans1!=ans2:
        print("error")