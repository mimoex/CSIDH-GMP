#include "CSIDH.hpp"

using namespace std;

int sign(const int& x) { return (x > 0) - (x < 0); }

void genCSIDHkey(seckey* K)
{
    random_device rd;
    default_random_engine eng(rd());
    uniform_int_distribution<int> distr(-5, 5);  //-5から5までの乱数を生成

    for (int i = 0; i < l; i++) K->e[i] = distr(eng);
}

//Algorithm 1: Verifying supersingularity.
//CSIDH論文 p15
bool validate(const Fp& a) {
    mpz_class k, d=1, mod;

    Fp one;
    one.buf[0] = 1;

    Fp x;
    Fp::random_fp(x);

    Point P, Q, temp_point, A_1;

    A_1.X = a;
    A_1.Z = one;

    P.X = x;
    P.Z = one;

    mpz_import(mod.get_mpz_t(), Fp::N, -1, 8, 0, 0, Fp::p.buf);

    bool isSupersingular = false;
    for (int i = 0; i < Fp::N; i++) {
        k = mod + 1;
        k /= primes[i];

        Q = xMUL(P, A_1, k);

        k = primes[i];
        temp_point = xMUL(Q, A_1, k);
        if (!mpn_zero_p(temp_point.Z.buf, Fp::N)) {
            isSupersingular = false;
            break;
        }

        if (!mpn_zero_p(Q.Z.buf, Fp::N))
            d *= primes[i];

        if (d > Fp::sqrt4) {
            isSupersingular = true;
            break;
        }
    }
    return isSupersingular;
}


Fp action(const Fp& A, const seckey& Key) {
    mpz_class mod, k, p_mul, q_mul, rhs_mpz;

    mpz_import(mod.get_mpz_t(), Fp::N, -1, 8, 0, 0, Fp::p.buf);

    Fp one;
    one.buf[0] = 1;

    Fp x, rhs;
    Point P, Q, R, temp1_point, temp2_point, A_point;

    int S[l], s;
    A_point.X = A;
    A_point.Z = one;

    int e[l]{};
    for (int i = 0; i < l; i++) {
        e[i] = Key.e[i];
    }


    //Evaluating the class group action.
    //A faster way to the CSIDH p4 https://eprint.iacr.org/2018/782.pdf
    bool flag = false;
    for (int i = 0; i < l; i++) {
        if (e[i] != 0) {
            flag = true;
            break;
        }
    }

    while (flag) {
        Fp::random_fp(x);
        rhs = calc_twist(A_point.X, x);
        mpz_import(rhs_mpz.get_mpz_t(), Fp::N, -1, 8, 0, 0, rhs.buf);

        s = mpz_kronecker(rhs_mpz.get_mpz_t(), mod.get_mpz_t());

        if (s == 0)
            continue;

        k = 1;
        memset(S, 0, sizeof S);
        for (int i = 0; i < l; i++) {
            if (sign(e[i]) == s) {
                S[i] = 1;
                k *= primes[i];
            }
        }
        if (k == 1)
            continue;


        p_mul = mod + 1;
        p_mul /= k;

        P.X = x;
        P.Z = one;

        Q = xMUL(P, A_point, p_mul);

        for (int i = 0; i < l; i++) {
            if (S[i] == 0)
                continue;

            q_mul = k / primes[i];

            //cout << "y^2 = x^3 + " << A_point.X << "*x^2 + x" << endl;
            R = xMUL(Q, A_point, q_mul);
            if(mpn_zero_p(R.Z.buf, Fp::N))
                continue;

            IsogenyCalc(A_point, Q, R, primes[i], &temp1_point, &temp2_point);

            Fp::div(A_point.X, temp1_point.X, temp1_point.Z);
            A_point.Z = one;

            Q = temp2_point;

            e[i] -= s;

        }



        flag = false;
        for (int i = 0; i < l; i++) {
            if (e[i] != 0) {
                flag = true;
                break;
            }
        }
    }

    return A_point.X;
}
