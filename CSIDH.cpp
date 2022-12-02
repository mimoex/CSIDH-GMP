#include "CSIDH.hpp"

int sign(int x) {
    if (x > 0) return 1;
    if (x == 0) return 0;
    return -1;
}

void genCSIDHkey(seckey* K)
{
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<int> distr(-5, 5);  //-5から5までの乱数を生成

    for (int i = 0; i < l; i++) K->e[i] = distr(eng);
}

//Algorithm 1: Verifying supersingularity.
//CSIDH論文 p15
bool validate(const mpz_class& a) {
    mpz_class k, d = 1;
    Point P, Q, temp_point, A_1;

    Fp::random_fp(P.X);
    P.Z = Fp::fpone;

    Fp::set_mpz(A_1.X, a);
    A_1.Z = Fp::fpone;

    bool isSupersingular = false;
    for (int i = 0; i < l; i++) {
        k = Fp::modmpz + 1;
        k /= primes[i];

        Q = xMUL(P, A_1, k);

        k = primes[i];
        temp_point = xMUL(Q, A_1, k);
        if (!mcl::bint::isZeroN(temp_point.Z.buf, Fp::N)) {
            isSupersingular = false;
            break;
        }

        if (!mcl::bint::isZeroN(Q.Z.buf, Fp::N))
            d *= primes[i];

        if (d > Fp::sqrt4) {
            isSupersingular = true;
            break;
        }
    }
    return isSupersingular;
}

bool isZero(const int* e, int n)
{
    for (int i = 0; i < n; i++) {
        if (e[i]) return false;
    }
    return true;
}


mpz_class action(const mpz_class& A, const seckey& Key) {
    mpz_class k, p_mul, q_mul, rhs_mpz, APointX_result;

    Fp rhs;
    Point P, Q, R, iso_A, iso_P, A_point;

    int s;
    Fp::set_mpz(A_point.X, A);
    //A_point.X = A;
    A_point.Z = Fp::fpone;

    int e[l];
    for (int i = 0; i < l; i++) {
        e[i] = Key.e[i];
    }

    //Evaluating the class group action.
    //A faster way to the CSIDH p4 https://eprint.iacr.org/2018/782.pdf
    while (!isZero(e, l)) {
        Fp::random_fp(P.X);

        rhs = calc_twist(A_point.X, P.X);

        s = Fp::isSquare(rhs);

        P.Z = Fp::fpone;

        k = 1;
        int S[l] = {};
        for (int i = 0; i < l; i++) {
            if (sign(e[i]) == s) {
                S[i] = 1;
                k *= primes[i];
            }
        }
        if (k == 1)
            continue;


        p_mul = Fp::modmpz + 1;
        p_mul /= k;

        Q = xMUL(P, A_point, p_mul);

        for (int i = 0; i < l; i++) {
            if (S[i] == 0)
                continue;

            q_mul = k / primes[i];



            R = xMUL(Q, A_point, q_mul);



            if (mcl::bint::isZeroN(R.Z.buf, Fp::N))
                continue;

            IsogenyCalc(A_point, Q, R, primes[i], &iso_A, &iso_P);

            Fp::div(A_point.X, iso_A.X, iso_A.Z);
            A_point.Z = Fp::fpone;

            Q = iso_P;

            e[i] -= s;

        }

        if (isZero(e, l)) break;
    }
    Fp::get_mpz(APointX_result, A_point.X);
    return APointX_result;
}