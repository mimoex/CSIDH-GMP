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
bool validate(const mpz_class& am) {
    mpz_class k, d = 1;
    Point Pm, Qm, temp_pointm, Am;

    Fp::RandomGen(Pm.X);
    //Pm.X = Fp::mrR2;
    Pm.Z = Fp::mrR2;

    Fp::set_mpz(Am.X, am);
    Am.Z = Fp::mrR2;

    bool isSupersingular = false;
    for (int i = 0; i < l; ++i) {
        k = Fp::modmpz + 1;
        k /= primes[i];

        Qm = xMULmon(Pm, Am, k);

        k = primes[i];
        temp_pointm = xMULmon(Qm, Am, k);
        if (!mcl::bint::isZeroN(temp_pointm.Z.buf, Fp::N)) {
            isSupersingular = false;
            break;
        }

        if (!mcl::bint::isZeroN(Qm.Z.buf, Fp::N))
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
    for (int i = 0; i < n; ++i) {
        if (e[i]) return false;
    }
    return true;
}


mpz_class action(const mpz_class& Am, const seckey& Key) {
    mpz_class k, mul, rhs_mpz, APointX_result;

    Fp rhs;
    Point Pm, Qm, Rm, A_pointm;
    int s;

    Fp::set_mpz(A_pointm.X, Am);
    A_pointm.Z = Fp::mrR2;

    int e[l];
    for (int i = 0; i < l; ++i) {
        e[i] = Key.e[i];
    }

    //Evaluating the class group action.
    //A faster way to the CSIDH p4 https://eprint.iacr.org/2018/782.pdf
    while (!isZero(e, l)) {
        Fp::RandomGen(Pm.X);

        calc_twist(rhs, A_pointm.X, Pm.X);
        //s = Fp::isSquare(rhs);    //これでも動く
        Fp::get_mpz(rhs_mpz, rhs);
        s = mpz_kronecker(rhs_mpz.get_mpz_t(), Fp::modmpz.get_mpz_t());     //こっちのほうが速い
        if (s == 0)
            continue;
        Pm.Z = Fp::mrR2;

        k = 1;
        int S[l] = {};
        for (int i = 0; i < l; ++i) {
            if (sign(e[i]) == s) {
                S[i] = 1;
                k *= primes[i];
            }
        }
        if (k == 1)
            continue;


        mul = Fp::modmpz + 1;
        mul /= k;

        Qm = xMULmon(Pm, A_pointm, mul);

        for (int i = 0; i < l; ++i) {
            if (S[i] == 0)
                continue;
            mul = k / primes[i];



            Rm = xMULmon(Qm, A_pointm, mul);



            if (mcl::bint::isZeroN(Rm.Z.buf, Fp::N))
                continue;

            IsogenyCalc(A_pointm, Qm, Rm, primes[i]);  //iso_A, iso_Point

            Fp::div(A_pointm.X, A_pointm.X, A_pointm.Z);
            A_pointm.Z = Fp::mrR2;

            e[i] -= s;

        }

        if (isZero(e, l)) break;
    }
    Fp::get_mpz(APointX_result, A_pointm.X);
    return APointX_result;
}