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

mpz_class action(const mpz_class& Am, const seckey& Key) {
    mpz_class rhs_mpz, APointX_result;

    mpz_class k[2];
    k[0] = 4;
    k[1] = 4;

    size_t e[2][l];

    for (size_t i = 0; i < l; ++i) {

        int t = Key.e[i];

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            k[1]*= primes[i];
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            k[0] *= primes[i];
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            k[0] *= primes[i];
            k[1] *= primes[i];
        }
    }

    Point A_pointm;
    Fp::set_mpz(A_pointm.X, Am);
    A_pointm.Z = Fp::mrR2;

    bool done[2] = { false, false };
    
    do {

        Point Pm, Qm;
        Fp::RandomGen(Pm.X);
        Pm.Z = Fp::mrR2;

        Fp rhs;
        calc_twist(rhs, A_pointm.X, Pm.X);

        bool s = Fp::isSquare(rhs);

        if (done[s])
            continue;

        Pm = xMULmon(Pm, A_pointm, k[s]);

        done[s] = true;

        for (size_t i = 0; i < l -1; ++i) {

            if (e[s][i] != 0) {
                mpz_class mul = 1;

                for (size_t j = i + 1; j < l; ++j) {
                    if (e[s][j] != 0) {
                        mul *= primes[j];
                    }
                }
                
                Point Rm;
                Rm = xMULmon(Pm, A_pointm, mul);

                if (!mcl::bint::isZeroN(Rm.Z.buf, Fp::N)) {
                    IsogenyCalc(&A_pointm, &Pm, Rm, primes[i]);  //iso_A, iso_Point

                    if (--e[s][i] == 0) {
                        k[s] *= primes[i];
                    }
                        
                }
            }
            done[s] &= (e[s][i] == 0);
            
        }
        Fp::div(A_pointm.X, A_pointm.X, A_pointm.Z);
        A_pointm.Z = Fp::mrR2;

    } while (!(done[0] && done[1]));

    Fp::get_mpz(APointX_result, A_pointm.X);
    return APointX_result;
}
