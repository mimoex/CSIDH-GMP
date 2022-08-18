#include "fp.h"

int sign(const int& x) { return (x > 0) - (x < 0); }

void genCSIDHkey(seckey* K)
{
    random_device rd;
    default_random_engine eng(rd());
    uniform_int_distribution<int> distr(-5, 5);  //-5から5までの乱数を生成

    for (int i = 0; i < N; i++) K->e[i] = distr(eng);
}

//Algorithm 1: Verifying supersingularity.
//CSIDH論文 p15
bool validate(const mpz_class& a, const mpz_class& mod) {
    mpz_class k, r2p, d=1;
	
    //4*(sqrt(mod))
    mpz_sqrt(r2p.get_mpz_t(), mod.get_mpz_t());
    r2p *= 4;

    mpz_class x;
    x=random_fp(mod);

    Point P, Q, temp_point, A_1;

    A_1.X = a;
    A_1.Z = 1;
	
    P.X = x;
    P.Z = 1;

    bool isSupersingular = false;
    for (int i = 0; i < N; i++) {
        k = mod + 1;
		k /= primes[i];

        temp_point =xMUL(P, A_1, k, mod);
        Q = temp_point;

        k=primes[i];
        temp_point=xMUL(Q, A_1, k, mod);
        if (!(temp_point.Z==0)) {
            isSupersingular = false;
            break;
        }

        if (!(Q.Z == 0))
            d *= primes[i];

        if (d > r2p) {
            isSupersingular = true;
            break;
        }
    }
    return isSupersingular;
}


mpz_class action(const mpz_class& A, const seckey& Key, const mpz_class& mod) {
    mpz_class x, rhs, k, p_mul, q_mul;

    Point P, Q, R, temp1_point, temp2_point, A_point;

    int S[N], s;
    A_point.X = A;
    A_point.Z = 1;

    int e[N];
    for (int i = 0; i < N; i++) {
        e[i] = Key.e[i];
    }


	//Evaluating the class group action.
	//A faster way to the CSIDH p4 https://eprint.iacr.org/2018/782.pdf
    bool flag = false;
    for (int i = 0; i < N; i++) {
        if (e[i] != 0) {
            flag = true;
            break;
        }
    }
    while (flag) {
        x=random_fp(mod);
        rhs=calc_twist(A_point.X, x, mod);

        s = mpz_kronecker(rhs.get_mpz_t(), mod.get_mpz_t());
		
        if (s == 0)
            continue;

        k=1;
        memset(S, 0, sizeof S);
        for (int i = 0; i < N; i++) {
            if (sign(e[i]) == s) {
                S[i] = 1;
                k *= primes[i];
            }
        }
        if (k == 1)
            continue;

        p_mul = mod + 1;
		p_mul /= k;

        P.X= x;
        P.Z = 1;
        Q=xMUL(P, A_point, p_mul, mod);

        for (int i = 0; i < N; i++) {
            if (S[i] == 0)
                continue;

            q_mul= k / primes[i];
            R=xMUL(Q, A_point, q_mul, mod);
            if (R.Z==0)
                continue;

            IsogenyCalc(A_point, Q, R, mod, primes[i], &temp1_point, &temp2_point);
			
            div_fp(temp1_point.X, temp1_point.Z, mod, &A_point.X);
            A_point.Z = 1;
			
            Q= temp2_point;

            e[i] -= s;

        }

        flag = false;
        for (int i = 0; i < N; i++) {
            if (e[i] != 0) {
                flag = true;
                break;
            }
        }
    }

    return A_point.X;
}