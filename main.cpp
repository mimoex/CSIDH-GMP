#include "CSIDH.hpp"
#include <time.h>

seckey testAkey = { -4, -3, 3, 1, 5, -5, 1, 5, -5, 1, -3, -2, 4, 3, 2, 5, -2, -1, -2, 3, 2, 1, -4, -4, -1, -3, -4, 2, -2, -4, -2, 5, 2, -1, 2, 0, 3, -2, 1, 1, -1, -5, -3, 2, 3, 2, 5, -4, 1, 2, -4, -3, 2, 4, 1, -1, -3, -5, 2, 4, -4, 3, -4, -4, 3, -3, -2, -2, 3, 0, 3, 5, -5, 3 };
seckey testBkey = { -2, -5, -1, -1, -5, 0, 5, 4, 2, 1, 4, 1, -5, -3, 2, 3, 3, 1, -2, 4, -3, -4, 3, -1, 5, -3, 1, -5, 5, -1, 0, 3, 5, 5, 3, -5, 1, -5, -4, -3, -2, -4, 4, -4, 1, 0, -5, -3, -2, -4, 1, -1, -5, 5, -3, 3, 4, 5, -4, 0, 4, 2, 4, -3, -2, 3, -5, -5, -2, 3, -4, 2, -5, -1 };


Fp Fp::p;
Fp Fp::inv4;
mpz_class Fp::sqrt4;
size_t Fp::nbit;
Fp Fp::R;
Fp Fp::R2;
Fp Fp::nr;
Fp Fp::mrR2;
Fp Fp::MR2;
Fp Fp::MR4;
mpz_class Fp::modmpz;
Fp Fp::fpone;
Fp Fp::fptwo;
mpz_class Fp::p_plus_1_quarte;
Fp Fp::MRinv4;

using namespace std;

bool csidh()
{
    mpz_class A;	//A=0 ( y^2 =x^3 +0*x^2 +x )
    clock_t t0, t1;

    seckey secA;

    t0 = clock();
    secA = testAkey;
    //genCSIDHkey(&secA);
    t1 = clock();
    cout << "Aさんの秘密鍵:" << endl;
    for (int i = 0; i < l; i++) cout << secA.e[i] << ", ";
    cout << endl;

    cout << "Aさんの秘密鍵生成時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n" << endl;

    //---
    mpz_class A_parm;

    t0 = clock();
    A_parm = action(A, secA);
    t1 = clock();

    //mpz_import(parm.get_mpz_t(), Fp::N, -1, 8, 0, 0, A_parm.buf);


    cout << "Aさんの公開情報:" << A_parm << endl;
    cout << "supersingular?: " << validate(A_parm) << endl;

    cout << "Aさんの公開鍵生成時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n" << endl;

    //---

    //Bさんのstep1
    seckey secB;

    t0 = clock();
    secB = testBkey;
    //genCSIDHkey(&secB);
    t1 = clock();

    cout << "Bさんの秘密鍵:" << endl;
    for (int i = 0; i < l; i++) cout << secB.e[i] << ", ";
    cout << endl;

    cout << "Bさんの秘密鍵生成時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n" << endl;


    mpz_class B_parm;

    t0 = clock();
    B_parm = action(A, secB);
    t1 = clock();

    //mpz_import(parm.get_mpz_t(), Fp::N, -1, 8, 0, 0, B_parm.buf);


    cout << "Bさんの公開情報:" << B_parm << endl;
    cout << "supersingular?: " << validate(B_parm) << endl;
    cout << endl;

    cout << "Bさんの公開鍵生成時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n" << endl;




    //Aさんのstep2
    mpz_class A_sec;
    mpz_class sec;

    t0 = clock();
    A_sec = action(B_parm, secA);
    t1 = clock();

   // mpz_import(sec.get_mpz_t(), Fp::N, -1, 8, 0, 0, A_sec);


    cout << "共有値(Aさん):" << A_sec << endl;
    cout << endl;

    cout << "Aさんの共有値を求める時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n" << endl;

    //Bさんのstep2
    mpz_class B_sec;

    t0 = clock();
    B_sec = action(A_parm, secB);
    t1 = clock();

    //mpz_import(sec.get_mpz_t(), Fp::N, -1, 8, 0, 0, B_sec.buf);


    cout << "共有値(Bさん):" << B_sec << endl;
    cout << endl;

    cout << "Bさんの共有値を求める時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n" << endl;


    //check
    if (A_sec == B_sec) return true;
    else return false;
}


int main()
{
    const_set();

    chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    bool check;
    check = csidh();

    end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "time(全体): " << elapsed << endl;

    if (check) cout << "OK" << endl;
    else cout << "NG" << endl;

}


