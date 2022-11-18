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

using namespace std;

bool csidh()
{
    mpz_class A;	//A=0 ( y^2 =x^3 +0*x^2 +x )
    clock_t t0, t1;

    seckey secA;

    t0 = clock();
    secA = testAkey;
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
    //genCSIDHkey(&secB);
    secB = testBkey;
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
    Fp::modmpz = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";

    mpz_export(&Fp::p.buf, NULL, -1, 8, 0, 0, Fp::modmpz.get_mpz_t());
    mpz_export(&Fp::p.inv4.buf, NULL, -1, 8, 0, 0, mpz_class("1331684699081905773686966904488651388517342873708180584403111660513502390006644134406723028256595313406156735410987361198165720310405343322235720072016415").get_mpz_t());
    Fp::p.sqrt4 = "72984510660328627668487265992516213834621203782665312994583680849552234545871";
    Fp::p.nbit = 511;
    mpz_export(&Fp::p.R.buf, NULL, -1, 8, 0, 0, mpz_class("13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095").get_mpz_t());
    mpz_export(&Fp::p.R2.buf, NULL, -1, 8, 0, 0, mpz_class("4129249336087068599567712096533367844627103091116178550771554439546768636879316215877587311580678557833245516716846632497161687522600395855704996009744164").get_mpz_t());
    mpz_export(&Fp::p.nr.buf, NULL, -1, 8, 0, 0, mpz_class("11352847703487789629542365874112686817491332780459616199453250137299746292403308035585887700708945576897648221145487706796396106115280843778209115643193677").get_mpz_t());

    mpz_export(&Fp::p.mrR2.buf, NULL, -1, 8, 0, 0, mpz_class("2754330337287350910078289762296635019340622830926948702498668159613744910020393901548090072114140920440777974898587161268428120328703823368547888429952778").get_mpz_t());
    mpz_export(&Fp::p.MR2.buf, NULL, -1, 8, 0, 0, mpz_class("181921878247078725408711906638664484611874167021175067384889677173480260014211265469288031201900587256929008153224877744193359415786273448152896571839897").get_mpz_t());
    mpz_export(&Fp::p.MR4.buf, NULL, -1, 8, 0, 0, mpz_class("363843756494157450817423813277328969223748334042350134769779354346960520028422530938576062403801174513858016306449755488386718831572546896305793143679794").get_mpz_t());


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


