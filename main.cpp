#include <time.h>
#include "fp.h"

bool csidh()
{
	mpz_class A;	//A=0 ( y^2 =x^3 +0*x^2 +x )
	clock_t t0, t1;

	//Aさんのstep1
	seckey secA;

	t0 = clock();
	genCSIDHkey(&secA);
	t1 = clock();
	cout << "Aさんの秘密鍵:" << endl;
	for(int i=0;i<N;i++) cout << secA.e[i] << ", ";
	cout << endl; 

	cout << "Aさんの秘密鍵生成時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n" << endl;

	//---

	mpz_class A_parm;

	t0 = clock();
	A_parm = action(A, secA);
	t1 = clock();

	cout << "Aさんの公開情報:" << A_parm << endl;
	cout << "supersingular?: " << validate(A_parm) << endl;

	cout << "Aさんの公開鍵生成時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n"<< endl;


	//Bさんのstep1
	seckey secB;

	t0 = clock();
	genCSIDHkey(&secB);
	t1 = clock();

	cout << "Bさんの秘密鍵:" << endl;
	for (int i = 0; i < N; i++) cout << secB.e[i] << ", ";
	cout << endl;

	cout << "Bさんの秘密鍵生成時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n" << endl;


	mpz_class B_parm;

	t0 = clock();
	B_parm = action(A, secB);
	t1 = clock();

	cout << "Bさんの公開情報:" << B_parm << endl;
	cout << "supersingular?: " << validate(B_parm) << endl;
	cout << endl;

	cout << "Bさんの公開鍵生成時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n" << endl;



	//Aさんのstep2
	mpz_class A_sec;

	t0 = clock();
	A_sec = action(B_parm, secA);
	t1 = clock();

	cout << "共有値(Aさん):" << A_sec << endl;
	cout << endl;

	cout << "Aさんの共有値を求める時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n" << endl;

	//Bさんのstep2
	mpz_class B_sec;

	t0 = clock();
	B_sec = action(A_parm, secB);
	t1 = clock();
	cout << "共有値(Bさん):" << B_sec << endl;
	cout << endl;

	cout << "Bさんの共有値を求める時間: " << 1000.000 * (t1 - t0) / CLOCKS_PER_SEC << "ms\n" << endl;


	//check
	if (A_sec == B_sec) return true;
	else return false;
}

int main()
{
	chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	bool check;
	check=csidh();

	end = std::chrono::system_clock::now();
	double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	cout << "time(全体): " << elapsed << endl;
	
	if (check) cout << "OK" << endl;
	else cout << "NG" << endl;


}