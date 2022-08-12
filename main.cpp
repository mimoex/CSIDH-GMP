#include "fp.h"

bool csidh()
{
	mpz_class mod, A;

	mod = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";

	//Aさんのstep1
	chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	seckey secA;
	genCSIDHkey(&secA);
	cout << "Aさんの秘密鍵:" << endl;
	for(int i=0;i<N;i++) cout << secA.e[i] << ", ";
	cout << endl; 

	end = std::chrono::system_clock::now();
	double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	cout << "Aさんの秘密鍵生成時間: " << elapsed << "\n" << endl;

	//---
	
	start = std::chrono::system_clock::now();

	mpz_class A_parm;

	A_parm = action(A, secA, mod);
	cout << "Aさんの公開情報:" << A_parm << endl;
	cout << "supersingular?: " << validate(A_parm, mod) << endl;

	end = std::chrono::system_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	cout << "Aさんの公開鍵生成時間: " << elapsed << "\n" << endl;


	//Bさんのstep1
	start = std::chrono::system_clock::now();

	seckey secB;
	genCSIDHkey(&secB);
	cout << "Bさんの秘密鍵:" << endl;
	for (int i = 0; i < N; i++) cout << secB.e[i] << ", ";
	cout << endl;

	end = std::chrono::system_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	cout << "Bさんの秘密鍵生成時間: " << elapsed << "\n" << endl;



	start = std::chrono::system_clock::now();

	mpz_class B_parm;

	B_parm = action(A, secB, mod);
	cout << "Bさんの公開情報:" << B_parm << endl;
	cout << "supersingular?: " << validate(B_parm, mod) << endl;
	cout << endl;

	end = std::chrono::system_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	cout << "Bさんの公開鍵生成時間: " << elapsed << "\n" << endl;



	//Aさんのstep2
	start = std::chrono::system_clock::now();
	mpz_class A_sec;
	A_sec = action(B_parm, secA, mod);
	cout << "共有値(Aさん):" << A_sec << endl;
	cout << endl;

	end = std::chrono::system_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	cout << "Aさんの共有値を求める時間: " << elapsed << "\n" << endl;

	//Bさんのstep2
	start = std::chrono::system_clock::now();
	mpz_class B_sec;
	B_sec = action(A_parm, secB, mod);
	cout << "共有値(Bさん):" << B_sec << endl;
	cout << endl;

	end = std::chrono::system_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	cout << "Bさんの共有値を求める時間: " << elapsed << "\n" << endl;


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