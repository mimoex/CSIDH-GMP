#include "fp.hpp"

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#diffadd-dadd-1987-m-3
Point Montgomery_ADDmon(const Point& Pm, const Point& Qm, const Point& Rm) {
	Fp A, B, C, D, DA, CB, temp;
	Point resultm;

	Fp::add(A, Qm.X, Qm.Z);
	Fp::sub(B, Qm.X, Qm.Z);
	Fp::add(C, Pm.X, Pm.Z);
	Fp::sub(D, Pm.X, Pm.Z);
	Fp::mul(DA, D, A);
	Fp::mul(CB, B, C);
	Fp::add(temp, DA, CB);
	Fp::sqr(temp, temp);
	Fp::mul(resultm.X, Rm.Z, temp);
	Fp::sub(temp, DA, CB);
	Fp::sqr(temp, temp);
	Fp::mul(resultm.Z, Rm.X, temp);
	return resultm;
}

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#doubling-dbl-1987-m-2
Point xDBLmon(const Point& Pm, const Point& Ap24m)
{
	Fp t0, t1, t2, t3;
	Point resultm;

	Fp::sub(t0, Pm.X, Pm.Z);
	Fp::add(t1, Pm.X, Pm.Z);
	Fp::sqr(t0, t0);
	Fp::sqr(t1, t1);
	Fp::mul(resultm.Z, Ap24m.Z, t0);
	Fp::mul(resultm.X, resultm.Z, t1);
	Fp::sub(t2, t1, t0);
	Fp::mul(t3, Ap24m.X, t2);
	Fp::add(resultm.Z, resultm.Z, t3);
	Fp::mul(resultm.Z, resultm.Z, t2);
	return resultm;
}

//Montgomery ladder
// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#ladder-ladd-1987-m-3
void xDBLADDmon(const Point& Pm, const Point& Qm, const Point& Rm, const Point& Ap24m,
	Point& DBLout, Point& ADDout) {
	Fp t0, t1, t2;

	Fp::add(t0, Pm.X, Pm.Z);
	Fp::sub(t1, Pm.X, Pm.Z);
	Fp::sqr(DBLout.X, t0);
	Fp::sub(t2, Qm.X, Qm.Z);

	Fp::add(ADDout.X, Qm.X, Qm.Z);
	Fp::mul(t0, t0, t2);
	Fp::sqr(DBLout.Z, t1);

	Fp::mul(t1, t1, ADDout.X);
	Fp::sub(t2, DBLout.X, DBLout.Z);
	Fp::mul(DBLout.X, DBLout.X, DBLout.Z);
	Fp::mul(ADDout.X, Ap24m.X, t2);
	Fp::sub(ADDout.Z, t0, t1);
	Fp::add(DBLout.Z, DBLout.Z, ADDout.X);
	Fp::add(ADDout.X, t0, t1);

	Fp::mul(DBLout.Z, DBLout.Z, t2);
	Fp::sqr(ADDout.Z, ADDout.Z);
	Fp::sqr(ADDout.X, ADDout.X);
	Fp::mul(ADDout.Z, ADDout.Z, Rm.X);
	Fp::mul(ADDout.X, ADDout.X, Rm.Z);
}

//モンゴメリ曲線のスカラー倍
Point xMUL(const Point& P, const Point& A_1, const mpz_class& n) {
	Point x0, Ap24;
	Point x0m, x1m, Pm, Ap24m;

	Fp::mul(Pm.X, P.X, Fp::p.R2);
	Fp::mul(Pm.Z, P.Z, Fp::p.R2);

	//x0 = P;
	x0m = Pm;

	//a24の計算	a24=(a+2)/4
	Fp::add(Ap24.X, A_1.X, Fp::fptwo);

	Fp::mul_test(Ap24.X, Ap24.X, Fp::p.inv4);
	Ap24.Z.buf[0] = 1;

	Fp::mul(Ap24m.X, Ap24.X, Fp::p.R2);
	Fp::mul(Ap24m.Z, Ap24.Z, Fp::p.R2);

	x1m = xDBLmon(Pm, Ap24m);

	size_t bit_size;
	bit_size = mpz_sizeinbase(n.get_mpz_t(), 2);

	for (int i = bit_size - 2; i >= 0; i--) {
		if (mpz_tstbit(n.get_mpz_t(), i) == 0) xDBLADDmon(x0m, x1m, Pm, Ap24m, x0m, x1m);
		else xDBLADDmon(x1m, x0m, Pm, Ap24m, x1m, x0m);
	}

	Fp::MR512(x0.X, x0m.X);
	Fp::MR512(x0.Z, x0m.Z);
	return x0;
}

//スカラー倍(Montgomeryで返す)
Point xMULmon(const Point& Pm, const Point& A_1, const mpz_class& n) {
	Point Ap24;
	Point x0m, x1m, Ap24m;

	x0m = Pm;

	//a24の計算	a24=(a+2)/4
	Fp::add(Ap24.X, A_1.X, Fp::fptwo);

	Fp::mul_test(Ap24.X, Ap24.X, Fp::p.inv4);
	Ap24.Z.buf[0] = 1;

	Fp::mul(Ap24m.X, Ap24.X, Fp::p.R2);
	Fp::mul(Ap24m.Z, Ap24.Z, Fp::p.R2);

	x1m = xDBLmon(Pm, Ap24m);

	size_t bit_size;
	bit_size = mpz_sizeinbase(n.get_mpz_t(), 2);

	for (int i = bit_size - 2; i >= 0; i--) {
		if (mpz_tstbit(n.get_mpz_t(), i) == 0) xDBLADDmon(x0m, x1m, Pm, Ap24m, x0m, x1m);
		else xDBLADDmon(x1m, x0m, Pm, Ap24m, x1m, x0m);
	}

	return x0m;
}

//モンゴメリ曲線右辺の計算
Fp calc_twist(const Fp& a, const Fp& x_mont) {
	Fp result, temp1, a_mont;

	Fp::mul(a_mont, a, Fp::p.R2);
	Fp::add(temp1, x_mont, a_mont);	// x+a
	Fp::mul(temp1, temp1, x_mont); // x^2 + ax
	Fp::add(temp1, temp1, Fp::p.mrR2); // x^2 + ax + 1
	Fp::mul(result, temp1, x_mont); // x^3 + ax^2 + x
	return result;
}

//Calc Isogeny
void IsogenyCalc(const Point& A, const Point& Pm, const Point& Km, const size_t& k,
	Point* Aout, Point* Pout) {

	Fp tmp0, tmp1, tmp2, Psum, Pdif;
	Point ed, prod, edt1;
	Point Am, Qm, Ap24_C24m, temp1;
	Point Aoutm, Poutm;

	Fp::mul(Am.X, A.X, Fp::p.R2);
	Fp::mul(Am.Z, A.Z, Fp::p.R2);

	Fp::add(Ap24_C24m.X, Am.X, Fp::MR2);
	Ap24_C24m.Z = Fp::p.MR4;

	Fp::add(Psum, Pm.X, Pm.Z);
	Fp::sub(Pdif, Pm.X, Pm.Z);

	Fp::sub(prod.X, Km.X, Km.Z);
	Fp::add(prod.Z, Km.X, Km.Z);

	Fp::mul(tmp1, prod.X, Psum);
	Fp::mul(tmp0, prod.Z, Pdif);
	Fp::add(Qm.X, tmp0, tmp1);
	Fp::sub(Qm.Z, tmp0, tmp1);

	Point M[3] = { Km };

	M[1] = xDBLmon(Km, Ap24_C24m);

	for (int i = 1; i < (k >> 1); ++i) {

		if (i >= 2)
			M[i % 3] = Montgomery_ADDmon(M[(i - 1) % 3], Km, M[(i - 2) % 3]);

		Fp::sub(tmp1, M[i % 3].X, M[i % 3].Z);
		Fp::add(tmp0, M[i % 3].X, M[i % 3].Z);
		Fp::mul(prod.X, prod.X, tmp1);
		Fp::mul(prod.Z, prod.Z, tmp0);
		Fp::mul(tmp1, tmp1, Psum);
		Fp::mul(tmp0, tmp0, Pdif);
		Fp::add(tmp2, tmp0, tmp1);
		Fp::mul(Qm.X, Qm.X, tmp2);
		Fp::sub(tmp2, tmp0, tmp1);
		Fp::mul(Qm.Z, Qm.Z, tmp2);
	}
	//Evaluation
	Fp::sqr(Qm.X, Qm.X);
	Fp::sqr(Qm.Z, Qm.Z);
	Fp::mul(Pout->X, Pm.X, Qm.X);
	Fp::mul(Pout->Z, Pm.Z, Qm.Z);

	//A faster way to the CSIDH p10くらい
	/* A = 2*(1+d)/(1-d) */
	Fp::add(edt1.Z, Am.Z, Am.Z);
	Fp::add(ed.X, Am.X, edt1.Z);
	Fp::sub(ed.Z, Am.X, edt1.Z);

	Fp::pow(ed.X, ed.X, k);
	Fp::pow(ed.Z, ed.Z, k);

	// compute prod.x^8, prod.z^8
	Fp::sqr(prod.X, prod.X);
	Fp::sqr(prod.X, prod.X);
	Fp::sqr(prod.X, prod.X);
	Fp::sqr(prod.Z, prod.Z);
	Fp::sqr(prod.Z, prod.Z);
	Fp::sqr(prod.Z, prod.Z);

	// compute image curve parameters
	Fp::mul(ed.Z, ed.Z, prod.X);
	Fp::mul(ed.X, ed.X, prod.Z);

	// compute Montgomery params
	Fp::add(temp1.X, ed.X, ed.Z);
	Fp::sub(Aout->Z, ed.X, ed.Z);
	Fp::add(Aout->X, temp1.X, temp1.X);
}