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
Point xDBLmon(const Point& Pm, const Point* Ap24m)
{
	Fp a, b, c;
	Point resultm;

	Fp::add(a, Pm.X, Pm.Z);
	Fp::sqr(a, a);
	Fp::sub(b, Pm.X, Pm.Z);
	Fp::sqr(b, b);
	Fp::sub(c, a, b);
	Fp::add(b, b, b); Fp::add(b, b, b);	// x4
	Fp::mul(b, b, Ap24m->Z);
	Fp::mul(resultm.X, a, b);
	Fp::add(a, Ap24m->Z, Ap24m->Z);
	Fp::add(a, a, Ap24m->X);
	Fp::mul(a, a, c);
	Fp::add(a, a, b);
	Fp::mul(resultm.Z, a, c);
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
	Fp::mul(DBLout.Z, DBLout.Z, Ap24m.Z);
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

//スカラー倍(Montgomeryで返す)
Point xMULmon(const Point& Pm, const Point& Am, const mpz_class& n) {
	Point x0m, x1m, Ap24m;

	x0m = Pm;	//R

	//a24の計算	a24=(a+2c:4c)
	Fp::add(Ap24m.X, Am.Z, Am.Z);
	Fp::add(Ap24m.Z, Ap24m.X, Ap24m.X);
	Fp::add(Ap24m.X, Ap24m.X, Am.X);

	x1m.X = Fp::mrR2;


	size_t bit_size;
	bit_size = mpz_sizeinbase(n.get_mpz_t(), 2) -1;

	do {
		bool bit = mpz_tstbit(n.get_mpz_t(), bit_size);

		if (bit) xDBLADDmon(x0m, x1m, Pm, Ap24m, x0m, x1m);
		else xDBLADDmon(x1m, x0m, Pm, Ap24m, x1m, x0m);
	} while (bit_size--);
	return x1m;
}


//Calc Isogeny
void IsogenyCalc(Point *Am, Point *Pm, const Point& Km, const size_t& k) {

	Fp tmp0, tmp1, tmp2, Psum, Pdif;
	Point ed, prod, edt1;
	Point Qm, Ap24_C24m, temp1;
	Point Aoutm, Poutm;

	//Fp::add(Ap24_C24m.X, Am->X, Fp::MR2);
	//Ap24_C24m.Z = Fp::p.MR4;

	Fp::add(Psum, Pm->X, Pm->Z);
	Fp::sub(Pdif, Pm->X, Pm->Z);

	Fp::sub(prod.X, Km.X, Km.Z);
	Fp::add(prod.Z, Km.X, Km.Z);

	Fp::mul(tmp1, prod.X, Psum);
	Fp::mul(tmp0, prod.Z, Pdif);
	Fp::add(Qm.X, tmp0, tmp1);
	Fp::sub(Qm.Z, tmp0, tmp1);

	Point M[3] = { Km };

	M[1] = xDBLmon(Km, Am);

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
	Fp::mul(Pm->X, Pm->X, Qm.X);
	Fp::mul(Pm->Z, Pm->Z, Qm.Z);

	//A faster way to the CSIDH p10くらい
	/* A = 2*(1+d)/(1-d) */
	Fp::add(ed.Z, Am->Z, Am->Z);
	Fp::add(ed.X, Am->X, ed.Z);
	Fp::sub(ed.Z, Am->X, ed.Z);

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
	Fp::sub(Am->Z, ed.X, ed.Z);
	Fp::add(Am->X, temp1.X, temp1.X);
}