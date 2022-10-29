#include "fp.hpp"

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#diffadd-dadd-1987-m-3
Point Montgomery_ADDmon(const Point& Pm, const Point& Qm, const Point& Rm) {
	Fp A, B, C, D, DA, CB, temp1, temp2, temp3, temp4;
	Fp::FpDbl multemp;
	Point resultm;

	Fp::add(A, Qm.X, Qm.Z);
	Fp::sub(B, Qm.X, Qm.Z);
	Fp::add(C, Pm.X, Pm.Z);
	Fp::sub(D, Pm.X, Pm.Z);
	Fp::mul(DA, D, A);
	Fp::mul(CB, B, C);
	Fp::add(temp1, DA, CB);
	Fp::sqr(temp2, temp1);
	Fp::mul(resultm.X, Rm.Z, temp2);
	Fp::sub(temp3, DA, CB);
	Fp::sqr(temp4, temp3);
	Fp::mul(resultm.Z, Rm.X, temp4);
	return resultm;
}

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#doubling-dbl-1987-m-2
Point xDBLmon(const Point& Pm, const Point& Ap24m)
{
	Fp t0, t1, t2, t3;
	Fp::FpDbl temp1, temp2, temp3, temp4, temp5, temp6;
	Point resultm;

	Fp::sub(t0, Pm.X, Pm.Z);
	Fp::add(t1, Pm.X, Pm.Z);
	mpn_sqr(temp1.buf, t0.buf, Fp::N);
	Fp::MR(t0, temp1);

	mpn_sqr(temp2.buf, t1.buf, Fp::N);
	Fp::MR(t1, temp2);

	mpn_mul_n(temp3.buf, Ap24m.Z.buf, t0.buf, Fp::N);
	Fp::MR(resultm.Z, temp3);
	mpn_mul_n(temp4.buf, resultm.Z.buf, t1.buf, Fp::N);
	Fp::MR(resultm.X, temp4);

	Fp::sub(t2, t1, t0);
	mpn_mul_n(temp5.buf, Ap24m.X.buf, t2.buf, Fp::N);
	Fp::MR(t3, temp5);

	Fp::add(resultm.Z, resultm.Z, t3);

	mpn_mul_n(temp6.buf, resultm.Z.buf, t2.buf, Fp::N);
	Fp::MR(resultm.Z, temp6);
	return resultm;
}

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
	Point x0, Ap24_1, Ap24_0;
	Point x0m, x1m, Pm, Ap24_1m;
	Fp::FpDbl multemp;
	Fp two;

	Fp::mul(Pm.X, P.X, Fp::p.R2);
	Fp::mul(Pm.Z, P.Z, Fp::p.R2);

	//x0 = P;
	x0m.X = Pm.X;
	x0m.Z = Pm.Z;

	//a24の計算
	Ap24_0.Z.buf[0] = 1;
	two.buf[0] = 2;
	Fp::add(Ap24_0.X, A_1.X, two);

	Fp::mul_test(Ap24_1.X, Ap24_0.X, Fp::p.inv4);
	Ap24_1.Z.buf[0] = 1;

	Fp::mul(Ap24_1m.X, Ap24_1.X, Fp::p.R2);
	Fp::mul(Ap24_1m.Z, Ap24_1.Z, Fp::p.R2);

	x1m = xDBLmon(Pm, Ap24_1m);

	size_t bit_size;
	bit_size = mpz_sizeinbase(n.get_mpz_t(), 2);

	for (int i = bit_size - 2; i >= 0; i--) {
		if (mpz_tstbit(n.get_mpz_t(), i) == 0) xDBLADDmon(x0m, x1m, Pm, Ap24_1m, x0m, x1m);
		else xDBLADDmon(x1m, x0m, Pm, Ap24_1m, x1m, x0m);
	}

	Fp::MR512(x0.X, x0m.X);
	Fp::MR512(x0.Z, x0m.Z);
	return x0;
}

//モンゴメリ曲線右辺の計算
Fp calc_twist(const Fp& a, const Fp& x) {
	Fp result, temp1, temp2, temp3, temp4;
	Fp a_mont, x_mont;

	mpz_class t1, t2, t3, t4, t5, t6;

	Fp::mul(a_mont, a, Fp::p.R2);
	Fp::mul(x_mont, x, Fp::p.R2);
	Fp::add(temp1, x_mont, a_mont);	// x+a
	Fp::mul(temp2, temp1, x_mont); // x^2 + ax
	Fp::add(temp3, temp2, Fp::p.mrR2); // x^2 + ax + 1
	Fp::mul(temp4, temp3, x_mont); // x^3 + ax^2 + x
	Fp::MR512(result, temp4);
	return result;
}

void Evaluation(const Point& Pm, const Point& Qm,
	Fp& Xout, Fp& Zout)
{
	Point temp;

	Fp::mul(temp.X, Qm.X, Qm.X);
	Fp::mul(temp.Z, Qm.Z, Qm.Z);
	Fp::mul(Xout, Pm.X, temp.X);
	Fp::mul(Zout, Pm.Z, temp.Z);
}

Point mont2ed(const Point& P)
{
	// compute twisted Edwards curve coefficients
	/* A = 2*(1+d)/(1-d) */
	Point t0, ed;
	Fp::add(t0.Z, P.Z, P.Z);
	Fp::add(ed.X, P.X, t0.Z);
	Fp::sub(ed.Z, P.X, t0.Z);
	return ed;
}

//Isogeny
void IsogenyCalc(const Point& A, const Point& P, const Point& K, const size_t& k,
	Point* Aout, Point* Pout) {

	Fp tmp0, tmp1, tmp2, Psum, Pdif;
	Point ed, prod;
	Point Am, Pm, Qm, Km, Ap24_C24m, temp1;
	Point Aoutm, Poutm;

	Fp::mul(Am.X, A.X, Fp::p.R2);
	Fp::mul(Am.Z, A.Z, Fp::p.R2);
	Fp::mul(Pm.X, P.X, Fp::p.R2);
	Fp::mul(Pm.Z, P.Z, Fp::p.R2);
	Fp::mul(Km.X, K.X, Fp::p.R2);
	Fp::mul(Km.Z, K.Z, Fp::p.R2);

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
	Evaluation(Pm, Qm, Pout->X, Pout->Z);	//in->montgomery, out->montgomery

	//A faster way to the CSIDH p10くらい
	ed = mont2ed(Am);	//確認済み

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
	Fp::sub(Aoutm.Z, ed.X, ed.Z);
	Fp::add(Aoutm.X, temp1.X, temp1.X);

	Fp::MR512(Aout->X, Aoutm.X);	//逆変換
	Fp::MR512(Aout->Z, Aoutm.Z);	//逆変換
}