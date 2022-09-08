#include "fp.h"

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#diffadd-dadd-1987-m-3
Point Montgomery_ADD(const Point& P, const Point& Q, const Point& R) {
	mpz_class A, B, C, D, DA, CB, temp1, temp2;
	Point result;
	Point Pm,Qm,Rm,resultm;

	Pm.X = MR(P.X * para.R2);
	Pm.Z = MR(P.Z * para.R2);
	Qm.X = MR(Q.X * para.R2);
	Qm.Z = MR(Q.Z * para.R2);
	Rm.X = MR(R.X * para.R2);
	Rm.Z = MR(R.Z * para.R2);

	add_fp(Qm.X, Qm.Z, &A);
	sub_fp(Qm.X, Qm.Z, &B);
	add_fp(Pm.X, Pm.Z, &C);
	sub_fp(Pm.X, Pm.Z, &D);
	DA = MR(D * A);
	CB = MR(B * C);
	add_fp(DA, CB, &temp1);
	temp1 = MR(temp1 * temp1);
	resultm.X = MR(Rm.Z * temp1);

	sub_fp(DA, CB, &temp2);
	temp2 = MR(temp2 * temp2);
	resultm.Z = MR(Rm.X * temp2);

	result.X = MR(resultm.X);
	result.Z = MR(resultm.Z);
	return result;
}

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#doubling-dbl-1987-m-2
Point xDBL(const Point& P, const Point& Ap24)
{
	Point result;
	mpz_class t0, t1, t2, t3;
	Point Pm, Ap24m,resultm;

	Pm.X = MR(P.X * para.R2);
	Pm.Z = MR(P.Z * para.R2);
	Ap24m.X = MR(Ap24.X * para.R2);
	Ap24m.Z = MR(Ap24.Z * para.R2);

	
	sub_fp(Pm.X, Pm.Z, &t0);
	add_fp(Pm.X, Pm.Z, &t1);
	t0 = MR(t0 * t0);

	t1 = MR(t1 * t1);
	resultm.Z = MR(Ap24m.Z * t0);
	resultm.X = MR(resultm.Z * t1);
	result.X = MR(resultm.X);		//逆変換

	sub_fp(t1, t0, &t1);
	t0 = MR(Ap24m.X * t1);
	add_fp(resultm.Z, t0, &resultm.Z);

	resultm.Z = MR(resultm.Z * t1);
	result.Z = MR(resultm.Z);		//逆変換
	return result;
}

void xDBLADD(const Point& P, const Point& Q, const Point& R, const Point& Ap24_1,
	Point& DBLout, Point& ADDout) {
	mpz_class t0, t1, t2, t3;
	add_fp(P.X, P.Z, &t0);
	sub_fp(P.X, P.Z, &t1);
	sqr_fp(t0, &DBLout.X);
	sub_fp(Q.X, Q.Z, &t2);
	add_fp(Q.X, Q.Z, &ADDout.X);
	mul_fp(t0, t2, &t0);
	sqr_fp(t1, &DBLout.Z);

	mul_fp(t1, ADDout.X, &t1);
	sub_fp(DBLout.X, DBLout.Z, &t2);
	mul_fp(DBLout.X, DBLout.Z, &DBLout.X);
	mul_fp(Ap24_1.X, t2, &ADDout.X);
	sub_fp(t0, t1, &ADDout.Z);
	add_fp(DBLout.Z, ADDout.X, &DBLout.Z);
	add_fp(t0, t1, &ADDout.X);

	mul_fp(DBLout.Z, t2, &DBLout.Z);
	sqr_fp(ADDout.Z, &ADDout.Z);
	sqr_fp(ADDout.X, &ADDout.X);
	mul_fp(ADDout.Z, R.X, &ADDout.Z);
	mul_fp(ADDout.X, R.Z, &ADDout.X);
}

//モンゴメリ曲線のスカラー倍
Point xMUL(const Point& P, const Point& A_1, const mpz_class& n) {
	Point x0, x1, Ap24_1, Point_temp1, Point_temp2;
	mpz_class div_ans;

	x0 = P;
	//a24の計算
	Ap24_1.Z = 1;
	add_fp(A_1.X, 2, &Ap24_1.X);

	mul_fp(Ap24_1.X, para.inv4, &div_ans);
	Ap24_1.X = div_ans;

	x1 = xDBL(P, Ap24_1);

	size_t bit_size;
	bit_size = mpz_sizeinbase(n.get_mpz_t(), 2);

	for (int i = bit_size - 2; i >= 0; i--) {
		if (mpz_tstbit(n.get_mpz_t(), i) == 0) xDBLADD(x0, x1, P, Ap24_1, x0, x1);
		else xDBLADD(x1, x0, P, Ap24_1, x1, x0);
	}
	return x0;
}

//モンゴメリ曲線右辺の計算
mpz_class calc_twist(const mpz_class& a, const mpz_class& x) {
	mpz_class result, temp1, temp2, temp3;
	mpz_class a_mont, x_mont;

	a_mont = MR(a * para.R2);
	x_mont = MR(x * para.R2);
	add_fp(x_mont, a_mont, &temp1);	// x+a
	temp2 = MR(temp1 * x_mont); // x^2 + ax
	add_fp(temp2, para.mrR2, &temp3); // x^2 + ax + 1
	result = MR(temp3 * x_mont); // x^3 + ax^2 + x
	return MR(result);
}

void Evaluation(const Point& P, const Point& Q,
	mpz_class& Xout, mpz_class& Zout)
{
	Point temp;
	Point Pm,Qm;
	Pm.X = MR(P.X);
	Pm.Z = MR(P.Z);
	Qm.X = MR(Q.X);
	Qm.Z = MR(Q.Z);

	//sqr_fp(Q.X, &temp.X);
	temp.X = MR(Qm.X * Qm.X);
	//sqr_fp(Q.Z, &temp.Z);
	temp.Z = MR(Qm.Z * Qm.Z);
	//mul_fp(P.X, temp.X, &Xout);
	Xout = MR(MR(Pm.X * temp.X));
	//mul_fp(P.Z, temp.Z, &Zout);
	Zout = MR(MR(Pm.Z * temp.Z));
}

Point mont2ed(const Point& P)
{
	// compute twisted Edwards curve coefficients
	/* A = 2*(1+d)/(1-d) */
	Point ed;
	add_fp(P.Z, P.Z, &ed.Z);
	add_fp(P.X, ed.Z, &ed.X);
	sub_fp(P.X, ed.Z, &ed.Z);
	return ed;
}


//Isogeny
void IsogenyCalc(const Point& A, const Point& P, const Point& K, const size_t& k,
	Point* Aout, Point* Pout) {

	mpz_class tmp0, tmp1, tmp2, Psum, Pdif;
	Point Q, ed, prod, Ap24_C24;

	add_fp(A.X, 2, &Ap24_C24.X);
	Ap24_C24.Z = 4;

	add_fp(P.X, P.Z, &Psum);
	sub_fp(P.X, P.Z, &Pdif);

	sub_fp(K.X, K.Z, &prod.X);
	add_fp(K.X, K.Z, &prod.Z);

	mul_fp(prod.X, Psum, &tmp1);
	mul_fp(prod.Z, Pdif, &tmp0);
	add_fp(tmp0, tmp1, &Q.X);
	sub_fp(tmp0, tmp1, &Q.Z);

	Point M[3] = { K };

	M[1] = xDBL(K, Ap24_C24);

	for (int i = 1; i < (k >> 1); ++i) {

		if (i >= 2)
			M[i % 3] = Montgomery_ADD(M[(i - 1) % 3], K, M[(i - 2) % 3]);

		sub_fp(M[i % 3].X, M[i % 3].Z, &tmp1);
		add_fp(M[i % 3].X, M[i % 3].Z, &tmp0);
		mul_fp(prod.X, tmp1, &prod.X);
		mul_fp(prod.Z, tmp0, &prod.Z);
		mul_fp(tmp1, Psum, &tmp1);
		mul_fp(tmp0, Pdif, &tmp0);
		add_fp(tmp0, tmp1, &tmp2);
		mul_fp(Q.X, tmp2, &Q.X);
		sub_fp(tmp0, tmp1, &tmp2);
		mul_fp(Q.Z, tmp2, &Q.Z);

	}

	Evaluation(P, Q, Pout->X, Pout->Z);

	//A faster way to the CSIDH p10くらい
	ed = mont2ed(A);

	pow_fp(ed.X, k, &ed.X);
	pow_fp(ed.Z, k, &ed.Z);

	// compute prod.x^8, prod.z^8
	pow_fp(prod.X, 8, &prod.X);
	pow_fp(prod.Z, 8, &prod.Z);

	// compute image curve parameters
	mul_fp(ed.Z, prod.X, &ed.Z);
	mul_fp(ed.X, prod.Z, &ed.X);

	// compute Montgomery params
	add_fp(ed.X, ed.Z, &Aout->X);
	sub_fp(ed.X, ed.Z, &Aout->Z);
	add_fp(Aout->X, Aout->X, &Aout->X);

}