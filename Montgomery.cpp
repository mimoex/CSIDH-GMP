#include "fp.h"

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#diffadd-dadd-1987-m-3
Point Montgomery_ADDmon(const Point& Pm, const Point& Qm, const Point& Rm) {
	mpz_class A, B, C, D, DA, CB, temp1, temp2;
	Point resultm;

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
	return resultm;
}

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#doubling-dbl-1987-m-2
Point xDBLmon(const Point& Pm, const Point& Ap24m)
{
	mpz_class t0, t1;
	Point resultm;

	sub_fp(Pm.X, Pm.Z, &t0);
	add_fp(Pm.X, Pm.Z, &t1);
	t0 = MR(t0 * t0);

	t1 = MR(t1 * t1);
	resultm.Z = MR(Ap24m.Z * t0);
	resultm.X = MR(resultm.Z * t1);

	sub_fp(t1, t0, &t1);
	t0 = MR(Ap24m.X * t1);
	add_fp(resultm.Z, t0, &resultm.Z);

	resultm.Z = MR(resultm.Z * t1);
	return resultm;
}


//動かない…
void xDBLADDmon(const Point& Pm, const Point& Qm, const Point& Rm, const Point& Ap24m,
	Point& DBLout, Point& ADDout) {
	mpz_class t0, t1, t2;
	//Point Pm, Qm, Rm, Ap24m;
	//
	//Pm.X = MR(P.X *para.R2);
	//Pm.Z = MR(P.Z * para.R2);
	//Qm.X = MR(Q.X * para.R2);
	//Qm.Z = MR(Q.Z * para.R2);
	//Rm.X = MR(R.X * para.R2);
	//Rm.Z = MR(R.Z * para.R2);
	//Ap24m.X = MR(Ap24.X * para.R2);
	
	add_fp(Pm.X, Pm.Z, &t0);
	sub_fp(Pm.X, Pm.Z, &t1);
	//sqr_fp(t0, &DBLout.X);
	mul_mon(t0, t0, &DBLout.X);
	sub_fp(Qm.X, Qm.Z, &t2);
	add_fp(Qm.X, Qm.Z, &ADDout.X);
	//mul_fp(t0, t2, &t0);
	mul_mon(t0, t2, &t0);
	//sqr_fp(t1, &DBLout.Z);
	mul_mon(t1, t1, &DBLout.Z);

	//OK

	//mul_fp(t1, ADDout.X, &t1);
	mul_mon(t1, ADDout.X, &t1);
	sub_fp(DBLout.X, DBLout.Z, &t2);
	//mul_fp(DBLout.X, DBLout.Z, &DBLout.X);
	mul_mon(DBLout.X, DBLout.Z, &DBLout.X);
	//mul_fp(Ap24.X, t2, &ADDout.X);
	mul_mon(Ap24m.X, t2, &ADDout.X);
	sub_fp(t0, t1, &ADDout.Z);
	add_fp(DBLout.Z, ADDout.X, &DBLout.Z);
	add_fp(t0, t1, &ADDout.X);

	//mul_fp(DBLout.Z, t2, &DBLout.Z);
	mul_mon(DBLout.Z, t2, &DBLout.Z);
	//sqr_fp(ADDout.Z, &ADDout.Z);
	mul_mon(ADDout.Z, ADDout.Z, &ADDout.Z);
	//sqr_fp(ADDout.X, &ADDout.X);
	mul_mon(ADDout.X, ADDout.X, &ADDout.X);
	//mul_fp(ADDout.Z, R.X, &ADDout.Z);
	mul_mon(ADDout.Z, Rm.X, &ADDout.Z);
	//mul_fp(ADDout.X, R.Z, &ADDout.X);
	mul_mon(ADDout.X, Rm.Z, &ADDout.X);
	
	//ADDout.X = MR(ADDout.X);
	//ADDout.Z = MR(ADDout.Z);
	//DBLout.X = MR(DBLout.X);
	//DBLout.Z = MR(DBLout.Z);
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
	Point x0, Ap24_1, Point_temp1, Point_temp2;
	mpz_class div_ans;
	Point x0m, x1m, Pm, Ap24_1m;

	Pm.X = MR(P.X * para.R2);
	Pm.Z = MR(P.Z * para.R2);

	//x0 = P;
	x0m.X = MR(P.X * para.R2);
	x0m.Z = MR(P.Z * para.R2);
	
	//a24の計算
	Ap24_1.Z = 1;
	add_fp(A_1.X, 2, &Ap24_1.X);
	
	mul_fp(Ap24_1.X, para.inv4, &div_ans);
	Ap24_1.X = div_ans;
	
	Ap24_1m.X = MR(Ap24_1.X * para.R2);
	Ap24_1m.Z = MR(Ap24_1.Z * para.R2);

	x1m = xDBLmon(Pm, Ap24_1m);

	size_t bit_size;
	bit_size = mpz_sizeinbase(n.get_mpz_t(), 2);

	for (int i = bit_size - 2; i >= 0; i--) {
		if (mpz_tstbit(n.get_mpz_t(), i) == 0) xDBLADDmon(x0m, x1m, Pm, Ap24_1m, x0m, x1m);
		else xDBLADDmon(x1m, x0m, Pm, Ap24_1m, x1m, x0m);
	}
	x0.X = MR(x0m.X);
	x0.Z = MR(x0m.Z);
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

void Evaluation(const Point& Pm, const Point& Qm,
	mpz_class& Xout, mpz_class& Zout)
{
	Point temp;
	
	temp.X = MR(Qm.X * Qm.X);
	temp.Z = MR(Qm.Z * Qm.Z);
	Xout = MR(Pm.X * temp.X);
	Zout = MR(Pm.Z * temp.Z);
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

//2* para.R2
const mpz_class MR2 = mpz_class("2754330337287350910078289762296635019340622830926948702498668159613744910020393901548090072114140920440777974898587161268428120328703823368547888429952778");
//4* para.R2
const mpz_class MR4 = mpz_class("181921878247078725408711906638664484611874167021175067384889677173480260014211265469288031201900587256929008153224877744193359415786273448152896571839897");

//Isogeny
void IsogenyCalc(const Point& A, const Point& P, const Point& K, const size_t& k,
	Point* Aout, Point* Pout) {

	mpz_class tmp0, tmp1, tmp2, Psum, Pdif;
	Point ed, prod;
	Point Am, Pm, Qm, Km, Ap24_C24m;
	Point Aoutm, Poutm;

	Am.X = MR(A.X * para.R2);
	Am.Z = MR(A.Z * para.R2);
	Pm.X = MR(P.X * para.R2);
	Pm.Z = MR(P.Z * para.R2);
	Km.X = MR(K.X * para.R2);
	Km.Z = MR(K.Z * para.R2);

	add_fp(Am.X, MR2, &Ap24_C24m.X);
	Ap24_C24m.Z = MR4;

	add_fp(Pm.X, Pm.Z, &Psum);
	sub_fp(Pm.X, Pm.Z, &Pdif);

	sub_fp(Km.X, Km.Z, &prod.X);
	add_fp(Km.X, Km.Z, &prod.Z);
	
	mul_mon(prod.X, Psum, &tmp1);
	mul_mon(prod.Z, Pdif, &tmp0);
	add_fp(tmp0, tmp1, &Qm.X);
	sub_fp(tmp0, tmp1, &Qm.Z);

	Point M[3] = { Km };

	M[1] = xDBLmon(Km, Ap24_C24m);

	for (int i = 1; i < (k >> 1); ++i) {

		if (i >= 2)
			M[i % 3] = Montgomery_ADDmon(M[(i - 1) % 3], Km, M[(i - 2) % 3]);

		sub_fp(M[i % 3].X, M[i % 3].Z, &tmp1);
		add_fp(M[i % 3].X, M[i % 3].Z, &tmp0);
		mul_mon(prod.X, tmp1, &prod.X);
		mul_mon(prod.Z, tmp0, &prod.Z);
		mul_mon(tmp1, Psum, &tmp1);
		mul_mon(tmp0, Pdif, &tmp0);
		add_fp(tmp0, tmp1, &tmp2);
		mul_mon(Qm.X, tmp2, &Qm.X);
		sub_fp(tmp0, tmp1, &tmp2);
		mul_mon(Qm.Z, tmp2, &Qm.Z);

	}

	Evaluation(Pm, Qm, Pout->X, Pout->Z);	//in->montgomery, out->montgomery

	//A faster way to the CSIDH p10くらい
	ed = mont2ed(A);

	pow_mon(ed.X, k, &ed.X);
	pow_mon(ed.Z, k, &ed.Z);

	// compute prod.x^8, prod.z^8
	mul_mon(prod.X, prod.X, &prod.X);
	mul_mon(prod.X, prod.X, &prod.X);
	mul_mon(prod.X, prod.X, &prod.X);
	mul_mon(prod.Z, prod.Z, &prod.Z);
	mul_mon(prod.Z, prod.Z, &prod.Z);
	mul_mon(prod.Z, prod.Z, &prod.Z);

	// compute image curve parameters
	mul_mon(ed.Z, prod.X, &ed.Z);
	mul_mon(ed.X, prod.Z, &ed.X);

	// compute Montgomery params
	add_fp(ed.X, ed.Z, &Aoutm.X);
	sub_fp(ed.X, ed.Z, &Aoutm.Z);
	add_fp(Aoutm.X, Aoutm.X, &Aoutm.X);

	Aout->X = MR(Aoutm.X);	//逆変換
	Aout->Z = MR(Aoutm.Z);	//逆変換
}