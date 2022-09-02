#include "fp.h"

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#diffadd-dadd-1987-m-3
Point Montgomery_ADD(const Point& P, const Point& Q, const Point& R) {
	mpz_class A, B, C, D, DA, CB, temp1, temp2;
	Point result;
	add_fp(Q.X, Q.Z, &A);
	sub_fp(Q.X, Q.Z, &B);
	add_fp(P.X, P.Z, &C);
	sub_fp(P.X, P.Z, &D);
	mul_fp(D, A, &DA);
	mul_fp(B, C, &CB);
	add_fp(DA, CB, &temp1);
	sqr_fp(temp1, &temp1);
	mul_fp(R.Z, temp1, &result.X);

	sub_fp(DA, CB, &temp2);
	sqr_fp(temp2, &temp2);
	mul_fp(R.X, temp2, &result.Z);

	return result;
}

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#doubling-dbl-1987-m-2
Point xDBL(const Point& P, const Point& Ap24_C24)
{
	Point result;
	mpz_class t0, t1, t2, t3;
	sub_fp(P.X, P.Z, &t0);
	add_fp(P.X, P.Z, &t1);
	sqr_fp(t0, &t0);

	sqr_fp(t1, &t1);
	mul_fp(Ap24_C24.Z, t0, &result.Z);
	mul_fp(result.Z, t1, &result.X);

	sub_fp(t1, t0, &t1);
	mul_fp(Ap24_C24.X, t1, &t0);
	add_fp(result.Z, t0, &result.Z);

	mul_fp(result.Z, t1, &result.Z);
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
	sqr_fp(t1,  &DBLout.Z);

	mul_fp(t1, ADDout.X, &t1);
	sub_fp(DBLout.X, DBLout.Z,  &t2);
	mul_fp(DBLout.X, DBLout.Z,  &DBLout.X);
	mul_fp(Ap24_1.X, t2, &ADDout.X);
	sub_fp(t0, t1, &ADDout.Z);
	add_fp(DBLout.Z, ADDout.X, &DBLout.Z);
	add_fp(t0, t1, &ADDout.X);

	mul_fp(DBLout.Z, t2, &DBLout.Z);
	sqr_fp(ADDout.Z ,  &ADDout.Z);
	sqr_fp(ADDout.X,&ADDout.X);
	mul_fp(ADDout.Z, R.X, &ADDout.Z);
	mul_fp(ADDout.X, R.Z, &ADDout.X);
}

//モンゴメリ曲線のスカラー倍
Point xMUL(const Point& P, const Point& A_1, const mpz_class& n) {
	Point x0, x1, Ap24_1, Point_temp1, Point_temp2;
	mpz_class div_ans;

	x0 = P;
	//a24の計算
	Ap24_1.Z=1;
	add_fp(A_1.X, 2, &Ap24_1.X);
	
	//div_fp(Ap24_1.X, 4, &div_ans);
	mul_fp(Ap24_1.X, para.inv4, &div_ans);
	Ap24_1.X = div_ans;

	x1=xDBL(P, Ap24_1);
	
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
	mpz_class result,temp1,temp2,temp3;
	add_fp(x, a, &temp1);	// x+a
	mul_fp(temp1, x, &temp2); // x^2 + ax
	add_fp(temp2, 1, &temp3); // x^2 + ax + 1
	mul_fp(temp3, x, &result); // x^3 + ax^2 + x
	return result;
}

void Evaluation(const Point& P, const Point& Q,
	mpz_class& Xout, mpz_class& Zout)
{
	Point temp;
	sqr_fp(Q.X, &temp.X);
	sqr_fp(Q.Z, &temp.Z);
	mul_fp(P.X, temp.X, &Xout);
	mul_fp(P.Z, temp.Z, &Zout);
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
void IsogenyCalc( const Point& A, const Point& P, const Point& K, const size_t& k,
	Point* Aout, Point* Pout) {

	mpz_class tmp0, tmp1, tmp2, Psum, Pdif;
	Point Q, ed, prod, Ap24_C24;

	add_fp(A.X, 2, &Ap24_C24.X);
	Ap24_C24.Z = 4;

	add_fp(P.X, P.Z, &Psum);
	sub_fp(P.X, P.Z, &Pdif);

	sub_fp(K.X, K.Z, &prod.X);
	add_fp(K.X, K.Z, &prod.Z);

	mul_fp(prod.X, Psum,  &tmp1);
	mul_fp(prod.Z, Pdif,  &tmp0);
	add_fp(tmp0, tmp1,  &Q.X);
	sub_fp(tmp0, tmp1, &Q.Z);

	Point M[3]={K};

	M[1]=xDBL(K, Ap24_C24);

	for (int i = 1; i < (k >>1); ++i) {

		if (i >= 2)
			M[i % 3]=Montgomery_ADD(M[(i - 1) % 3], K, M[(i - 2) % 3]);

		sub_fp(M[i % 3].X, M[i % 3].Z, &tmp1);
		add_fp(M[i % 3].X, M[i % 3].Z, &tmp0);
		mul_fp(prod.X,tmp1,  &prod.X);
		mul_fp(prod.Z, tmp0, &prod.Z);
		mul_fp(tmp1, Psum, &tmp1);
		mul_fp(tmp0, Pdif, &tmp0);
		add_fp(tmp0, tmp1,  &tmp2);
		mul_fp(Q.X, tmp2,  &Q.X);
		sub_fp(tmp0, tmp1,  &tmp2);
		mul_fp(Q.Z, tmp2,  &Q.Z);

	}

	Evaluation(P, Q, Pout->X, Pout->Z);

	//A faster way to the CSIDH p10くらい
	ed=mont2ed(A);
	
	pow_fp(ed.X, k, &ed.X);
	pow_fp(ed.Z, k, &ed.Z);

	// compute prod.x^8, prod.z^8
	pow_fp(prod.X, 8, &prod.X);
	pow_fp(prod.Z, 8, &prod.Z);

	// compute image curve parameters
	mul_fp(ed.Z, prod.X,  &ed.Z);
	mul_fp(ed.X, prod.Z,  &ed.X);

	// compute Montgomery params
	add_fp(ed.X, ed.Z,  &Aout->X);
	sub_fp(ed.X, ed.Z,  &Aout->Z);
	add_fp(Aout->X, Aout->X, &Aout->X);
	
}