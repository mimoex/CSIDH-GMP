#include "fp.h"

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#diffadd-dadd-1987-m-3
Point Montgomery_ADD(const Point& P, const Point& Q, const Point& R, const mpz_class& mod) {
	mpz_class A, B, C, D, DA, CB, temp1, temp2;
	Point result;
	add_fp(Q.X, Q.Z, mod, &A);
	sub_fp(Q.X, Q.Z, mod, &B);
	add_fp(P.X, P.Z, mod, &C);
	sub_fp(P.X, P.Z, mod, &D);
	mul_fp(D, A, mod, &DA);
	mul_fp(B, C, mod, &CB);
	add_fp(DA, CB, mod, &temp1);
	sqr_fp(temp1, mod, &temp1);
	mul_fp(R.Z, temp1, mod, &result.X);

	sub_fp(DA, CB, mod, &temp2);
	sqr_fp(temp2, mod, &temp2);
	mul_fp(R.X, temp2, mod, &result.Z);

	return result;
}

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#doubling-dbl-1987-m-2
Point xDBL(const Point& P, const Point& Ap24_C24, const mpz_class& mod)
{
	Point result;
	mpz_class t0, t1, t2, t3;
	sub_fp(P.X, P.Z,mod,&t0);
	add_fp(P.X, P.Z,mod,&t1);
	sqr_fp(t0,mod,&t0);

	sqr_fp(t1, mod, &t1);
	mul_fp(Ap24_C24.Z, t0, mod, &result.Z);
	mul_fp(result.Z, t1, mod, &result.X);

	sub_fp(t1, t0,mod,&t1);
	mul_fp(Ap24_C24.X, t1, mod, &t0);
	add_fp(result.Z, t0, mod, &result.Z);

	mul_fp(result.Z, t1, mod, &result.Z);
	return result;
}

void xDBLADD(const Point& P, const Point& Q, const Point& R, const Point& Ap24_1, const mpz_class& mod,
	Point& DBLout, Point& ADDout) {
	mpz_class t0, t1, t2, t3;
	add_fp(P.X, P.Z, mod, &t0);
	sub_fp(P.X, P.Z, mod, &t1);
	sqr_fp(t0, mod, &DBLout.X);
	sub_fp(Q.X, Q.Z, mod, &t2);
	add_fp(Q.X, Q.Z, mod, &ADDout.X);
	mul_fp(t0, t2, mod, &t0);
	sqr_fp(t1,mod, &DBLout.Z);

	mul_fp(t1, ADDout.X, mod, &t1);
	sub_fp(DBLout.X, DBLout.Z,mod, &t2);
	mul_fp(DBLout.X, DBLout.Z,mod, &DBLout.X);
	mul_fp(Ap24_1.X, t2, mod, &ADDout.X);
	sub_fp(t0, t1, mod, &ADDout.Z);
	add_fp(DBLout.Z, ADDout.X,mod,&DBLout.Z);
	add_fp(t0, t1, mod, &ADDout.X);

	mul_fp(DBLout.Z, t2, mod, &DBLout.Z);
	sqr_fp(ADDout.Z ,mod, &ADDout.Z);
	sqr_fp(ADDout.X, mod,&ADDout.X);
	mul_fp(ADDout.Z, R.X, mod, &ADDout.Z);
	mul_fp(ADDout.X, R.Z, mod, &ADDout.X);
}

//モンゴメリ曲線のスカラー倍
Point xMUL(const Point& P, const Point& A_1, const mpz_class& n, const mpz_class& mod) {
	Point x0, x1, Ap24_1, Point_temp1, Point_temp2;
	mpz_class div_ans;

	x0 = P;
	//a24の計算
	Ap24_1.Z=1;
	add_fp(A_1.X, 2, mod, &Ap24_1.X);
	div_fp(Ap24_1.X, 4,mod,&div_ans);
	Ap24_1.X = div_ans;

	x1=xDBL(P, Ap24_1,mod);

	string bit;
	size_t bit_size;
	bit = n.get_str(2);
	bit_size = bit.size();

	for (int i = bit_size - 2; i >= 0; i--) {
		if (mpz_tstbit(n.get_mpz_t(), i) == 0) xDBLADD(x0, x1, P, Ap24_1, mod, x0, x1);
	else xDBLADD(x1, x0, P, Ap24_1, mod, x1, x0);
	}
	return x0;
}

//モンゴメリ曲線右辺の計算
mpz_class calc_twist(const mpz_class& a, const mpz_class& x,const mpz_class& mod) {
	mpz_class result,temp1,temp2,temp3;
	add_fp(x, a, mod, &temp1);	// x+a
	mul_fp(temp1, x, mod, &temp2); // x^2 + ax
	add_fp(temp2, 1, mod, &temp3); // x^2 + ax + 1
	mul_fp(temp3, x, mod, &result); // x^3 + ax^2 + x
	return result;
}

void Evaluation(const Point& P, const Point& Q, const mpz_class& mod,
	mpz_class& Xout, mpz_class& Zout)
{
	Point temp;
	sqr_fp(Q.X, mod, &temp.X);
	sqr_fp(Q.Z, mod, &temp.Z);
	mul_fp(P.X, temp.X, mod, &Xout);
	mul_fp(P.Z, temp.Z, mod, &Zout);
}

Point mont2ed(const Point& P, const mpz_class& mod)
{
	// compute twisted Edwards curve coefficients
	/* A = 2*(1+d)/(1-d) */
	Point ed;
	add_fp(P.Z, P.Z, mod, &ed.Z);
	add_fp(P.X, ed.Z, mod, &ed.X);
	sub_fp(P.X, ed.Z, mod, &ed.Z);
	return ed;
}


//Isogeny
void IsogenyCalc( const Point& A, const Point& P, const Point& K, const mpz_class& mod, const size_t& k, Point* Aout, Point* Pout) {

	mpz_class tmp0, tmp1, tmp2, Psum, Pdif;
	Point Q, ed, prod, Ap24_C24;

	add_fp(A.X, 2, mod, &Ap24_C24.X);
	Ap24_C24.Z = 4;

	add_fp(P.X, P.Z, mod, &Psum);
	sub_fp(P.X, P.Z, mod, &Pdif);

	sub_fp(K.X, K.Z, mod, &prod.X);
	add_fp(K.X, K.Z, mod, &prod.Z);

	mul_fp(prod.X, Psum,mod, &tmp1);
	mul_fp(prod.Z, Pdif,mod, &tmp0);
	add_fp(tmp0, tmp1,mod, &Q.X);
	sub_fp(tmp0, tmp1, mod, &Q.Z);

	Point M[3]={K};

	M[1]=xDBL(K, Ap24_C24,mod);

	for (int i = 1; i < (k >>1); ++i) {

		if (i >= 2)
			M[i % 3]=Montgomery_ADD(M[(i - 1) % 3], K, M[(i - 2) % 3], mod);

		sub_fp(M[i % 3].X, M[i % 3].Z, mod, &tmp1);
		add_fp(M[i % 3].X, M[i % 3].Z, mod, &tmp0);
		mul_fp(prod.X,tmp1,mod, &prod.X);
		mul_fp(prod.Z, tmp0, mod, &prod.Z);
		mul_fp(tmp1, Psum, mod, &tmp1);
		mul_fp(tmp0, Pdif, mod, &tmp0);
		add_fp(tmp0, tmp1,mod, &tmp2);
		mul_fp(Q.X, tmp2,mod, &Q.X);
		sub_fp(tmp0, tmp1,mod, &tmp2);
		mul_fp(Q.Z, tmp2,mod, &Q.Z);

	}

	Evaluation(P, Q, mod, Pout->X, Pout->Z);

	//A faster way to the CSIDH p10くらい
	ed=mont2ed(A, mod);

	pow_fp(ed.X, k, mod, &ed.X);
	pow_fp(ed.Z, k, mod, &ed.Z);

	// compute prod.x^8, prod.z^8
	pow_fp(prod.X, 8, mod, &prod.X);
	pow_fp(prod.Z, 8, mod, &prod.Z);

	// compute image curve parameters
	mul_fp(ed.Z, prod.X,mod, &ed.Z);
	mul_fp(ed.X, prod.Z,mod, &ed.X);

	// compute Montgomery params
	add_fp(ed.X, ed.Z,mod, &Aout->X);
	sub_fp(ed.X, ed.Z,mod, &Aout->Z);
	add_fp(Aout->X, Aout->X, mod, &Aout->X);
	
}