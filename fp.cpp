#include "fp.h"

void add_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	*c = a + b;
	if (*c >= p) *c -= p;
}

void sub_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	if (a >= b) *c = a - b;
	else *c = a + p - b;
}


void mul_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	*c = a * b;
	*c %= p;
}


//Montgomery multiplicatuon
struct mon para = {
	511,
	mpz_class("6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042047") ,
	mpz_class("1032312334021767149891928024133341961156775772779044637692888609886692159219829053969396827895169639458311379179211658124290421880650098963926249002436041"),
	mpz_class("4648943738516491079755353375009763753751649870163419510591469415438864277366534547184950551625493863052632292052244681369519164709307558804992291140151629")

};

mpz_class MR(const mpz_class& t, const mpz_class& mod)
{
	mpz_class c;
	c = t * para.nr;
	c &= para.R;
	c *= mod;
	c += t;
	c >>= para.nbit;
	if (c >= mod) c -= mod;
	
	return c;
}

void mon_mul_fp(const mpz_class& a, const mpz_class& b, const mpz_class& mod,
	mpz_class* c)
{
	//mon para;
	mpz_class test;
	test = MR((a * b), mod);
	test *= para.R2;
	*c= MR(test, mod);
}


void sqr_fp(const mpz_class& a, const mpz_class& p, mpz_class* c)
{
	*c = a * a;
	*c %= p;
}


void pow_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	mpz_powm(c->get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), p.get_mpz_t());
}

//逆元
void inv_fp(const mpz_class& a, const mpz_class& p, mpz_class* c)
{
	if (a == 0) {
		throw std::range_error("Divided by zero.");

	}
	else {
		pow_fp(a, p - 2, p, c);
	}

}


void div_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	inv_fp(b, p, c);
	mul_fp(a, *c, p, c);
}

//有限体pから乱数生成
mpz_class random_fp(const mpz_class& mod)
{
	mpz_class x, cnt;
	size_t n;
	string bit;

	bit = mod.get_str(2);
	n = bit.size();

	random_device rnd;
	gmp_randclass r(gmp_randinit_default);
	r.seed(rnd());

	while (1) {
		x = r.get_z_bits(n);
		cnt++;

		if (x < mod) {
			return x;
		}
	}
}
