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
