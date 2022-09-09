#include "fp.h"

//Montgomery multiplicatuon
mon para = {
	mpz_class("5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659"),
	mpz_class("1331684699081905773686966904488651388517342873708180584403111660513502390006644134406723028256595313406156735410987361198165720310405343322235720072016415"),
	mpz_class("72984510660328627668487265992516213834621203782665312994583680849552234545871"),
	511,
	mpz_class("6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042047"),
	mpz_class("1032312334021767149891928024133341961156775772779044637692888609886692159219829053969396827895169639458311379179211658124290421880650098963926249002436041"),
	mpz_class("4648943738516491079755353375009763753751649870163419510591469415438864277366534547184950551625493863052632292052244681369519164709307558804992291140151629"),
	mpz_class("1377165168643675455039144881148317509670311415463474351249334079806872455010196950774045036057070460220388987449293580634214060164351911684273944214976389")

};


void add_fp(const mpz_class& a, const mpz_class& b, mpz_class* c)
{
	*c = a + b;
	if (*c >= para.mod) *c -= para.mod;
}

void sub_fp(const mpz_class& a, const mpz_class& b, mpz_class* c)
{
	if (a >= b) *c = a - b;
	else *c = a + para.mod - b;
}


mpz_class MR(const mpz_class& t)
{
	mpz_class c;
	c = (t & para.R) * para.nr;
	c &= para.R;
	c *= para.mod;
	c += t;
	c >>= para.nbit;
	if (c >= para.mod) c -= para.mod;
	
	return c;
}

void mul_fp(const mpz_class& a, const mpz_class& b,
	mpz_class* c)
{
	mpz_class test;
	test = MR(a * b) * para.R2;
	*c= MR(test);
}


void sqr_fp(const mpz_class& a, mpz_class* c)
{
	mpz_class test;
	test = MR(a * a) * para.R2;
	*c = MR(test);
}


void pow_fp(const mpz_class& a, const mpz_class& b, mpz_class* c)
{
	mpz_class x;
	x = MR(a * para.R2);

	mpz_class result = MR(para.R2);

	for (int i = para.nbit - 1; i >= 0; i--) {
		result = MR(result * result);
		if (mpz_tstbit(b.get_mpz_t(), i) == 1) result = MR(result * x);
	}
	*c = MR(result);
}

void mul_mon(const mpz_class& a, const mpz_class& b,
	mpz_class* c)
{
	*c = MR(a * b);
}

void pow_mon(const mpz_class& a, const mpz_class& b, mpz_class* c)
{
	mpz_class result = MR(para.R2);

	for (int i = para.nbit - 1; i >= 0; i--) {
		result = MR(result * result);
		if (mpz_tstbit(b.get_mpz_t(), i) == 1) result = MR(result * a);
	}
	*c = result;
}

//逆元
void inv_fp(const mpz_class& a, mpz_class* c)
{
	if (a == 0) {
		throw std::range_error("Divided by zero.");

	}
	else {
		pow_fp(a, para.mod - 2, c);
	}

}


void div_fp(const mpz_class& a, const mpz_class& b, mpz_class* c)
{
	inv_fp(b, c);
	mul_fp(a, *c, c);
}

//有限体pから乱数生成
mpz_class random_fp()
{
	mpz_class x;

	random_device rnd;
	gmp_randclass r(gmp_randinit_default);
	r.seed(rnd());

	while (1) {
		x = r.get_z_bits(para.nbit);

		if (x < para.mod) {
			return x;
		}
	}
}
