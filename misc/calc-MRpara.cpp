#include <iostream>
#include <gmpxx.h>
#include <chrono>

using namespace std;

struct mon {
	size_t nbit;
	mpz_class R, R2, nr;
};

void mont(const mpz_class& mod, mon& para)
{
	mpz_class t, vi;
	para.nr = 0;
	t = 0;
	vi = 1;

	para.nbit = mpz_sizeinbase(mod.get_mpz_t(), 2);
	//para.nbit = 256;

	para.R2 = (mpz_class(1) << (para.nbit * 2)) % mod;
	para.R = (mpz_class(1) << para.nbit) - 1;

	for (int i = 0; i < para.nbit; i++) {
		if ((t & 1) == 0) {
			t += mod;
			para.nr += vi;
		}
		t >>= 1;
		vi <<= 1;
	}
}

mpz_class MR(const mpz_class& t, const mpz_class& mod, const mon& para)
{
	mpz_class c;
	c = t * para.nr;
	c &= para.R;
	c *= mod;
	cout << "c:" << c << endl;
	uint64_t buf[16];
	mpz_export(&buf, NULL, -1, 16, 0, 0, c.get_mpz_t());
	for (int i = 0; i < 16; i++) {
		std::cout << i << ":" << buf[i] << std::endl;
	}
	c += t;
	c >>= para.nbit;

	if (c >= mod) c -= mod;
	return c;
}

mpz_class MRmul(const mpz_class& a, const mpz_class& b, const mpz_class& mod, const mon& para)
{
	mpz_class test;
	test = MR((a * b), mod, para);
	test *= para.R2;
	return MR(test, mod, para);
}

int main()
{
	mpz_class a, b, mod, result;
	mon para;

	a = "0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc";
	b = "0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b";
	mod = "0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed";

	mont(mod, para);
	//result = MRmul(a, b, mod, para);
	//cout << "result:" << result << endl;

	cout << "nbit:"<< para.nbit << "," << endl;
	cout << "nr:"<< para.nr << "," << endl;
	cout << "R:"<< para.R << "," << endl;
	cout << "R2:"<< para.R2 << "," << endl;

	//cout << "2*R2:"<< MR(2 * para.R2, mod, para) << endl;
	//cout << MR(3*para.R2, mod, para) << endl;
	//cout << "4*R2:"<< MR(4 * para.R2, mod, para) << endl;
	cout << "MR2:"<< MR(para.R2, mod, para) << endl;

}