#pragma once
#include <iostream>
#include <gmpxx.h>
#include <random>
#include <chrono>

using namespace std;

const size_t N = 74;	//log(11^74)/log(2) ≒ 255.9979

const size_t primes[N] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
						67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131,
						137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
						199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,
						277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
						359, 367, 373, 587 };

struct Point {
	mpz_class X, Y, Z;	//inf is (0, 0, 1)
};

struct seckey {
	int e[N];
};

//montgomery剰余乗算
struct mon {
	const size_t nbit;
	const mpz_class R;
	const mpz_class R2;
	const mpz_class nr;
};

//有限体の加算	c=a+b %p
void add_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//有限体の減算	c=a-b %p
void sub_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//有限体の乗算	c=a*b %p
void mul_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//有限体のべき乗	c=a^b %p
void pow_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//有限体の逆元	c=a^-1
void inv_fp(const mpz_class& a, const mpz_class& p, mpz_class* c);
//有限体の除算	c=a/b %p
void div_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//有限体の2乗	c=a^2 %p
void sqr_fp(const mpz_class& a, const mpz_class& p, mpz_class* c);

//p以下の乱数を生成
mpz_class random_fp(const mpz_class& p);


Point xMUL(const Point& P, const Point& A_1, const mpz_class& k, const mpz_class& mod);

mpz_class calc_twist(const mpz_class& a, const mpz_class& x, const mpz_class& mod);

void IsogenyCalc(const Point& A, const Point& P, const Point& K, const mpz_class& mod, const size_t& k, Point* Aout, Point* Pout);

void genCSIDHkey(seckey* K);


bool validate(const mpz_class& a, const mpz_class& mod);

mpz_class action(const mpz_class& A, const seckey& Key, const mpz_class& mod);