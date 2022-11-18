#pragma once
#include "fp.hpp"

const size_t l = 74;

const size_t primes[l] = {
	3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131,
137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,
277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
359, 367, 373, 587 };

struct seckey {
	int e[l];
};

//MontgomeryCurve関連
Point xMUL(const Point& P, const Point& A_1, const mpz_class& n);
Fp calc_twist(const Fp& a, const Fp& x);
void IsogenyCalc(const Point& A, const Point& P, const Point& K, const size_t& k, Point* Aout, Point* Pout);

//CSIDH関連
bool validate(const mpz_class& a);
mpz_class action(const mpz_class& A, const seckey& Key);
