#ifndef UNLINEAR
#define UNLINEAR
#include "function.hpp"
#include <vector>
#include "rand.hpp"
#include <iostream>
using namespace std;

class Solver{
private:
	double (*f) (int, Matrix&);
	int m, n; // n为变量数目, m为方程数目.
	double eps;
public:
	Solver(int M, int N, double(*func) (int, Matrix&)) :m(M), n(N), f(func), eps(1e-8) {}

	Solver(const Solver &b) {
		f = b.f, n = b.n, m = b.m;
	}

	Solver& operator =(const Solver &b) {
		if (this == (&b)) return *this;
		f = b.f, n = b.n, m = b.m;
		return *this;
	}
	
	inline double partial(int i, int j, Matrix x) {
		double f1, f2;
		x[j] += eps;
		f1 = f(i, x);
		x[j] -= 2. * eps;
		f2 = f(i, x);
		return (f1 - f2) / (2. * eps);
	}
	
	Matrix gradient(int i, const Matrix& x) {
		Matrix c(n, 1);
		for (int j = 1; j <= n; j++) {
			c[j] = partial(i, j, x);
		}
		return c;
	}

	Matrix jacobian(const Matrix& x) {
		Matrix c(m, n);
		for (int i = 1; i <= m; i++) {
			Matrix g = gradient(i, x);
			for (int j = 1; j <= n; j++) c(i, j) = g[j];
		}
		return c;
	}

	Matrix eval(Matrix& x) {
		Matrix c(m, 1);
		for (int i = 1; i <= m; i++) c[i] = f(i, x);
		return c;
	}

	Matrix newton() {
		if (n != m) throw Exception("I can't do it!");
		Randmat rm1;
		Matrix x(rm1.unirand(n, 1));
		Matrix y(rm1.unirand(n, 1));
		while (y.norm1() > 1e-8) {
			y = eval(x);
			x -= (jacobian(x).inv() * y);
		}
		return x;
	 }

	Matrix broydens() {
		if (n != m) throw Exception("I can't do it!");
		Randmat rm1;
		Matrix x(rm1.unirand(n, 1)), s(n, 1), y(n, 1);
		Matrix B(rm1.unirand(n, n));
		while ((y = eval(x)).norm1() > 1e-8) {
			s = Matrix(n, 1) - B / y;
			x += s;
			y = eval(x) - y;
			B += ((y - B * s) * s.transpose()) / (s.transpose() * s)(1, 1);
		}
		return x;
	}
};

#endif