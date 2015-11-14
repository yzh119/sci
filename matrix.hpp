#ifndef MATRIX
#define MATRIX
#include <fstream>
#include "exception.hpp"
#include "complex.hpp"
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>

const double eps = 1e-9;

class Matrix{
private:
	double **p;
	int n, m;

	inline double sqr(double t) const{
		return t * t;
	}

	inline double sgn(double t) const{
		return t > eps ? 1.: (t < eps ? -1.: 0.);
	}

	void caleig(Matrix& H, std::vector<comp>& eig) {
		if ((int)eig.size() != n) throw Exception("Initialization error!");
		for (int i = 1; i <= n; i++) {
			if (i == n) {eig[n - 1] = (comp)H(i, i); break;}
			if (fabs(H(i + 1, i)) < eps) {eig[i - 1] = (comp)H(i, i); continue;}
			std::pair <comp, comp> eigpair = caleig22(H(i, i), H(i, i + 1), H(i + 1, i), H(i + 1, i + 1));
			eig[i - 1] = eigpair.first;
			eig[i++] = eigpair.second;
		}
	}

public:
	Matrix(int N, int M): n(N), m(M) {
		if (N <= 0 || M <= 0) throw ("Out of Range!");
		p = new double* [n];
		for (int i = 0; i < n; i++) {
			p[i] = new double [m];
			for (int j = 0; j < m; j++)
				p[i][j] = 0.;
		}
	}

	Matrix(int N, char c = 'N'): n(N) {
		if (N <= 0) throw ("Out of Range!");
		m = n;
		p = new double* [n];
		for (int i = 0; i < n; i++) {
			p[i] = new double [n];
			for (int j = 0; j < n; j++)
				p[i][j] = 0.;
		}
		if (c == 'I') {
			for (int i = 0; i < n; i++) p[i][i] = 1.;
		}
	}

	Matrix(const Matrix &b) {
		n = b.n; m = b.m;
		p = new double* [n];
		for (int i = 0; i < n; i++) {
			p[i] = new double [m];
			for (int j = 0; j < m; j++)
				p[i][j] = b.p[i][j];
		}
	}

	Matrix& operator =(const Matrix &b) {
		if (this == &b) return *this;
		for (int i = 0; i < n; i++) delete [] p[i];
		delete [] p;
		n = b.n; m = b.m;
		p = new double* [n];
		for (int i = 0; i < n; i++) {
			p[i] = new double [m];
			for (int j = 0; j < m; j++)
				p[i][j] = b.p[i][j];
		}
		return *this;
	}

	inline double& operator ()(int i, int j) {
		if (i > n || i < 1 || j > m || j < 1) throw Exception("Out of Range!");
		return p[i - 1][j - 1];
	}

	inline Matrix operator ()(int i) {
		if (i < 1 || i > n) throw Exception("Out of Range!");
		Matrix c(1, m);
		for (int j = 0; j < m; j++) c.p[0][j] = p[i - 1][j];
		return c;
	}

	inline Matrix operator ()(char star, int j) {
		if (star != '*') throw Exception("star must be a *!");
		if (j < 1 || j > m) throw Exception("Out of Range!");
		Matrix c(n, 1);
		for (int i = 0; i < n; i++) c.p[i][0] = p[i][j - 1];
		return c;
	}

	Matrix sub(int i1, int j1, int i2, int j2) const{
		if (i2 < i1 || j2 < j1 || i1 < 1 || i2 > n || j1 < 1 || j2 > m) throw Exception("Out of Range!");
		Matrix c(i2 - i1 + 1, j2 - j1 + 1);
		for (int i = 0; i < i2 - i1 + 1; i++)
			for (int j = 0; j < j2 - j1 + 1; j++)
				c.p[i][j] = p[i + i1 - 1][j + j1 - 1];
		return c;
	}

	friend Matrix operator *(const Matrix&, const Matrix&);
	friend Matrix operator *(double, const Matrix&);
	friend Matrix operator *(const Matrix&, double);
	friend Matrix operator +(const Matrix&, const Matrix&);
	friend Matrix operator -(const Matrix&, const Matrix&);
	friend Matrix operator ^(Matrix, int);
	friend Matrix operator /(Matrix, Matrix);
	friend Matrix operator /(const Matrix&, double);
	friend Matrix& operator +=(Matrix&, const Matrix&);
	friend Matrix& operator -=(Matrix&, const Matrix&);
	friend Matrix& operator *=(Matrix&, double);
	friend Matrix& operator *=(Matrix&, const Matrix&);
	friend Matrix& operator /=(Matrix&, const Matrix&);
	friend Matrix& operator /=(Matrix&, double);
	friend std::ostream& operator<< (std::ostream&, const Matrix&);

	inline Matrix transpose() {
		Matrix c(m, n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				c.p[j][i] = p[i][j];
		return c;
	}

	Matrix ladder() {
		Matrix c(*this);
		int k = -1;
		for (int i = 0; i < n; i++) {
			bool flag = true;
			do {
				flag = true;
				k++;
				for (int j = i; j < n; j++) {
					if (fabs(c.p[j][k]) > eps) {flag = false; break;}
				}
			} while (k < m && flag);

			if (k >= m) break;

			int cmp = i;
			for (int j = i + 1; j < n; j++) {
				if (fabs(c.p[j][k]) > fabs(c.p[cmp][k])) cmp = j;
			}
			for (int l = k; l < m; l++) std::swap(c.p[i][l], c.p[cmp][l]);

			for (int j = i + 1; j < n; j++) {
				double coef = -c.p[j][k] / c.p[i][k];
				for (int l = k; l < m; l++) {
					c.p[j][l] += coef * c.p[i][l];
				}
			}

		}
		return c;
	}

	Matrix pseudoinv() {
		int r = rank();
		Matrix AT = transpose();
		if (r == n) return AT * ((*this) * AT).inv();
		if (r == m) return (AT * (*this)).inv() * AT;
		std::pair <Matrix, std::pair<Matrix, Matrix> > SUV = SVD();
		Matrix S1(SUV.first.transpose());
		for (int i = 0; i < std::min(S1.n, S1.m); i++) {
			if (fabs(S1.p[i][i]) > eps)
				S1.p[i][i] = 1. / S1.p[i][i];
		}
		return SUV.second.second * S1 * SUV.second.first.transpose();
	}

	Matrix inv() {
		if (n != m) throw("No inverse exists");
		Matrix A(*this), c(n, 'I');
		for (int i = 0; i < n; i++) {
			if (fabs(A.p[i][i]) < eps) throw("No inverse exists");
			double coef = 1. / A.p[i][i];
			for (int j = i; j < n; j++)
			 	A.p[i][j] *= coef;
			for (int j = 0; j < n; j++)
				c.p[i][j] *= coef;
			for (int k = 0; k < n; k++) {
				if (k == i) continue;
				coef = -A.p[k][i];
				for (int j = i; j < n; j++) {
					A.p[k][j] += coef * A.p[i][j];
				}
				for (int j = 0; j < n; j++) {
					c.p[k][j] += coef * c.p[i][j];
				}
			}
		}
		return c;
	}

// LU returns LU decomposition of A, in which L * U = A;
	std::pair <Matrix, Matrix> LU() const{
		if (n != m) throw Exception("Matrix is not a square!");
		Matrix L(n, 'I'), U(*this);
		for (int i = 1; i <= n; i++) {
			for (int j = i + 1; j <= n; j++) {
				double coef = -U(j, i) / U(i, i);
				L(j, i) = -coef;
				for (int k = i; k <= n; k++)
					U(j, k) += coef * U(i, k);
			}
		}
		return std::make_pair(L, U);
	}

// LUP returns make_pair(P, make_pair(L, U)), in which P * A = L * U;
	std::pair <Matrix, std::pair<Matrix, Matrix> > LUP() const{
		if (n != m) throw Exception("Matrix is not a square!");
		Matrix L(n, 'I'), U(*this), P(n, 'I');



		return std::make_pair(P, std::make_pair(L, U));
	}

/* QR returns QR decomposition of A,
 * make_pair(Q, R), in which Q * R = A;
 */
	std::pair <Matrix, Matrix> QR() const{
		if (n < m) throw Exception("n should be greater than m!");
		Matrix Q(n, 'I'), R(*this);
		for (int k = 1; k <= n; k++) {
			double alpha = -sgn(R(k, k)) * (R.sub(k, k, m, k).norm2());
			Matrix v(m, 1); v(k, 1) -= alpha;
			for (int j = k; j <= m; j++) v(j, 1) += R(j, k);
			double beta = (v.transpose() * v)(1, 1);
			R -= 2. * v * (v.transpose() * R)/ beta;
			Q -= 2. * (Q * v) * v.transpose()/ beta;
		}
		return std::make_pair(Q, R);
	}

/* SVD returns singular value decomposition of A,
 * make_pair(Sigma, make_pair(U, V)), in which U * Sigma * V' = A;
 */
	std::pair <Matrix, std::pair<Matrix, Matrix> > SVD() const{



	}

// arnoldi returns make_pair(Q, H)
	std::pair<Matrix, Matrix> arnoldi() {
		if (n != m) throw Exception("Matrix size doesn't match!");
		Matrix Q(n, n), H(n, n);
		Matrix x(n, 1), u(n, 1);
		for (int i = 1; i <= n; i++) x(i, 1) = 1. * rand() / RAND_MAX;
		double norm  = x.norm2();
		for (int i = 1; i <= n; i++) Q(i, 1) = x(i, 1) / norm;

		for (int k = 1; k <= n; k++) {
			u = (*this) * Q('*', k);
    	for (int j = 1; j <= k; j++) {
				H(j, k) = (Q('*', j).transpose() * u)(1, 1);
				u -= H(j, k) * Q('*', j);
			}
			if (k != n) {
				H(k + 1, k) = u.norm2();
				for (int j = 1; j <= n; j++) {
					Q(j, k + 1) = u(j, 1) / H(k + 1, k);
				}
			}
		}
		return std::make_pair(Q, H);
	}

	std::pair<Matrix, Matrix> lanczos() {
		if (n != m) throw Exception("Matrix size doesn't match!");
		Matrix Q(n, n), H(n, n);
		Matrix x(n, 1), u(n, 1);
		for (int i = 1; i <= n; i++) x(i, 1) = 1. * rand() / RAND_MAX;
		double norm = x.norm2();
		for (int i = 1; i <= n; i++) Q(i, 1) = x(i, 1) / norm;

		for (int k = 1; k <= n; k++) {
			u = (*this) * Q('*', k);
			H(k, k) = (Q('*', k).transpose() * u)(1, 1);
			if (k > 1) u -= (H(k - 1, k) * Q('*', k - 1) + H(k, k) * Q('*', k));
				else u -= H(k, k) * Q('*', k);
			if (k != n) {
				H(k + 1, k) = H(k, k + 1) = u.norm2();
				for (int j = 1; j <= n; j++) {
					Q(j, k + 1) = u(j, 1) / H(k + 1, k);
				}
			}
		}
		return std::make_pair(Q, H);
	}

	std::vector<comp> eigenvalue() {
		if (n != m) throw Exception("Matrix size doesn't match!");
		std::pair <Matrix, Matrix> QH = arnoldi();
		Matrix Q1(QH.first), H(QH.second);
		std::vector <double> s, c;
		std::vector <comp> Eig(n), Eig1(n);
		caleig(H, Eig);

		while (1) {
			double err = 0;
			for (int i = 0; i < n; i++) {
				if ((Eig[i] - Eig1[i]).modulo() > err) err = (Eig[i] - Eig1[i]).modulo();
			}
			if (err < eps) break;
			s.clear(); c.clear();

			for (int i = 2; i <= n; i++) {
					double r = sqrt(H(i, i - 1) * H(i, i - 1) + H(i - 1, i - 1) * H(i - 1, i - 1));
					c.push_back(H(i - 1, i - 1) / r);
					s.push_back(H(i, i - 1) / r);
					for (int j = i - 1; j <= n; j++) {
						double hi = H(i, j), hi_1 = H(i - 1, j);
						H(i, j) = -s.back() * hi_1 + c.back() * hi;
						H(i - 1, j) = c.back() * hi_1 + s.back() * hi;
					}
			}

			for (int j = 2; j <= n; j++) {
					for (int i = 1; i <= n; i++) {
						double hj = H(i, j), hj_1 = H(i, j - 1);
						H(i, j) = -s[j - 2] * hj_1 + c[j - 2] * hj;
						H(i, j - 1) = c[j - 2] * hj_1 + s[j - 2] * hj;
					}
			}

			Eig1 = Eig;
			caleig(H, Eig);
		}

		return Eig;
	}

	inline std::pair <int, int> sz() const {
		return std::make_pair(n, m);
	}

	inline bool isSquare() {
		return n == m;
	}

	inline double elem(int x) const {
		if (x < 1 || x > m) throw Exception("Out of Range!");
		return p[0][x - 1];
	}

	inline double norm1() {
		double ret = 0.;
		if (n == 1) {
			for (int i = 0; i < m; i++) ret += fabs(p[0][i]);
			return ret;
		}
		if (m == 1) {
			for (int i = 0; i < n; i++) ret += fabs(p[i][0]);
			return ret;
		}
		for (int j = 0; j < m; j++) {
			double ans = 0.;
			for (int i = 0; i < n; i++)
				ans += fabs(p[i][j]);
			if (ans > ret) ret = ans;
		}
		return ret;
	}

	inline double norm2() {
		double ret = 0.;
		if (n == 1) {
			for (int i = 0; i < m; i++) ret += p[0][i] * p[0][i];
			return sqrt(ret);
		}
		if (m == 1) {
			for (int i = 0; i < n; i++) ret += p[i][0] * p[i][0];
			return sqrt(ret);
		}
		return SVD().first.p[0][0];
	}

	inline double fnorm() {
		double ret = 0.;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				ret += p[i][j] * p[i][j];
		return sqrt(ret);
	}

	inline double det() {
		if (n != m) throw Exception("Matrix size doesn't match");
		Matrix l(ladder());
		double ans = 1.;
		for (int i = 0; i < n; i++) ans *= l.p[i][i];
		return ans;
	}

	inline int rank() {
		Matrix l(ladder());
		int r = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (fabs(l.p[i][j]) > eps) {r = i; break;}
			}
		}
		return r + 1;
	}

	~Matrix() {
		for (int i = 0; i < n; i++)	delete [] p[i];
		delete [] p;
	}
};

Matrix operator +(const Matrix& a, const Matrix& b) {
	if (a.n != b.n || a.m != b.m) throw Exception("Matrix size doesn't match!");
	Matrix c(a.n, a.m);
	for (int i = 0; i < a.n; i++)
		for (int j = 0; j < a.m; j++)
			c.p[i][j] = a.p[i][j] + b.p[i][j];
	return c;
}

Matrix operator -(const Matrix& a, const Matrix& b) {
	if (a.n != b.n || a.m != b.m) throw Exception("Matrix size doesn't match!");
	Matrix c(a.n, a.m);
	for (int i = 0; i < a.n; i++)
		for (int j = 0; j < a.m; j++)
			c.p[i][j] = a.p[i][j] - b.p[i][j];
	return c;
}

Matrix operator *(const Matrix& a, const Matrix& b) {
	if (a.m != b.n) throw Exception("Matrix size doesn't match!");
	Matrix c(a.n, b.m);
	for (int k = 0; k < a.m; k++)
		for (int i = 0; i < a.n; i++)
			for (int j = 0; j < b.m; j++)
				c.p[i][j] += a.p[i][k] * b.p[k][j];
	return c;
};

Matrix operator *(double a, const Matrix& b) {
	Matrix c(b.n, b.m);
	for (int i = 0; i < b.n; i++)
		for (int j = 0; j < b.m; j++)
			c.p[i][j] = a * b.p[i][j];
	return c;
}

Matrix operator *(const Matrix& b, double a) {
	return a * b;
}

Matrix operator ^(Matrix a, int m) {
	if (!a.isSquare()) throw Exception("Matrix size doesn't match!");
	int n = a.sz().first;
	if (m < 0) return a.inv() ^ (-m);
	Matrix ans(n, 'I');
	while (m > 0) {
		if (m & 1) ans = ans * a;
		a = a * a;
		m >>= 1;
	}
	return ans;
}

/* If A is nonsingular, solve A * x = b;
 * else solve A * x ~= b, in which x has the minimum 2-norm;
 */
Matrix operator /(Matrix a, Matrix b) {
	if (a.n != b.n) throw Exception("Matrix size doesn't match!");
	if (b.m != 1) throw Exception("Doesn't make sense!");
	if (a.n < a.m) throw Exception("Condition is not enough!");
	Matrix c(a.n, 1);
	if (a.n == a.m) {
		for (int i = 0; i < a.n; i++) {
			int cmp = i;
			for (int j = i + 1; j < a.n; j++)
				if (a.p[j][i] > a.p[cmp][i]) cmp = j;
			std::swap(b.p[i][0], b.p[cmp][0]);
			for (int k = i; k < a.m; k++) {
				std::swap(a.p[i][k], a.p[cmp][k]);
			}
			if (fabs(a.p[i][i]) < eps) throw("Matrix is not full rank");
			for (int j = 0; j < a.n; j++) {
				if (j == i) continue;
				double coef = -a.p[j][i] / a.p[i][i];
				for (int k = i; k < a.m; k++) {
					a.p[j][k] += coef * a.p[i][k];
				}
				b.p[j][0] += coef * b.p[i][0];
			}
		}
		for (int i = 0; i < a.n; i++) {
			c.p[i][0] = b.p[i][0] / a.p[i][i];
		}
	} else {
		for (int i = 0; i < a.m; i++) {




		}
	}
	return c;
}

Matrix operator /(const Matrix&a, double b) {
	if (fabs(b) == 0) throw Exception("Matrix can't be divided by zero!");
	return a * (1. / b);
}

Matrix& operator +=(Matrix& a, const Matrix& b) {
	if (a.n != b.n || a.m != b.m) throw Exception("Matrix size doesn't match!");
	for (int i = 0; i < a.n; i++)
		for (int j = 0; j < a.m; j++)
			a.p[i][j] += b.p[i][j];
	return a;
}

Matrix& operator -=(Matrix& a, const Matrix& b) {
	if (a.n != b.n || a.m != b.m) throw Exception("Matrix size doesn't match!");
	for (int i = 0; i < a.n; i++)
		for (int j = 0; j < a.m; j++)
			a.p[i][j] -= b.p[i][j];
	return a;
}

Matrix& operator *=(Matrix& a, double b) {
	for (int i = 0; i < a.n; i++)
		for (int j = 0; j < a.m; j++)
			a.p[i][j] *= b;
	return a;
}

Matrix& operator /=(Matrix& a, double b) {
	if (fabs(b) == 0) throw Exception("Matrix can't be divided by zero!");
	for (int i = 0; i < a.n; i++)
		for (int j = 0; j < a.m; j++)
			a.p[i][j] /= b;
	return a;
}

Matrix& operator *=(Matrix& a, const Matrix& b) {
	a = a * b;
	return a;
}

Matrix& operator /=(Matrix& a, const Matrix& b) {
	a = a / b;
	return a;
}

std::ostream& operator << (std::ostream& out, const Matrix& b) {
	for (int i = 0; i < b.n; i++) {
		for (int j = 0; j < b.m; j++)
			out << std::setw(10) << std::fixed << std::setprecision(15) << b.p[i][j] << " ";
		out << std::endl;
	}
	return out;
}

#endif
