#ifndef SCFUN
#define SCFUN
#include "matrix.hpp"
#include <cmath>

class Function {
public:
  double eps;
  Function() {eps = 1e-8;}
  Function(const Function &b) {
    eps = b.eps;
  }

  Function& operator = (const Function &b) {
    if (&b == this) return *this;
    eps = b.eps;
    return *this;
  }

  inline void modifyprecision(double Neps) {
    eps = Neps;
  }
  virtual inline int numofArguments() const { return 0; }
  virtual inline double eval(const Matrix&) { return 0.; }
};

class funx: public Function{
private:
  double (* f) (const Matrix &x);
  int n;
public:
  double eps;

  funx(double (* fun)(const Matrix &x), int N): f(fun), n(N) {}

  ~funx() {}

  inline int numofArguments() const{
    return n;
  }

  inline double eval(const Matrix &x) const{
    return f(x);
  }

  funx(const funx &b): f(b.f), n(b.n){}

  funx& operator = (const funx &b) {
    if (this == (&b)) return *this;
    n = b.n;
    f = b.f;
    return *this;
  }

  double partial(Matrix x, int k) const{
    double f1, f2;
    if (k < 1 || k > n) throw Exception("Augument out of range!");
    x[k] += eps;
    f1 = f(x);
    x[k] -= 2. * eps;
    f2 = f(x);
    return (f1 - f2) / (2. * eps);
  }

  double secondpartial(Matrix x, int k1, int k2) const{
    double eps1 = sqrt(eps);
    double f11, f12, f21, f22;
    x[k1] -= eps1;
    x[k2] -= eps1;
    f11 = f(x);
    x[k2] += 2. * eps1;
    f12 = f(x);
    x[k1] += 2. * eps1;
    f22 = f(x);
    x[k2] -= 2. * eps1;
    f21 = f(x);
    return ((f22 - f12) / (2. * eps1) - (f21 - f11) / (2. * eps1)) / (2. * eps1);
  }

  Matrix hessian(const Matrix &x) {
    Matrix c(n, n);
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        if (i <= j) c(i, j) = secondpartial(x, i, j);
          else c(i, j) = c(j, i);
      }
    }
    return c;
  }

  Matrix gradient(Matrix x) {
    Matrix c(n, 1);
    for (int i = 1; i <= n; i++) {
      c[i] = partial(x, i);
    }
    return c;
  }

};

class fun1:public Function{
private:
  double (* f) (double x);
public:
  fun1() {f = NULL;}
  fun1(double (* fun)(double )): f(fun) {}
  ~fun1() {}
  inline int numofArguments() {return 1;}

  virtual inline double eval(double x) const{
    if (f == NULL) return 0.;
    return f(x);
  }

  inline double eval(Matrix& M) const{
    return eval(M(1, 1));
  }

  inline double derivative(double x) const{
    return (eval(x + eps) - eval(x - eps)) / (2. * eps);
  }

  inline double twoorder(double x) const {
	  double eps1 = 1e-4;
	  return (eval(x + eps1) + eval(x - eps1) - 2. * eval(x)) / (eps1 * eps1);
  }
};
#endif
