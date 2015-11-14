#ifndef SCFUN
#define SCFUN
#include "matrix.hpp"
#include <cmath>

class Function{
private:
  double (* f) (const Matrix &x);
  int n;
  double eps;
public:
  Function(double (* fun)(const Matrix &x), int N): f(fun), n(N) {
    eps = 1e-8;
  }

  inline int numofArguments() const{
    return n;
  }

  inline double eval(const Matrix &x) const{
    return f(x);
  }

  Function(const Function &b): eps(b.eps), f(b.f), n(b.n){}

  Function& operator = (const Function &b) {
    if (this == (&b)) return *this;
    n = b.n;
    f = b.f;
    eps = b.eps;
    return *this;
  }

  double partial(Matrix x, int k) const{
    double f1, f2;
    if (k < 1 || k > n) throw Exception("Augument out of range!");
    x(1, k) += eps;
    f1 = f(x);
    x(1, k) -= 2. * eps;
    f2 = f(x);
    return (f1 - f2) / (2. * eps);
  }

  double secondpartial(Matrix x, int k1, int k2) const{
    double eps1 = sqrt(eps);
    double f11, f12, f21, f22;
    x(1, k1) -= eps1;
    x(1, k2) -= eps1;
    f11 = f(x);
    x(1, k2) += 2. * eps1;
    f12 = f(x);
    x(1, k1) += 2. * eps1;
    f22 = f(x);
    x(1, k2) -= 2. * eps1;
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
    Matrix c(1, n);
    for (int i = 1; i <= n; i++) {
      c(1, i) = partial(x, i);
    }
    return c;
  }

  void modifyprecision(double Neps) {
    eps = Neps;
  }

  ~Function() {}
};

#endif
