#ifndef UNLINEAR
#define UNLINEAR
#include "function.hpp"
#include <vector>
#include "rand.hpp"

class Solver{
private:
  std::vector <Function> fs;
  int n;
public:
  Solver(const std::vector <Function>& func): fs(func) {
    if (func.empty()) throw Exception("There is no function here!");
    n = fs[0].numofArguments();
  }

  Solver(const Solver &b) {
    fs = b.fs;
    n = b.n;
  }

  Solver& operator =(const Solver &b) {
    if (this == (&b)) return *this;
    fs = b.fs;
    n = b.n;
    return *this;
  }

  Matrix jacobian(const Matrix& x) {
    Matrix c((int)fs.size(), n);
    for (int i = 1; i <= (int)fs.size(); i++) {
      Matrix g = fs[i - 1].gradient(x);
      for (int j = 1; j <= n; j++) c.get(i, j) = g.elem(j);
    }
    return c;
  }

  Matrix eval(const Matrix& x) {
    Matrix c((int)fs.size(), 1);
    for (int i = 1; i <= (int)fs.size(); i++) {
      c.get(i, 1) = fs[i - 1].eval(x);
    }
    return c;
  }

  Matrix newton() {
    if (n != (int)fs.size()) throw Exception("Number of variables and equations doesn't match!");
    Randmat rm1;
    Matrix x(rm1.unirand(1, n));
    int cnt = 0; Matrix y(rm1.unirand(n, 1));
    while (y.norm1() > 1e-10) {
      y = eval(x);
      x -= (jacobian(x).inv() * y).transpose();
      cnt++;
    }
    return x.transpose();
  }

  Matrix broydens() {
    if (n != (int)fs.size()) throw Exception("Number of variables and equations doesn't match!");
    Randmat rm1;
    Matrix x(rm1.unirand(1, n));
    Matrix B(rm1.unirand(n, n)), S(n, 1), y(rm1.unirand(n, 1));
    while (y.norm1() > 1e-10) {
      y = eval(x);
      S = B / (Matrix(n, 1) - y);
      x += S.transpose();
      y = eval(x) - y;
      B += ((y - B * S) * S.transpose()) / (S.transpose() * S).get(1, 1);
    }
    return x.transpose();
  }
};

class optimization{

};

#endif
