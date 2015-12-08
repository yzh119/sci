#ifndef POLYNOMIAL
#define POLYNOMIAL
#include "function.hpp"
#include "matrix.hpp"
#include <vector>
class Poly: public fun1{
private:
  std::vector <double> a;
public:
  Poly(const std::vector <double> &A) {
    a = A;
  }

  ~Poly() {}

  Poly(const Poly& b) {
    a = b.a;
  }

  Poly& operator = (const Poly& b) {
    if (&b == this) return *this;
    a = b.a;
    return *this;
  }

  inline double eval(double x) const{
    double ret = 0.;
    for (int i = (int)a.size() - 1; i >= 0; i--) {
      ret = ret * x + a[i];
    }
    return ret;
  }
};
#endif
