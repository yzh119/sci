#ifndef COMP
#define COMP
#include <cmath>
#include "exception.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <utility>

class comp{
private:
  double a, b;
public:
  comp(double t = 0.): a(t), b(0.) {}
  comp(double a, double b): a(a), b(b) {}
  ~comp() {}

  comp(const comp &y) {
    a = y.a;
    b = y.b;
  }

  comp& operator = (const comp &y) {
    if (this == (&y)) return *this;
    a = y.a;
    b = y.b;
    return *this;
  }
  inline double re() {return a;}
  inline double im() {return b;}
  inline bool isreal() {
    return (b < 1e-10);
  }

  inline double modulo() {
    return sqrt(a * a + b * b);
  }

  friend comp operator + (const comp&, const comp&);
  friend comp operator - (const comp&, const comp&);
  friend comp operator * (const comp&, const comp&);
  friend comp operator * (const comp&, double);
  friend comp operator * (double, const comp&);
  friend comp operator / (const comp&, const comp&);
  friend comp operator / (const comp&, double);
  friend std::ostream& operator << (std::ostream&, const comp&);
};

comp operator + (const comp& x, const comp& y) {
  return comp(x.a + y.a, x.b - y.b);
}

comp operator - (const comp& x, const comp& y) {
  return comp(x.a - y.a, x.b - y.b);
}

comp operator * (const comp& x, const comp& y) {
  return comp(x.a * y.a - x.b * y.b, x.a * y.b + x.b * y.a);
}

comp operator * (double x, const comp& y) {
  return comp(x * y.a, x * y.b);
}

comp operator * (const comp& x, double y) {
  return comp(x.a * y, x.b * y);
}

comp operator / (const comp& x, const comp& y) {
  double r2 = y.a * y.a + y.b * y.b;
  return comp((x.a * y.a + x.b * y.b) / r2, (x.b * y.a - x.a * y.b) / r2);
}

comp operator / (const comp& x, double y) {
  return comp(x.a / y, x.b / y);
}

std::ostream& operator << (std::ostream& out, const comp& x) {
	out << std::fixed << std::setprecision(15) << x.a << " " << ((x.b < 0) ? '-': '+') << " " << ((x.b < 0) ? -x.b: x.b) << 'i';
	return out;
}

std::pair<comp, comp> caleig22(double x11, double x12, double x21, double x22) {
  double b = x11 + x22, delta = (x11 + x22) * (x11 + x22) - 4 * (x11 * x22 - x12 * x21);
  if (delta > 0) return std::make_pair(comp((b - sqrt(delta)) / 2., 0.), comp((b + sqrt(delta)) / 2., 0.));
    else return std::make_pair(comp(b / 2., -sqrt(-delta) / 2.), comp(b / 2., sqrt(-delta) / 2.));
}
#endif
