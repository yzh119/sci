#ifndef INTERPOLATION
#define INTERPOLATION
#include "polynomial.hpp"
#include "exception.hpp"
#include <vector>
#include <algorithm>
#include <iostream>
using std::vector;
using std::swap;
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::lower_bound;

class Inter: public fun1{
private:
  vector <double> b, x, y, d;
public:
  Inter() {}
  ~Inter() {}

// Newton method
  inline void interpolate(double x_i, double y_i) {
    if ((int)b.size() == 0) {
      b.push_back(y_i); d.push_back(y_i);
    } else {
      double temp = d[0];
      d[0] = y_i;
      int n = d.size();
      for (int i = 1; i < n; i++) {
        temp = (d[i - 1] - temp) / (x_i - x[n - i]);
        swap(d[i], temp);
      }
      d.push_back((d[n - 1] - temp) / (x_i - x[0]));
      b.push_back(d.back());
    }
    x.push_back(x_i); y.push_back(y_i);
  }

  Inter(vector <double>& X, vector <double>& Y) {
    if (X.size() != Y.size()) throw("Interpolation points' size doesn't match!");
    for (int i = 0; i < (int)X.size(); i++)
      interpolate(X[i], Y[i]);
  }

  inline double eval(double X) const{
    double ret = 0.;
    int n = (int)x.size();
    for (int i = n - 1; i >= 0; i--) {
      ret = b[i] + ret * (X - x[i]);
    }
    return ret;
  }
};

class Cspline: public fun1{
private:
  vector <Poly> spline;
  vector <double> pts;
  int n;
public:
  Cspline(vector<double> x, vector<double> y) {
    if (x.size() != y.size()) throw Exception("Interpolation points' size doesn't match!");
    for (int i = 0; i < (int)x.size(); i++)
      for (int j = i + 1; j < (int)x.size() - 1; j++) {
        if (x[i] > x[j]) swap(x[i], x[j]), swap(y[i], y[j]);
      }

    pts = x;

    n = (int)x.size();
    Matrix M(4 * (n - 1)), V(4 * (n - 1), 1);
    int i = 1, j = 1;
    for (int t = 0; t < n - 1; t++, j += 4, i += 2) {
      M(i, j) = 1.; M(i, j + 1) = x[t]; M(i, j + 2) = x[t] * x[t]; M(i, j + 3) = x[t] * x[t] * x[t];
      M(i + 1, j) = 1.; M(i + 1, j + 1) = x[t + 1]; M(i + 1, j + 2) = x[t + 1] * x[t + 1]; M(i + 1, j + 3) = x[t + 1] * x[t + 1] * x[t + 1];
      V(i, 1) = y[t];
      V(i + 1, 1) = y[t + 1];
    }
    j = 1;
    for (int t = 1; t < n - 1; t++, i++) {
      M(i, j + 1) = 1.; M(i, j + 2) = 2 * x[t]; M(i, j + 3) = 3 * x[t] * x[t]; j += 4;
      M(i, j + 1) = -1.; M(i, j + 2) = -2 * x[t]; M(i, j + 3) = -3 * x[t] * x[t];
    }
    j = 1;
    for (int t = 1; t < n - 1; t++, i++) {
      M(i, j + 2) = 2.; M(i, j + 3) = 6 * x[t]; j += 4;
      M(i, j + 2) = -2.; M(i, j + 3) = -6 * x[t];
    }
    M(i, 3) = 2.; M(i, 4) = 6 * x[0]; i++;
    M(i, 4 * n - 5) = 2.; M(i, 4 * n - 4) = 6 * x[n - 1];

    M /= V;
    vector <double> coef(4);
    for (int i = 1; i <= 4 * (n - 1); i += 4) {
      for (int j = i; j < i + 4; j++) coef[j - i] = M(j, 1);
      spline.push_back(Poly(coef));
    }
  }
  ~Cspline() {}

  inline double eval(double x) const{
    if (x <= pts.front() || x >= pts.back()) return 0.;
    return spline[lower_bound(pts.begin(), pts.end(), x) - pts.begin() - 1].eval(x);
  }
};

class bfun: public fun1{
private:
  vector <double> x;
public:
  bfun(const vector <double>& X): x(X) {}
  ~bfun() {}

  double cal(int l, int r, double t) const{
    if (t < x[l]) return 0.;
    if (t >= x[r]) return 0.;
    if (l >= r) return 0.;
    if (r == l + 1) return 1.;
    return (t - x[l]) * cal(l, r - 1, t) / (x[r - 1] - x[l]) + (x[r] - t) * cal(l + 1, r, t) / (x[r] - x[l + 1]);
  }

  inline double eval(double t) const{
    return cal(0, (int)x.size() - 1, t);
  }
};

class Bspline: public fun1{
private:
  vector <double> pts;
  vector <bfun> spline;
  vector <double> a;
  int k, n, l, r;
  inline double pt(int x) {
    double base = 0., step = (pts.back() - pts.front()) * n / (n - 1.);
    while (x < 0) x += n, base -= step;
    while (x >= n) x -= n, base += step;
    return pts[x] + base;
  }
public:
  Bspline(vector <double> x, vector <double> y, int K = 3): k(K){
    if (x.size() != y.size()) throw Exception("Interpolation points' size doesn't match!");
    for (int i = 0; i < (int)x.size(); i++)
      for (int j = i + 1; j < (int)x.size() - 1; j++) {
        if (x[i] > x[j]) swap(x[i], x[j]), swap(y[i], y[j]);
      }

    pts = x;
    n = x.size();
    l = -(k + 1)/ 2, r = k + 1 + l;

    vector <double> t(k + 2);
    for (int i = 0; i < n; i++) {
      for (int j = i + l; j <= i + r; j++) {
        t[j - i - l] = pt(j);
      }
      spline.push_back(bfun(t));
    }

    Matrix M(n), V(n, 1);
    for (int i = 1; i <= n; i++) {
      V(i, 1) = y[i - 1];
      for (int j = max(1, i - r - 2); j <= min(n, i - l); j++)
        M(i, j) = spline[j - 1].eval(x[i - 1]);
    }

    M /= V;

    for (int i = 1; i <= n; i++)
      a.push_back(M(i, 1));
  }

  ~Bspline() {}

  inline double eval(double x) const{
    double ret = 0.;
    int i = lower_bound(pts.begin(), pts.end(), x) - pts.begin();
    for (int j = max(i - r - 1, 0); j < min(n, i - l); j++) {
      ret += spline[j].eval(x) * a[j];
    }
    return ret;
  }
};
#endif
