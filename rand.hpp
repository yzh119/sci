#ifndef RANDOMMATRIX
#define RANDOMMATRIX
#include <cstdlib>
#include <ctime>
#include <climits>
#include "matrix.hpp"

class Randmat {
private:
// Return a real number between -1 and 1.
  double uni() {
    return 2. * rand() / RAND_MAX - 1.;
  }

public:
  Randmat() {
    srand((unsigned)time(NULL));
  }

  ~Randmat() {}
  Matrix unirand(int n, int m) {
    if (n < 1 || m < 1) throw Exception("Matrix size doesn't make sense!");
    Matrix c(n, m);
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= m; j++) {
        c(i, j) = uni();
      }
    return c;
  }

  Matrix unirand(int n) {
    if (n < 1) throw Exception("Matrix size doesn't make sense!");
    Matrix c(n);
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= n; j++) {
        c(i, j) = uni();
      }
    return c;
  }

  Matrix uniposdef(int n) {
    Matrix x = unirand(n);
    return x.transpose() * x;
  }
};
#endif
