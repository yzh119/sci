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
    Matrix c(n, m);
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= m; j++) {
        c.get(i, j) = uni();
      }
    return c;
  }

};
#endif
