#ifndef PLOT2D
#define PLOT2D
#include <string>
#include <cstring>
#include <cstdlib>
#include <utility>
#include <cstdio>
#include <cmath>
#include <vector>
#include <ctime>
#include "exception.hpp"
const int DIM = 2048;
using std::vector;
using std::pair;
using std::min;
using std::max;
using std::make_pair;

class Pl2d{
private:
  double eps;
  unsigned char color[DIM][DIM][3];
  enum Type{light, bold};

  vector <double> ptx;
  vector <double> pty;
  vector <int> clr;
  vector <int> r, g, b;
  vector <Type> type;
  double minx, maxx, miny, maxy;
  double a1, b1, a2, b2;
  int ccnt;

  inline double sgn(double x) {
      return (x < -eps) ? -1.: (x > eps? 1.: 0);
  }

  inline pair<double, double> map(double x, double y) {
    return make_pair(int(a1 * y + b1), int(a2 * x + b2));
  }

  void newcolor() {
    ccnt++;
    r.push_back(rand() & 255);
    g.push_back(rand() & 255);
    b.push_back(rand() & 255);
  }

public:
  Pl2d() {
    ccnt = -1;
    eps = 1e-10;
    minx = 1e18, maxx = -1e18;
    miny = 1e18; maxy = -1e18;
    srand((unsigned)time(NULL));
    memset(color, 255, sizeof(color));
  }

  ~Pl2d() {}

  void adjustrange() {
    if (fabs(maxy - miny) > eps) {
      a1 = 2000. / (miny - maxy);
      b1 = (2022. * maxy - 22. * miny) / (maxy - miny);
    } else {a1 = 0; b1 = DIM / 2.;}
    a2 = 2000. / (maxx - minx);
    b2 = (22. * maxx - 2022. * minx) / (maxx - minx);
  }

  void plot(const vector <double>& x, const vector <double>&y, int weight = 1) {
    int n = (int)x.size();
    newcolor();
    for (int i = 0; i < n; i++) {
      minx = min(x[i], minx);
      maxx = max(x[i], maxx);
      miny = min(y[i], miny);
      maxy = max(y[i], maxy);
      ptx.push_back(x[i]);
      pty.push_back(y[i]);
      clr.push_back(ccnt);
      type.push_back((weight == 1) ? bold: light);
    }
    adjustrange();
  }

  void plot(double x, double y, int weight = 1) {
    minx = min(x, minx);
    maxx = max(x, maxx);
    miny = min(y, miny);
    maxy = max(y, maxy);
    newcolor();
    ptx.push_back(x);
    pty.push_back(y);
    clr.push_back(ccnt);
    type.push_back((weight == 1) ? bold: light);
    adjustrange();
  }

  void plot(const fun1 &f, double l, double r) {
    if (l > r) throw Exception("WTF!, l > r!");
    minx = min(l, minx);
    maxx = max(r, maxx);
    double dx = (r - l) / 2000., dy = dx / 10.;
    double f1 = f.eval(l), f2;
    miny = min(f1, miny);
    maxy = max(f1, maxy);
    newcolor();
    for (double x = l + dx; x <= r; x += dx, f1 = f2) {
      f2 = f.eval(x);
      miny = min(f2, miny);
      maxy = max(f2, maxy);
      dy = max(dy, (maxy - miny) / 2000.);
      for (double y = min(f1, f2); y <= max(f1, f2); y += dy) {
        if (fabs(f2 - f1) > eps) {
          ptx.push_back(x + (y - f2) * dx / (f2 - f1));
        } else ptx.push_back(x);
        pty.push_back(y);
        clr.push_back(ccnt);
        type.push_back(light);
      }
    }
    adjustrange();
  }

  void generate(const string& filename) {
    FILE *fp = fopen((filename + ".ppm").c_str(), "wb");
    fprintf(fp, "P6\n%d %d\n255\n", DIM, DIM);

    int n = (int)ptx.size();
    for (int i = 0; i < n; i++) {
      pair <int, int> cord = map(ptx[i], pty[i]);
      if (type[i] == light) {
        for (int xx = cord.first - 1; xx <= cord.first + 1; xx++)
          for (int yy = cord.second - 1; yy <= cord.second + 1; yy++) {
            color[xx][yy][0] = r[clr[i]];
            color[xx][yy][1] = g[clr[i]];
            color[xx][yy][2] = b[clr[i]];
          }
      } else {
        for (int xx = cord.first - 10; xx <= cord.first + 10; xx++)
          for (int yy = cord.second - 10; yy <= cord.second + 10; yy++) {
            if ((xx - cord.first) * (xx - cord.first) + (yy - cord.second) * (yy - cord.second) <= 100) {
              color[xx][yy][0] = r[clr[i]];
              color[xx][yy][1] = g[clr[i]];
              color[xx][yy][2] = b[clr[i]];
            }
          }
      }
    }

    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
        fwrite(color[i][j], 1, 3, fp);
    fclose(fp);
  }
};

#endif
