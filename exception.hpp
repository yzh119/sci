#ifndef NUMEXCEPTION
#define NUMEXCEPTION
#include <string>

using std::string;

class Exception {
private:
  string sym;
public:
  Exception(const string &s): sym(s) {}
  string toString() {
    return sym;
  }
};

#endif
