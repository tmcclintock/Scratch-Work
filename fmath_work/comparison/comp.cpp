#include <cmath>
#include "/home/tom/code/fmath/fmath.hpp"
#include "comp.h"

extern "C" float fexpf(float x){
  return fmath::exp(x);
}

extern "C" float cexpf(float x){
  return std::exp(x);
}

extern "C" double fexpd(double x){
  return fmath::expd(x);
}

extern "C" double cexpd(double x){
  return std::exp(x);
}

