#include <cmath>
#include "/home/tom/code/fmath/fmath.hpp"
#include "comp.h"

extern "C"{
  float fexp(float x){
    return fmath::exp(x);
  }

  float cexp(float x){
    return std::exp(x);
  }
}
