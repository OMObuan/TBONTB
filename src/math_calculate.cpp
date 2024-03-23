// Copyright [2024] <OMObuan>

#include <math_calculate.h>

#include <cmath>

namespace MATH_CALCULATE {

bool doubleEqual(double lft, double rgt, double eps) {
    if (std::fabs(lft - rgt) <= eps) {
        return true;
    } else {
        return false;
    }
}

}  // namespace MATH_CALCULATE
