// Copyright [2024] <OMObuan>

#include <math_calculate.h>

#include <cmath>

namespace MATH_CALCULATE {

bool f64Equal(f64 lft, f64 rgt, f64 eps) {
    if (std::fabs(lft - rgt) > eps) {
        return false;
    } else {
        return true;
    }
}

bool f64Less(f64 lft, f64 rgt, f64 eps) {
    if (std::fabs(lft - rgt) > eps) {
        return lft < rgt;
    } else {
        return false;
    }
}

bool f64Greater(f64 lft, f64 rgt, f64 eps) {
    if (std::fabs(lft - rgt) > eps) {
        return lft > rgt;
    } else {
        return false;
    }
}

bool f64LessEqual(f64 lft, f64 rgt, f64 eps) {
    if (std::fabs(lft - rgt) > eps) {
        return lft <= rgt;
    } else {
        return true;
    }
}

bool f64GreaterEqual(f64 lft, f64 rgt, f64 eps) {
    if (std::fabs(lft - rgt) > eps) {
        return lft >= rgt;
    } else {
        return true;
    }
}

}  // namespace MATH_CALCULATE
