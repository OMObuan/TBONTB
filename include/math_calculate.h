// Copyright [2024] <OMObuan>

#pragma once

#include <rename_type.h>

#include <type_traits>

namespace MATH_CALCULATE {

bool doubleEqual(double, double, double);

template <typename _T1, typename _T2, typename _Result>
    requires((std::is_same_v<_T2, i32> || std::is_same_v<_T2, usize>) &&
             requires(_Result result) { result* result; })
constexpr auto quickNatureNumberPow(_T1 const& baseNumber,
                                    _T2 const& indexNumber) {
    _Result result{1}, nowBaseNumber{baseNumber};
    while (indexNumber) {
        if (indexNumber % 2 == 0) {
            nowBaseNumber = nowBaseNumber * nowBaseNumber;
            indexNumber /= 2;
        } else {
            result *= nowBaseNumber;
            --indexNumber;
        }
    }
    return result;
}

}  // namespace MATH_CALCULATE
