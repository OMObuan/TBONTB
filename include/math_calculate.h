// Copyright [2024] <OMObuan>

#pragma once

#include <rename_type.h>

#include <type_traits>

namespace MATH_CALCULATE {

bool f64Equal(f64, f64, f64);

bool f64Less(f64, f64, f64);

bool f64Greater(f64, f64, f64);

bool f64LessEqual(f64, f64, f64);

bool f64GreaterEqual(f64, f64, f64);

template <typename _T1, typename _T2,
          typename _Result = std::decay_t<decltype(true ? std::declval<_T1>()
                                                        : std::declval<_T2>())>>
    requires((std::is_same_v<_T2, i32> || std::is_same_v<_T2, usize> ||
              std::is_same_v<_T2, i16>) &&
             requires(_Result result) { result* result; })
constexpr _Result quickNatureNumberPow(_T1 const& baseNumber, _T2 indexNumber) {
    _Result result{static_cast<_Result>(1)},
        nowBaseNumber{static_cast<_Result>(baseNumber)};
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
