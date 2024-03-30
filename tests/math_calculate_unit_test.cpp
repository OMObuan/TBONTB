// Copyright [2024] <OMObuan>

#include <TBONTB/math_calculate.h>
#include <gtest/gtest.h>

namespace TEST_FOR_F64_EQUAL {

TEST(MathCalculateUnitTest, F64EqualUnitTest) {
    ASSERT_TRUE(MATH_CALCULATE::f64Equal(1.0, 0.99, 0.1));
    ASSERT_FALSE(MATH_CALCULATE::f64Equal(1. / 3., 0.33, 0.00001));
}

}  // namespace TEST_FOR_F64_EQUAL

namespace TEST_FOR_QUICK_NATURE_NUMBER_POW {

TEST(MathCalculateUnitTest, QuickNatureNumberPowUnitTest) {
    ASSERT_EQ(MATH_CALCULATE::quickNatureNumberPow(10, 3), 1000);
    ASSERT_TRUE(MATH_CALCULATE::f64Equal(
        MATH_CALCULATE::quickNatureNumberPow(2.2, 10), 2655.99227914, 0.0001));
}

}  // namespace TEST_FOR_QUICK_NATURE_NUMBER_POW
