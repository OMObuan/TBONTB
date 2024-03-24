// Copyright [2024] <OMObuan>

#include <gtest/gtest.h>

#include <GeneticAlgorithm.hpp>
#include <cstring>

using i32 = int;

bool test{};

i32 main(i32 argc, char *argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
