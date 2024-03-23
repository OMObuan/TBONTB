// Copyright [2024] <OMObuan>

#include <gtest/gtest.h>

#include <GeneticAlogorithm.hpp>
#include <cstring>

using i32 = int;

bool test{};

i32 main(i32 argc, char *argv[]) {
    if (argc == 2) {
        if (strcmp(argv[1], "test") == 0) {
            test = true;
        }
    }
    if (test) {
        testing::InitGoogleTest(&argc, argv);
    }
    if (test) {
        return RUN_ALL_TESTS();
    }
    return 0;
}
