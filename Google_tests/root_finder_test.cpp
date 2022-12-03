//
// Created by Moin on 12/2/2022.
//
#include <cmath>

#include "gtest/gtest.h"
#include "../root_finder.h"

double func_root(double x){
    return (x-1) * (1 + std::pow((x-1),2));
}

TEST(RootFinderTest, Brent) {
    EXPECT_EQ(Brent(func_root, 3.0, 0.0, 1e-5, 90), 1) << "Brent not working.";
}
