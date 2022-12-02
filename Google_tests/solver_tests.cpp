//
// Created by Moin on 12/2/2022.
//
#include <iostream>
#include "gtest/gtest.h"
#include "../odesolvers.h"

double func(double x, double y){
    return 4 * std::exp(0.8*x) - 0.5*y;
}

double x_prev=0, y_init = 2, step_size = 1;
Eigen::ArrayXd x_array = Eigen::ArrayXd::LinSpaced(5, 0, 4);

TEST(odeSolverTest, Euler) {
    EXPECT_EQ(Euler(x_prev, y_init, step_size, func), 5) << "Euler Not Working.";

    Eigen::ArrayXd real_result(5);
    real_result << 2.0000,5.0000,11.40216,25.51321,56.84931;
    Eigen::ArrayXd res_array = Euler_vec(x_array, y_init, step_size, func);
    EXPECT_TRUE(res_array.isApprox(real_result, 0.05)) << "Euler_vec Not Working.";
}

TEST(odeSolverTest, rk4) {
    EXPECT_NEAR(rk4(x_prev, y_init, step_size, func), 6.20104, 1e-4) << "rk4 Not Working.";

    Eigen::ArrayXd realResults(5);
    realResults << 2.0000, 6.20103707, 14.86248359, 33.72134801, 75.43917199;
    Eigen::ArrayXd resArray = rk4_vec(x_array, y_init, step_size, func);
    EXPECT_TRUE(resArray.isApprox(resArray, 1e-6)) << "rk4_vec Not Working";
}