//
// Created by Moin on 12/2/2022.
//
#include <iostream>

#include "gtest/gtest.h"
#include "../odesolvers.h"


double func(double x, double y);
double x_prev=0, y_init = 2, step_size = 1;
Eigen::ArrayXd x_array = Eigen::ArrayXd::LinSpaced(5, 0, 4);
Eigen::ArrayXd func_sys_ode(double x, Eigen::ArrayXd y);

TEST(odeSolverTest, Euler) {
    EXPECT_EQ(Euler(x_prev, y_init, step_size, func), 5) << "Euler Not Working.";

    Eigen::ArrayXd real_result(5);
    real_result << 2.0000,5.0000,11.40216,25.51321,56.84931;
    Eigen::ArrayXd res_array = Euler_vec(x_array, y_init, step_size, func);
    EXPECT_TRUE(res_array.isApprox(real_result, 0.05)) << "Euler_vec Not Working.";

    Eigen::ArrayXd x_array = Eigen::ArrayXd::LinSpaced(6,0,10);
    Eigen::ArrayXd y_init = Eigen::ArrayXd::Zero(2);
    y_init << 0,0;
    Eigen::MatrixXd realResultsMatrix(x_array.size(), y_init.size());
    realResultsMatrix << 0,0,
            0 ,19.62,
            39.24,36.4137,
            112.0674,46.2983,
            204.6640,50.1802,
            305.0244, 51.3123;
    auto y_matrix = sys_Euler_vec(x_array, y_init, func_sys_ode);
    EXPECT_TRUE(y_matrix.isApprox(realResultsMatrix, 0.001)) << "sys_rk4_vec Not Working.";
}

TEST(odeSolverTest, rk4) {
    EXPECT_NEAR(rk4(x_prev, y_init, step_size, func), 6.20104, 1e-4) << "rk4 Not Working.";

    Eigen::ArrayXd realResults(5);
    realResults << 2.0000, 6.20103707, 14.86248359, 33.72134801, 75.43917199;
    Eigen::ArrayXd resArray = rk4_vec(x_array, y_init, step_size, func);
    EXPECT_TRUE(resArray.isApprox(resArray, 1e-6)) << "rk4_vec Not Working";

    Eigen::ArrayXd x_array = Eigen::ArrayXd::LinSpaced(6,0,10);
    Eigen::ArrayXd y_init = Eigen::ArrayXd::Zero(2);
    y_init << 0,0;
    Eigen::MatrixXd realResultsMatrix(x_array.size(), y_init.size());
    realResultsMatrix << 0,0,
    19.1656 ,18.7256,
    71.9311,33.0995,
    147.9521, 42.0547,
    237.5104,46.9345,
    334.1626, 49.4027;
    auto y_matrix = sys_rk4_vec(x_array, y_init, func_sys_ode);
    EXPECT_TRUE(y_matrix.isApprox(realResultsMatrix, 0.001)) << "sys_rk4_vec Not Working.";
}

double func(double x, double y){
    return 4 * std::exp(0.8*x) - 0.5*y;
}

Eigen::ArrayXd func_sys_ode(double x, Eigen::ArrayXd y) {
    double y1 = y[0];
    double y2 = y[1];

    Eigen::ArrayXd funcRes = Eigen::ArrayXd::Zero(y.size());
    double f1 = y2;
    double f2 = 9.81 - (0.25/68.1) * std::pow(y2,2);
    funcRes << f1,f2;
    return funcRes;
}