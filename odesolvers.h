//
// Created by Moin on 12/2/2022.
//

#ifndef ODE_SOLVERS_ODESOLVERS_H
#define ODE_SOLVERS_ODESOLVERS_H

#include <string>

#include "lib/Eigen/Dense"

double Euler(const double x_prev, const double y_prev, const double step_size, double (*func)(double, double));
double Euler(const double y_prev, const double step_size, double func_value);
Eigen::ArrayXd Euler_vec(Eigen::ArrayXd x_array, double y_init, double step_size, double (*func)(double, double));
Eigen::MatrixXd sys_Euler_vec(Eigen::ArrayXd x_array, Eigen::ArrayXd y_init_array, Eigen::ArrayXd (*func)(double, Eigen::ArrayXd));

double rk4(const double x_prev, const double y_prev, const double step_size, double (*func)(double x, double y));
Eigen::ArrayXd rk4_vec(Eigen::ArrayXd x_array, double y_init, double step_size, double (*func)(double, double));
Eigen::MatrixXd sys_rk4_vec(Eigen::ArrayXd x_array, Eigen::ArrayXd y_init_array, Eigen::ArrayXd (*func)(double, Eigen::ArrayXd));

void saveArray(std::string filename, Eigen::ArrayXd array);

#endif //ODE_SOLVERS_ODESOLVERS_H
