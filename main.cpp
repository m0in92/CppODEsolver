#include <iostream>
#include <cmath>

#include "odesolvers.h"
#include "lib/Eigen/Dense"
#include "root_finder.h"

double func_root(double);
double func(double x, double y);
double x_prev=0, y_init = 2, step_size = 1;
Eigen::ArrayXd func_sys_ode(double x, Eigen::ArrayXd y);

// Example usage below:
int main() {
    Eigen::ArrayXd x_array = Eigen::ArrayXd::LinSpaced(6,0,10);
    Eigen::ArrayXd y_init = Eigen::ArrayXd::Zero(2);
    y_init << 0,0;
//    std::cout << func_sys_ode(x_array[0],y_init) << std::endl;
    auto y_matrix = sys_rk4_vec(x_array, y_init, func_sys_ode);
    std::cout << y_matrix << std::endl;

//    //Usage of the odesolver
//    double ynew = rk4(x_prev, y_init, step_size, func); //single data-point
//    std::cout << ynew << std::endl;
//
//    Eigen::ArrayXd x_array = Eigen::ArrayXd::LinSpaced(5, 0, 4); // an array of data-points
//    Eigen::ArrayXd res_array = rk4_vec(x_array, y_init, step_size, func);
//    std::cout << res_array << std::endl;
//    saveArray("results.csv", res_array);
//
//    // Usage of the root finder library
//    double root = Brent(func_root, 3.0, 0.0, 1e-5, 90);
//    std::cout << root << std::endl;

    return 0;
}

double func_root(double x){
    return (x-1) * (1 + std::pow((x-1),2));
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