//
// Created by Moin on 12/2/2022.
//
#include <iostream>
#include <fstream>

//#include "lib/Eigen/Dense"
#include "odesolvers.h"

/*
 * These functions solve the initial value ODEs of the form:
 * dy/dx = func(x,y).
 *
 * When using the below functions, formulate the input func such that its first argument and second arguments are x and
 * y, respectively. As an example, consider double func(double x, double y).
 * */

double Euler(const double x_prev, const double y_prev, const double step_size, double (*func)(double x, double y)){
    return y_prev + func(x_prev, y_prev) * step_size;
}

Eigen::ArrayXd Euler_vec(Eigen::ArrayXd x_array, double y_init, double step_size, double (*func)(double, double)){
    Eigen::ArrayXd resArray = Eigen::ArrayXd::Zero(x_array.size()); // initialize the results vector to zeros.
    for (int index = 1; index < x_array.size(); index++){
        resArray[index] = Euler(x_array[index-1], resArray[index-1], step_size, func);
    }
    return resArray;
}

double rk4(const double x_prev, const double y_prev, const double step_size, double (*func)(double x, double y)){
    double k1 = func(x_prev, y_prev);
    double k2 = func(x_prev + 0.5*step_size, y_prev + 0.5*k1*step_size);
    double k3 = func(x_prev + 0.5*step_size, y_prev + 0.5*k2*step_size);
    double k4 = func(x_prev + step_size, y_prev + k3*step_size);
    return y_prev + (k1 + 2*k2 + 2*k3 + k4) * (step_size/6.0);
}

Eigen::ArrayXd rk4_vec(Eigen::ArrayXd x_array, double y_init, double step_size, double (*func)(double, double)){
    Eigen::ArrayXd res = Eigen::ArrayXd::Zero(x_array.size());
    res[0] = y_init;
    for (int i=1; i < res.size(); i++){
        res[i] = rk4(x_array[i-1], res[i-1], step_size, func);
    }
    return res;
}

void saveArray(std::string fileName, Eigen::ArrayXd array) {
//    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file;
    file.open(fileName, std::ios::out | std::ios::app);
    file << array;
    file.close();
}
