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

Eigen::ArrayXd matrix_row_to_array(const Eigen::MatrixXd& matrix, const int rowIndex) {
    Eigen::ArrayXd rowArray = Eigen::ArrayXd::Zero(matrix.cols());
    for (int colIndex=0; colIndex <  rowArray.size(); colIndex++){
        rowArray[colIndex] = matrix(rowIndex, colIndex);
    }
    return rowArray;
}

double Euler(const double x_prev, const double y_prev, const double step_size, double (*func)(double x, double y)){
    return y_prev + func(x_prev, y_prev) * step_size;
}

double Euler(const double y_prev, const double step_size, double func_value){
    return y_prev + func_value * step_size;
}

Eigen::ArrayXd Euler_vec(Eigen::ArrayXd x_array, double y_init, double step_size, double (*func)(double, double)){
    Eigen::ArrayXd resArray = Eigen::ArrayXd::Zero(x_array.size()); // initialize the results vector to zeros.
    resArray[0] = y_init;
    for (int index = 1; index < x_array.size(); index++){
        resArray[index] = Euler(x_array[index-1], resArray[index-1], step_size, func);
    }
    return resArray;
}

Eigen::MatrixXd sys_Euler_vec(Eigen::ArrayXd x_array, Eigen::ArrayXd y_init_array, Eigen::ArrayXd (*func)(double, Eigen::ArrayXd)){
    Eigen::MatrixXd resMatrix = Eigen::MatrixXd::Zero(x_array.size(), y_init_array.size());
    Eigen::ArrayXd y_prev = Eigen::ArrayXd::Zero(y_init_array.size());
    Eigen::ArrayXd res_prev = Eigen::ArrayXd::Zero(y_init_array.size());
    double dx = 0;
    for (int i = 0; i < resMatrix.rows(); i++){
        if (i==0) {
            for (int j=0; j<resMatrix.cols(); j++) {
                resMatrix(i,j) = y_init_array[j];
            }
        }
        else {
            dx = x_array[i] - x_array[i-1];
            y_prev = matrix_row_to_array(resMatrix, i-1);
            res_prev = func(x_array[i-1], y_prev);
            for (int j =0; j < resMatrix.cols(); j++){
                resMatrix(i,j) = Euler(y_prev[j], dx, res_prev[j]);
//                resMatrix(i,j) = y_prev[j] + res_prev[j] * dx;
            }
        }
    }
    return resMatrix;
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

Eigen::MatrixXd sys_rk4_vec(Eigen::ArrayXd x_array, Eigen::ArrayXd y_init_array, Eigen::ArrayXd (*func)(double, Eigen::ArrayXd)){
    Eigen::MatrixXd resMatrix = Eigen::MatrixXd::Zero(x_array.size(), y_init_array.size());
    Eigen::ArrayXd y_prev = Eigen::ArrayXd::Zero(y_init_array.size());
    Eigen::ArrayXd res_prev = Eigen::ArrayXd::Zero(y_init_array.size());
    double dx = 0;
    for (int i = 0; i < resMatrix.rows(); i++){
        if (i==0) {
            for (int j=0; j<resMatrix.cols(); j++) {
                resMatrix(i,j) = y_init_array[j];
            }
        }
        else {
            dx = x_array[i] - x_array[i-1];
            y_prev = matrix_row_to_array(resMatrix, i-1);
            res_prev = func(x_array[i-1], y_prev);
            // calc k1,k2,k3,k4
            auto k1 = func(x_array[i-1], y_prev);
            auto k2 = func(x_array[i-1] + 0.5*dx, y_prev + 0.5*k1*dx);
            auto k3 = func(x_array[i-1] + 0.5*dx, y_prev + 0.5*k2*dx);
            auto k4 = func(x_array[i-1] + dx, y_prev + k3*dx);
            for (int j =0; j < resMatrix.cols(); j++){
                resMatrix(i,j) = y_prev[j] + (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]) * (dx/6.0);
            }
        }
    }
    return resMatrix;
}

void saveArray(std::string fileName, Eigen::ArrayXd array) {
//    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file;
    file.open(fileName, std::ios::out | std::ios::app);
    file << array;
    file.close();
}
