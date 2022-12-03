//
// Created by Moin on 12/2/2022.
//
#include <cmath>
#include "root_finder.h"

int sign(double x){
    if (x>=0)
        return 1;
    else
        return -1;
}

double Brent(double a, double b, double (*func)(double), double tol){
    double c = a;
    double m = (c-b)/2.0;
    double root = b;
    double b__, b_, b_new;

    while ((func(a) * func(b)) > 0) {
        if ((std::abs(b)<=tol) | (std::abs(m)<=tol)){
            root = b;
            break;
        }
        double i = b - func(b) * (b-a) / (func(b) - func(a));
        if ((i >= b) | (i <= (b+m)))
            b__ = i;
        else
            b__ = b+m;
        if (std::abs(b - b__) > tol)
            b_ = b__;
        else
            b_ = b__ + tol*sign(m);
        b_new = b_;
        a = b;
        if (sign(func(b_new)) != sign(func(b))) {
            c = b;
        }
        b = b_new;
    }
    return root;
}
