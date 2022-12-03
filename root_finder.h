//
// Created by Moin on 12/2/2022.
//

#ifndef ODE_SOLVERS_ROOT_FINDER_H
#define ODE_SOLVERS_ROOT_FINDER_H

int sign(double x);
double Brent(double a, double b, double (*func)(double), double tol);

#endif //ODE_SOLVERS_ROOT_FINDER_H
