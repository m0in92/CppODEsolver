//
// Created by Moin on 12/2/2022.
//

#ifndef ODE_SOLVERS_ROOT_FINDER_H
#define ODE_SOLVERS_ROOT_FINDER_H

double Brent(double (*f)(double), double lower_bound, double upper_bound, double TOL, double MAX_ITER);

#endif //ODE_SOLVERS_ROOT_FINDER_H
