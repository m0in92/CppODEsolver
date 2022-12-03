//
// Created by Moin on 12/2/2022.
//
#include <iostream>
#include <cmath>
#include "root_finder.h"


// below code is from the internet:
double Brent(double (*f)(double), double lower_bound, double upper_bound, double TOL, double MAX_ITER)
{
    double a = lower_bound;
    double b = upper_bound;
    double fa = f(a);   // calculated now to save function calls
    double fb = f(b);   // calculated now to save function calls
    double fs = 0;      // initialize
    double s = b;

    if (!(fa * fb < 0))
    {
        std::cout << "Signs of f(lower_bound) and f(upper_bound) must be opposites" << std::endl; // throws exception if root isn't bracketed
        return s;
    }

    if (std::abs(fa) < std::abs(b)) // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
    {
        std::swap(a,b);
        std::swap(fa,fb);
    }

    double c = a;           // c now equals the largest magnitude of the lower and upper bounds
    double fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
    bool mflag = true;      // boolean flag used to evaluate if statement later on
//    double s = 0;           // Our Root that will be returned
    double d = 0;           // Only used if mflag is unset (mflag == false)

    for (unsigned int iter = 1; iter < MAX_ITER; ++iter)
    {
        // stop if converged on root or error is less than tolerance
        if (std::abs(b-a) < TOL)
        {
            std::cout << "After " << iter << " iterations the root is: " << s << std::endl;
            return s;
        } // end if

        if (fa != fc && fb != fc)
        {
            // use inverse quadratic interopolation
            s =   ( a * fb * fc / ((fa - fb) * (fa - fc)) )
                  + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
                  + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
        }
        else
        {
            // secant method
            s = b - fb * (b - a) / (fb - fa);
        }

        /*
            Crazy condition statement!:
            -------------------------------------------------------
            (condition 1) s is not between  (3a+b)/4  and b or
            (condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
            (condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
            (condition 4) (mflag is set and |b−c| < |TOL|) or
            (condition 5) (mflag is false and |c−d| < |TOL|)
        */
        if (    ( (s < (3 * a + b) * 0.25) || (s > b) ) ||
                ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
                ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
                ( mflag && (std::abs(b-c) < TOL) ) ||
                ( !mflag && (std::abs(c-d) < TOL))  )
        {
            // bisection method
            s = (a+b)*0.5;

            mflag = true;
        }
        else
        {
            mflag = false;
        }

        fs = f(s);  // calculate fs
        d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
        c = b;      // set c equal to upper bound
        fc = fb;    // set f(c) = f(b)

        if ( fa * fs < 0)   // fa and fs have opposite signs
        {
            b = s;
            fb = fs;    // set f(b) = f(s)
        }
        else
        {
            a = s;
            fa = fs;    // set f(a) = f(s)
        }

        if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
        {
            std::swap(a,b);     // swap a and b
            std::swap(fa,fb);   // make sure f(a) and f(b) are correct after swap
        }

    } // end for

    return s;

} // end brent_fun
