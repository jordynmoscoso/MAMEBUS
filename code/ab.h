//
//  ab.h
//  The header file for ab.c -- the time integration using the
//  Adams-Bashforth methods with a variable time step
//
//  Created by Jordyn Moscoso on 9/24/18.
//

#ifndef ab_h
#define ab_h

#include "nsode.h"
#include <stdio.h>
#include <string.h>

void ab1 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f);

void ab2 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            real *                  dxdt_1,
            const real              h,
            const real              h1,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f);

void ab3 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            real *                  dxdt_1,
            real *                  dxdt_2,
            const real              h,
            const real              h1,
            const real              h2,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f);

 /*void ab4 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            real *                  dxdt_1,
            real *                  dxdt_2,
            real *                  dxdt_3,
            real *                  h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f);*/



#endif /* ab_h */
