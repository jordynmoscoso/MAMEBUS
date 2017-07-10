// Header file for rk.c
#ifndef _RKTVD_H_
#define _RKTVD_H_

#include "nsode.h"

real rktvd1 ( real *                  t,
              real *                  x,
              real *                  xout,
              const real              cfl,
              const uint              numvars,
              DERIVATIVE_FUNCTION_CFL f);

real rktvd2 ( real *                  t,
              real *                  x,
              real *                  xout,
              real *                  k,
              const real              cfl,
              const uint              numvars,
              DERIVATIVE_FUNCTION_CFL f);

real rktvd3 ( real *                  t,
              real *                  x,
              real *                  xout,
              real *              	  k,
              const real              cfl,
              const uint              numvars,
              DERIVATIVE_FUNCTION_CFL f);

#endif
