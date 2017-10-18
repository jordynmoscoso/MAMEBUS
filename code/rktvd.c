/****************************************************************************************
**
**  File        ::      rk.c
**
**  Author      ::      Andrew Stewart
**
**  Description ::      Provides a series of total variation diminishing (TVD)
**                      Runge-Kutta methods of various orders, for numerical integration.
**
****************************************************************************************/
#include "rktvd.h"

/****************************************************************************************
**
**  Function    ::      rktvd1
**
**  Purpose     ::      Performs a single iteration of the standard Runge-Kutta first
**                      order method, incrementing t by the step-size h and approximating
**                      the next set of dependent variable values (x), which are updated
**                      with their new values.
**
**                      The arrays x and xout MUST each have a length at least
**                      as great as numvars
**
**  Input       ::      PARAMETERS:
**
**                      t
**                          Pointer to the independent variable value.
**                      x
**                          Array of values of dependent variables
**                      xout
**                          Array into which the incremented x-values will be put
**                      cfl
**                          Specifies the CFL number to use. The derivative function f
**                          must return an estimate of the minimum time dt_w for a wave
**                          to cross a grid cell, or otherwise reach a point that
**                          invalidates the method. The time step will be chosen such
**                          that dt = cfl*dt_w.
**                      numvars
**                          Specifies the size of all input arrays
**                      f
**                          Derivative function
**
**  Output      ::      The value of the independent variable (t) will be incremented by h,
**                      and the new values of the dependent variables will be put in
**                      the xout array.
**
**                      Returns the time step used, based on the supplied CFL number.
**
****************************************************************************************/
real rktvd1 ( real *                  t,
              real *                  x,
              real *                  xout,
              const real              cfl,
              const uint              numvars,
              DERIVATIVE_FUNCTION_CFL f)
{
    // For looping
    uint i;

    // Max time step, based on estimate from derivative function
    real dt_w;

    // Actual time step
    real h;

    // Calculate the derivatives of x at t
    dt_w = (*f)(*t, x, xout, numvars);

    // Calculate time step
    h = cfl*dt_w;

    // Increment the dependent variable
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = xout[i] * h + x[i];
    }

    // Increment t
    *t += h;

    return h;
}

/****************************************************************************************
 **
 **  Function    ::      rkt1
 **
 **  Purpose     ::      Performs a single iteration of the standard Runge-Kutta first
 **                      order method (Euler integration), with a prescribed time step
 **
 **                      The arrays x and xout MUST each have a length at least
 **                      as great as numvars
 **
 **  Input       ::      PARAMETERS:
 **
 **                      t
 **                          Pointer to the independent variable value.
 **                      x
 **                          Array of values of dependent variables
 **                      xout
 **                          Array into which the incremented x-values will be put
 **                      cfl
 **                          Specifies the CFL number to use. The derivative function f
 **                          must return an estimate of the minimum time (dt_w = 1) for a wave
 **                          to cross a grid cell, or otherwise reach a point that
 **                          invalidates the method. The time step will be chosen such
 **                          that dt = cfl*dt_w.
 **                      numvars
 **                          Specifies the size of all input arrays
 **                      f
 **                          Derivative function
 **
 **  Output      ::      The value of the independent variable (t) will be incremented by h,
 **                      and the new values of the dependent variables will be put in
 **                      the xout array.
 **
 **                      Returns the time step used, based on the supplied CFL number.
 **
 ****************************************************************************************/

real rk1   ( real *                  t,
             real *                  x,
             real *                  xout,
             const real              cfl,
             const uint              numvars,
             DERIVATIVE_FUNCTION_CFL f)
{
    // For looping
    uint i;
    
    // Max time step, based on estimate from derivative function
    real dt_w;
    
    // Actual time step
    real h;
    
    // Calculate the derivatives of x at t
    dt_w = (*f)(*t, x, xout, numvars);
    
    // Calculate time step
    h = cfl*dt_w;
    
    // Increment the dependent variable
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = xout[i] * h + x[i];
    }
    
    // Increment t
    *t += h;
    
    return h;
}

/****************************************************************************************
**
**  Function    ::      rktvd2
**
**  Purpose     ::      Performs a single iteration of the TVD Runge-Kutta second
**                      order method, incrementing t by the step-size h and approximating
**                      the next set of dependent variable values (x), which are updated
**                      with their new values.
**
**                      The arrays x, xout, and k MUST each have a length at least
**                      as great as numvars, and MUST all be distinct arrays that are not
**                      modified as the method proceeds. If these conditions are not met,
**                      the results of this procedure are undefined.
**
**                      k is a temporary storage buffer used
**                      by the Runge-Kutta algorithm. It is passed in as a parameter
**                      so that memory allocation need not happen every time rktvd2 is
**                      called - the calling code can supply the same buffer
**                      each time, greatly reducing memory allocations.
**
**  Input       ::      PARAMETERS:
**
**                      t
**                          Pointer to the independent variable value
**                      x
**                          Array of values of dependent variables
**                      xout
**                          Array into which the incremented x-values will be put
**                      k
**                          Temporary storage buffer for intermediate values
**                      cfl
**                          Specifies the CFL number to use. The derivative function f
**                          must return an estimate of the minimum time dt_w for a wave
**                          to cross a grid cell, or otherwise reach a point that
**                          invalidates the method. The time step will be chosen such
**                          that dt = cfl*dt_w.
**                      numvars
**                          Specifies the size of all input arrays
**                      f
**                          Derivative function
**
**  Output      ::      The value of the independent variable (t) will be incremented by h,
**                      and the new values of the dependent variables will be put in
**                      the xout array.
**
**                      Returns the time step used, based on the supplied CFL number.
**
****************************************************************************************/
real rktvd2 ( real *                  t,
              real *                  x,
              real *                  xout,
              real *                  k,
              const real              cfl,
              const uint              numvars,
              DERIVATIVE_FUNCTION_CFL f)
{
    // For looping
    uint i;

    // Max time step, based on estimate from derivative function
    real dt_w;

    // Actual time step
    real h;

    // Calculate time derivative and obtain max time step
    dt_w = (*f)(*t, x, k, numvars);

    // Calculate time step
    h = cfl*dt_w;

    // Calculate x1 = x+h*f(x,t)
    for (i = 0; i < numvars; i ++)
    {
        k[i] = k[i] * h + x[i];
    }

    // Needed to calculate x2
    (*f)(*t + h, k, xout, numvars);

    // Calculate the result, xout = 0.5*(x + x2), x2 = x1 + h*f(x1,t+dt)
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = 0.5*(x[i] + k[i] + h*xout[i]);
    }

    // Increment t
    *t += h;

    return h;
}

/****************************************************************************************
**
**  Function    ::      rktvd3
**
**  Purpose     ::      Performs a single iteration of the standard Runge-Kutta fourth
**                      order method, incrementing t by the step-size h and approximating
**                      the next set of dependent variable values (x), which are updated
**                      with their new values.
**
**                      The arrays x, xout, and k MUST each have a length at least
**                      as great as numvars, and MUST all be distinct arrays that are not
**                      modified as the method proceeds. If these conditions are not met,
**                      the results of this procedure are undefined.
**
**                      k is a temporary storage buffer used
**                      by the Runge-Kutta algorithm. It is passed in as a parameter
**                      so that memory allocation need not happen every time rktvd3 is
**                      called - the calling code can supply the same buffer
**                      each time, greatly reducing memory allocations.
**
**  Input       ::      PARAMETERS:
**
**                      t
**                          Pointer to the independent variable value
**                      x
**                          Array of values of dependent variables
**                      xout
**                          Array into which the incremented x-values will be put
**                      k
**                          Temporary storage buffer for intermediate values
**                      cfl
**                          Specifies the CFL number to use. The derivative function f
**                          must return an estimate of the minimum time dt_w for a wave
**                          to cross a grid cell, or otherwise reach a point that
**                          invalidates the method. The time step will be chosen such
**                          that dt = cfl*dt_w.
**                      numvars
**                          Specifies the size of all input arrays
**                      f
**                          Derivative function
**
**  Output      ::      The value of the independent variable (t) will be incremented by h,
**                      and the new values of the dependent variables will be put in
**                      the xout array.
**
**                      Returns the time step used, based on the supplied CFL number.
**
****************************************************************************************/
real rktvd3 ( real *                  t,
              real *                  x,
              real *                  xout,
              real *                  k,
              const real              cfl,
              const uint              numvars,
              DERIVATIVE_FUNCTION_CFL f)
{
    // For looping
    uint i;

    // Max time step, based on estimate from derivative function
    real dt_w;

    // Actual time step
    real h;

    // Rename pointers for convenience. Note that we use the xout array to store both
    // x1 and x3, as there is no overlap in the calculation.
    real * x1 = xout;
    real * x2 = k;
    real * x3 = xout;

    // Calculate time derivative and obtain max time step
    dt_w  = (*f)(*t, x, x1, numvars);

    // Calculate time step
    h = cfl*dt_w;

    // Calculate x1 = x + h*f(x)
    for (i = 0; i < numvars; i ++)
    {
        x1[i] = x[i] + h*x1[i];
    }

    // Calculate x2 = (3/4)*x + (1/4)*(x1 + h*f(x1))
    (*f)(*t+h, x1, x2, numvars);
    for (i = 0; i < numvars; i ++)
    {
        x2[i] = 0.75*x[i] + 0.25*(x1[i] + h*x2[i]);
    }

    // Calculate x3 = (1/4)*x + (3/4)*(x2 + h*f(x2))
    (*f)(*t+2*h, x2, x3, numvars);
    for (i = 0; i < numvars; i ++)
    {
        x3[i] = 0.25*x[i] + 0.75*(x2[i] + h*x3[i]);
    }

    // xout == x3, so we're finished with x

    // Increment t
    *t += h;

    return h;
}
