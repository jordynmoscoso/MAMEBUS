/****************************************************************************************
 **
 **  File        ::      ab.c
 **
 **  Author      ::      Jordyn Moscoso
 **
 **  Description ::      Provides a series of Adams-Bashforth methods of various orders
 **                      with variable time step sizes for use in numerical integration
 **                      problems.
 **
 ****************************************************************************************/
#include "ab.h"


/****************************************************************************************
 **
 **  Function    ::      ab1
 **
 **  Purpose     ::      Performs a single iteration of the standard Adams-Bashforth first
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
 **                      dxdt
 **                          Storage buffer for the time derivative of x
 **                      h
 **                          Step size to take
 **                      numvars
 **                          Specifies the size of all input arrays
 **                      f
 **                          Derivative function
 **
 **  Output      ::      The value of the independent variable (t) will be incremented by h,
 **                      and the new values of the dependent variables will be put in
 **                      the xout array.
 **
 ****************************************************************************************/
real ab1 (  real *                      t,
            real *                      x,
            real *                      xout,
            real *                      dxdt,
            const real                  cfl,
            const uint                  numvars,
            DERIVATIVE_FUNCTION_CFL     f)
{
    // For looping
    uint i;
    
    // Actual time step
    real h;
    real dt_w;
    
    // Calculate the derivatives of x at t
    dt_w = (*f)(*t, x, dxdt, numvars);
    
    // Calculate the timestep
    h = cfl*dt_w;
#pragma parallel
    
    // Increment the dependent variable
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = dxdt[i]*h + x[i];
    }
    
    // Increment t
    *t += h;
    
    return h;
}







/****************************************************************************************
 **
 **  Function    ::      ab2
 **
 **  Purpose     ::      Performs a single iteration of the standard Adams-Bashforth
 **                      second order method, incrementing t by the variable step-size h
 **                      and h1, and approximating the next set of dependent variable
 **                      values (x), which are updated with their new values.
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
 **                      dxdt
 **                          Storage buffer for the time derivative of x
 **                      dxdt_1
 **                          Array containing the time derivative of x from the
 **                          previous iteration
 **                      h
 **                          The stepsize for the current iteration
 **                      h1
 **                          The stepsize for the previous iteration
 **                      numvars
 **                          Specifies the size of all input arrays
 **                      f
 **                          Derivative function
 **
 **  Output      ::      The value of the independent variable (t) will be incremented by h,
 **                      and the new values of the dependent variables will be put in
 **                      the xout array.
 **
 ****************************************************************************************/
real ab2 (  real *                      t,
            real *                      x,
            real *                      xout,
            real *                      dxdt,
            real *                      dxdt_1,
            const real                  cfl,
            const real                  h1,
            const uint                  numvars,
            DERIVATIVE_FUNCTION_CFL     f)
{
    // For looping
    uint i;
    
    // Real cfl condition
    real h;
    real dt_w;
    
    // Calculate dx/dt = f(x,t)
    dt_w = (*f)(*t, x, dxdt, numvars);
    
    // Calculate the current time step
    h = cfl*dt_w;
    
    
    // To avoid multiple calculation
    real c0 = h/(2*h1);
    
#pragma parallel
    
    // Calculate the result
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = x[i] + c0*( 2*h1*dxdt[i] + h*dxdt[i] - h*dxdt_1[i]);
    }
    
    // The new time derivative becomes the old time derivative
    memcpy(dxdt_1,dxdt,numvars*sizeof(real));
    
    // Increment t
    *t += h;
    
    return h;
}





/****************************************************************************************
 **
 **  Function    ::      ab3
 **
 **  Purpose     ::      Performs a single iteration of the standard Adams-Bashforth
 **                      third order method, incrementing t by the variable step-size h,
 **                      h1, and h2, and approximating the next set of dependent variable
 **                      values (x), which are updated with their new values.
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
 **                      dxdt
 **                          Storage buffer for the time derivative of x
 **                      dxdt_1
 **                          Array containing the time derivative of x from the
 **                          previous iteration
 **                      dxdt_2
 **                          Array containing the time derivative of x from the
 **                          second previous iteration
 **                      h
 **                          The stepsize for the current iteration
 **                      h1
 **                          The stepsize for the previous iteration
 **                      h2
 **                          The stepsize for the second previous iteration
 **                      numvars
 **                          Specifies the size of all input arrays
 **                      f
 **                          Derivative function
 **
 **  Output      ::      The value of the independent variable (t) will be incremented by h,
 **                      and the new values of the dependent variables will be put in
 **                      the xout array.
 **
 ****************************************************************************************/
real ab3 (  real *                      t,
            real *                      x,
            real *                      xout,
            real *                      dxdt,
            real *                      dxdt_1,
            real *                      dxdt_2,
            const real                  cfl,
            const real                  h1,
            const real                  h2,
            const uint                  numvars,
            DERIVATIVE_FUNCTION_CFL     f)
{
    
    // For looping
    uint i;
    
    //Real cfl condition
    real h;
    real dt_w;
    
    // Calculate dx/dt = f(x,t)
    dt_w = (*f)(*t, x, dxdt, numvars);
    
    // Calculate the current timestep
    h = dt_w*cfl;
    
    // To avoid multiple calculation
    real c2 = ( h*h )*( 2*h + 3*h1 )/( 6*h2*( h1+h2 ) );
    real c1 = ( h*h )*( 2*h + 3*h1 + 3*h1 )/( 6*h1*h2 );
    real c0 = h*( 2*h*h + 6*h*h1 + 3*h*h2 + 6*h1*h1 + 6*h1*h2 )/( 6*( h1 + h2)*h1 );
    
#pragma parallel
    
    // Calculate the result
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = x[i] + c0*dxdt[i] - c1*dxdt_1[i] + c2*dxdt_2[i];
    }
    
    // The new time derivative becomes the old time derivative
    memcpy(dxdt_2,dxdt_1,numvars*sizeof(real));
    memcpy(dxdt_1,dxdt,numvars*sizeof(real));
    
    // Increment t
    *t += h;
    
    return h;
}




// TO DO: CHECK THE VALIDITY
/****************************************************************************************
 **
 **  Function    ::      ab4
 **
 **  Purpose     ::      Performs a single iteration of the standard Adams-Bashforth
 **                      fourth order method, incrementing t by the variable step-size h,
 **                      h1, h2, and h3, and approximating the next set of dependent
 **                      variable values (x), which are updated with their new values.
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
 **                      dxdt
 **                          Storage buffer for the time derivative of x
 **                      dxdt_1
 **                          Array containing the time derivative of x from the
 **                          previous iteration
 **                      dxdt_2
 **                          Array containing the time derivative of x from the
 **                          second previous iteration
 **                      dxdt_3
 **                          Array containing the time derivative of x from the
 **                          third previous iteration
 **                      h
 **                          The stepsize for the current iteration
 **                      h1
 **                          The stepsize for the previous iteration
 **                      h2
 **                          The stepsize for the second previous iteration
 **                      h3
 **                          The stepsize for the third previous iteration
 **                      numvars
 **                          Specifies the size of all input arrays
 **                      f
 **                          Derivative function
 **
 **  Output      ::      The value of the independent variable (t) will be incremented by h,
 **                      and the new values of the dependent variables will be put in
 **                      the xout array.
 **
 ****************************************************************************************/
/*
void ab4 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            real *                  dxdt_1,
            real *                  dxdt_2,
            real *                  dxdt_3,
            const real              h,
            const real              h1,
            const real              h2,
            const real              h3,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f)
{
    
    // For looping
    uint i;
    
    // To avoid multiple calculation
    c3 = ( h*h )*( 3*h*h + 8*h*h1 + 4*h*h2 + 6*h1*h1 + 6*h1*h2 )/( 12*( h1 + h2 + h3 )*( h2 +h3 )*h3 );
    c2 = ( h*h )*( 8*h*h1 + 4*h*h2 + 4*h*h3 + 6*h1*h2 + 6*h1*h3 + 3*h*h + 6*h1*h1 )/( 12*h3*h2( h2 + h3 ) );
    c1 = ( h*h )*( 3*h*h + 8*h*h1 + 8*h*h2 + 4*h*h3 + 6*h1*h1 + 12*h1*h2 + 6*h1*h3 + 6*h2*h2 + 6*h2*h3 )/( 12*( h2 + h3 )*h1*h2 );
    c0 = h*( 3*h*h*h + 12*h*h*h1 + 8*h*h*h2 + 4*h*h*h3 + 18*h1*h1*h + 24*h*h1*h2 + 12*h*h1*h3 + 6*h*h2*h3 + 12*h1*h1*h1 + 24*h1*h1*h2 + 12*h1*h1*h3 + 12*h1*h2*h2 + 12*h1*h2*h3 )/ ( 12*( h1 + h2 + h3 )*( h1 + h2 )*h1 );
    
    // Calculate dx/dt = f(x,t)
    (*f)(*t, x, dxdt, numvars);
    
#pragma parallel
    
    // Calculate the result
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = x[i] + c0*dxdt[i] - c1*dxdt_1[i] + c2*dxdt_2[i] - c3*dxdt_3[i];
    }
    
    // The new time derivative becomes the old time derivative
    memcpy(dxdt_3,dxdt_2,numvars*sizeof(real));
    memcpy(dxdt_2,dxdt_1,numvars*sizeof(real));
    memcpy(dxdt_1,dxdt,numvars*sizeof(real));
    
    // Increment t
    *t += h;
}
*/
