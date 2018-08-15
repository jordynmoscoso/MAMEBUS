%%%
%%% paramTypes.m
%%%
%%% Defines integers used to specify different input parameter types for
%%% MAMEBUS.
%%%

PARM_INT = 1; %%% Integer number
PARM_REALF = 2; %%% Real number, decimal form
PARM_REALE = 3; %%% Real number, floating-point form
PARM_STR = 4; %%% String literal

%%% Time-stepping scheme identifiers
TIMESTEPPING_RKTVD1 = 0;  %%% First-order total-variation-diminishing Runge-Kutta
TIMESTEPPING_RKTVD2 = 1;  %%% Second-order total-variation-diminishing Runge-Kutta
TIMESTEPPING_RKTVD3 = 2;  %%% Third-order total-variation-diminishing Runge-Kutta

%%% Biogeochemical model identifiers
BGC_NONE = 0;             %%% No biogeochemistry
BGC_NITRATEONLY = 1;      %%% Single nitrate model
BGC_NPZD = 1;             %%% Size structured NPZD model

%%% Advection scheme identifiers
ADVECTION_CENTERED = 0;   %%% Centered differencing of tracers
ADVECTION_KT00  = 1;      %%% Kurganov-Tadmor scheme for conservation laws (see K&T 2000)

%%% Momentum scheme identifiers
MOMENTUM_NONE = 0;        %%% No explicit momentum time stepping, mean velocities have prescribed boundary layer structure only
MOMENTUM_TTTW = 1;        %%% Momentum evolves under time-dependend turbulent thermal wind approximation