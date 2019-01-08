%%%
%%% paramTypes.m
%%%
%%% Defines integers used to specify different parameter types for
%%% MAMEBUS.
%%%

%%% Input parameter types
PARM_INT = 1;             %%% Integer number
PARM_REALF = 2;           %%% Real number, decimal form
PARM_REALE = 3;           %%% Real number, floating-point form
PARM_STR = 4;             %%% String literal

%%% Time-stepping scheme identifiers
TIMESTEPPING_RKTVD1 = 0;  %%% First-order total-variation-diminishing Runge-Kutta
TIMESTEPPING_RKTVD2 = 1;  %%% Second-order total-variation-diminishing Runge-Kutta
TIMESTEPPING_RKTVD3 = 2;  %%% Third-order total-variation-diminishing Runge-Kutta

%%% Biogeochemical model identifiers
BGC_NONE = 0;             %%% No biogeochemistry
BGC_NITRATEONLY = 1;      %%% Single nitrate model
BGC_NPZD = 2;             %%% NPZD model
BGC_SSEM = 3;             %%% Size structured NPZD model

%%% Advection scheme identifiers
ADVECTION_CENTERED = 0;   %%% Centered differencing of tracers
ADVECTION_KT00  = 1;      %%% Kurganov-Tadmor scheme for conservation laws (see K&T 2000)

%%% Momentum scheme identifiers
MOMENTUM_NONE = 0;        %%% No explicit momentum time stepping, mean velocities have prescribed boundary layer structure only
MOMENTUM_TTTW = 1;        %%% Momentum evolves under time-dependend turbulent thermal wind approximation

%%% Pressure scheme identifiers
PRESSURE_LINEAR = 0;
PRESSURE_CUBIC = 1;

%%% Fixed tracer array indices
IDX_UVEL = 1;             %%% Zonal velocity
IDX_VVEL = 2;             %%% Meridional velocity
IDX_BUOY = 3;             %%% Buoyancy
IDX_NITRATE = 4;          %%% Nitrate

