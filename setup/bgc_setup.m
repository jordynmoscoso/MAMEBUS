%%% 
%%% bgc_setup.m
%%%
%%% Sets biogeochemical parameters to be used in MAMEBUS in the NPZD Model
%%%
%%% The initial values will be output as a matrix, and biogeochemical
%%% parameters are output as a vector to be read into the model. 
%%%


function [params, bgc_init, nbgc] = bgc_setup(ZZ_tr,Nx,Nz)

% length parameters
lp = 5; % micrometers
lz = 10; 

% light and temperature parameters
qsw = 340; % W/m^2
kw = 0.04; % light attenuation of water
kc = 0.01; % light attenuation of phytoplankton (nitrogen units)
Tref = 20; % deg C
r = 0.05;  % temperature dependence

% phytoplankton parameters
au = 2.6; %1/d
bu = -0.45;
umax = au*lp^bu; % maximum uptake
kn = 0.1; % mmol N/m^3
mp = 0.2; % mortality as a fraction of uptake

% zooplankton parameters
ag = 26;
bg = -0.4;
gmax = ag*lz^bg; % maximum grazing rate
lambda = 1.0/3.0; % grazing efficiency
kp = 3; % mmol N/m^3
delta_x = 0.25; % width of grazing profile

ap = 0.65;
bp = 0.56;
preyopt = ap*lz^bp; %optimal predator prey length

% detrital parameters
r_remin = 0.04; % 1/d
wsink = 10; % m/d


%%% save all parameters
params = [lp, lz, qsw, kw, kc, Tref, r, ...
            umax, kn, mp, ...
            gmax, lambda, kp, delta_x, preyopt, ...
            r_remin, wsink];
        
nbgc = length(params);


%%% Create initial conditions
Pcline = 200;
Pmax = 0.1; % mmol/m3
Dmax = 0;

euph_init = Pmax*(tanh(ZZ_tr/Pcline))+0.1;

%%% Initial nitrate profile (Hyperbolic)
Nmax = 30; %%% Maximum concentration of nutrient at the ocean bed
Ncline = 80; % Approximate guess of the depth of the nutracline


bgc_init(:,:,1) = -Nmax*tanh(ZZ_tr/Ncline);
bgc_init(:,:,2) = euph_init;
bgc_init(:,:,3) = 0.1*euph_init;
bgc_init(:,:,4) = Dmax*zeros(Nx,Nz);
end
