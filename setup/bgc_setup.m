%%% 
%%% bgc_setup.m
%%%
%%% Sets biogeochemical parameters to be used in MAMEBUS in the NPZD Model
%%%
%%% The initial values will be output as a matrix, and biogeochemical
%%% parameters are output as a vector to be read into the model. 
%%%


function [params, bgc_init, bgc_restore, nbgc, lp, lz, NP, NZ, bgcRates, nallo, idxAllo] = bgc_setup(XX_tr,ZZ_tr,Lx,Nx,Nz,modeltype,MP,MZ,data_dir,dx,nexp,nsml,nmax)

% light and temperature parameters
qsw = 340; % W/m^2
kw = 0.04; % light attenuation of water
kc = 0.01; % light attenuation of phytoplankton (nitrogen units)
Tref = 18; % deg C (same as tmax)
r = 0.05;  % temperature dependence

% phytoplankton parameters
au = 2.6; %1/d
bu = -0.45;
kn = 0.1; % mmol N/m^3
mp = 0.02; % mortality as a fraction of uptake

% zooplankton parameters
ag = 26;
bg = -0.4;
eff = 1.0/3.0; % grazing efficiency
seff = 1.0/3.0; % self grazing efficiency
kp = 3; % mmol N/m^3
bl = 0.56;
al = 0.65;
delta_x = dx; % width of grazing profile
idxAllo = 4;
pdia = 1e-9;

%%% Depending on the model type, return size grid and the data grid. 
[lp, lz, NP, NZ, bgcRates, idxUptake, idxSat, idxGrazing, idxPredPrey, zeta] ...
    = buildSizeGrid(modeltype,MP,MZ,delta_x,data_dir, au, bu, kn, ag, bg, al, bl ,eff);

nallo = max(NP,NZ);

% detrital parameters
r_remin = 0.04; % 1/d
wsink = 10; % m/d


%%% save all parameters
params = [qsw, kw, kc, Tref, r, ... % physical parameters
            mp, pdia, ... % singular phytoplankton parameters
            eff, seff, kp, delta_x, zeta ... % singular zooplankton parameters
            r_remin, wsink ... % detritus parameters
            idxUptake, idxSat, idxGrazing, idxPredPrey]; % indices for allometric parameters
        
nbgc = length(params);


%%% Create initial conditions
Pcline = 180;
Pmax = 2.5/MP; % mmol/m3
Dmax = 0.001;

%%% Initial Nitrate
% Hexp = 120;
Hexp = nexp;
Nmin = 2*XX_tr/Lx;
Nmax = nmax;
Hsml = nsml;
Ninit =  Nmin-Nmax*tanh((ZZ_tr + Hsml)./Hexp);
Ninit(Ninit < 0) = 0;

euph_init = Pmax*(tanh(ZZ_tr/Pcline))+Pmax;


ind = 1;
% nutrients
bgc_init(:,:,ind) = Ninit; ind = ind+1;

% phytoplankton
for ii = 1:NP
    bgc_init(:,:,ind) = euph_init; ind = ind+1;
end

%zooplankton
for ii = 1:NZ
    bgc_init(:,:,ind) = 0.1*euph_init; ind = ind +1;
end

% detritus
bgc_init(:,:,ind) = Dmax*ones(Nx,Nz);


% restoring Nitrate
% Hexp = nexp;
% Nmin = 2*XX_tr/Lx;
% Nmax = 30;
% Hsml = 25;
Nrestore =  Nmin-Nmax*tanh((ZZ_tr + Hsml)./Hexp);
Nrestore(Nrestore < 0) = 0;

bgc_restore(:,:,1) = Nrestore;
end
