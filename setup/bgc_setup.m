%%% 
%%% bgc_setup.m
%%%
%%% Sets biogeochemical parameters to be used in MAMEBUS. model_type
%%% indicates the complexity of the biogeochemical model to be used. Based
%%% on the type of model indicated, the output will create a vector with
%%% the appropriate uptake values, and intial conditions.
%%%
%%% The initial values will be output as either a vector (for the nitrate
%%% only model of a matrix of initial values for phytoplankton and
%%% zooplankton).
%%%


function [params, bgc_init,nparams] = bgc_setup(model_type,NP,NZ,ND,XX_tr,ZZ_tr,Nx,Nz)
%%% Note, all biogeochemical parameters are calculated in seconds
t1day = 24*60*60; %%% Seconds in 1 day

Pcline = 200;
Pmax = 0.1; % mmol/m3
Dmax = 0;

euph_init = Pmax*(tanh(ZZ_tr/Pcline))+0.1;


%%% Initial nitrate profile (Hyperbolic)
Nmax = 30; %%% Maximum concentration of nutrient at the ocean bed
Ncline = 80; % Approximate guess of the depth of the nutracline
bgc_init(:,:,1) = -Nmax*tanh(ZZ_tr/Ncline);

        figure(300)
        pcolor(XX_tr,ZZ_tr,bgc_init(:,:,1))
        hold on
        [C, h] = contour(XX_tr,ZZ_tr,bgc_init(:,:,1),[1 6.5 12.1 17.6 23.1],'-k');
        clabel(C,h)
        hold off
        title('Initial Nitrate Concentration')
        colorbar
        shading interp


bgc_init(:,:,2) = euph_init;
bgc_init(:,:,3) = 0.1*euph_init;
bgc_init(:,:,4) = Dmax*zeros(Nx,Nz);

params = 0; 
nparams = 0;


end
