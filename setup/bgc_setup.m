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

switch (model_type)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% Nitrate %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 1 %%% Nitrate Only Model
        umax = 2.6/t1day; % maximum uptake rate
        f = 0.1;          % fraction of exported material
        
        params = [umax,f];
        
        %%% Initial nitrate profile (Hyperbolic)
        Nmax = 40; %%% Maximum concentration of nutrient at the ocean bed
        Ncline = 75; % Approximate guess of the depth of the nutracline
        bgc_init = -Nmax*tanh(ZZ_tr./Ncline);
    
        nparams = length(params);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% NPZ %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 2 %%% NPZD Model

        %%% Here's another model. The HADOCC model.
        % can be updated to keep track of dic and alk
        
        % nitrate parameters
        cp = 6.625;
        cz = 5.626;
        cd = 7.5;
        
        % phytoplankton parameters
        umax  = 0.8/t1day;     % 1/d max photosynthetic rate
        pi_a  =  0.0055/t1day; % m^2 /(W d) initial slope of PI curve
        q10   = 2.2;             % temp-dependent growth
        kn    = 0.1;           % mmol/m^3 half saturation
        mu1p  = 0.02/t1day;    % 1/d linear mortality
        mu2p  = 0.05/t1day;    % m^3/mmol d quadratic mortality
        
        % zooplantkon parameters
        gmax  = 1.0/t1day;     % 1/d max grazing rate
        kg    = 0.75;          % half saturation
        gthr  = 0.1/t1day;     % 1/d grazing threshold
        bp    = 0.3;           % assimilation efficiency
        bd    = 0.5;
        mu1z  = 0.05/t1day;    % 1/d linear mortality
        mu2z  = 0.2/t1day;     % m^3/mmol d density dependent mortality
        
        % detritus
        rs    = 0.075/t1day;     % remin shallow
        rd    = 0.02/t1day;    % remin deep
        wsink = 10/t1day;      % m/d sinking speed
        
        
        params = [cp, cz, cd, ...
                    umax, pi_a, q10, kn, mu1p, mu2p, ...
                    gmax, kg, gthr, bp, bd, mu1z, mu2z, ...
                    rs, rd, wsink];
        nparams = length(params);
        
        spec_tot = NP + NZ + ND + 1; % total species available,
        bgc_init = zeros(Nx,Nz,spec_tot); %%% Add one for nitrate
        
        Pcline = 200;
        Pmax = 0.1; % mmol/m3
        Dmax = 0;
        
%         figure(301)
%         pcolor(XX_tr,ZZ_tr,bgc_init(:,:,1))
%         title('Initial Nitrate Concentration')
%         colorbar
%         shading interp
        
        euph_init = Pmax*(tanh(ZZ_tr/Pcline))+0.1;

        
        %%% Initial nitrate profile (Hyperbolic)
        Nmax = 25; %%% Maximum concentration of nutrient at the ocean bed
        Ncline = 100; % Approximate guess of the depth of the nutracline
        bgc_init(:,:,1) = -Nmax*tanh(ZZ_tr/Ncline);
%         bgc_init(:,:,1) = Nmax*ones(Nx,Nz);
        
%         figure(300)
%         pcolor(XX_tr,ZZ_tr,bgc_init(:,:,1))
%         hold on
%         [C, h] = contour(XX_tr,ZZ_tr,bgc_init(:,:,1),[1 6.5 12.1 17.6 23.1],'-k');
%         clabel(C,h)
%         hold off
%         title('Initial Nitrate Concentration')
%         colorbar
%         shading interp
        
        
        bgc_init(:,:,2) = euph_init;
        bgc_init(:,:,3) = 0.1*euph_init;
        bgc_init(:,:,4) = Dmax*zeros(Nx,Nz);
        
        
%         figure(302)
%         pcolor(XX_tr,ZZ_tr,euph_init)
%         shading interp
%         colorbar
% %         caxis([0 Pmax])
%         title('Initial Phytoplankton Concentration')
%         
%         
%         figure(302)
%         pcolor(XX_tr,ZZ_tr,bgc_init(:,:,3))
%         shading interp
%         colorbar
% %         caxis([0 Pmax])
%         title('Initial Zooplankton Concentration')

    case 3 %%% Size Structured NPZ model based on Banas
        nparams = 0; 
        
        % Zooplankton Spectra
        minp = 1; % micrometers
        minz = exp(log(minp/0.65)/0.56); % micrometers
        maxz = 460;
        lz = logspace(log10(minz),log10(maxz),NZ);   
        lz = lz.';
    
        % Phytoplankton Spectra
        maxp = floor(0.65*(maxz)^0.56);
        lp = logspace(log10(minp),log10(maxp),NP);
        lp = lp.';
        
        %%%%%%%%%%%% Nutrient Parameters %%%%%%%%%%%%%%%
        % maximum uptake rate
        vmax = zeros(NP,1); 
        av = 2.6/t1day; % s^-1
        bv = -0.45;
    
        % nutrient limiting Monod constant
        monod = zeros(NP,1);
        am = 0.1;
    
%         %%%%%%%%% Phytoplankton Parameters %%%%%%%%%%%%
%         as = 0;
%         bs = 0.39;
    
        %%%%%%%%% Zooplankton Parameters %%%%%%%%%%%%%%
        maxgraze = zeros(NZ,1);
        amg = 26/t1day; % s^-1
        bmg = -0.4;
        
        preyopt = zeros(NZ,1);
        ap = 0.65;
        bp = 0.56;
                
        for ii = 1:NZ
            maxgraze(ii) = amg*lz(ii)^bmg;
            preyopt(ii)  = ap*lz(ii)^bp;
        end
    
        for ii = 1:NP
            vmax(ii) = av*lp(ii)^bv;
            monod(ii) = am*lp(ii);
       
%             sink(ii) = as*lp(ii)^bs;
        end
        
        kp = 3*ones(NZ,1);
        
        %%% Determine of NP or NZ is larger, 
       Ntot = 0;
       Ncase = -1;
       if (NP > NZ)
           Ntot = NP; % More phytoplankton size classes
           Ncase = 0;
       else
           Ntot = NZ; % Equal size classes or more zooplankton
           Ncase = 1;
       end
       
       %%% Build the params matrix
       params = zeros(Ntot,nparams);
       
       if (NP == NZ)
           params = [lp;lz;vmax;monod;maxgraze;preyopt;kp];
       else
           if (Ncase == 0) % must add appropriate number of zeros to the zooplankton parameters to create a matrix
               addZ = NP - NZ;
               append = zeros(addZ,1);
               maxgraze = [maxgraze ; append];
               preyopt = [preyopt ; append];
               kp = [kp ; append];
               lp = [lp; append];
               lz = [lz; append];
               if length(kp) ~= length(vmax)
                   disp('Parameter vectors are not the same length')
                   return
               end
           else
               addP = NZ - NP;
               append = zeros(addP,1);
               vmax = [vmax ; append];
               monod = [monod; append];
               lp = [lp; append];
               lz = [lz; append];
               if length(kp) ~= length(vmax)
                   disp('Parameter vectors are not the same length')
                   return
               end
           end
           params = [lp; lz; vmax; monod; maxgraze; preyopt; kp];
       end
       
       
       
       
       % Initialize bgc to have a maximum amount of phytoplankton at the
       % top, use the initial value of 1 mmol /m^3 of phytoplankton
       % concentration at the surface for every size class and 0.1 mmol/m^3
       % of zooplankton
       
       spec_tot = NP + NZ + ND + 1; % total species available,
       bgc_init = zeros(40,40,spec_tot); %%% Add one for nitrate
       
       
       %%% Initial values of phytoplankton and zooplankton concentration
       %%% These worked in the matlab model, for now.
       pmax = 2;
       zmax = 1;
       dmax = 5;
       bgc_cline = 50; %m
       
       %%% Initial nitrate profile (Hyperbolic)
       Nmax = 15; %%% Maximum concentration of nutrient at the ocean bed
       
       bgc_init(:,:,1) = Nmax*ones(size(ZZ_tr));
       
       for ii = 1:NP
           bgc_init(:,:,ii+1) = pmax*ones(size(ZZ_tr));
       end
       
       for ii = 1:NZ
           bgc_init(:,:,ii+1+NP) = zmax*ones(size(ZZ_tr));
       end
       
       for ii = 1:ND
           bgc_init(:,:,ii+1+NP+NZ) = dmax*(1-tanh(ZZ_tr./bgc_cline));
       end
       
       nparams = length(params);

end   
end
