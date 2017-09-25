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


function [params, bgc_init,lp,lz] = bgc_setup(model_type,NP,NZ,ND,XX_tr,ZZ_tr)

fignum = 100;

%%% Note, all biogeochemical parameters are calculated in seconds
t1day = 86400; %%% Seconds in 1 day

switch (model_type)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% Nitrate %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 0 %%% Nitrate Only Model
        a_temp = 0.6/t1day;
        b_temp = 1.066;
        c_temp = 1;
        alpha = 0.025/t1day; % s^-1 (W/m^2)^-1
        monod = 0;
        
        params = [a_temp;b_temp;c_temp;alpha;monod];
        
        %%% Initial nitrate profile (Hyperbolic)
        Nmax = 30; %%% Maximum concentration of nutrient at the ocean bed
        Ncline = 150; % Approximate guess of the depth of the nutracline
        bgc_init = -Nmax*tanh(ZZ_tr./Ncline);
    
        figure(fignum);
        fignum = fignum+1;
        pcolor(XX_tr,ZZ_tr,bgc_init);
%      plot(N_init(1,:),ZZ_tr(1,:))
        shading interp
        title('Initial Nitrate Profile')
        colorbar
        caxis([0 Nmax])
  
        figure(fignum);
        fignum = fignum+1;
        plot(bgc_init(1,:),ZZ_tr(1,:))
        title('Initial Nutracline Profile')
        
        lp = 0; lz = 0;
        nparams = length(params);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% NPZ %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 1 %%% Size Structured NPZ model based on Banas
        nparams = 0; 
        
        % Zooplankton Spectra
        minp = 1; % micrometers
        minz = exp(log(minp/0.65)/0.56); % micrometers
        maxz = 460;
        lz = logspace(log10(minz),log10(maxz),NZ);   
        lz = lz.';
        dz = [lz(1) ; lz(2:end) - lz(1:end-1)];
    
        % Phytoplankton Spectra
        maxp = floor(0.65*(maxz)^0.56);
        lp = logspace(log10(minp),log10(maxp),NP);
        lp = lp.';
        dp = [lp(1) ; lp(2:end) - lp(1:end-1)];
        
        %%%%%%%%%%%% Nutrient Parameters %%%%%%%%%%%%%%%
        % maximum uptake rate
        vmax = zeros(NP,1); 
        av = 2.6/t1day; % s^-1
        bv = -0.45;
        nparams = nparams + 1;
    
        % nutrient limiting Monod constant
        monod = zeros(NP,1);
        am = 0.1;
        nparams = nparams + 1;
    
%         %%%%%%%%% Phytoplankton Parameters %%%%%%%%%%%%
%         sink = zeros(NP,1);
%         as = 0;
%         bs = 0.39;
    
        %%%%%%%%% Zooplankton Parameters %%%%%%%%%%%%%%
        maxgraze = zeros(NZ,1);
        amg = 26/t1day; % s^-1
        bmg = -0.4;
        nparams = nparams + 1;
        
        preyopt = zeros(NZ,1);
        ap = 0.65;
        bp = 0.56;
        nparams = nparams + 1;
                
        for ii = 1:NZ
            maxgraze(ii) = amg*lz(ii)^bmg;
            preyopt(ii)  = ap*lz(ii)^bp;
        end
    
        for ii = 1:NP
            size = lp(ii);
        
            vmax(ii) = av*size^bv;
            monod(ii) = am*size;
       
%             sink(ii) = as*size^bs;
        end
        
        kp = 3*ones(NZ,1);
        nparams = nparams + 1;
        
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
           params = [vmax,monod,maxgraze,preyopt,kp];
       else
           if (Ncase == 0) % must add appropriate number of zeros to the zooplankton parameters to create a matrix
               addZ = NP - NZ;
               append = zeros(addZ,1);
               maxgraze = [maxgraze ; append];
               preyopt = [preyopt ; append];
               kp = [kp ; append];
               if length(kp) ~= length(vmax)
                   disp('Parameter vectors are not the same length')
                   return
               end
           else
               addP = NZ - NP;
               append = zeros(addP,1);
               vmax = [vmax ; append];
               monod = [monod; append];
               if length(kp) ~= length(vmax)
                   disp('Parameter vectors are not the same length')
                   return
               end
           end
           params = [vmax,monod,maxgraze,preyopt,kp];
       end
       
       
       
       %%%%%%% Detritus Parameter %%%%%%%
       dparam = zeros(Ntot,1);
       dparam(2) = 10/t1day;  % Ward model m/s ~10m/day
       
       params = [params, dparam];
       
       % Initialize bgc to have a maximum amount of phytoplankton at the
       % top, use the initial value of 1 mmol /m^3 of phytoplankton
       % concentration at the surface for every size class and 0.1 mmol/m^3
       % of zooplankton
       
       spec_tot = NP + NZ + ND + 1; % total species available,
       bgc_init = zeros(40,40,spec_tot); %%% Add one for nitrate
       
       
       %%% Initial values of phytoplankton and zooplankton concentration
       
       pmax = 1;
       zmax = 0.1;
       dmax = 0.3;
       bgc_cline = 50; %m
       
       %%% Initial nitrate profile (Hyperbolic)
       Nmax = 30; %%% Maximum concentration of nutrient at the ocean bed
       Ncline = 150; % Approximate guess of the depth of the nutracline
       
       bgc_init(:,:,1) = -Nmax*tanh(ZZ_tr./Ncline);
       
       for ii = 2:1:NP
           bgc_init(:,:,ii) = pmax*(1-tanh(ZZ_tr./bgc_cline));
       end
       
       for ii = NP+2:NZ+1
           bgc_init(:,:,ii) = zmax*(1-tanh(ZZ_tr./bgc_cline));
       end
       
       for ii = NZ+2:ND+1
           bgc_init(:,:,ii) = dmax*(1-tanh(ZZ_tr./bgc_cline));
       end
       

       
end   
end