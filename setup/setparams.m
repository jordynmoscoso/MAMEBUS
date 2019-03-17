%%%
%%% setparams.m
%%%
%%% Sets parameters for MAMEBUS, the Meridionally-Averaged Model of
%%% Eastern Boundary Upwelling Systems.
%%%
%%% run_name specifies the name of the simulation. All simulation files
%%% except the executable will be placed in a folder of this name within
%%% the 'runs' directory.
%%%
%%% local_home_dir specifies the directory in the local system into which
%%% the run files will be written. N.B. a directory called 'run_name' will
%%% be created within local_home_dir to house the files.
%%%
function setparams (local_home_dir,run_name)  

%%% TODO move depth tracer and dye tracer to end of tracer matrix as
%%% examples of how to prescribe arbitrary tracer inputs

  %%% Convenience scripts used in this function
  addpath ../utils;
  
  %%% Load globally-defined constants
  paramTypes;
  
  
%%% TODO this definitely ought to appear later in this script


  
  %%% The number of biogeochemical classes are entered here. 
  modeltype = BGC_NONE; %%% This automatically defaults so that the model runs without biogeochemistry

  switch (modeltype)
    case BGC_SSEM
      MN = 2; %%% The number of nutrients in the model (must be 2) one active one dye.
      MP = 5;
      MZ = 5;
      MD = 2; %%% Currently this variable is not set to change, and more than two size classes are not resolved.
      bio = MP+MZ+MD;
      disp(['Number of: (Phytoplankton, Zooplankton) = (',num2str(MP),', ', num2str(MZ),')']);
    case BGC_NPZD
      MN = 2; %%% one active dye
      MP = 1;
      MZ = 1;
      MD = 1;
      bio = MP+MZ+MD;
      disp(['Number of: (Phytoplankton, Zooplankton) = (',num2str(MP),', ', num2str(MZ),')']);
    case BGC_NITRATEONLY
      MN = 2; %%% The number of nutrients in the model (must be 2) one active one dye.
      MP = 0;
      MZ = 0;
      MD = 0;
      bio = 0;
    case BGC_NONE
      MN = 0;
      MP = 0;
      MZ = 0;
      MD = 0;
      bio = 0;
  end
  Nbgc = bio + MN; 

  %%% Check to see if a valid model type is indicated for biogeochemistry,
  %%% if not use the default single nitrate model (modeltype = 0)
%   if (MP < 1 || MZ < 1)
%       modeltype = 0;
%   end
  

%   pressureScheme = PRESSURE_CUBIC;
  pressureScheme = PRESSURE_LINEAR;
        
  
  %%% For plotting figures of setup
  fignum = 1;

  %%% If set true, set up this run for the cluster
  use_cluster = false;
  walltime = 24;
  cluster_home_dir = '/data1/astewart/MAMEBUS/runs';
  cluster_username = 'astewart';
  cluster_address = 'ardbeg.atmos.ucla.edu';
  
  %%% Run directory
  run_name = strtrim(run_name);  
  exec_name = 'mamebus.exe';      
  local_run_dir = fullfile(local_home_dir,run_name);
  pfname = fullfile(local_run_dir,[run_name,'_in']);   
  mkdir(local_run_dir);

  %%% To store parameters  
  PARAMS = {};
  
  %%% Time parameters
  t1day = 86400; %%% Seconds in 1 day
  t1year = 365*t1day; %%% Seconds in 1 year
  endTime = 50*t1year;
  restart = false;
  startIdx = 15;
  outputFreq = 0.01*t1day;
    
  %%% Domain dimensions
  m1km = 1000; %%% Meters in 1 km    
  H = 3*m1km; %%% Depth, excluding the mixed layer
  Lx = 300*m1km; %%% Computational domain width
  
  %%% Scalar parameter definitions 
  tau0 = -1e-1; %%% Northward wind stress (N m^{-2})
  rho0 = 1025; %%% Reference density
  f0 = 1e-4; %%% Coriolis parameter (CCS)
  Kgm0 = 500; %%% Reference GM diffusivity
  Kiso0 = 2000; %%% Reference surface isopycnal diffusivity m^2/s
  Kiso_hb = 200; %%% Reference interior isopycnal diffusivity
  Kdia0 = 1e-5; %%% Reference diapycnal diffusivity
  Cp = 4e3; %%% Heat capacity
  g = 9.81; %%% Gravity
  s0 = tau0/rho0/f0/Kgm0; %%% Theoretical isopycnal slope    
  Hsml = 60; %%% Surface mixed layer thickness
  Hbbl = 60; %%% Bottom boundary layer thickness
  r_bbl = 1e-3; %%% Bottom boundary layer drag coefficient
  
  %%% Biogeochemical Parameters

  %%% Grid parameters
  h_c = 300; %%% Sigma coordinate surface layer thickness parameter (must be > 0)
  theta_s = 6; %%% Sigma coordinate surface stretching parameter (must be in [0,10])
  theta_b = 4; %%% Sigma coordinage bottom stretching parameter (must be in [0,4])
  
  %%% Grid parameters (no stretching)
%   h_c = 1e16; %%% Sigma coordinate surface layer thickness parameter (must be > 0)
%   theta_s = 0; %%% Sigma coordinate surface stretching parameter (must be in [0,10])
%   theta_b = 0; %%% Sigma coordinage bottom stretching parameter (must be in [0,4])
   
  %%% Grids  
  Nphys = 3; %%% Number of physical tracers (u-velocity, v-velocity and buoyancy)
  Ntracs = Nphys + 1 + Nbgc; %%% Number of tracers (physical plus bgc plus any other user-defined tracers)
  Nx = 40; %%% Number of latitudinal grid points 
  Nz = 40; %%% Number of vertical grid points
  dx = Lx/Nx; %%% Latitudinal grid spacing (in meters)
  xx_psi = 0:dx:Lx; %%% Streamfunction latitudinal grid point locations
  xx_tr = dx/2:dx:Lx-dx/2; %%% Tracer latitudinal grid point locations  
  xx_topog = [-dx/2 xx_tr Lx+dx/2]; %%% Topography needs "ghost" points to define bottom slope
  
  %%% Create tanh-shaped topography
  shelfdepth = 3000;
  disp(['Shelf Depth: ', num2str(shelfdepth)])
  if shelfdepth < 50
      disp('Shelf is smaller than sml and bbl')
      return
  end
  
  Xtopog = 200*m1km;
  Ltopog = 25*m1km;
  Htopog = H-shelfdepth;  
  hb = H - Htopog*0.5*(1+tanh((xx_topog-Xtopog)/(Ltopog)));
%   hb = H*ones(size(hb));
  hb_psi = 0.5*(hb(1:end-1)+hb(2:end));  
  hb_tr = hb(2:end-1);
  
  % calculate shelf slope:
  hb_slope = (hb_psi(2:end) - hb_psi(1:end-1))./(dx);
  min(hb_slope)
  
  %%% Generate full sigma-coordinate grids
  [XX_tr,ZZ_tr,XX_psi,ZZ_psi,XX_u,ZZ_u,XX_w,ZZ_w] ...
                    = genGrids(Nx,Nz,Lx,h_c,theta_s,theta_b,hb_tr,hb_psi);  % Full output [XX_tr,ZZ_tr,XX_psi,ZZ_psi,XX_u,ZZ_u,XX_w,ZZ_w]
  slopeidx = max((hb_psi>Htopog/2));
  disp(['slopeidx = ',num2str(slopeidx)])
  disp(['Vertical grid spacing at (',num2str(XX_psi(1,1)),',',num2str(ZZ_psi(1,1)),'): ',num2str(ZZ_psi(1,2)-ZZ_psi(1,1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(1,end)),',',num2str(ZZ_psi(1,end)),'): ',num2str(ZZ_psi(1,end)-ZZ_psi(1,end-1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(end,1)),',',num2str(ZZ_psi(end,1)),'): ',num2str(ZZ_psi(end,2)-ZZ_psi(end,1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(end,end)),',',num2str(ZZ_psi(end,end)),'): ',num2str(ZZ_psi(end,end)-ZZ_psi(end,end-1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(slopeidx,1)),',',num2str(ZZ_psi(slopeidx,1)),'): ',num2str(ZZ_psi(slopeidx,2)-ZZ_psi(slopeidx,1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(slopeidx,end)),',',num2str(ZZ_psi(slopeidx,end)),'): ',num2str(ZZ_psi(slopeidx,end)-ZZ_psi(slopeidx,end-1))])
  disp(['Max Topographic Slope : ' num2str(min(hb_slope))])
  
  %%% ZZ_tr size: 40 40 (centers)
  %%% ZZ_psi size: 41 41 (edges) n = 0 is base, n = N is top
  
  %%% Calculate grid stiffness  
  rx1 = abs(diff(0.5*(ZZ_psi(:,1:Nz)+ZZ_psi(:,2:Nz+1)),1,1) ./ diff(0.5*(ZZ_psi(1:Nx,:)+ZZ_psi(2:Nx+1,:)),1,2) );
  disp(['Grid stiffness: ' num2str(max(max(rx1)))]  )
  
  %%% Define parameter
  PARAMS = addParameter(PARAMS,'Ntracs',Ntracs,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nx',Nx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nz',Nz,PARM_INT);  
  PARAMS = addParameter(PARAMS,'Lx',Lx,PARM_REALF);
  PARAMS = addParameter(PARAMS,'Lz',H,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'cflFrac',0.1,PARM_REALF);
  PARAMS = addParameter(PARAMS,'endTime',endTime,PARM_REALF);
  PARAMS = addParameter(PARAMS,'monitorFrequency',outputFreq,PARM_REALF);
  PARAMS = addParameter(PARAMS,'restart',restart,PARM_INT);
  PARAMS = addParameter(PARAMS,'startIdx',startIdx,PARM_INT);
  PARAMS = addParameter(PARAMS,'rho0',rho0,PARM_REALF);
  PARAMS = addParameter(PARAMS,'f0',f0,PARM_REALF);    
  PARAMS = addParameter(PARAMS,'h_c',h_c,PARM_REALE);    
  PARAMS = addParameter(PARAMS,'theta_s',theta_s,PARM_REALF);    
  PARAMS = addParameter(PARAMS,'theta_b',theta_b,PARM_REALF);    
  PARAMS = addParameter(PARAMS,'Hsml',Hsml,PARM_REALF);    
  PARAMS = addParameter(PARAMS,'Hbbl',Hbbl,PARM_REALF);
  PARAMS = addParameter(PARAMS,'r_bbl',r_bbl,PARM_REALF);
  
  %%% Indicate number of phytoplankton, zooplankton and detrital pools
  PARAMS = addParameter(PARAMS,'bgcModel',modeltype,PARM_INT);
  PARAMS = addParameter(PARAMS,'pressureScheme',pressureScheme,PARM_INT);
  PARAMS = addParameter(PARAMS,'MP',MP,PARM_INT);
  PARAMS = addParameter(PARAMS,'MZ',MZ,PARM_INT);
  
  %%% Save biogeochemical parameters in vector form call bgc_setup function
  %%% TODO do we really need a switch here?
  switch(modeltype)
      case BGC_NONE
        bgc_params = [];
        nbgc = 0;
        bgc_init = [];
      case BGC_NITRATEONLY
        [bgc_params, bgc_init,nbgc] = bgc_setup(modeltype,MP,MZ,MD,XX_tr,ZZ_tr,Nx,Nz);
        
        disp('Nitrate only')
      case BGC_NPZD 
         [bgc_params, bgc_init,nbgc] = bgc_setup(modeltype,MP,MZ,MD,XX_tr,ZZ_tr,Nx,Nz);
          
         disp('NPZD Model')
      case BGC_SSEM
        [bgc_params, bgc_init, nbgc] = bgc_setup(modeltype,MP,MZ,MD,XX_tr,ZZ_tr,Nx,Nz);
        disp('NPZD')
        %%% Store phytoplankton size and zooplankton size to determine what size
        %%% pool of detritus they go into (large or small) when passed into
        %%% mamebus.c code.
  end
  
  bgcFile = 'bgcFile.dat';
  writeDataFile(fullfile(local_run_dir,bgcFile),bgc_params);
  PARAMS = addParameter(PARAMS,'bgcFile',bgcFile,PARM_STR);
  PARAMS = addParameter(PARAMS,'nbgc',nbgc,PARM_INT);
  
  
 
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Tracer initial conditions %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% To store all tracers
  phi_init = zeros(Ntracs,Nx,Nz);
  
  %%% Initial velocities
  uvel_init = zeros(Nx,Nz);
  vvel_init = zeros(Nx,Nz);
  
  %%% Initial buoyancy
  Hexp = 500;
  Tmax = 20 - 5*XX_tr/Lx;
  Tmin = 0;
%   buoy_init = Tmin + (Tmax-Tmin).*(exp(ZZ_tr/Hexp+1)-exp(-H/Hexp+1))./(exp(1)-exp(-H/Hexp+1));
  buoy_init = 20*(1-tanh(-ZZ_tr./Hexp));
  
  
  %%% Initial depth tracer
  dtr_init = ZZ_tr;
  
  %%% Store physical tracers in 3D matrix
  phi_init(IDX_UVEL,:,:) = reshape(uvel_init,[1 Nx Nz]);
  phi_init(IDX_VVEL,:,:) = reshape(vvel_init,[1 Nx Nz]);
  phi_init(IDX_BUOY,:,:) = reshape(buoy_init,[1 Nx Nz]);  
  
  %%% Count number of bgc tracers
  %%% TODO NEEDS TO BE UPDATED
  switch (modeltype)
      case BGC_NITRATEONLY
          phi_init(IDX_NITRATE,:,:) = reshape(bgc_init,[1 Nx Nz]);
          phi_init(Nphys+Nbgc,:,:) = reshape(bgc_init,[1,Nx,Nz]);
      case BGC_NPZD
          bgc_tracs = MP + MZ + MD + 1;
          for ii = 1:bgc_tracs
            phi_init(IDX_BUOY + ii,:,:) = reshape(bgc_init(:,:,ii),[1 Nx Nz]);
          end
          % dye
          phi_init(Nphys+Nbgc,:,:) = reshape(bgc_init(:,:,1),[1,Nx,Nz]);
      case BGC_SSEM
          bgc_tracs = MP + MZ + MD + 1;
          for ii = 1:bgc_tracs
              phi_init(IDX_BUOY+ii,:,:) = reshape(bgc_init(:,:,ii),[1 Nx Nz]); 
          end
  end

  %%% Additional tracers
  phi_init(Nphys+Nbgc+1,:,:) = reshape(dtr_init,[1 Nx Nz]);
  Nphys
  Nbgc
  Ntracs
  %%% Write to data file
  initFile = 'initFile.dat';  
  writeDataFile(fullfile(local_run_dir,initFile),phi_init);
  PARAMS = addParameter(PARAMS,'initFile',initFile,PARM_STR);  
    
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% Topography %%%%%
  %%%%%%%%%%%%%%%%%%%%%%
 
  topogFile = 'topog.dat';  
  writeDataFile(fullfile(local_run_dir,topogFile),hb);
  PARAMS = addParameter(PARAMS,'topogFile',topogFile,PARM_STR);  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Surface wind stress %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
  %%% Load in the surface wind stress.
%   [tau,tlength] = sfc_wind_stress(tau0,Lx,xx_psi);
%   tau = tau0*tanh(((Lx)-xx_psi)/(Lx/4));
%   Lmax = 300*m1km;
%   tau = tau0*cos(pi/2*(xx_psi-Lmax)/Lmax).^2;
%   tau(xx_psi < 200) = 0;
  tau = tau0*ones(size(xx_psi));
  tau(xx_psi < 200) = 0;
  tau(xx_psi > Lx - 200) = 0;
  tlength = length(tau);
  
  tauFile = 'tau.dat';  
  writeDataFile(fullfile(local_run_dir,tauFile),tau);
  PARAMS = addParameter(PARAMS,'tlength',tlength,PARM_INT);
  PARAMS = addParameter(PARAMS,'tauFile',tauFile,PARM_STR); 

%  
  
 
  
  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Tracer relaxation concentrations and timescales %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Buoyancy relaxation parameters
  L_relax = 50*m1km;  
  T_relax_max = 30*t1day; %%% Fastest relaxation time

  %%% Relax to initial buoyancy at the western boundary
  uvel_relax = zeros(size(uvel_init));
  vvel_relax = zeros(size(vvel_init));
  buoy_relax = buoy_init;
  T_relax_buoy = -ones(Nx,Nz);
  T_relax_buoy(XX_tr<L_relax) = 1 ./ (1/T_relax_max * (1 - XX_tr(XX_tr<L_relax) / L_relax));
  T_relax_buoy(XX_tr>=L_relax) = -1;
  
  T_relax_veloc = -ones(Nx,Nz);
  T_relax_veloc(XX_tr<L_relax) = 1./T_relax_max;
  T_relax_veloc(XX_tr>=L_relax) = -1;
  
  %%% Add relaxation to an atmospheric temperature profile
  buoy_surf_max = 20;
  buoy_surf_min = 15;
  buoy_surf = buoy_surf_max + (buoy_surf_min-buoy_surf_max)*xx_tr/Lx;
  buoy_relax((xx_tr>=L_relax),Nz) = buoy_surf((xx_tr>=L_relax)); 
  T_relax_buoy((xx_tr>=L_relax),Nz) = 1*t1day; 
  
  %%% Depth tracer relaxation  
  dtr_relax = dtr_init;
  T_relax_dtr = 5*t1year * ones(Nx,Nz);
  
  %%% Relax nitrate to initial conditions
  bgc_relax = bgc_init;
 
  %%% Store tracer relaxation data in 3D matrices
  phi_relax_all = zeros(Ntracs,Nx,Nz);
  phi_relax_all(IDX_UVEL,:,:) = reshape(uvel_relax,[1 Nx Nz]);
  phi_relax_all(IDX_VVEL,:,:) = reshape(vvel_relax,[1 Nx Nz]);
  phi_relax_all(IDX_BUOY,:,:) = reshape(buoy_relax,[1 Nx Nz]);  
  %%% TODO NEEDS TO BE UPDATED
  switch (modeltype)
      case BGC_NITRATEONLY
          phi_relax_all(IDX_NITRATE,:,:) = reshape(bgc_relax,[1 Nx Nz]);
          phi_relax_all(Nphys+Nbgc,:,:) = reshape(bgc_relax,[1 Nx Nz]);
      case BGC_NPZD
          phi_relax_all(IDX_NITRATE,:,:) = reshape(bgc_relax(:,:,1),[1 Nx Nz]);
          phi_relax_all(Nphys+Nbgc,:,:) = reshape(bgc_relax(:,:,1),[1 Nx Nz]);
      case BGC_SSEM
          phi_relax_all(3:end,:,:) = reshape(bgc_relax,[Nbgc Nx Nz]);
  end
  phi_relax_all(Nphys+Nbgc+1,:,:) = reshape(dtr_relax,[1 Nx Nz]);
  
  %%% Store tracer relaxation timescales
  T_relax_all = zeros(Ntracs,Nx,Nz);
  T_relax_all(IDX_UVEL,:,:) = reshape(T_relax_veloc,[1 Nx Nz]);
  T_relax_all(IDX_VVEL,:,:) = reshape(T_relax_veloc,[1 Nx Nz]);
  
  
  T_relax_all(IDX_BUOY,:,:) = reshape(T_relax_buoy,[1 Nx Nz]);  
  %%% TODO NEEDS TO BE UPDATED
  switch (modeltype)
      case BGC_NITRATEONLY
          T_relax_all(IDX_NITRATE,:,:) = 100*t1day*ones(1,Nx,Nz); % Nitrate restored at 100 days conserved
          T_relax_all(Nphys+Nbgc,:,:) = -ones(1,Nx,Nz); % Total dye conserved
      case BGC_NPZD
          T_relax_all(IDX_NITRATE,:,:) = 100*t1day*ones(1,Nx,Nz); % Nitrate restored at 100 days conserved
          T_relax_all(IDX_NITRATE+1:IDX_NITRATE+3,:,:) = -ones(3,Nx,Nz); % Conserve all other parts of npzd model
          T_relax_all(Nphys+Nbgc,:,:) = -ones(1,Nx,Nz); % Total dye conserved
      case BGC_SSEM
          T_relax_all(IDX_NITRATE:end,:,:) = -ones(Nbgc,Nx,Nz);
  end
  T_relax_all(Nphys+Nbgc+1,:,:) = reshape(T_relax_dtr,[1 Nx Nz]);

  relaxTracerFile = 'relaxTracer.dat';
  relaxTimeFile = 'relaxTime.dat';
  writeDataFile(fullfile(local_run_dir,relaxTracerFile),phi_relax_all);
  writeDataFile(fullfile(local_run_dir,relaxTimeFile),T_relax_all);
  PARAMS = addParameter(PARAMS,'relaxTracerFile',relaxTracerFile,PARM_STR);     
  PARAMS = addParameter(PARAMS,'relaxTimeFile',relaxTimeFile,PARM_STR);
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Buoyancy diffusivity %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Uniform diffusivity
  Kgm = Kgm0*ones(Nx+1,Nz+1);             
  KgmFile = 'Kgm.dat';
  writeDataFile(fullfile(local_run_dir,KgmFile),Kgm);
  PARAMS = addParameter(PARAMS,'KgmFile',KgmFile,PARM_STR);
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Isopycnal diffusivity %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Uniform diffusivity
%   Kiso = Kiso0*ones(Nx+1,Nz+1);        
  %%% First guess is a linearly decreasing profile with depth from Kiso0 to
  %%% Kiso_int, with respect to depth. 
%   Kiso = (((Kiso0 - Kiso_hb)/H).*ZZ_psi) + Kiso0;

  
  %%% Another guess is a hyperbolic profile decreasing to 200 at the lower
  %%% boundary. 
  Kefold = 1000;
  Kiso = Kiso0 + (Kiso0-Kiso_hb)*tanh(ZZ_psi./Kefold);
  
%   Kiso = Kgm;
  KisoFile = 'Kiso.dat';
  writeDataFile(fullfile(local_run_dir,KisoFile),Kiso);
  PARAMS = addParameter(PARAMS,'KisoFile',KisoFile,PARM_STR);
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Diapycnal diffusivity %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

  %%% Uniform diffusivity
  Kdia = Kdia0*ones(Nx+1,Nz+1);  
  Ksml = 1e-1;
  Kbbl = 1e-1;
  HB_psi = repmat(reshape(hb_psi,[Nx+1 1]),[1 Nz+1]);
  
  %%% Check if sml and bbl overlap and add the profiles and create crude
  %%% mixed layers.
  kvec = Kdia0*ones(Nz+1,1);
  for ii = 1:Nx+1
      if ZZ_psi(ii,1) > -(Hsml + Hbbl) % Overlapping boundary layers
          H = -ZZ_psi(ii,1);
          for jj = 1:Nz+1
              Kdia(ii,jj) = Kdia0 + 2*(Ksml * -4*(ZZ_psi(ii,jj)/H).*(ZZ_psi(ii,jj)/H+1));
          end
      else
          for jj = 1:Nz+1
              if (ZZ_psi(ii,jj) > -Hsml) % Builds profile when sml and bbl don't overlap
                  Kdia(ii,jj) = Kdia(ii,jj) + (Ksml * -4*(ZZ_psi(ii,jj)/Hsml).*(ZZ_psi(ii,jj)/Hsml+1));
              elseif (ZZ_psi(ii,jj) < -hb_psi(ii)+Hbbl)
                  Kdia(ii,jj) = Kdia(ii,jj) + Kbbl * -4*((ZZ_psi(ii,jj)+HB_psi(ii,jj))/Hbbl).*((ZZ_psi(ii,jj)+HB_psi(ii,jj))/Hbbl-1);
              end
          end
      end
  end
  
  
  %%% Write to file
  KdiaFile = 'Kdia.dat';
  writeDataFile(fullfile(local_run_dir,KdiaFile),Kdia);
  PARAMS = addParameter(PARAMS,'KdiaFile',KdiaFile,PARM_STR); 
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Create scripts %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Create a run script
  createRunScript (local_home_dir,run_name,exec_name, ...
                   use_cluster,cluster_username,cluster_address, ...
                   cluster_home_dir,walltime)

  %%% Create the input parameter file
  writeParamFile(pfname,PARAMS);    
  
  
  
  
  
  
  
  
  
  
  
  %%% 
  %%% The following is for visualization purposes and can be commented out.
  %%% Plot some figures to show some initial values
  %%%
  
  % Wind Stress Profile
  
  set(0,'DefaultAxesFontSize',14)
  
  figure(fignum);
  fignum = fignum+1;
  plot(xx_psi,tau(1,:),'k','LineWidth',2)
  title('Surface Wind Stress')
  axis tight
  xticks([])
  yticks([-0.09 -0.05 0])
  
  % Wind Stress Curl
  figure(fignum)
  fignum = fignum + 1;
  wsc = (tau(1,2:end) - tau(1,1:end-1))/(dx);  
  plot(xx_psi(1:end-1),wsc,'k','LineWidth',2)
  title('Windstress Curl')
  
  %%% Plot Buoyancy relaxation profile

  figure(fignum)
  fignum = fignum + 1;
  pcolor(XX_tr(1:10,:),ZZ_tr(1:10,:),buoy_relax(1:10,:))
%   title('Buoyancy Restoring Profile')
  colorbar('westoutside')
  shading interp
  yticks([])
  xticks([])
  
  %%% Plot diapynal and isopycnal diffusivities together.
  % Diapycnal Diffusivities
  figure(fignum)
  subplot(1,2,1)
  pcolor(XX_psi,ZZ_psi,Kdia)
  title('Diapycnal diffusivity, \kappa_{\nu}')
  shading interp
  colorbar
%   caxis([-7 -2])
%   c = [0 1e-4 1e-3 1e-2];
%   ctick = [-7 -4 -3 -2];
%   colorbar('FontSize',14,'Ytick', ctick, 'YTickLabel',c);
%   caxis([0,1e-3])
  ylabel('Depth (m)')
  xlabel('Distance (km)')
  
  % Isopycnal diffusivity
  figure(fignum);
  subplot(1,2,2)
  pcolor(XX_psi,ZZ_psi,Kiso)
  title('Isopycnal Diffusivity, \kappa_{iso}')
  shading interp
  colorbar
  xlabel('Distance (km)')
  fignum = fignum+1;
  
  % Initial buoyancy
  figure(fignum);
  fignum = fignum+1;
  pcolor(XX_tr,ZZ_tr,buoy_init);
  title('Initial Buoyancy with Grid')
  colorbar
  
  % Initial buoyancy
  figure(fignum);
  fignum = fignum+1;
  contourf(XX_tr,ZZ_tr,buoy_init,10);
  title('Initial Buoyancy')
  colorbar
  
  figure(fignum)
  fignum = fignum+1;
  contourf(XX_tr,ZZ_tr,dtr_init,-(0:200:H));
  colorbar
  title('Depth Tracer, t = 0')
  xlabel('Distance (km)')
  ylabel('Depth (m)')
  
%   % Slope
%   figure(fignum)
%   fignum = fignum+1;
%   plot(xx_tr,hb_slope)
%   title('Topographic Slope')
  
end
