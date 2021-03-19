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
  %%% Convenience scripts used in this function
  addpath ../utils; % cmocean.m should be downloaded here
  model_code_dir = '../code';
 
  %%% Load globally-defined constants
  paramTypes;
  
  % User defined choices:
  modeltype = BGC_NONE; % BGC_NONE (physical model only) or BGC_NPZD (NPZD model)
  if (modeltype > 1)
    modeltype = BGC_NONE; %%% This automatically defaults so that the model runs without biogeochemistry
  end
  
  %%% Plot the setup figures
  plotfigs = true;
  
  %%% If set true, set up this run for the cluster
  %%% This may need tuning for individual systems
  use_cluster = false;
  use_intel = false;
  use_pbs = use_cluster;
  cluster_home_dir = '/data3/jmoscoso/MAMEBUS/runs';
  uname = 'jmoscoso';
  cluster_addr = 'caolila.atmos.ucla.edu';
  
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
  endTime = 20*t1year;
  restart = false;
  startIdx = 0;
  outputFreq = 10*t1day;
    
  %%% Domain dimensions
  m1km = 1000; %%% Meters in 1 km    
  Lx = 1000*m1km; %%% Computational domain width
  H = 2500; %%% Depth of domain including mixed layers
  
  %%% Scalar parameter definitions 
  tau0 = 0.05;%-0.025e-1; %%% Northward wind stress (N m^{-2})
  shelfdepth = 300; %%% Depth of shelf on western boundary
  rho0 = 1025; %%% Reference density
  f0 = 1e-4; %%% Coriolis parameter (CCS)
  Kgm0 = 1000; %%% Reference GM diffusivity  
  Hsml = 100; %%% Surface mixed layer thickness
  Hbbl = 40; %%% Bottom boundary layer thickness
  r_bbl = 1e-3; %%% Bottom boundary layer drag coefficient

  %%% Grid parameters
  h_c = 250; %%% Sigma coordinate surface layer thickness parameter (must be > 0)
  theta_s = 0; %%% Sigma coordinate surface stretching parameter (must be in [0,10])
  theta_b = 2; %%% Sigma coordinage bottom stretching parameter (must be in [0,4])
  
  %%% Grids  
  Nx = 96; %%% Number of latitudinal grid points 
  Nz = 32; %%% Number of vertical grid points
  dx = Lx/Nx; %%% Latitudinal grid spacing (in meters)
  xx_psi = 0:dx:Lx; %%% Streamfunction latitudinal grid point locations
  xx_tr = dx/2:dx:Lx-dx/2; %%% Tracer latitudinal grid point locations  
  xx_topog = [-dx/2 xx_tr Lx+dx/2]; %%% Topography needs "ghost" points to define bottom slope

  %%% For plotting figures of setup
  fignum = 1;

  
  disp(['Shelf Depth: ', num2str(shelfdepth)])
  if shelfdepth < Hsml + Hbbl
      disp('Full depth mixed layer on the shelf')
  end

  %%% Define the topography
  Xtopog_w = 150*m1km;
  Xtopog_e = 850*m1km;
  Ltopog = 25*m1km;
  Htopog = H-shelfdepth;
  hb = H - Htopog*0.5*(1+tanh((xx_topog-Xtopog_e)/(Ltopog))) - Htopog*0.5*(1-tanh((xx_topog-Xtopog_w)/(Ltopog)));
  hb_psi = 0.5*(hb(1:end-1)+hb(2:end));  
  hb_tr = hb(2:end-1);
  
  % calculate shelf slope:
  hb_slope = (hb_psi(2:end) - hb_psi(1:end-1))./(dx);
  
  %%% Generate full sigma-coordinate grids
  [XX_tr,ZZ_tr,XX_psi,ZZ_psi,XX_u,ZZ_u,XX_w,ZZ_w] ...
                    = genGrids(Nx,Nz,Lx,h_c,theta_s,theta_b,hb_tr,hb_psi); 
  slopeidx = max((hb_psi>Htopog/2));
  disp(['slopeidx = ',num2str(slopeidx)])
  disp(['Vertical grid spacing at (',num2str(XX_psi(1,1)),',',num2str(ZZ_psi(1,1)),'): ',num2str(ZZ_psi(1,2)-ZZ_psi(1,1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(1,end)),',',num2str(ZZ_psi(1,end)),'): ',num2str(ZZ_psi(1,end)-ZZ_psi(1,end-1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(end,1)),',',num2str(ZZ_psi(end,1)),'): ',num2str(ZZ_psi(end,2)-ZZ_psi(end,1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(end,end)),',',num2str(ZZ_psi(end,end)),'): ',num2str(ZZ_psi(end,end)-ZZ_psi(end,end-1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(slopeidx,1)),',',num2str(ZZ_psi(slopeidx,1)),'): ',num2str(ZZ_psi(slopeidx,2)-ZZ_psi(slopeidx,1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(slopeidx,end)),',',num2str(ZZ_psi(slopeidx,end)),'): ',num2str(ZZ_psi(slopeidx,end)-ZZ_psi(slopeidx,end-1))])
  disp(['Max Topographic Slope : ' num2str(min(hb_slope))])
  
  
  %%% Calculate grid stiffness  
  rx1 = abs(diff(0.5*(ZZ_psi(:,1:Nz)+ZZ_psi(:,2:Nz+1)),1,1) ./ diff(0.5*(ZZ_psi(1:Nx,:)+ZZ_psi(2:Nx+1,:)),1,2) );
  disp(['Grid stiffness: ' num2str(max(max(rx1)))]  )
  
  
  %%% Define parameter
  PARAMS = addParameter(PARAMS,'Nx',Nx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nz',Nz,PARM_INT);  
  PARAMS = addParameter(PARAMS,'Lx',Lx,PARM_REALF);
  PARAMS = addParameter(PARAMS,'Lz',H,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'cflFrac',0.25,PARM_REALF);
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
  PARAMS = addParameter(PARAMS,'pressureScheme',PRESSURE_CUBIC,PARM_INT);
  
   
  % Calculate the number of tracers
  Nphys = 3; %%% Number of physical tracers (u-velocity, v-velocity and buoyancy)
  Ntracs = Nphys;
  PARAMS = addParameter(PARAMS,'nbgc',0,PARM_INT);  
  PARAMS = addParameter(PARAMS,'MP',0,PARM_INT);
  PARAMS = addParameter(PARAMS,'MZ',0,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ntracs',Ntracs,PARM_INT);
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Tracer initial conditions %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% To store all tracers
  phi_init = zeros(Ntracs,Nx,Nz);
  
  %%% Initial velocities
  uvel_init = zeros(Nx,Nz);
  vvel_init = zeros(Nx,Nz);
  
  %%% Initial temperature
  Hexp = 500; 
  Tmax = 10;
  Tmin = 0;
  buoy_init = Tmin + (Tmin + (Tmax-Tmin)*(XX_tr/Lx) - Tmin).*(exp(ZZ_tr/Hexp+1)-exp(-H/Hexp+1))./(exp(1)-exp(-H/Hexp+1));  
  
  
  %%% Store physical tracers in 3D matrix
  phi_init(IDX_UVEL,:,:) = reshape(uvel_init,[1 Nx Nz]);
  phi_init(IDX_VVEL,:,:) = reshape(vvel_init,[1 Nx Nz]);
  phi_init(IDX_BUOY,:,:) = reshape(buoy_init,[1 Nx Nz]);  
   
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
 
  %%% Create the profile of wind stress
%   tau = tau0*( tanh(((Lx)-xx_psi)/(Lx/4)) );
  tau = tau0*ones(size(xx_psi));
  figure(10)
  plot(xx_psi,tau);
  tlength = length(tau);
  
  tauFile = 'tau.dat';  
  writeDataFile(fullfile(local_run_dir,tauFile),tau);
  PARAMS = addParameter(PARAMS,'tlength',tlength,PARM_INT);
  PARAMS = addParameter(PARAMS,'tauFile',tauFile,PARM_STR); 
  

  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Tracer relaxation concentrations and timescales %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  %%% Relax to initial buoyancy at the western boundary
  uvel_relax = -ones(size(uvel_init));
  vvel_relax = -ones(size(vvel_init));
  T_relax_veloc = -ones(Nx,Nz);
  

  
  %%% Add relaxation to an atmospheric temperature profile
  buoy_surf_max = Tmax;
  buoy_surf_min = Tmin;
  buoy_surf = buoy_surf_min + (buoy_surf_max-buoy_surf_min)*xx_tr/Lx;
  buoy_relax = repmat(reshape(buoy_surf,[Nx 1]),[1 Nz]); 
  figure(100);   
  plot(xx_tr,buoy_relax(:,Nz));
  T_relax_buoy = -ones(Nx,Nz);
  T_relax_buoy(:,Nz) = 30*t1day; 
  
 
  %%% Store tracer relaxation data in 3D matrices
  phi_relax_all = zeros(Ntracs,Nx,Nz); 
  phi_relax_all(IDX_UVEL,:,:) = reshape(uvel_relax,[1 Nx Nz]);
  phi_relax_all(IDX_VVEL,:,:) = reshape(vvel_relax,[1 Nx Nz]);
  phi_relax_all(IDX_BUOY,:,:) = reshape(buoy_relax,[1 Nx Nz]); 
  
  
  %%% Store tracer relaxation timescales
  T_relax_all = -ones(Ntracs,Nx,Nz);
  T_relax_all(IDX_UVEL,:,:) = reshape(T_relax_veloc,[1 Nx Nz]);
  T_relax_all(IDX_VVEL,:,:) = reshape(T_relax_veloc,[1 Nx Nz]);
  T_relax_all(IDX_BUOY,:,:) = reshape(T_relax_buoy,[1 Nx Nz]);  
  
  %%% Store the files
  relaxTracerFile = 'relaxTracer.dat';
  relaxTimeFile = 'relaxTime.dat';
  writeDataFile(fullfile(local_run_dir,relaxTracerFile),phi_relax_all);
  writeDataFile(fullfile(local_run_dir,relaxTimeFile),T_relax_all);
  PARAMS = addParameter(PARAMS,'relaxTracerFile',relaxTracerFile,PARM_STR);     
  PARAMS = addParameter(PARAMS,'relaxTimeFile',relaxTimeFile,PARM_STR);
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Meridional Profiles of Tracers %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Northern Profiles
  phi_north = zeros(Ntracs,Nx,Nz);
  
  buoy_north = buoy_init;
%   phi_north(IDX_BUOY,:,:) = buoy_north - ones(Nx,Nz);
  phi_north(IDX_BUOY,:,:) = zeros(Nx,Nz);
  
  % Store the files
  northTracerFile = 'northTracer.dat';
  writeDataFile(fullfile(local_run_dir,northTracerFile),phi_north);
  PARAMS = addParameter(PARAMS,'northTracerFile',northTracerFile,PARM_STR);
  
  %%% Southern Profiles
  phi_south = zeros(Ntracs,Nx,Nz);
  
  buoy_south = buoy_init;
%   phi_south(IDX_BUOY,:,:) = buoy_south + ones(Nx,Nz);
  phi_south(IDX_BUOY,:,:) = ones(Nx,Nz);
  
  % Store the files
  southTracerFile = 'southTracer.dat';
  writeDataFile(fullfile(local_run_dir,southTracerFile),phi_south);
  PARAMS = addParameter(PARAMS,'southTracerFile',southTracerFile,PARM_STR);
  
  %%% Length of domain
  Ly = 1000*m1km; % negative values indicate an infinite domain
  PARAMS = addParameter(PARAMS,'Ly',Ly,PARM_REALF);
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Profiles of Surface Fluxes %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % surface fluxes must have units of (units of [c] m/s)
  phi_flux = zeros(Ntracs,Nx);
  
  fluxFile = 'fluxFile.dat';
  writeDataFile(fullfile(local_run_dir,fluxFile),phi_flux);
  PARAMS = addParameter(PARAMS,'fluxFile',fluxFile,PARM_STR);
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Buoyancy diffusivity %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Initalize the profile
  Kgm = Kgm0*ones(Nx+1,Nz+1);
%   H_int = hb_psi - (Hbbl + Hsml); 
%   H_int(H_int < 0) = 0;
%   Hmax = max(H_int);
%   
%   lambda = 0.25; % tuning parameter for KGM diffusivity
%   for jj = 1:Nx+1
%       for kk = 1:Nz+1 
%         Kgm(jj,kk) = Kgm(jj,kk)*H_int(jj)*exp(ZZ_psi(jj,kk)/(lambda*Hmax))/Hmax;
% %         Kgm(jj,kk) = Kgm(jj,kk)*exp(ZZ_psi(jj,kk)/(lambda*Hmax));
%       end
%   end
  
  KgmFile = 'Kgm.dat';
  writeDataFile(fullfile(local_run_dir,KgmFile),Kgm);
  PARAMS = addParameter(PARAMS,'KgmFile',KgmFile,PARM_STR);
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Isopycnal diffusivity %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Kiso = 0*Kgm;
  KisoFile = 'Kiso.dat';
  writeDataFile(fullfile(local_run_dir,KisoFile),Kiso);
  PARAMS = addParameter(PARAMS,'KisoFile',KisoFile,PARM_STR);
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Diapycnal diffusivity %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

  Kdia0 = 1e-5; % interior diffusivity
  Ksml = 1e-1;  % SML diffusivity
  Kbbl = 1e-1;  % BBL diffusivity

  %%% Write to file
  PARAMS = addParameter(PARAMS,'Kdia0',Kdia0,PARM_REALF);
  PARAMS = addParameter(PARAMS,'Ksml',Ksml,PARM_REALF);
  PARAMS = addParameter(PARAMS,'Kbbl',Kbbl,PARM_REALF);
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Create scripts %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Create a run script
  createRunScript (  local_home_dir, run_name, model_code_dir, ...
                     exec_name, use_intel, use_pbs, use_cluster, ...
                     uname, cluster_addr, cluster_home_dir)

  %%% Create the input parameter file
  writeParamFile(pfname,PARAMS);    
  
  

  
  
  
  %%% 
  %%% The following is for visualization purposes and can be commented out.
  %%% Plot some figures to show some initial values
  %%%
  
  % Wind Stress Profile
  

  
  if(plotfigs)
  
  set(0,'DefaultAxesFontSize',14)
  
  %%% Plot diapynal and isopycnal diffusivities together.
  difmap = cmocean('deep');
  
  % Isopycnal diffusivity
  figure(fignum);
  subplot(1,2,1)
  pcolor(XX_psi,ZZ_psi,Kiso)
  title('Isopycnal Diffusivity, \kappa_{iso}')
  colormap(difmap)
  shading interp
  colorbar
  xlabel('Distance (km)')
  
  % Buoyancy diffusivity
  figure(fignum);
  subplot(1,2,2)
  pcolor(XX_psi,ZZ_psi,Kgm)
  title('Gent/McWilliams Diffusivity, \kappa_{gm}')
  colormap(difmap)
  shading interp
  colorbar
  xlabel('Distance (km)')
  fignum = fignum+1;
  
  % Initial buoyancy
  bmap = cmocean('thermal');
  figure(fignum);
  [C, h] = contourf(XX_tr,ZZ_tr,buoy_init,0:1:10);
  clabel(C,h)
  colormap(bmap)
  title('Initial Buoyancy')
  colorbar
  end
  
  license('inuse')
  
end
