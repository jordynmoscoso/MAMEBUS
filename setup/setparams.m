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
  addpath ../utils;
  
  %%% Load globally-defined constants
  paramTypes;
  
  % User defined choices:
  modeltype = BGC_NPZD; % BGC_NONE (physical model only) or BGC_NPZD (NPZD model)
  if (modeltype > 1)
    modeltype = BGC_NONE; %%% This automatically defaults so that the model runs without biogeochemistry
  end
  
  %%% Plot the setup figures
  plotfigs = true;
  
  %%% If set true, set up this run for the cluster
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
  endTime = 30*t1year;
  restart = false;
  startIdx = 15;
  outputFreq = 1*t1day;
    
  %%% Domain dimensions
  m1km = 1000; %%% Meters in 1 km    
  Lx = 400*m1km; %%% Computational domain width
  H = 3000; %%% Depth of domain including mixed layers
  
  %%% Scalar parameter definitions 
  tau0 = -0.025e-1; %%% Northward wind stress (N m^{-2})
  shelfdepth = 50; %%% Depth of shelf on western boundary
  rho0 = 1025; %%% Reference density
  f0 = 1e-4; %%% Coriolis parameter (CCS)
  Kgm0 = 1200; %%% Reference GM diffusivity
  Kdia0 = 1e-5; %%% Reference diapycnal diffusivity  
  Hsml = 40; %%% Surface mixed layer thickness
  Hbbl = 40; %%% Bottom boundary layer thickness
  r_bbl = 1e-3; %%% Bottom boundary layer drag coefficient

  %%% Grid parameters
  h_c = 300; %%% Sigma coordinate surface layer thickness parameter (must be > 0)
  theta_s = 9; %%% Sigma coordinate surface stretching parameter (must be in [0,10])
  theta_b = 4; %%% Sigma coordinage bottom stretching parameter (must be in [0,4])
  
  %%% Grids  
  Nx = 32; %%% Number of latitudinal grid points 
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
  Xtopog = 290*m1km;
  Ltopog = 30*m1km;
  shelfdepth = 50;
  Htopog = H-shelfdepth;
  hb = H - Htopog*0.5*(1+tanh((xx_topog-Xtopog)/(Ltopog)));
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
  
  
  
  %%%%%%%%%%%%%%%%%%%%%
  %%%%% BGC Model %%%%%
  %%%%%%%%%%%%%%%%%%%%%
  
  if modeltype
      MP = 1; MZ = 1; MD = 1; MN = 1;
      Nbgc = MP + MZ + MD + MN;
  else
      MP = 0; MZ = 0; MD = 0;
      Nbgc = 0;
  end
  
  PARAMS = addParameter(PARAMS,'bgcModel',modeltype,PARM_INT);
  PARAMS = addParameter(PARAMS,'pressureScheme',PRESSURE_CUBIC,PARM_INT);
  PARAMS = addParameter(PARAMS,'MP',MP,PARM_INT);
  PARAMS = addParameter(PARAMS,'MZ',MZ,PARM_INT);
  
  % Calculate the number of tracers
  Nphys = 3; %%% Number of physical tracers (u-velocity, v-velocity and buoyancy)
  Ntracs = Nphys + Nbgc; %%% Number of tracers (physical plus bgc plus any other user-defined tracers)
  
  % Get the BGC parameters and initial conditions
  [bgc_params, bgc_init, nbgc] = bgc_setup(ZZ_tr,Nx,Nz);
  
  % Write the biogeochemistry to file
  bgcFile = 'bgcFile.dat';
  writeDataFile(fullfile(local_run_dir,bgcFile),bgc_params);
  PARAMS = addParameter(PARAMS,'bgcFile',bgcFile,PARM_STR);
  PARAMS = addParameter(PARAMS,'Ntracs',Ntracs,PARM_INT);
  PARAMS = addParameter(PARAMS,'nbgc',nbgc,PARM_INT);
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Tracer initial conditions %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% To store all tracers
  phi_init = zeros(Ntracs,Nx,Nz);
  
  %%% Initial velocities
  uvel_init = zeros(Nx,Nz);
  vvel_init = zeros(Nx,Nz);
  
  %%% Initial temperature
  Hexp = 150; 
  Tmax = 22 - 5*XX_tr/Lx;
  Tmin = 4;
  buoy_init = Tmin + (Tmax-Tmin).*(exp(ZZ_tr/Hexp+1)-exp(-H/Hexp+1))./(exp(1)-exp(-H/Hexp+1));  
  
  
  %%% Store physical tracers in 3D matrix
  phi_init(IDX_UVEL,:,:) = reshape(uvel_init,[1 Nx Nz]);
  phi_init(IDX_VVEL,:,:) = reshape(vvel_init,[1 Nx Nz]);
  phi_init(IDX_BUOY,:,:) = reshape(buoy_init,[1 Nx Nz]);  
  
  if (modeltype == 1)
      phi_init(Nphys+1:Nphys+Nbgc,:,:) = reshape(bgc_init,[Nbgc Nx Nz]);
  end
  

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
 
  %%% Create the profile of wind stress
  tau = tau0*( tanh(((Lx)-xx_psi)/(Lx/4)) );
  tlength = length(tau);
  
  tauFile = 'tau.dat';  
  writeDataFile(fullfile(local_run_dir,tauFile),tau);
  PARAMS = addParameter(PARAMS,'tlength',tlength,PARM_INT);
  PARAMS = addParameter(PARAMS,'tauFile',tauFile,PARM_STR); 
  

  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Tracer relaxation concentrations and timescales %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Buoyancy relaxation parameters
  L_relax = 50*m1km;  
  T_relax_max = 30*t1day; %%% Fastest relaxation time

  %%% Relax to initial buoyancy at the western boundary
  uvel_relax = -ones(size(uvel_init));
  vvel_relax = -ones(size(vvel_init));
  T_relax_veloc = -ones(Nx,Nz);
  
  buoy_relax = buoy_init;
  T_relax_buoy = -ones(Nx,Nz);
  T_relax_buoy(XX_tr<L_relax) = 1 ./ (1/T_relax_max * (1 - XX_tr(XX_tr<L_relax) / L_relax));
  T_relax_buoy(XX_tr>=L_relax) = -1;

  
  %%% Add relaxation to an atmospheric temperature profile
  buoy_surf_max = 22;
  buoy_surf_min = 18;
  buoy_surf = buoy_surf_max + (buoy_surf_min-buoy_surf_max)*xx_tr/Lx;
  buoy_relax((xx_tr>=L_relax),Nz) = buoy_surf((xx_tr>=L_relax)); 
  T_relax_buoy((xx_tr>=L_relax),Nz) = 1*t1day; 
  
 
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
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Buoyancy diffusivity %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Initalize the profile
  Kgm = Kgm0*ones(Nx+1,Nz+1);
  H_int = hb_psi - (Hbbl + Hsml); 
  H_int(H_int < 0) = 0;
  Hmax = max(H_int);
  
  lambda = 0.25; % tuning parameter for KGM diffusivity
  for jj = 1:Nx+1
      for kk = 1:Nz+1 
        Kgm(jj,kk) = Kgm(jj,kk)*H_int(jj)*exp(ZZ_psi(jj,kk)/(lambda*Hmax))/Hmax;
      end
  end
  
  KgmFile = 'Kgm.dat';
  writeDataFile(fullfile(local_run_dir,KgmFile),Kgm);
  PARAMS = addParameter(PARAMS,'KgmFile',KgmFile,PARM_STR);
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Isopycnal diffusivity %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Kiso = 2*Kgm;
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
  
  % Isopycnal diffusivity
  figure(fignum);
  subplot(1,2,1)
  pcolor(XX_psi,ZZ_psi,Kiso)
  title('Isopycnal Diffusivity, \kappa_{iso}')
  shading interp
  colorbar
  xlabel('Distance (km)')
  
  % Buoyancy diffusivity
  figure(fignum);
  subplot(1,2,2)
  pcolor(XX_psi,ZZ_psi,Kgm)
  title('Gent/McWilliams Diffusivity, \kappa_{gm}')
  shading interp
  colorbar
  xlabel('Distance (km)')
  fignum = fignum+1;
  
  % Initial buoyancy
  figure(fignum);
  [C, h] = contourf(XX_tr,ZZ_tr,buoy_init,0:2:20);
  clabel(C,h)
  title('Initial Buoyancy')
  colorbar

  end
  
end
