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
function setparams (run_name)    

  %%% Convenience scripts used in this function
  addpath ../utils;

  %%% If set true, set up this run for the cluster
  use_cluster = false;
  walltime = 24;
  
  %%% Run directory
  run_name = strtrim(run_name);  
  exec_name = 'mamebus.exe';
  local_home_dir = '../runs';    
  local_run_dir = fullfile(local_home_dir,run_name);
  pfname = fullfile(local_run_dir,[run_name,'_in']);   
  mkdir(local_run_dir);

  %%% To store parameters
  paramTypes;
  PARAMS = {};
  
  %%% Domain dimensions
  m1km = 1000; %%% Meters in 1 km
  t1day = 86400; %%% Seconds in 1 day
  t1year = 365*t1day; %%% Seconds in 1 year
  H = 4*m1km; %%% Depth, excluding the mixed layer
  Lx = 500*m1km; %%% Computational domain width
  
  %%% Scalar parameter definitions 
  tau0 = -1e-1; %%% Northward wind stress
  rho0 = 1e3; %%% Reference density
  f0 = 1e-4; %%% Coriolis parameter (Southern Ocean)
  Kgm0 = 500; %%% Reference GM diffusivity
  Kiso0 = 500; %%% Reference isopycnal diffusivity
  Kdia0 = 1e-5; %%% Reference diapycnal diffusivity
  Cp = 4e3; %%% Heat capacity
  g = 9.81; %%% Gravity
  s0 = tau0/rho0/f0/Kgm0; %%% Theoretical isopycnal slope  

  %%% Grids  
  Ntracs = 2; %%% Number of tracers
  Nx = 50; %%% Number of latitudinal grid points 
  Nz = 40; %%% Number of vertical grid points
  dx = Lx/Nx; %%% Latitudinal grid spacing (in meters)
  dz = H/Nz; %%% Vertical grid spacing (in meters)
  xx_psi = 0:dx:Lx; %%% Streamfunction latitudinal grid point locations
  zz_psi = -H:dz:0; %%% Streamfunction vertical grid point locations
  xx_tr = dx/2:dx:Lx-dx/2; %%% Tracer latitudinal grid point locations
  zz_tr = -H+dz/2:dz:0-dz/2; %%% Tracer vertical grid point locations
  [ZZ_psi XX_psi] = meshgrid(zz_psi,xx_psi); %%% Streamfunction meshgrid
  [ZZ_tr XX_tr] = meshgrid(zz_tr,xx_tr); %%% Density meshgrid        
  xx_topog = [-dx/2 xx_tr Lx+dx/2]; %%% Topography needs "ghost" points to define bottom slope
  
  %%% Create tanh-shaped topography
  Xtopog = 400*m1km;
  Ltopog = 25*m1km;
  Htopog = H-200;  
  hb = H - Htopog*0.5*(1+tanh((xx_topog-Xtopog)/Ltopog));
  hb_psi = 0.5*(hb(1:end-1)+hb(2:end));  
  hb_phi = hb(2:end-1);  
  
  %%% Modify vertical grids to account for topography   
  for j=1:Nx
    ZZ_tr(j,:) = ZZ_tr(j,:) * hb_phi(j)/H;
  end
  for j=1:Nx+1
    ZZ_psi(j,:) = ZZ_psi(j,:) * hb_psi(j)/H;
  end
  
  %%% Calculate grid stiffness  
  rx1 = abs(diff(0.5*(ZZ_psi(:,1:Nz)+ZZ_psi(:,2:Nz+1)),1,1) ./ diff(0.5*(ZZ_psi(1:Nx,:)+ZZ_psi(2:Nx+1,:)),1,2) ) ;
  ['Grid stiffness: ' num2str(max(max(rx1)))]  
  
  %%% Define parameter
  PARAMS = addParameter(PARAMS,'Ntracs',Ntracs,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nx',Nx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nz',Nz,PARM_INT);
  PARAMS = addParameter(PARAMS,'Lx',Lx,PARM_REALF);
  PARAMS = addParameter(PARAMS,'Lz',H,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'cflFrac',1,PARM_REALF);
  PARAMS = addParameter(PARAMS,'monitorFrequency',1*t1year,PARM_REALF);
  PARAMS = addParameter(PARAMS,'maxTime',50*t1year,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'rho0',rho0,PARM_REALF);
  PARAMS = addParameter(PARAMS,'f0',f0,PARM_REALF);     
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Target residuals %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
  
  targetRes = 1e-16 * ones(Ntracs,1);
  targetResFile = 'targetRes.dat';  
  writeDataFile(fullfile(local_run_dir,targetResFile),targetRes);
  PARAMS = addParameter(PARAMS,'targetResFile',targetResFile,PARM_STR); 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Tracer initial conditions %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% To store all tracers
  phi_init = zeros(Ntracs,Nx,Nz);
  
  %%% Initial buoyancy
  Hexp = 1000;
  Tmax = 20;
  Tmin = 0;
  buoy_init = Tmin + (Tmax-Tmin)*(exp(ZZ_tr/Hexp+1)-exp(-H/Hexp+1))./(exp(1)-exp(-H/Hexp+1));
  
  %%% Initial depth tracer
  dtr_init = ZZ_tr;
  
  %%% Store tracers in 3D matrix
  phi_init(1,:,:) = reshape(buoy_init,[1 Nx Nz]);
  phi_init(2,:,:) = reshape(dtr_init,[1 Nx Nz]);
  
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
 
  
  tau = tau0*cos(pi*xx_psi/(2*Lx));
  tauFile = 'tau.dat';  
  writeDataFile(fullfile(local_run_dir,tauFile),tau);
  PARAMS = addParameter(PARAMS,'tauFile',tauFile,PARM_STR);  
    
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
  Kiso = Kiso0*ones(Nx+1,Nz+1);           
  KisoFile = 'Kiso.dat';
  writeDataFile(fullfile(local_run_dir,KisoFile),Kiso);
  PARAMS = addParameter(PARAMS,'KisoFile',KisoFile,PARM_STR);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Diapycnal diffusivity %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

  %%% Uniform diffusivity
  Kdia = Kdia0*ones(Nx+1,Nz+1);
  
%   %%% Crude mixed layers
%   Kdia(ZZ_psi>-300) = 1e0;
  
  %%% Write to file
  KdiaFile = 'Kdia.dat';
  writeDataFile(fullfile(local_run_dir,KdiaFile),Kdia);
  PARAMS = addParameter(PARAMS,'KdiaFile',KdiaFile,PARM_STR); 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Tracer relaxation concentrations and timescales %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Buoyancy relaxation parameters
  L_relax = 50*m1km;  
  T_relax_max = 30*t1day; %%% Fastest relaxation time

  %%% Relax to initial buoyancy at the western boundary
  buoy_relax = buoy_init;
  T_relax_buoy = -ones(Nx,Nz);
  T_relax_buoy(XX_tr<L_relax) = 1 ./ (1/T_relax_max * (1 - XX_tr(XX_tr<L_relax) / L_relax));
  T_relax_buoy(XX_tr>=L_relax) = -1;
  
  %%% Add relaxation to an atmospheric temperature profile
%   buoy_surf_max = 20;
%   buoy_surf_min = 15;
%   buoy_surf = buoy_surf_max + (buoy_surf_min-buoy_surf_max)*xx_tr/Lx;
%   buoy_relax(find(xx_tr>=L_relax),Nz) = buoy_surf(find(xx_tr>=L_relax)); 
%   T_relax_buoy(find(xx_tr>=L_relax),Nz) = 10*t1day; 
  
  %%% Depth tracer relaxation  
  dtr_relax = dtr_init;
  T_relax_dtr = 5*t1year * ones(Nx,Nz);
  
  %%% Store tracer relaxation data in 3D matrices
  phi_relax_all = zeros(Ntracs,Nx,Nz);
  phi_relax_all(1,:,:) = reshape(buoy_relax,[1 Nx Nz]);
  phi_relax_all(2,:,:) = reshape(dtr_relax,[1 Nx Nz]);
  T_relax_all = zeros(Ntracs,Nx,Nz);  
  T_relax_all(1,:,:) = reshape(T_relax_buoy,[1 Nx Nz]);
  T_relax_all(2,:,:) = reshape(T_relax_dtr,[1 Nx Nz]);
  
  relaxTracerFile = 'relaxTracer.dat';
  relaxTimeFile = 'relaxTime.dat';
  writeDataFile(fullfile(local_run_dir,relaxTracerFile),phi_relax_all);
  writeDataFile(fullfile(local_run_dir,relaxTimeFile),T_relax_all);
  PARAMS = addParameter(PARAMS,'relaxTracerFile',relaxTracerFile,PARM_STR);     
  PARAMS = addParameter(PARAMS,'relaxTimeFile',relaxTimeFile,PARM_STR);
    
  
  
  %%% Create a run script
  createRunScript('.',run_name,exec_name,use_cluster,walltime);

  %%% Create the input parameter file
  writeParamFile(pfname,PARAMS);    

  
end
