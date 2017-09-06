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
  paramTypes;
  PARAMS = {};
  
  %%% Domain dimensions
  m1km = 1000; %%% Meters in 1 km
  t1day = 86400; %%% Seconds in 1 day
  t1year = 365*t1day; %%% Seconds in 1 year
  H = 3*m1km; %%% Depth, excluding the mixed layer
  Lx = 300*m1km; %%% Computational domain width
  
  %%% Scalar parameter definitions 
  tau0 = -1e-1; %%% Northward wind stress (N m^{-2})
  rho0 = 1e3; %%% Reference density
  f0 = 1e-4; %%% Coriolis parameter (CCS)
  Kgm0 = 500; %%% Reference GM diffusivity
  Kiso0 = 2000; %%% Reference surface isopycnal diffusivity m^2/s
  Kiso_hb = 200; %%% Reference interior isopycnal diffusivity
  
  Kdia0 = 1e-5; %%% Reference diapycnal diffusivity
  Cp = 4e3; %%% Heat capacity
  g = 9.81; %%% Gravity
  s0 = tau0/rho0/f0/Kgm0; %%% Theoretical isopycnal slope    
  Hsml = 50; %%% Surface mixed layer thickness
  Hbbl = 50; %%% Bottom boundary layer thickness
  
  %%% Biogeochemical Parameters
  a_temp = 0.6/t1day;
  b_temp = 1.066;
  c_temp = 1;
  alpha = 0.025/t1day; % d^-1 (W/m^2)^-1
  monod = 0;
  
  %%% Grid parameters
  h_c = 300; %%% Sigma coordinate surface layer thickness parameter (must be > 0)
  theta_s = 6; %%% Sigma coordinate surface stretching parameter (must be in [0,10])
  theta_b = 4; %%% Sigma coordinage bottom stretching parameter (must be in [0,4])
  
  %%% Grid parameters (no stretching)
%   h_c = 1e16; %%% Sigma coordinate surface layer thickness parameter (must be > 0)
%   theta_s = 0; %%% Sigma coordinate surface stretching parameter (must be in [0,10])
%   theta_b = 0; %%% Sigma coordinage bottom stretching parameter (must be in [0,4])
   
  %%% Grids  
  Ntracs = 3; %%% Number of tracers
  Nx = 40; %%% Number of latitudinal grid points 
  Nz = 40; %%% Number of vertical grid points
  dx = Lx/Nx; %%% Latitudinal grid spacing (in meters)
  xx_psi = 0:dx:Lx; %%% Streamfunction latitudinal grid point locations
  xx_tr = dx/2:dx:Lx-dx/2; %%% Tracer latitudinal grid point locations  
  xx_topog = [-dx/2 xx_tr Lx+dx/2]; %%% Topography needs "ghost" points to define bottom slope
  
  %%% Create tanh-shaped topography
  shelfdepth = 105;
  if shelfdepth < 50
      disp('Shelf is smaller than sml and bbl')
      return
  end
  
  Xtopog = 200*m1km;
  Ltopog = 25*m1km;
  Htopog = H-shelfdepth;  
  hb = H - Htopog*0.5*(1+tanh((xx_topog-Xtopog)/(Ltopog)));
  hb_psi = 0.5*(hb(1:end-1)+hb(2:end));  
  hb_tr = hb(2:end-1);
  
  %%% Generate full sigma-coordinate grids
  [XX_tr,ZZ_tr,XX_psi,ZZ_psi,~,~,~,~] ...
                    = genGrids(Nx,Nz,Lx,h_c,theta_s,theta_b,hb_tr,hb_psi);  % Full output [XX_tr,ZZ_tr,XX_psi,ZZ_psi,XX_u,ZZ_u,XX_w,ZZ_w]
  slopeidx = max((hb_psi>Htopog/2));
  disp(['slopeidx = ',num2str(slopeidx)])
  disp(['Vertical grid spacing at (',num2str(XX_psi(1,1)),',',num2str(ZZ_psi(1,1)),'): ',num2str(ZZ_psi(1,2)-ZZ_psi(1,1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(1,end)),',',num2str(ZZ_psi(1,end)),'): ',num2str(ZZ_psi(1,end)-ZZ_psi(1,end-1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(end,1)),',',num2str(ZZ_psi(end,1)),'): ',num2str(ZZ_psi(end,2)-ZZ_psi(end,1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(end,end)),',',num2str(ZZ_psi(end,end)),'): ',num2str(ZZ_psi(end,end)-ZZ_psi(end,end-1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(slopeidx,1)),',',num2str(ZZ_psi(slopeidx,1)),'): ',num2str(ZZ_psi(slopeidx,2)-ZZ_psi(slopeidx,1))])
  disp(['Vertical grid spacing at (',num2str(XX_psi(slopeidx,end)),',',num2str(ZZ_psi(slopeidx,end)),'): ',num2str(ZZ_psi(slopeidx,end)-ZZ_psi(slopeidx,end-1))])
  
  %%% ZZ_tr size: 40 40
  %%% ZZ_psi size: 41 41
  
  %%% Calculate grid stiffness  
  rx1 = abs(diff(0.5*(ZZ_psi(:,1:Nz)+ZZ_psi(:,2:Nz+1)),1,1) ./ diff(0.5*(ZZ_psi(1:Nx,:)+ZZ_psi(2:Nx+1,:)),1,2) );
  disp(['Grid stiffness: ' num2str(max(max(rx1)))]  )
  
  %%% Define parameter
  PARAMS = addParameter(PARAMS,'Ntracs',Ntracs,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nx',Nx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nz',Nz,PARM_INT);
  PARAMS = addParameter(PARAMS,'H',H,PARM_REALF);
  PARAMS = addParameter(PARAMS,'Lx',Lx,PARM_REALF);
  PARAMS = addParameter(PARAMS,'Lz',H,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'cflFrac',0.5,PARM_REALF);
  PARAMS = addParameter(PARAMS,'maxTime',100*t1year,PARM_REALF);
  PARAMS = addParameter(PARAMS,'monitorFrequency',0.1*t1year,PARM_REALF);
  PARAMS = addParameter(PARAMS,'rho0',rho0,PARM_REALF);
  PARAMS = addParameter(PARAMS,'f0',f0,PARM_REALF);    
  PARAMS = addParameter(PARAMS,'h_c',h_c,PARM_REALE);    
  PARAMS = addParameter(PARAMS,'theta_s',theta_s,PARM_REALF);    
  PARAMS = addParameter(PARAMS,'theta_b',theta_b,PARM_REALF);    
  PARAMS = addParameter(PARAMS,'Hsml',Hsml,PARM_REALF);    
  PARAMS = addParameter(PARAMS,'Hbbl',Hbbl,PARM_REALF);
  
  %%% Save biogeochemical parameters
  PARAMS = addParameter(PARAMS,'a_temp',a_temp,PARM_REALF);
  PARAMS = addParameter(PARAMS,'b_temp',b_temp,PARM_REALF);
  PARAMS = addParameter(PARAMS,'c_temp',c_temp,PARM_REALF);
  PARAMS = addParameter(PARAMS,'alpha',alpha,PARM_REALF);
  PARAMS = addParameter(PARAMS,'monod',monod,PARM_REALF);
  
  
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
  Hexp = 500;
  Tmax = 20 - 5*XX_tr/Lx;
  Tmin = 0;
  buoy_init = Tmin + (Tmax-Tmin).*(exp(ZZ_tr/Hexp+1)-exp(-H/Hexp+1))./(exp(1)-exp(-H/Hexp+1));

  %%% Plot initial buoyancy
  figure(fignum);
  fignum = fignum+1;
  pcolor(XX_tr,ZZ_tr,buoy_init);
  title('Initial Buoyancy')
  colorbar
  
  %%% Initial depth tracer
  dtr_init = ZZ_tr;
  
  %%% Initial nitrate profile (Hyperbolic)
  Nmax = 30; %%% Maximum concentration of nutrient at the ocean bed
  Ncline = 150; % Approximate guess of the depth of the nutracline
  N_init = -Nmax*tanh(ZZ_tr./Ncline);
    
  figure(fignum);
  fignum = fignum+1;
  pcolor(XX_tr,ZZ_tr,N_init);
%   plot(N_init(1,:),ZZ_tr(1,:))
  shading interp
  title('Initial Nitrate Profile')
  colorbar
  caxis([0 Nmax])
  
  
  figure(fignum);
  fignum = fignum+1;
  plot(N_init(1,:),ZZ_tr(1,:))
  title('Initial Nutracline Profile')
  
  %%% Store tracers in 3D matrix
  phi_init(1,:,:) = reshape(buoy_init,[1 Nx Nz]);
  phi_init(2,:,:) = reshape(dtr_init,[1 Nx Nz]);
  phi_init(3,:,:) = reshape(N_init,[1 Nx Nz]);  
  
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
 
  
%  tau = tau0*cos(pi*xx_psi/(2*Lx));
   temp = tau0*tanh(((Lx)-xx_psi)/(Lx/16));
  
   amp = 0.7846/4;                % Scaling amplitude for seasonal forcing
   per = 2*pi/52;               % Period for seasonal forcing of one year in seconds
   peak = 17;                   % Peak wind stress at the end of April (Haack, et al 2005).
   bb = 1.0392;                 % Shift so that the max wind stress is at 1.6 (April 30)
  
   %Use weekly averaged wind forcing (if this value is changed, it must be
   %changed in the mamebus.c code as well in the windInterp function.
  tyear = 0:1:52;
%   fcing = amp*(bb + cos((tyear-peak)*per));
  fcing = ones(size(tyear));            % Constant forcing to determine upwelling. 
  tlength = length(fcing);                        % Determine the number of points of wind stress data
  tau = zeros(length(fcing),length(xx_psi));
  
  for ii = 1:1:53
      tau(ii,:) = fcing(ii)*temp;
  end
  
  tauFile = 'tau.dat';  
  writeDataFile(fullfile(local_run_dir,tauFile),tau);
  PARAMS = addParameter(PARAMS,'tlength',tlength,PARM_INT);
  PARAMS = addParameter(PARAMS,'tauFile',tauFile,PARM_STR); 
  
  figure(fignum);
  fignum = fignum+1;
  plot(xx_psi,tau(1,:))
%   surf(tau)
  shading interp
  title('Surface Wind Stress')
  view(2)
  axis tight
%   colorbar
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Irradiance Profile %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  irr0 = 340;                   %%% W/m^2 (Can eventually include a seasonal amplitude)
  irr_scaleheight = 30;          %%% m
  irradiance = irr0.*exp(ZZ_tr./irr_scaleheight);
  irFile = 'irradiance.dat';
  writeDataFile(fullfile(local_run_dir,irFile),irradiance);
  PARAMS = addParameter(PARAMS,'irFile',irFile,PARM_STR);
  
%   figure(fignum)
%   fignum = fignum+1;
%   surf(XX_tr,ZZ_tr,irradiance)
%   title('Irradiance Profile')
%   shading interp
%   view(2)
%   colorbar
  
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
%   figure(fignum);
%   fignum = fignum+1;
%   pcolor(XX_psi,ZZ_psi,Kiso)
%   title('Isopycnal Diffusivity')
%   shading interp
%   colorbar
  
  %%% Another guess is a hyperbolic profile decreasing to 200 at the lower
  %%% boundary. 
  Kefold = 1000;
  
  Kiso = Kiso0 + (Kiso0-Kiso_hb)*tanh(ZZ_psi./Kefold);
  figure(fignum);
  fignum = fignum+1;
  pcolor(XX_psi,ZZ_psi,Kiso)
  title('Isopycnal Diffusivity')
  shading interp
  colorbar
  
  Kiso(:,1)
  
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

  %%% Crude mixed layers
  idx_sml = ZZ_psi>-Hsml;
  Kdia(idx_sml) = Kdia(idx_sml) + Ksml * -4*(ZZ_psi(idx_sml)/Hsml).*(ZZ_psi(idx_sml)/Hsml+1);  
  idx_bbl = ZZ_psi<-HB_psi+Hbbl;
  Kdia(idx_bbl) = Kdia(idx_bbl) + Kbbl * -4*((ZZ_psi(idx_bbl)+HB_psi(idx_bbl))/Hbbl).*((ZZ_psi(idx_bbl)+HB_psi(idx_bbl))/Hbbl-1); 
  
  % Check if sml and bbl overlap 
  depth = hb_psi;
  for jj = 1:length(depth)
      if abs(depth(jj)) < (Hsml + Hbbl)
        for kk = 1:Nz
        	Kdia(jj,kk) = Kdia0 + (Ksml+Kbbl)* 0.5*(-4*(ZZ_psi(jj,kk)/depth(jj))*(ZZ_psi(jj,kk)/depth(jj)+1) ...
                -4*((ZZ_psi(jj,kk)+HB_psi(jj,kk))/depth(jj))*((ZZ_psi(jj,kk)+HB_psi(jj,kk))/depth(jj)-1));
        end
      end
  end
  
  %%% Plot open-ocean diapycnal diffusivity
  figure(fignum);
  fignum = fignum+1;
  plot(Kdia(1,:),ZZ_psi(1,:));
  title('Open ocean diapycnal diffusivity')
  
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
  buoy_surf_max = 20;
  buoy_surf_min = 15;
  buoy_surf = buoy_surf_max + (buoy_surf_min-buoy_surf_max)*xx_tr/Lx;
  buoy_relax((xx_tr>=L_relax),Nz) = buoy_surf((xx_tr>=L_relax)); 
  T_relax_buoy((xx_tr>=L_relax),Nz) = 10*t1day; 
  
  figure(100)
  plot(xx_tr,buoy_surf)
  
  %%% Depth tracer relaxation  
  dtr_relax = dtr_init;
  T_relax_dtr = 5*t1year * ones(Nx,Nz);
  
  %%% Relax nitrate to initial conditions
  N_relax = N_init;
 
  %%% Store tracer relaxation data in 3D matrices
  phi_relax_all = zeros(Ntracs,Nx,Nz);
  phi_relax_all(1,:,:) = reshape(buoy_relax,[1 Nx Nz]);
  phi_relax_all(2,:,:) = reshape(dtr_relax,[1 Nx Nz]);
  phi_relax_all(3,:,:) = reshape(N_relax,[1 Nx Nz]);
  T_relax_all = zeros(Ntracs,Nx,Nz);
  T_relax_all(1,:,:) = reshape(T_relax_buoy,[1 Nx Nz]);
  T_relax_all(2,:,:) = reshape(T_relax_dtr,[1 Nx Nz]);
  T_relax_all(3,:,:) = -ones(size(T_relax_all(1,:,:))); % Nutrients are conserved

  relaxTracerFile = 'relaxTracer.dat';
  relaxTimeFile = 'relaxTime.dat';
  writeDataFile(fullfile(local_run_dir,relaxTracerFile),phi_relax_all);
  writeDataFile(fullfile(local_run_dir,relaxTimeFile),T_relax_all);
  PARAMS = addParameter(PARAMS,'relaxTracerFile',relaxTracerFile,PARM_STR);     
  PARAMS = addParameter(PARAMS,'relaxTimeFile',relaxTimeFile,PARM_STR);
  
  %%% Create a run script
  createRunScript (local_home_dir,run_name,exec_name, ...
                   use_cluster,cluster_username,cluster_address, ...
                   cluster_home_dir,walltime)

  %%% Create the input parameter file
  writeParamFile(pfname,PARAMS);    

  
end
