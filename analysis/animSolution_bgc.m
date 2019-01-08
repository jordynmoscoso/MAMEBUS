%%%
%%% Reads in data from the output of 'Overturning' and makes a movie 
%%% of the solution.
%%%
%%% local_home_dir specifies the directory in the local system in which
%%% run files are stored. N.B. this function will search within a 
%%% subdirectory called 'run_name' to find the run's output files.
%%%
%%% run_name specifies the name of the run.
%%%
%%% If plot_strat is true then the density will be plotted, whereas if it
%%% is false then the streamfunction will be plotted.
%%%
%%% var_id Specifies the tracer number to plot (if plot_trac is true) or
%%% the streamfunction to plot (if plot_trac is false).
%%%
function M = animSolution_bgc (local_home_dir,run_name,plot_trac,var_id,plot_NPP)
 
  %%% Load convenience functions
  addpath ../utils;
  addpath ./redblue
  
  mov_on = 0;
  mov_name = strcat(run_name,'_name');
  
  if var_id < 3
      var_id = 3;
      disp('Defaulting to show Nitrate')
  end

  
  %%%%%%%%%%%%%%%%%%%%%
  %%%%% VARIABLES %%%%%
  %%%%%%%%%%%%%%%%%%%%% 

  %%% Parameter and data file names
  run_name = strtrim(run_name);
  dirpath = fullfile(local_home_dir,run_name);
  params_file = fullfile(dirpath,[run_name,'_in']);  

  %%% Plotting grid
  [Nx Nx_found] = readparam(params_file,'Nx','%u');
  [Nz Nz_found] = readparam(params_file,'Nz','%u');
  [Lx Lx_found] = readparam(params_file,'Lx','%lf');
  [H H_found] = readparam(params_file,'Lz','%lf');
  if ((~Nx_found) || (~Nz_found) || (~Lx_found) || (~H_found))
    error('Could not read grid parameters');
  end    
  
  %%% Read grid parameters
  [h_c h_c_found] = readparam(params_file,'h_c','%le');
  [theta_s theta_s_found] = readparam(params_file,'theta_s','%lf');
  [theta_b theta_b_found] = readparam(params_file,'theta_b','%lf');
  
  %%% Read bottom topography
  hb = readDataFile (params_file,dirpath,'topogFile',Nx+2,1,H*ones(Nx+2,1));
  hb_psi = 0.5*(hb(1:end-1)+hb(2:end));  
  hb_tr = hb(2:end-1); %%% Remove "ghost" points
  
  %%% Parameters related to number of iterations
  [dt_s dt_s_found] = readparam(params_file,'monitorFrequency','%lf');
  [startTime startTime_found] = readparam(params_file,'startTime','%lf');
  [endTime endTime_found] = readparam(params_file,'endTime','%lf'); 
  [restart restart_found] = readparam(params_file,'restart','%d');
  [n0 n0_found] = readparam(params_file,'startIdx','%u');

  %%% Default is that we're not picking up from a previous simulation
  if (~restart_found)
    restart = false;
  end

  %%% Default is that we pickup from the first output file
  if ((~restart) || (~n0_found))
    n0 = 0;
  end
  
  %%% If the start time isn't specified then it may be specified implicitly
  %%% by the pickup file
  if (~startTime_found)
    if (restart && dt_s_found)
      startTime = n0*dt_s;
    else
      startTime = 0;
    end
  end
  
  %%% For convenience
  t1year = 365*86400; %%% Seconds in one year
  t1day = 86400;
  
%%% Load grids from model output
%   dx = (Lx/Nx);
%   dz = (H/Nz);
%   [xx_phi zz_phi XX_tr ZZ_tr] = createmesh(0.5*dx,Lx-0.5*dx,Nx,-H+0.5*dz,-0.5*dz,Nz);  
%   [xx_psi zz_psi XX_psi ZZ_psi] = createmesh(0,Lx+dx,Nx+1,-H,0,Nz+1);  
%   fid = fopen(fullfile(dirpath,['ZZ_PHI.dat']),'r');
%   if (fid == -1)
%     error(['Could not open ',paramFile]);
%   end
%   ZZ_tr = fscanf(fid,'%f',[Nx Nz]);
%   fclose(fid);
%   fid = fopen(fullfile(dirpath,['ZZ_PSI.dat']),'r');
%   if (fid == -1)
%     error(['Could not open ',paramFile]);
%   end
%   ZZ_psi = fscanf(fid,'%f',[Nx Nz]);
%   fclose(fid);  

    Nstore = [];
    tsave = [];
    
      dx = Lx/Nx; %%% Latitudinal grid spacing (in meters)
  xx_psi = 0:dx:Lx;
    
  %%% Generate full sigma-coordinate grids
  [XX_tr,ZZ_tr,XX_psi,ZZ_psi,XX_u,ZZ_u,XX_w,ZZ_w] ...
                        = genGrids(Nx,Nz,Lx,h_c,theta_s,theta_b,hb_tr,hb_psi);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% PLOTTING LOOP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%
  surf_layer = 75;
  sl_ind = find(ZZ_tr > - surf_layer);
  denom = size(sl_ind);
  
  if (var_id ~= 0 && var_id ~= 1)
     n_var = var_id; 
  end
  
  %%% Max number of output files is the number of whole time steps that can
  %%% fit into tmax, plus one initial save, plus one final save
  Noutput = floor(endTime/dt_s) + 2;
  
  %%% Make a movie of the data - allocate a movie array large enough to
  %%% hold the largest possible number of frames
  nfig = 107;
  figure(nfig);
  clf;
  axes('FontSize',18);
  M = moviein(Noutput);  
  
  %%% Tracks whether we should still read data
  stillReading = true;
  counter = 1;
  n = n0;
  
  % Determines whether or not to make a movie and writes a new file to
  % store the data.
  if mov_on
      vidObj = VideoWriter(mov_name);
      vidObj.FrameRate = 10;
      open(vidObj)
  end
  
  
  
  
  % Since we are focusing on plotting biogeochemical tracers, calculate the
  % temperature and irradiance profiles
  qsw = 340;
  hsml = 50;
  kpar = 0.04;
  r = 0.05;
  T0 = 20;
  kn = 0.1;
  
  % Light dependent uptake
  I0 = 0.45*qsw;
  IR = I0*exp(kpar*ZZ_tr);
  
  II = IR./sqrt(IR.^2 + I0.^2);
  
  figure(100)
  surf(XX_tr,ZZ_tr,II)
  shading interp
  view(2)
  colorbar
  
  f = 0.01;
    a_temp = 0.6;
    b_temp = 1.066;
    c_temp = 1;
  alpha = 0.025; % s^-1 (W/m^2)^-1
  %%% At each time iteration...
  
  
      data_file = fullfile(dirpath,['TRAC',num2str(var_id),'_n=',num2str(1795),'.dat']);
      phi = readOutputFile(data_file,Nx,Nz);
      
      buoy_file = fullfile(dirpath,['TRAC',num2str(2),'_n=',num2str(1795),'.dat']);
      buoy = readOutputFile(buoy_file,Nx,Nz);
  
      figure(200)
          pcolor(XX_tr,ZZ_tr,phi)
          shading interp
          colorbar
          colormap jet;
          h=colorbar;
          hold on
          plot(xx_psi,-hb_psi,'k')
          hold off
          
          figure(201)
          contourf(XX_tr,ZZ_tr,phi,[0:5:10 10:1:24 24:2:30])
          shading interp
          colorbar
          colormap jet;
          h=colorbar;
          hold on
          plot(xx_psi,-hb_psi,'k')
          hold off
          
          figure(202)
          contourf(XX_tr,ZZ_tr,buoy,10)
          shading interp
          colorbar
          colormap default;
          h=colorbar;
          hold on
          plot(xx_psi,-hb_psi,'k')
          hold off
          pause
          
  while (stillReading)
      
    %%% Get the time value 
    t = startTime + (n-n0)*dt_s;
    
    %%% If plot_trac==true, load the tracer specified by trac_num and plot
    %%% it
    if (plot_trac)    
        
        titlestr = ['t = ', num2str(t/t1year), ' yr'];
        
      buoy_file = fullfile(dirpath,['TRAC',num2str(2),'_n=',num2str(n),'.dat']);
      buoy = readOutputFile(buoy_file,Nx,Nz);
      Tlim = exp(r*(buoy-T0));

      %%% Data file name
      data_file = fullfile(dirpath,['TRAC',num2str(var_id),'_n=',num2str(n),'.dat']);
      phi = readOutputFile(data_file,Nx,Nz);
      
      R = zeros(Nx,Nz);

      
      if (plot_NPP)
          NPP = Tlim.*II.*phi;
            
          %%% Plot NPP     
          figure(107)
          pcolor(XX_tr,ZZ_tr,NPP)
          shading interp
          colorbar
          title(titlestr)
          colormap jet;
          h=colorbar;
          title(titlestr)
          hold on
          plot(xx_psi,-hb_psi,'k')
          hold off
      else
          %%% Plot the tracer     
          pcolor(XX_tr,ZZ_tr,phi)
          shading interp
          colorbar
          title(titlestr)
          colormap jet;
          h=colorbar;
          title(titlestr)
          hold on
          plot(xx_psi,-hb_psi,'k')
          hold off
      end

          
%       caxis([0 20]);
      set(h,'FontSize',18);
%       axis([0 1 -1 0]);
% axis([0 Lx -500 0])


    %%% If plot_trac==false, plot the streamfunction
    else    
    
      %%% Load different streamfunctions      
      switch (var_id)
        case 0 %%% Residual streamfunction
          data_file = fullfile(dirpath,['PSIR_n=',num2str(n),'.dat']);
        case 1 %%% Mean streamfunction
          data_file = fullfile(dirpath,['PSIM_n=',num2str(n),'.dat']);
        case 2 %%% Eddy streamfunction
          data_file = fullfile(dirpath,['PSIE_n=',num2str(n),'.dat']);
      end
                
      %%% Get the psi values on the gridpoints
      psi = readOutputFile (data_file,Nx+1,Nz+1);       

      %%% Plot the streamfunction
      psi_r_lim = psi;
      limval = 2;
      psi_r_lim = min(psi_r_lim,limval);
      psi_r_lim = max(psi_r_lim,-limval);
%       [C h] = contourf(XX_psi,ZZ_psi,psi_r_lim,-limval:limval/40:limval,'EdgeColor','k');  
      figure(1)
      pcolor(XX_psi,ZZ_psi,psi_r_lim);
      shading interp;     
      colormap redblue;
      h=colorbar;        
      caxis([-limval limval]);
      set(h,'FontSize',18);

    end
      
    
    %%% Store the image in the movie buffer  
    xlabel('Distance (km)');    
    ylabel('Depth (m)');        
%     axis tight;
%     set(gca,'XTick',0:Lx/5:Lx);
%     set(gca,'YTick',-H:H/5:0);       
    
    nextframe = getframe(gcf);    
    M(counter) = nextframe; 
    
    
    % Stores the frame to the video file.
    if mov_on
        writeVideo(vidObj,nextframe);
        
    end
    
    counter = counter+1;   
    n = n + 1;
        
  end    
  
        figure(101)
      surf(XX_tr,ZZ_tr,Tlim)
      shading interp
      view(2)
      colorbar
  
  
  if (var_id ~= 0 && var_id ~= 1 && plot_trac)
      titlestr = ['Final time = ', num2str(t/t1year), ' yr'];
      
      figure(nfig+1)
      size(tsave)
      size(Nstore)
      plot(tsave/t1year,Nstore)
      title('D (small) Surface Layer Averaged (mmol/m^3)')
      xlabel(titlestr)
  end
  
  % Close the movie writer
  if mov_on
      close(vidObj);
  end
  
end