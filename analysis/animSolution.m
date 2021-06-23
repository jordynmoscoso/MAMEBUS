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
function M = animSolution (local_home_dir,run_name,plot_trac,var_id,...
                            mov_on,mov_name)
 
  %%% Load convenience functions
  addpath ../utils;

  
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
  if (~restart || ~n0_found)
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
  
  %%% If user is plotting nitrate, indicate what to plot
%   if (var_id == 2 && plot_trac)
%     prompt = 'Please indicate display \n 0 for Primary Productivity (default) \n 1 for Nitrate Concentration \n';
%     ncase = input(prompt);
%     %%% Check for valid input, if not choose default
%     if (isempty(ncase) || ncase < 0 || ncase > 1)
%         ncase = 0;
%     end
    Nstore = [];
    tsave = [];
%   end

ncase = 1; 
  
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
  
  %%% At each time iteration...
  while (stillReading)
      
    %%% Get the time value 
    t = startTime + (n-n0)*dt_s;
    
    %%% If plot_trac==true, load the tracer specified by trac_num and plot
    %%% it
    if (plot_trac)    

      %%% Data file name
      data_file = fullfile(dirpath,['TRAC',num2str(var_id),'_n=',num2str(n),'.dat']);
      phi = readOutputFile(data_file,Nx,Nz);    
      
      %%% Plot the tracer     
      switch (var_id)
        case 0 %%% u-velocity
          cmap = cmocean('balance');
          pcolor(XX_tr,ZZ_tr,phi)
%             shading interp
          h = colorbar; 
          colormap(cmap)
          h.TickLabelInterpreter = 'latex';
          maxspeed = 0.01;          
          caxis([-maxspeed maxspeed]);
          shading interp
        case 1 %%% v-velocity
          cmap = cmocean('balance');
          pcolor(XX_tr,ZZ_tr,phi)
%             shading interp
          h = colorbar; 
          colormap(cmap)
          h.TickLabelInterpreter = 'latex';
          maxspeed = 0.3;
          caxis([-maxspeed maxspeed]);
          shading interp
        case 2 %%% Buoyancy (temperature)
          size(XX_tr)
          size(phi)
%           [C h] = contourf(XX_tr/1000,ZZ_tr,phi,0:0.5:10);
          pcolor(XX_tr/1000,ZZ_tr,phi);
          shading interp
          set(gca, 'CLim', [0, 20]);
            colormap(cmocean('thermal',20));
          h=colorbar;
          title('Potential temp. (^oC)');          
%       caxis([0 20]);
%       axis([0 1 -1 0]);

          otherwise 
%           [C h] = contourf(XX_tr/1000,ZZ_tr,phi,0:0.5:10);
          pcolor(XX_tr/1000,ZZ_tr,phi);
          shading interp
%           set(gca, 'CLim', [0, 1]);
            colormap(cmocean('speed',20));
          h=colorbar;
          title(['BGC, t=', num2str(t/t1year)]);          
%       caxis([0 20]);
          axis([min(min(XX_tr/1000)) max(max(XX_tr/1000)) -300 0]);
      end
%       clabel(C,h,'Color','w');  
%       set(h,'ShowText','on'); 
%       pcolor(XX_phi,ZZ_phi,phi);
            
    
      set(h,'FontSize',18);
      


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
    xlabel('x (km)');    
    ylabel('z (km)','Rotation',0);        
%     set(gca,'XTick',(0:Lx/5:Lx)/1000);
%     set(gca,'YTick',-H:H/5:0);
%     title(strcat(['t=',num2str(round(t/t1day)),' days (',num2str(round(t/t1year)),' yr)', ' mmol N/m^3']));           
    
    nextframe = getframe(gcf);    
    M(counter) = nextframe; 
    
    
    % Stores the frame to the video file.
    if mov_on
        writeVideo(vidObj,nextframe);
        
    end
    
    counter = counter+1;   
    n = n + 1;
        
  end    
  
  if (var_id ~= 0 && var_id ~= 1 && plot_trac)
      titlestr = ['Final time = ', num2str(t/t1year), ' yr'];
      
      figure(nfig+1)
      size(tsave)
      size(Nstore)
      plot(tsave/t1year,Nstore)
      xlabel(titlestr)
  end
  
  % Close the movie writer
  if mov_on
      close(vidObj);
  end
  
end