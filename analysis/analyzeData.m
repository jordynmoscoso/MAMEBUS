%%%
%%% Reads in data from the output of 'Overturning' and makes a movie 
%%% of the velocity field
%%%
%%% local_home_dir specifies the directory in the local system in which
%%% run files are stored. N.B. this function will search within a 
%%% subdirectory called 'run_name' to find the run's output files.
%%%
%%% run_name specifies the name of the run.
%%%
%%% This program is used to analyze the data from MAMEBUS runs.
%%%
function M = analyzeData(local_home_dir,run_name,var_id,velocity)
    %%% Calcuate the velocity field from the residual streamfunction
    %%% Load convenience functions
    addpath ../utils;
    addpath ./redblue

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
    dt_s = readparam(params_file,'monitorFrequency','%lf');
    tmax = readparam(params_file,'maxTime','%lf');
  
    %%% For convenience
    t1year = 365*86400; %%% Seconds in one year
    t1day = 86400;
    
    %%% Generate full sigma-coordinate grids
    [XX_tr,ZZ_tr,XX_psi,ZZ_psi,XX_u,ZZ_u,XX_w,ZZ_w] ...
                        = genGrids(Nx,Nz,Lx,h_c,theta_s,theta_b,hb_tr,hb_psi);
                    
    %%% Allocate vectors to hold the data
    sx = size(XX_u);
    dx = zeros(size(XX_tr));
    dz = dx;
    u = dx;
    w = dz;
    for ii = 1:sx(2)
        dx(:,ii) = XX_u(2:end,ii) - XX_u(1:end-1,ii);
        dz(ii,:) = ZZ_w(ii,2:end) - ZZ_w(ii,1:end-1);
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% PLOTTING LOOP %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%% Max number of time steps is the number of whole time steps that can
    %%% fit into tmax, plus one initial save, plus one final save
    Nt = floor(tmax/dt_s) + 2;
  
    %%% Make a movie of the data - allocate a movie array large enough to
    %%% hold the largest possible number of frames
    figure(1);
    clf;
    axes('FontSize',18);
    M = moviein(Nt);  
  
    %%% Tracks whether we should still read data
    stillReading = true;
    counter = 1;
    n = 0;
  
    % Determines whether or not to make a movie and writes a new file to
    % store the data.
%     if mov_on
%     	vidObj = VideoWriter(mov_name);
%         vidObj.FrameRate = 10;
%         open(vidObj)
%     end
    
     %%% At each time iteration...
    while (stillReading)
      
      %%% Get the time value 
      t = n*dt_s;
    
      %%% Load different streamfunctions      
      switch (var_id)
        case 0 %%% Residual streamfunction
          data_file = fullfile(dirpath,['PSIR_n=',num2str(n),'.dat']);
        case 1 %%% Mean streamfunction
          data_file = fullfile(dirpath,['PSIM_n=',num2str(n),'.dat']);
        case 2 %%% Eddy streamfunction
          data_file = fullfile(dirpath,['PSIE_n=',num2str(n),'.dat']);
      end

      %%% Open the output file for reading    
      dfid = fopen(data_file,'r');
      if (dfid == -1)
        stillReading = false;
        continue;
        %%% Ignore any missing data - might just be the end of the
        %%% computation
      end

      %%% Get the phi values on the gridpoints
      psi = fscanf(dfid,'%le',[Nx+1,Nz+1]);           
      if (size(psi,1)~=Nx+1 || size(psi,2)~=Nz+1)
        error(['ERROR: Could not find data file: ',data_file]);
      end    
      
      %%% Close data file
      fclose(dfid);
      
      %%% Calculate the zonal and vertial velocities with respect to the
      %%% stream function.
      for ii = 1:sx(2)
        u(:,ii) = (psi(2:end,ii) - psi(1:end-1,ii))./dx(:,ii);
        w(ii,:) = (psi(ii,2:end) - psi(ii,1:end-1))./dz(ii,:);
      end
      
      %%% Determine which velocity to plot
      if strcmp(velocity,'zonal')
        %%% Plot the zonal velocity
        psi_r_lim = psi;
        limval = 10e-4;
        psi_r_lim = min(psi_r_lim,limval);
        psi_r_lim = max(psi_r_lim,-limval);
%       [C h] = contourf(XX_psi,ZZ_psi,psi_r_lim,-limval:limval/40:limval,'EdgeColor','k');  
        figure(1)
        pcolor(XX_tr,ZZ_tr,u);
        shading interp;     
        colormap redblue;
        h=colorbar;        
        caxis([-limval limval]);
        set(h,'FontSize',18);
        hold on
        plot(XX_tr,-hb_tr,'k','LineWidth',2)
        hold off
      
        %%% Store the image in the movie buffer  
        xlabel('x');    
        ylabel('z','Rotation',0);        
        axis tight;
        set(gca,'XTick',0:Lx/5:Lx);
        set(gca,'YTick',-H:H/5:0);
        title(strcat(['t=',num2str(round(t/t1day)),' days (',num2str(round(t/t1year)),' yr)']));   
      else
        %%% Plot the vertical velocity otherwise
        psi_r_lim = psi;
        limval = 0.25;
        psi_r_lim = min(psi_r_lim,-limval);
        psi_r_lim = max(psi_r_lim,limval);
%       [C h] = contourf(XX_psi,ZZ_psi,psi_r_lim,-limval:limval/40:limval,'EdgeColor','k');  
        figure(1)
        pcolor(XX_tr,ZZ_tr,w);
        shading interp;     
        colormap redblue;
        h=colorbar;        
        caxis([-limval limval]);
        set(h,'FontSize',18);
        hold on
        plot(XX_tr,-hb_tr,'k','LineWidth',1)
        hold off
      
        %%% Store the image in the movie buffer  
        xlabel('x');    
        ylabel('z','Rotation',0);        
        axis tight;
        set(gca,'XTick',0:Lx/5:Lx);
        set(gca,'YTick',-H:H/5:0);
        title(strcat(['t=',num2str(round(t/t1day)),' days (',num2str(round(t/t1year)),' yr)']));  
      end
    
      nextframe = getframe(gcf);    
     M(counter) = nextframe; 
        
      counter = counter+1;   
      n = n + 1;
    
%     if mov_on
%         close(vidObj)
%     end
    end
end