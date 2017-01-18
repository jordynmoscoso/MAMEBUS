%%%
%%% Reads in data from the output of 'Overturning' and makes a movie 
%%% of the solution.
%%%
%%% run_name specifies the name of the run.
%%%
%%% If plot_strat is true then the density will be plotted, whereas if it
%%% is false then the streamfunction will be plotted.
%%%
%%% trac specifies the tracer number to plot
%%%
function M = animSolution (run_name,plot_trac,trac_num)
 
  %%% Load convenience functions
  addpath ../utils;
  addpath ./redblue

  %%%%%%%%%%%%%%%%%%%%%
  %%%%% VARIABLES %%%%%
  %%%%%%%%%%%%%%%%%%%%% 

  %%% Parameter and data file names
  run_name = strtrim(run_name);
  dirpath = fullfile('../runs',run_name);
  params_file = fullfile(dirpath,[run_name,'_in']);  

  %%% Plotting grid
  [Nx Nx_found] = readparam(params_file,'Nx','%u');
  [Nz Nz_found] = readparam(params_file,'Nz','%u');
  [Lx Lx_found] = readparam(params_file,'Lx','%lf');
  [H H_found] = readparam(params_file,'Lz','%lf');
  if ((~Nx_found) || (~Nz_found) || (~Lx_found) || (~H_found))
    error('Could not read grid parameters');
  end  
  dx = (Lx/Nx);
  dz = (H/Nz);
  [xx_phi zz_phi XX_phi ZZ_phi] = createmesh(0.5*dx,Lx-0.5*dx,Nx,-H+0.5*dz,-0.5*dz,Nz);  
  [xx_psi zz_psi XX_psi ZZ_psi] = createmesh(0,Lx+dx,Nx+1,-H,0,Nz+1);  
  
  %%% Read bottom topography
  hb = readDataFile (params_file,dirpath,'topogFile',Nx+2,1,H*ones(Nx+2,1));
  hb_psi = 0.5*(hb(1:end-1)+hb(2:end));  
  hb_phi = hb(2:end-1); %%% Remove "ghost" points
  
  %%% Parameters related to number of iterations
  dt_s = readparam(params_file,'monitorFrequency','%lf');
  tmax = readparam(params_file,'maxTime','%lf');
  
  %%% For convenience
  t1year = 365*86400; %%% Seconds in one year
  
  %%% Modify plotting grids to account for topography
  for j=1:Nx
    ZZ_phi(j,:) = ZZ_phi(j,:) * hb_phi(j)/H;
  end
  for j=1:Nx+1
    ZZ_psi(j,:) = ZZ_psi(j,:) * hb_psi(j)/H;
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
  
  %%% At each time iteration...
  while (stillReading)
      
    %%% Get the time value 
    t = n*dt_s;
    
    %%% If plot_trac==true, load the tracer specified by trac_num and plot
    %%% it
    if (plot_trac)    

      %%% Data file name
      data_file = fullfile(dirpath,['TRAC',num2str(trac_num),'_n=',num2str(n),'.dat']);

      %%% Open the output file for reading    
      dfid = fopen(data_file,'r');
      if (dfid == -1)
        stillReading = false;
        continue;
        %%% Ignore any missing data - might just be the end of the
        %%% computation
      end

      %%% Get the phi values on the gridpoints
      phi = fscanf(dfid,'%le',[Nx,Nz]);            
      if (size(phi,1)~=Nx || size(phi,2)~=Nz)
        error(['ERROR: Could not find data file: ',data_file]);
      end          

      %%% Close data file
      fclose(dfid);
      
      %%% Plot the tracer     
      switch (trac_num)
        case 0 %%% Buoyancy (temperature)
          [C h] = contourf(XX_phi,ZZ_phi,phi,0:1:20);
        case 1 %%% Depth tracer
          [C h] = contourf(XX_phi,ZZ_phi,phi,-(0:200:H));
      end
      clabel(C,h,'Color','w');      
      set(h,'ShowText','on'); 
%       pcolor(XX_phi,ZZ_phi,phi);
            
      colormap jet;
      h=colorbar;
%       caxis([0 20]);
      set(h,'FontSize',18);
%       axis([0 1 -1 0]);
      
    %%% If plot_trac==false, plot the residual streamfunction
    else    
    
      %%% Data file name
      data_file = fullfile(dirpath,['PSIR_n=',num2str(n),'.dat']);
      
      %%% Open the data file for reading    
      dfid = fopen(data_file,'r');
      if (dfid == -1)
        stillReading = false;
        continue;
        %%% Ignore any missing data - might just be the end of the
        %%% computation
      end             
      
      %%% Get the psi values on the gridpoints
      psi_r = fscanf(dfid,'%le',[Nx+1,Nz+1]);           
      if (size(psi_r,1)~=Nx+1 || size(psi_r,2)~=Nz+1);
        error(['ERROR: Could not find data file: ',data_file]);
      end    
      
      %%% Close data file
      fclose(dfid);

      %%% Plot the streamfunction
      psi_r_lim = psi_r;
      limval = 2;
      psi_r_lim = min(psi_r_lim,limval);
      psi_r_lim = max(psi_r_lim,-limval);
      [C h] = contourf(XX_psi,ZZ_psi,psi_r_lim,-limval:limval/40:limval,'EdgeColor','k');                 
      colormap redblue;
      h=colorbar;        
      caxis([-limval limval]);
      set(h,'FontSize',18);

    end
      
    
    %%% Store the image in the movie buffer          
    xlabel('x');    
    ylabel('z','Rotation',0);        
    axis tight;
    set(gca,'XTick',0:Lx/5:Lx);
    set(gca,'YTick',-H:H/5:0);
    title(strcat(['t=',num2str(round(t/t1year)),' yr']));           \
    
    nextframe = getframe(gcf);    
    M(counter) = nextframe; 
    
    counter = counter+1;   
    n = n + 1;
    
  end    
  
end