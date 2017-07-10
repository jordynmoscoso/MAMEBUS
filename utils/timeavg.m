%%%
%%% Reads in data from the output of 'Overturning' and takes the
%%% time-average of the output over the specified interval.
%%%
function [phi_avg,psi_avg] = timeavg (run_name,tmin,tmax)
  
  %%%%%%%%%%%%%%%%%%%%%
  %%%%% VARIABLES %%%%%
  %%%%%%%%%%%%%%%%%%%%%

  %%% Parameter and data file names
  run_name = strtrim(run_name);
  dirpath = fullfile('../runs',run_name);
  params_file = fullfile(dirpath,[run_name,'_in']);
  data_file = fullfile(dirpath,[run_name,'_out']);

  %%% Plotting grid
  [Ny Ny_found] = readparam(params_file,'Ny','%u');
  [Nz Nz_found] = readparam(params_file,'Nz','%u');
  if ((~Ny_found) || (~Nz_found))
    error('Could not read grid parameters');
  end   
  
  %%% Parameters related to number of iterations  
  maxTime = readparam(params_file,'maxTime','%lf');
  
  %%% Perform some sanity checks
  if (tmin > tmax)
    error('tmin must be <= tmax');
  end
  if (tmax > maxTime)
    error('tmax must be <= maxTime');
  end
  if (tmin < 0)
    error('tmin must be >= 0');
  end
  
  %%% To store averages
  phi_avg = zeros(Ny,Nz);
  psi_avg = zeros(Ny+1,Nz+1);
     
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% AVERAGING LOOP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  %%% Open the output file for reading      
  dfid = fopen(data_file,'r');
  if (dfid == -1)
    error('ERROR: Could not find data file.');
  end
  
  %%% Tracks whether we should still read data
  stillReading = true;
  counter = 0; 
  tprev = 0;
  
  %%% At each time iteration...
  while (stillReading)
      
    %%% Get the time value 
    t = fscanf(dfid,'%le',1);
    if (size(t)==0)      
      error(['timeavg ran out of data at t=',num2str(tprev)]);%%% computation      
    end
        
    %%% Get the phi values on the gridpoints
    phi = fscanf(dfid,'%le',[Ny,Nz]);            
    if (size(phi,1)~=Ny || size(phi,2)~=Nz)      
      error(['timeavg ran out of data at t=',num2str(t)]);
    end          

    %%% Residual streamfunction only available after first time step
    if (t > 0)      
      %%% Get the psi values on the gridpoints
      psi = fscanf(dfid,'%le',[Ny+1,Nz+1]);            
      if (size(psi,1)~=Ny+1 || size(psi,2)~=Nz+1)        
        error(['timeavg ran out of data at t=',num2str(t)]);
      end     
    end  
    
    %%% If we're in the right time bracket, add the data to the average
    if (t >= tmin)
      if (t <= tmax)        
        phi_avg = phi_avg + phi;
        psi_avg = psi_avg + psi;        
        counter = counter+1;           
      end
      if (t >= tmax)
        stillReading = false;
      end
    end             
    
    %%% Store last time iteration
    tprev = t;
    
  end
  
  %%% Equally-spaced data saves, so the time-average may be 
  %%% calculated as follows
  phi_avg = phi_avg / counter;
  psi_avg = psi_avg / counter;
    
  %%% Close the data file
  fclose(dfid); 

end