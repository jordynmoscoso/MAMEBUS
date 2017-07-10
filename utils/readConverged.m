%%%
%%% readConverged.m
%%%
%%% Reads in data from the output of 'Overturning' and returns the final
%%% data dump, i.e. the converged solution.
%%%
function [phi psi tend] = readConverged (run_name)
  
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
  
  %%% Default return values
  phi = zeros(Ny,Nz);
  psi = zeros(Ny+1,Nz+1);  
     
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
  
  %%% At each time iteration...
  while (stillReading)
      
    %%% Get the time value 
    t = fscanf(dfid,'%le',1);
    if (size(t)==0)      
      break;
    end
    
    %%% Store last time save found
    tend = t;
        
    %%% Get the phi values on the gridpoints
    phi = fscanf(dfid,'%le',[Ny,Nz]);            
    if (size(phi,1)~=Ny || size(phi,2)~=Nz)      
      error(['readConverged ran out of data at t=',num2str(t)]);
    end          

    %%% Residual streamfunction only available after first time step
    if (t > 0)      
      %%% Get the psi values on the gridpoints
      psi = fscanf(dfid,'%le',[Ny+1,Nz+1]);            
      if (size(psi,1)~=Ny+1 || size(psi,2)~=Nz+1)        
        error(['readConverged ran out of data at t=',num2str(t)]);
      end     
    end  
            
  end
    
  %%% Close the data file
  fclose(dfid); 

end