%%%
%%% readTracer.m
%%%
%%% Reads in data from the output of 'Tracer.exe' and returns the final
%%% data dump, i.e. the converged solution.
%%%
function [cc tend] = readTracer (run_name)
  
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
    cc = fscanf(dfid,'%le',[Ny,Nz]);            
    if (size(cc,1)~=Ny || size(cc,2)~=Nz)      
      error(['readTracer ran out of data at t=',num2str(t)]);
    end          
            
  end
    
  %%% Close the data file
  fclose(dfid); 

end