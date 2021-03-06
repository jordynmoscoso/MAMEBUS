%%%
%%% readOutputFile.m
%%%
%%% Reads a 2D array of real model output data values for the parameter specified
%%% by 'paramName'. The specified file is opened, and an Nx by Nz array of data is
%%% read from it. If the output file cannot be found then an empty matrix
%%% is returned.
%%%
function data = readOutputFile (outputFileName,Nx,Nz)
  
  fid = fopen(outputFileName,'r','b');
  if (fid == -1)
    data = [];
    return;
  end
  data = fread(fid,[Nx Nz],'real*8','ieee-le');
  fclose(fid);
    
  if ((size(data,1) ~= Nx) || (size(data,2) ~= Nz))
    error(['Insufficient data found in ',outputFileName]);
  end   

end

