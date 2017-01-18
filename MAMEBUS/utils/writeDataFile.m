%%%
%%% writeDataFile.m
%%%
%%% Writes the array of real values specified by 'data' to the file 
%%% specified by the string file name 'fname'.
%%%
function writeDataFile (fname,data)
  
  fid = fopen(fname,'w');
  if (fid == -1)
    error(['Could not open ',fname]);
  end
  fprintf(fid,'%e ',data);
  fclose(fid);

end

