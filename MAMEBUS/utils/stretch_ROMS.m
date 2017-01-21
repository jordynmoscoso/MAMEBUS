%%%
%%% stretch_ROMS.m
%%%
%%% Implements the topographic stretching function of A. Shchepetkin (2010). 
%%%
%%% sigma is the fractional vertical stretching coordinate, ranging from -1 to 0.
%%% h_c is a (positive) depth parameter controlling the range of depths over which 
%%%    the coordinates are approximately aligned with geopotentials.
%%% theta_s is the surface stretching parameter, defined "meaningfully" between 0 and 10
%%% theta_b is the bottom stretching parameter, defined "meaningfully" between 0 and 4
%%% h_b is the thickness of the fluid column
%%%
%%%
function z = stretch_ROMS (sigma, h_c, theta_s, theta_b, h_b)

  %%% Surface refinement function
  if (theta_s > 0)  
    C = (1 - cosh(theta_s*sigma)) / (cosh(theta_s) - 1);  
  else  
    C = - sigma.*sigma;
  end
  
  %%% Augment surface refinement function to include bottom refinement
  if (theta_b > 0)  
    C = (exp(theta_b*C) - 1) / (1 - exp(-theta_b));
  end
  
  %%% Calculate actual depth
  z = h_b .* (h_c*sigma + h_b.*C) ./ (h_c + h_b);
  
end