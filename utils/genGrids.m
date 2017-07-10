%%%
%%% genGrids.m
%%%
%%% Generates sigma-coordinate grids using the stretching function of A.
%%% Shchepetkin (2010).
%%%
%%% Nx - number of gridpoints in x-direction
%%% Nz - number of gridpoints in z-direction
%%% Lx - Domain length in x-direction
%%% h_c - a (positive) depth parameter controlling the range of depths over which 
%%%    the coordinates are approximately aligned with geopotentials.
%%% theta_s - surface stretching parameter, defined "meaningfully" between 0 and 10
%%% theta_b - bottom stretching parameter, defined "meaningfully" between 0 and 4
%%% hb_tr - Nx x 1 array of water column thicknesses on tracer points
%%% hb_psi - Nx+1 x 1 array of water column thicknesses on streamfunction points
%%%
function [XX_tr,ZZ_tr,XX_psi,ZZ_psi,XX_u,ZZ_u,XX_w,ZZ_w] ...
                     = genGrids(Nx,Nz,Lx,h_c,theta_s,theta_b,hb_tr,hb_psi)
  
  %%% Grid spacings
  dx = Lx/Nx; %%% Latitudinal grid spacing (in meters)
  ds = 1/Nz; %%% Sigma coordinate grid spacing
  
  %%% One-dimensional grids in y and sigma space
  xx_tr = dx/2:dx:Lx-dx/2; %%% Tracer latitudinal grid point locations
  ss_tr = -1+ds/2:ds:0-ds/2; %%% Tracer sigma grid point locations  
  xx_psi = 0:dx:Lx; %%% Streamfunction latitudinal grid point locations
  ss_psi = -1:ds:0; %%% Streamfunction sigma grid point locations
    
  %%% Meshgrids in y and sigma space
  [SS_tr,XX_tr] = meshgrid(ss_tr,xx_tr); %%% Tracer meshgrids       
  [SS_psi,XX_psi] = meshgrid(ss_psi,xx_psi); %%% Streamfunction meshgrids  
  [SS_u,XX_u] = meshgrid(ss_tr,xx_psi); %%% Zonal velocity meshgrids
  [SS_w,XX_w] = meshgrid(ss_psi,xx_tr); %%% Vertical velocity meshgrids          
  
  %%% Ocean depth meshgrids
  HB_tr = repmat(reshape(hb_tr,[Nx 1]),[1 Nz]);
  HB_psi = repmat(reshape(hb_psi,[Nx+1 1]),[1 Nz+1]);
  HB_u = repmat(reshape(hb_psi,[Nx+1 1]),[1 Nz]);
  HB_w = repmat(reshape(hb_tr,[Nx 1]),[1 Nz+1]);
  
  %%% Convert sigma grids to z grids
  ZZ_tr = stretch_ROMS(SS_tr,h_c,theta_s,theta_b,HB_tr);
  ZZ_psi = stretch_ROMS(SS_psi,h_c,theta_s,theta_b,HB_psi);
  ZZ_u = stretch_ROMS(SS_u,h_c,theta_s,theta_b,HB_u);
  ZZ_w = stretch_ROMS(SS_w,h_c,theta_s,theta_b,HB_w);
 
end

