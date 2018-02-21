%%% 
%%% sfc_wind_stress.m
%%%
%%% Sets the surface wind stress profile to be read in by MAMEBUS. This
%%% separate file allows for users to choose an idealized surface wind
%%% stress or load in a wind stress profile. Note that, the "Seasonal" 
%%% wind stress profile is also an idealized profile.
%%% 
%%% MAMEBUS interpolates data at every time step so that weekly/monthly 
%%% data can be used. 
%%%

function [tau,tlength] = sfc_wind_stress(tau0,Lx,xx_psi)
   %%% Determine the type of wind forcing used by MAMEBUS:
   IDEALIZED = 1;
   SEASONAL = 0;
   DATA = 0;
       
   
   %%% Determine the averaging on the wind stress (daily, monthly, etc.)
   %%% How many data points per year?
   tlength = 52;
   tyear = 0:1:tlength-1;

   %%% Build the wind stress profile based on the wind forcing above
   if (DATA)
       %%% To be updated when data is present
       
       %%% Must include a match so that the data "sits" correctly on the
       %%% grid, along with have the correct number of data points on the
       %%% time average. 

   elseif (SEASONAL)         
       %%% Idealized structure of the wind stress profile "sharp" profile.
       temp = tau0*tanh(((Lx)-xx_psi)/(Lx/16));
  
       amp = 0.7846/4;              % Scaling amplitude for seasonal forcing
       per = 2*pi/52;               % Period for seasonal forcing of one year in seconds
       peak = 17;                   % Peak wind stress at the end of April (Haack, et al 2005).
       bb = 1.0392;                 % Shift so that the max wind stress is at 1.6 (April 30)
  
       fcing = amp*(bb + cos((tyear-peak)*per));
       tau = zeros(length(fcing),length(xx_psi));
  
       for ii = 1:1:53
           tau(ii,:) = fcing(ii)*temp;
       end
   elseif (IDEALIZED)
       temp = tau0*tanh(((Lx)-xx_psi)/(Lx/4));
%        temp = tau0*ones(1,length(xx_psi));

%        m1km = 1000; %%% Meters in 1 km    
%        Lmax = 200*m1km;
%        temp = tau0*cos(pi/2*(xx_psi-Lmax)/Lmax).^2;
%        temp(xx_psi>2*Lmax) = 0;
       
       
       fcing = ones(size(tyear));   % Constant forcing to determine upwelling. 
       tau = zeros(length(fcing),length(xx_psi));
  
       for ii = 1:1:tlength
           tau(ii,:) = fcing(ii)*temp;
       end
   else
       disp('No wind stress profile is determined')
       return
   end

  
end