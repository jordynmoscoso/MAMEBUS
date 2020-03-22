%%%
%%% Reads in data from the output of 'Overturning' and makes a movie 
%%% of the solution.
%%%
%%% local_home_dir specifies the directory in the local system in which
%%% run files are stored. N.B. this function will search within a 
%%% subdirectory called 'run_name' to find the run's output files.
%%%
%%% run_name specifies the name of the run.
%%%
%%% If plot_strat is true then the density will be plotted, whereas if it
%%% is false then the streamfunction will be plotted.
%%%
%%% var_id Specifies the tracer number to plot (if plot_trac is true) or
%%% the streamfunction to plot (if plot_trac is false).
%%%
function [XX_tr,ZZ_tr,XX_psi,ZZ_psi,avgVals] = plotSolution (local_home_dir,run_name,plot_trac,var_id,avgTime)

    %%% Load convenience functions
    addpath ../utils;
    addpath ./redblue
    
    fs = 24; % control the size of the font here.
    fsplt = 24;
    m1km = 1000;

    mov_on = 0;
    mov_name = strcat(run_name,'_name');

    if var_id > 3
      disp('Plotting biogeochemistry solutions')
    end
    
    if (var_id > 3 && ~plot_trac)
        plot_trac = true;
        disp('Check "plot_trac" and switch to "true" for this variable')
    end
    
    show_w = true; % plots vertical velocities
    plot_uptake = false;
    
    plot_density = false;
    T0 = 20; % reference temperature of water
    alpha = 1e-4; % thermal coefficient
    rho0 = 1000;
    rhosw = 1026;
    if var_id == 2
        plot_density = false;
        if plot_density
            disp('Plotting density')
        end
    end

    %%%%%%%%%%%%%%%%%%%%%
    %%%%% VARIABLES %%%%%
    %%%%%%%%%%%%%%%%%%%%% 

    %%% Parameter and data file names
    run_name = strtrim(run_name);
    dirpath = fullfile(local_home_dir,run_name);
    params_file = fullfile(dirpath,[run_name,'_in']);  
    
    [dt_s dt_s_found] = readparam(params_file,'monitorFrequency','%f');
    [modeltype modeltype_found] = readparam(params_file,'bgcModel','%u');
    
    %%% Plotting grid
    [Nx Nx_found] = readparam(params_file,'Nx','%u');
    [Nz Nz_found] = readparam(params_file,'Nz','%u');
    [Lx Lx_found] = readparam(params_file,'Lx','%lf');
    [H H_found] = readparam(params_file,'Lz','%lf');
    if ((~Nx_found) || (~Nz_found) || (~Lx_found) || (~H_found))
    error('Could not read grid parameters');
    end    

    dx = Lx/Nx;
    %%% Read grid parameters
    [h_c h_c_found] = readparam(params_file,'h_c','%le');
    [theta_s theta_s_found] = readparam(params_file,'theta_s','%lf');
    [theta_b theta_b_found] = readparam(params_file,'theta_b','%lf');

    %%% Read bottom topography
    hb = readDataFile (params_file,dirpath,'topogFile',Nx+2,1,H*ones(Nx+2,1));
    hb_psi = 0.5*(hb(1:end-1)+hb(2:end));  
    hb_tr = hb(2:end-1); %%% Remove "ghost" points

    %%% Parameters related to number of iterations
    
    [startTime startTime_found] = readparam(params_file,'startTime','%lf');
    [endTime endTime_found] = readparam(params_file,'endTime','%lf'); 
    [restart restart_found] = readparam(params_file,'restart','%d');
    [n0 n0_found] = readparam(params_file,'startIdx','%u');

    %%% Default is that we're not picking up from a previous simulation
    if (~restart_found)
    restart = false;
    end

    %%% Default is that we pickup from the first output file
    if ((~restart) || (~n0_found))
    n0 = 0;
    end

    %%% If the start time isn't specified then it may be specified implicitly
    %%% by the pickup file
    if (~startTime_found)
    if (restart && dt_s_found)
      startTime = n0*dt_s;
    else
      startTime = 0;
    end
    end

    %%% For convenience
    t1year = 365*86400; %%% Seconds in one year
    t1day = 86400;
    
    CtoChl = 100/12; %mmolC /mg Chl
    NtoChl = 0.8; %mg Chl/mmol N
    

    dx = Lx/Nx; %%% Latitudinal grid spacing (in meters)
    xx_psi = 0:dx:Lx;

    %%% Generate full sigma-coordinate grids
    [XX_tr,ZZ_tr,XX_psi,ZZ_psi,XX_u,ZZ_u,XX_w,ZZ_w] ...
                        = genGrids(Nx,Nz,Lx,h_c,theta_s,theta_b,hb_tr,hb_psi);

    if (var_id ~= 0 && var_id ~= 1)
     n_var = var_id; 
    end
  
    %%% Max number of output files is the number of whole time steps that can
    %%% fit into tmax, plus one initial save, plus one final save
    Noutput = floor(endTime/dt_s) + 2;

    
    %%% Upload the names of the files in order to pick out the max N value
    A = dir(dirpath);
    M = size(A);
    M = floor(max(M)/6); % there are six output variables saved
                         % so this cuts down on the for loop.
    lastVal = 0;
    
    for ii = 1:M
        temp = A(ii).name;
        temp = strsplit(temp,'=');
        if max(size(temp) > 1)
            temp2 = char(temp(2));
            compVal = strsplit(temp2,'.');
            compVal = str2num(compVal{1});

            if compVal > lastVal
                lastVal = compVal;
            end
        end
    end
    
    
    %%% Pick out the last N value in the saved files. 
    outputFrac = dt_s/t1year;
    yearLength = 1/outputFrac; % Count the number of N values in one year.
    
    
%     avgTime = 10; % Choose the amount of time to average over (in years)
    LL = yearLength*avgTime;
    avgStart = round(lastVal - LL); % Calculate where we start the average.
    
    
    if avgStart < 1
        avgStart = 0;
    end
    
    ndt = 0;
    
    if (plot_trac)
        if (var_id == 0)
            title_name = 'Zonal Velocity';
        elseif (var_id == 1)
            title_name = 'Meridional Velocity';
        elseif (var_id == 2)
            title_name = 'Buoyancy';
            if plot_density
                title_name = 'Density (kg/m$^3$)';
            end
        else
            % title names
            switch (modeltype)
                case 0
                    title_name = 'Depth Tracer';
                case 1 % nitrate only
                    if (var_id == 3)
                        title_name = 'Nitrate (mmol/m^3)';
                    else
                        title_name = 'Depth Tracer';
                    end
                case 2 % npzd
                    if (var_id == 3)
                        title_name = 'Nitrate (mmol/m$^3$)';
                    elseif (var_id == 4)
                        title_name = 'Phytoplankton (mg Chl/m$^3$)';
                    elseif (var_id == 5)
                        title_name = 'Zooplankton (mmol N/m$^3$)';
                    elseif (var_id == 6)
                        title_name = 'Detritus (mmol N/m$^3$)';
                    elseif (var_id == 7)
                        title_name = 'Passive Tracer';
                    else
                        title_name = 'Depth Tracer';
                    end
                otherwise
                    title_name = 'Depth Tracer';
            end
        end
    else
        if (var_id == 0)
            title_name = 'Residual Streamfunction';
        elseif (var_id == 1)
            title_name = 'Mean Streamfunction';
        else 
            title_name = 'Eddy Streamfunction';
        end
    end
        
            
    
    
    disp(['Start time is ' num2str(avgStart*dt_s/t1year) ' yrs']);
    disp(['End time is ' num2str(lastVal*dt_s/t1year) ' yrs']);
    
    %%% Create a placeholder for avgVals
    if (plot_trac)
        avgVals = zeros(Nx,Nz);
        
        if (modeltype == 1 || modeltype == 2)
            avgUptake = zeros(Nx,Nz);
            
            %%% BGC parameters
            qsw = 340;
            kw = 0.04;
            I0 = 0.45*qsw;
            T0 = 20;
            r = 0.05;
            umax = 2.6;
            kpar = kw;
            
            %%% Light Profile.
            IR = I0*exp(kpar.*ZZ_tr);
            IR = IR./(sqrt(IR.^2 + I0.^2));
        end
    else
        avgVals = zeros(Nx+1,Nz+1);
    end
    
%     lastVal = 1200;
%     avgStart = lastVal - LL;
    
    %%% Averaging loop
    ndt = lastVal - avgStart + 1; % add one to account for Matlab Indexing
    for n = avgStart:lastVal
        if (plot_trac)
            data_file = fullfile(dirpath,['TRAC',num2str(var_id),'_n=',num2str(n),'.dat']);
            phi = readOutputFile(data_file,Nx,Nz);
            
            
            
            %%% calculate the uptake offline for the single nitrate bgc model
            if (modeltype == 1 && var_id == 3)
                
                data_file = fullfile(dirpath,['TRAC',num2str(var_id),'_n=',num2str(n),'.dat']);
                N = readOutputFile(data_file,Nx,Nz);
                
                buoy_id = var_id-1;
                buoy_file = fullfile(dirpath,['TRAC',num2str(buoy_id),'_n=',num2str(n),'.dat']);
                buoy = readOutputFile(buoy_file,Nx,Nz);
                
                t_uptake = exp(r.*(buoy - T0));

                uptake = umax.*t_uptake.*IR.*N;
                
                avgUptake = avgUptake + log10(uptake);
%                 avgUptake = avgUptake + uptake;
                plot_uptake = true;
            end
            
            

        else
            %%% Load different streamfunctions      
            switch (var_id)
                case 0 %%% Residual streamfunction
                    data_file = fullfile(dirpath,['PSIR_n=',num2str(n),'.dat']);
                case 1 %%% Mean streamfunction
                    data_file = fullfile(dirpath,['PSIM_n=',num2str(n),'.dat']);
                case 2 %%% Eddy streamfunction
                    data_file = fullfile(dirpath,['PSIE_n=',num2str(n),'.dat']);
            end
            
          %%% Get the psi values on the gridpoints
          phi = readOutputFile (data_file,Nx+1,Nz+1);
        end
        

        avgVals = avgVals + phi;
    end
    avgVals = avgVals/ndt;
    
%     if (modeltype == 1 || modeltype == 2 && plot_trac == true)
%         avgUptake = avgUptake/ndt;
%     end
    
    
    Nmax = 30; %%% Maximum concentration of nutrient at the ocean bed
    Ncline = 250; % Approximate guess of the depth of the nutracline
    Ninit(:,:,1) = -Nmax*tanh(ZZ_tr./Ncline);

    
    if plot_density
        tempvals = avgVals;
        avgVals = 23.5 + rho0*(alpha*(T0 - tempvals));
    end

%     v = avgVals;
%     save('meridveloc.mat','v')
    
    %%% Plot the average values %%%
    timelengthstr = lastVal*dt_s/t1year - avgStart*dt_s/t1year;
    titlestr = title_name;
    figure
    if (plot_trac)
        if (var_id == 0 || var_id == 1)
            cmap = cmocean('balance');
            pcolor(XX_tr,ZZ_tr,avgVals)
%             shading interp
            h = colorbar; 
            colormap(cmap)
            h.TickLabelInterpreter = 'latex';
            maxspeed = 0.01;
            minval = max(max(max(avgVals)),maxspeed);
            minval = abs(min(min(min(avgVals)),-minval));
            caxis([-minval minval]);
            title(titlestr,'interpreter','latex')
            shading interp
        elseif (var_id == 2) % plot buoyancy or density
            if plot_density
                cmap = cmocean('dense');
%                 cvec = round(linspace(23,max(max(avgVals)),10),2);
                cvec = [23.9 24.6 25.2 25.8 26.4];
                c = colorbar; 
                c.TickLabelInterpreter = 'latex';
                pcolor(XX_tr,ZZ_tr,avgVals)
                shading interp
                hold on
                [C h] = contour(XX_tr,ZZ_tr,avgVals,cvec,'-w');
                hold off
                axis([min(min(XX_tr)) max(max(XX_tr)) -180 0])
                caxis([23.4 27])
                set(gca,'TickLabelInterpreter','latex')
                set(gca,'FontSize',fsplt)
                set(gcf,'Position',[100 100 724 654])
                ylabel('Depth (m) ','FontSize',fs,'interpreter','latex')
                xticks([0:50e3:350e3])
                xticklabels({'$350$' '300' '250' '200' '150' '100' '50' '0'})
                xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
                yticks(-180:20:0)
                yticklabels({'-180' '-160' '-140' '-120' '-100' '-80' '-60' '-40' '-20' '0'})
                colorbar off
            else
                cmap = cmocean('thermal');
                pcolor(XX_tr,ZZ_tr,avgVals)
                cvec = 0:2:18;
                shading interp
                hold on
                [C h] = contour(XX_tr,ZZ_tr,avgVals+1,cvec,'-w');
                hold off
                axis([min(min(XX_tr)) max(max(XX_tr)) -180 0])
            end
            clabel(C,h,'Color','w')
%             pcolor(XX_tr,ZZ_tr,avgVals)
%             shading interp;
            colorbar; 
            colormap(cmap);
            title(titlestr,'interpreter','latex')
%             axis([min(min(XX_tr)) max(max(XX_tr)) -200 0])
%             clabel(C,h,'Color','w');  
%             set(h,'ShowText','on'); 
        else % plot biogeochemistry
%             [C h] = contourf(XX_tr,ZZ_tr,avgVals,[0:2:30]);
%             clabel(C,h,'Color','w');   
%             set(h,'ShowText','on'); 
%             set(h,'FontSize',18);
            if (var_id == 4 && modeltype == 2)
                avgVals = avgVals*NtoChl;
            end
            clr = 'w';
            % chose the appropriate colormap
            if (var_id == 3)
                cmap = cmocean('matter');
                cvec = [1 6.5 12.1 17.6 23.1];
                
            elseif (var_id == 4)
                cmap = cmocean('speed');
                cvec = [0.2 0.4 0.9 1.8];
                clr = 'k';
            elseif (var_id == 5)
                cmap = cmocean('amp');
                cvec = 5;
            elseif (var_id == 6)
                cmap = cmocean('turbid');
                
                cvec = 5;
            else
                disp('Not explicitly chosen yet')
                cmap = cmocean('deep');
                cvec = 5;
            end
                
            pcolor(XX_tr,ZZ_tr,avgVals);
            hold on
            [C h] = contour(XX_tr,ZZ_tr,avgVals,cvec,clr);
            clabel(C,h,'Color',clr)
            hold off
            shading interp
            c = colorbar; 
            colormap(cmap);
            c.TickLabelInterpreter = 'latex';
            set(gca,'TickLabelInterpreter','latex')
            set(gca,'FontSize',fsplt)
            ylabel('Depth (m) ','FontSize',fs,'interpreter','latex')
            xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
            set(gcf,'Position',[100 100 724 654])
            xticks(0:50e3:350e3)
            xticklabels({'$350$' '300' '250' '200' '150' '100' '50' '0'})
            yticks(-180:20:0)
            yticklabels({'-180' '-160' '-140' '-120' '-100' '-80' '-60' '-40' '-20' '0'})
%             caxis([min(min(avgVals)) max(max(avgVals))])
%             caxis([0 0.11])
            axis([min(min(XX_tr)) max(max(XX_tr)) -180 0])
            title(titlestr,'interpreter','latex')
            plot_uptake = 0;
            
            if (var_id == 3 && modeltype == 1 && plot_uptake == 1)
                titlestr = [ 'Uptake log(mmol N/m^3)/d averaged over the final ', num2str(timelengthstr) , ' year(s)'];
                figure
                pcolor(XX_tr,ZZ_tr,avgUptake);
                shading interp
                colorbar; 
                colormap jet;
%                 caxis([min(min(avgVals)) max(max(avgVals))])
                caxis([-1 1])
                title(titlestr)
                
            elseif (var_id == 4 && modeltype == 2 && plot_uptake == 1)
                titlestr = [ 'Uptake (mmol N/m^3/d) averaged over the final ', num2str(timelengthstr) , ' year(s)'];
                figure
                pcolor(XX_tr,ZZ_tr,avgUptake);
                shading interp
                colorbar; 
                colormap jet;
%                 caxis([min(min(avgVals)) max(max(avgVals))])
%                 caxis([0 5])
                title(titlestr)
            end
            
            
            
            
        end
        
%         set(h,'FontSize',18);
        
        % Use a divergent colorscheme if we are looking at positive or
        % negative values
        if (var_id == 0 || var_id == 1)
            colormap redblue;
        end
    else
        psi_r_lim = avgVals;
        limval = 2;
        psi_r_lim = min(psi_r_lim,limval);
        psi_r_lim = max(psi_r_lim,-limval);
%       	contourf(XX_psi,ZZ_psi,psi_r_lim,[-2:0.5:2]); 
        pcolor(XX_psi,ZZ_psi,psi_r_lim)
        shading interp
        h = colorbar; 
        colormap redblue;
        caxis([-limval limval]);
        title(titlestr)
    end
    
%     psie = psi_r_lim;
    save('XX_psi_hi.mat','XX_psi')
    save('XX_tr_hi.mat','XX_tr')
    save('ZZ_psi_hi.mat','ZZ_psi')
    save('ZZ_tr_hi.mat','ZZ_tr')
    save('hb_psi_hi.mat','hb_psi')
    save('x_psi_hi.mat','xx_psi')
    
    
    

    
%         title(titlestr)
        hold on
        plot(xx_psi,-hb_psi,'k')
        hold off
        

%     if (show_w)
%         avgVals = zeros(Nx+1,Nz+1);
%         wvel = zeros(Nx,Nz+1);
%         dx = Lx/Nx;
%         
%         for n = avgStart:lastVal
%             data_file = fullfile(dirpath,['PSIR_n=',num2str(n),'.dat']);
%             phi = readOutputFile (data_file,Nx+1,Nz+1);
%             for j = 1:Nx
%                 wvel(j,:) = (avgVals(j+1,:) - avgVals(j,:))/dx;
%             end
%             avgVals = avgVals + phi;
%         end
%         
%         % calculate the average of the residual stream function
%         avgVals = avgVals/ndt;
%         
%         
%         for j = 1:Nx
%             wvel(j,:) = (avgVals(j+1,:) - avgVals(j,:))/dx;
%         end
%         
%         figure
%         pcolor(XX_w,ZZ_w,wvel)
%         shading interp
%         colormap redblue;
%         colorbar
%         
%     end
        
        
        
    %%% TO DO: UPDATE FUNCTION TO INCLUDE VERTICAL VELOCITIES
end