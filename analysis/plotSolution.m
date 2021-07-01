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
function [XX_tr,ZZ_tr,XX_psi,ZZ_psi,avgVals,hb_psi,xx_psi] = plotSolution (local_home_dir,run_name,plot_trac,var_id,avgTime)

    %%% Load convenience functions
    addpath ../utils;
    
    fs = 22; % control the size of the font here.
    fsplt = fs;

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

    dx = Lx/Nx; %%% Latitudinal grid spacing (in meters)
    xx_psi = 0:dx:Lx;

    %%% Generate full sigma-coordinate grids
    [XX_tr,ZZ_tr,XX_psi,ZZ_psi,XX_u,ZZ_u,XX_w,ZZ_w] ...
                        = genGrids(Nx,Nz,Lx,h_c,theta_s,theta_b,hb_tr,hb_psi);

    
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
    
    LL = yearLength*avgTime;
    avgStart = round(lastVal - LL); % Calculate where we start the average.
    
    
    if avgStart < 1
        avgStart = 0;
    end
    
    
    if (plot_trac)
        if (var_id == 0)
            title_name = 'Zonal Velocity';
        elseif (var_id == 1)
            title_name = 'Meridional Velocity';
        elseif (var_id == 2)
            title_name = 'Temperature ($^o$C)';
        else
            % title names
            switch (modeltype)
                case 0
                    title_name = 'Depth Tracer';
                case 1 % npzd
                    if (var_id == 3)
                        title_name = 'Nitrate (mmol/m$^3$)';
                    elseif (var_id == 4)
                        title_name = 'Phytoplankton (mg Chl/m$^3$)';
                    elseif (var_id == 5)
                        title_name = 'Zooplankton (mmol N/m$^3$)';
                    elseif (var_id == 6)
                        title_name = 'Detritus (mmol N/m$^3$)';
                    else
                        title_name = 'Depth Tracer';
                    end
                otherwise
                    disp('No Modeltype Found')
                    title_name = 'Phytoplankton (mg Chl/m$^3$)';
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
    else
        avgVals = zeros(Nx+1,Nz+1);
    end
    
    
    %%% Averaging loop
    ndt = lastVal - avgStart + 1; % add one to account for Matlab Indexing
    for n = avgStart:lastVal
        if (plot_trac)
            data_file = fullfile(dirpath,['TRAC',num2str(var_id),'_n=',num2str(n),'.dat']);
            phi = readOutputFile(data_file,Nx,Nz);
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
    
    
    Nmax = 30; %%% Maximum concentration of nutrient at the ocean bed
    Ncline = 250; % Approximate guess of the depth of the nutracline
    Ninit(:,:,1) = -Nmax*tanh(ZZ_tr./Ncline);

    
    
    %%% Plot the average values %%%
    timelengthstr = lastVal*dt_s/t1year - avgStart*dt_s/t1year;
    titlestr = title_name;
    figure(7)
    if (plot_trac)
        if (var_id == 0 || var_id == 1)
            cmap = cmocean('balance');
            pcolor(XX_tr,ZZ_tr,avgVals)
%             shading interp
            c = colorbar; 
            colormap(cmap);
            c.TickLabelInterpreter = 'latex';
            set(gca,'TickLabelInterpreter','latex')
            set(gca,'FontSize',fsplt)
            ylabel('Depth (m) ','FontSize',fs,'interpreter','latex')
            xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
%             set(gcf,'Position',[100 100 724 654])
            xticks(0:50e3:400e3)
            xticklabels({'400' '350' '300' '250' '200' '150' '100' '50' '0'})
            maxspeed = 0.01;
            minval = max(max(max(avgVals)),maxspeed);
            minval = abs(min(min(min(avgVals)),-minval));
%             yticks(-180:20:0)
%             yticklabels({'-180' '-160' '-140' '-120' '-100' '-80' '-60' '-40' '-20' '0'})
            axis([min(min(XX_tr))+50e3 max(max(XX_tr)) -H 0])
%             minval = 0.35;
            caxis([-minval minval]);
%             caxis([-0.8 0.8])
            title(titlestr,'interpreter','latex')
            shading interp
        else % plot tracers 
            % chose the appropriate colormap and contour vector
            if (var_id == 2)
                cmap = cmocean('thermal');
                cvec = 0:4:22;
                clr = 'w';
            elseif (var_id == 3) % nitrogen
                cmap = cmocean('matter');
                cvec = [1 6.5 12.1 17.6 23.1];
                clr = 'w';
            elseif (var_id == 4) % phytoplankton
                cmap = cmocean('speed');
                cvec = [0.05 0.2 0.4 0.9 1.8];
                clr = 'k';
            elseif (var_id == 5) % zooplankton
                cmap = cmocean('amp');
                cvec = 5;
                clr = 'w';
            elseif (var_id == 6) % detritus
                cmap = cmocean('turbid');
                cvec = 5;
                clr = 'k';
            else
                disp('Not explicitly chosen yet')
                cmap = cmocean('speed');
                cvec = [0.2 0.4 0.9 1.8];
                clr = 'k';
            end
                
            pcolor(XX_tr,ZZ_tr,avgVals);
            hold on
            [C h] = contour(XX_tr,ZZ_tr,avgVals,cvec,clr);
            clabel(C,h,'Color','w')
            hold off
            shading interp
            c = colorbar; 
            colormap(cmap);
            c.TickLabelInterpreter = 'latex';
            set(gca,'TickLabelInterpreter','latex')
            set(gca,'FontSize',fsplt)
            ylabel('Depth (m) ','FontSize',fs,'interpreter','latex')
            xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
%             set(gcf,'Position',[100 100 724 654])
            xticks(0:50e3:400e3)
            xticklabels({'400' '350' '300' '250' '200' '150' '100' '50' '0'})
            yticks(-180:20:0)
            yticklabels({'-180' '-160' '-140' '-120' '-100' '-80' '-60' '-40' '-20' '0'})
            axis([min(min(XX_tr))+50e3 max(max(XX_tr)) -180 0])
            title(titlestr,'interpreter','latex')
        end
    else
        cmap = cmocean('balance');
        psi_r_lim = avgVals;
        limval = 300;
        psi_r_lim = min(psi_r_lim,limval);
        psi_r_lim = max(psi_r_lim,-limval);
        pcolor(XX_psi,ZZ_psi,psi_r_lim)
        shading interp
        c = colorbar; 
        colormap(cmap);
        c.TickLabelInterpreter = 'latex';
        set(gca,'TickLabelInterpreter','latex')
        set(gca,'FontSize',fsplt)
        ylabel('Depth (m) ','FontSize',fs,'interpreter','latex')
        xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
        xticks(0:50e3:400e3)
        xticklabels({'400' '350' '300' '250' '200' '150' '100' '50' '0'})
        title(titlestr,'interpreter','latex')
%         caxis([min(min(psi_r_lim)) max(max(psi_r_lim))])
        caxis([-limval limval])
        axis([min(min(XX_tr))+50e3 max(max(XX_tr)) -H 0])
    end
   
        hold on
        plot(xx_psi,-hb_psi,'k')
        hold off
end