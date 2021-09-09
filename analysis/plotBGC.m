%
% plotBGC.m
%
%
%
%

function [XX_tr,ZZ_tr,XX_psi,ZZ_psi,avgVals,hb_psi,xx_psi] = plotBGC (local_home_dir,run_name,bgc_id,avgTime,plot_type,plot_roem,overlay)
    %%% Load convenience functions
    addpath ../utils;
    
    fs = 22; % control the size of the font here.
    fsplt = fs;
    
    % index values and plot types
    idx_nitrate = 0;        % (default)
    idx_phyto = 1;
    idx_zoo = 2;
    idx_det = 3;
    
    idx_totconc = 0;    % plots the total concentration (default)
    idx_size = 1;       % plots concentration weighted size
    idx_nums = 2;       % plots number of size classes present above a certain threshold value, tol
    
    tol = 1e-2;         % tolerance for the concentration for number of size classes. 
    
    overlay_trac = 0;
    overlay_b = 1;
    overlay_n = 2;
    

    %%%%%%%%%%%%%%%%%%%%%
    %%%%% VARIABLES %%%%%
    %%%%%%%%%%%%%%%%%%%%% 

    %%% Parameter and data file names
    run_name = strtrim(run_name);
    dirpath = fullfile(local_home_dir,run_name);
    params_file = fullfile(dirpath,[run_name,'_in']);  
    
    [dt_s dt_s_found] = readparam(params_file,'monitorFrequency','%f');
    [modeltype modeltype_found] = readparam(params_file,'bgcModel','%u');
    [NP NP_found] = readparam(params_file,'MP','%u');
    [NZ NZ_found] = readparam(params_file,'MZ','%u');
    
    disp([num2str(NP) ' phytoplankton classes'])
    disp([num2str(NZ) ' zooplankton classes'])
    
    
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
   
    
    lpFile = fullfile(dirpath,'lpFile.dat');
    lp = readOutputFile(lpFile,NP,1);
    
    
    lzFile = fullfile(dirpath,'lzFile.dat');
    lz = readOutputFile(lzFile,NZ,1);
    
    %%%%%%%%% CODE TO COMPARE ROEM TO A SIMILAR RUN IN SSEM
    if plot_roem
        % assumes a standard labeling for ROEM and SSEM
        vals = split(run_name,'_');

%         ref_run = char(strcat('~/Desktop/PhysTests-Su21/0802/runs/ssem_dx_',(vals(3)),'_',vals(4),'_',(vals(5)))); 
%         run_name2 = char(strcat('ssem_dx_',(vals(3)),'_',vals(4),'_',(vals(5))));
        ref_run = char('~/Desktop/PhysTests-Su21/0802/runs/ssem_dx_0.3_tau_0.075');
        run_name2 = char('ssem_dx_0.3_tau_0.075');
        params_file2 = fullfile(ref_run,[run_name2,'_in']); 
        [NPs NP_found] = readparam(params_file2,'MP','%u');
        [NZs NZ_found] = readparam(params_file2,'MZ','%u');

        lpFile = fullfile(ref_run,'lpFile.dat');
        lpssem = readOutputFile(lpFile,NPs,1);
        
        lzFile = fullfile(ref_run,'lzFile.dat');
        lzssem = readOutputFile(lzFile,NZs,1);
        
        %%% Upload the names of the files in order to pick out the max N value
        A = dir(ref_run);
        M = size(A);
        M = floor(max(M)/6); % there are six output variables saved
                             % so this cuts down on the for loop.
        lastVal2 = 0;
        for ii = 1:M
            temp = A(ii).name;
            temp = strsplit(temp,'=');
            if max(size(temp) > 1)
                temp2 = char(temp(2));
                compVal2 = strsplit(temp2,'.');
                compVal2 = str2num(compVal2{1});

                if compVal2 > lastVal2
                    lastVal2 = compVal2;
                end
            end
        end
        
        if (bgc_id == idx_nitrate)
            title_name = 'Nitrate (mmol/m$^3$)';
            inds2 = 1;
            disp('Plotting Nitrate')
            cmap = cmocean('matter');
        elseif (bgc_id == idx_phyto)
            title_name = 'Phytoplankton (mg Chl/m$^3$)';
            inds2 = 2:2+NPs-1;
            disp('Plotting Phytoplankton')
            cmap = cmocean('speed');
            Ntracs = NP;
        elseif (bgc_id == idx_zoo)
            disp('Plotting Zooplankton')
            title_name = 'Zooplankton (mmol N/m$^3$)';
            tol = 1e-6; % uptate tolerance if zooplankton
            cmap = cmocean('dense');
            inds2 = 2+NPs:2+NPs+NZs-1;
            Ntracs = NZ;
        elseif (bgc_id == idx_det)
            title_name = 'Detritus (mmol N/m$^3$)';
            inds2 = 2+NPs+NZs;
            disp('Plotting Detritus')
            cmap = cmocean('turbid');
        else
            disp('Default to Nutrient')
            title_name = 'Nitrate (mmol/m$^3$)';
            inds2 = 1;
        end
        inds2 = inds2+2;
    end
    

    
    
    % based on the bgc_id, determine the plotting tile and calculate the correct indices
    if (bgc_id == idx_nitrate)
        title_name = 'Nitrate (mmol/m$^3$)';
        inds = 1;
        disp('Plotting Nitrate')
        cmap = cmocean('matter');
    elseif (bgc_id == idx_phyto)
        title_name = 'Phytoplankton (mg Chl/m$^3$)';
        inds = 2:2+NP-1;
        disp('Plotting Phytoplankton')
        cmap = cmocean('speed');
    elseif (bgc_id == idx_zoo)
        disp('Plotting Zooplankton')
        title_name = 'Zooplankton (mmol N/m$^3$)';
        tol = 1e-6; % uptate tolerance if zooplankton
        cmap = cmocean('dense');
        inds = 2+NP:2+NP+NZ-1;
    elseif (bgc_id == idx_det)
        title_name = 'Detritus (mmol N/m$^3$)';
        inds = 2+NP+NZ;
        disp('Plotting Detritus')
        cmap = cmocean('turbid');
    else
        disp('Default to Nutrient')
        title_name = 'Nitrate (mmol/m$^3$)';
        inds = 1;
    end
    
    % add to the index to include physical tracers
    inds = inds+2;
    Ntracs = length(inds);
    
    disp(['Start time is ' num2str(avgStart*dt_s/t1year) ' yrs']);
    disp(['End time is ' num2str(lastVal*dt_s/t1year) ' yrs']);
    
    
    t1year = 365*86400; %%% Seconds in one year
    %%% Create a placeholder for avgVals
    avgVals = zeros(Nx,Nz,Ntracs);
    ssemVals = zeros(Nx,Nz,length(inds));
    
    PMAT = zeros(Nx,Nz,Ntracs);
    lxvec = zeros(Nx,1);
    PLTMAT = zeros(Nx,Ntracs);
    
%     lastVal = min(lastVal2,lastVal);
    
    if plot_roem
        PLTSSEM = zeros(Nx,NPs);
        for var_id = inds2
            data_file = fullfile(ref_run,['TRAC',num2str(var_id),'_n=',num2str(lastVal2),'.dat']);
            phi_ref = readOutputFile(data_file,Nx,Nz);  

            ssemVals(:,:,var_id+1-min(inds2)) = phi_ref;
        end
    end
    %%% Averaging loop
%     ndt = lastVal - avgStart + 1; % add one to account for Matlab Indexing
%     for n = 1:2*365
    n = lastVal;
        ndt = 1;
            for var_id = inds
                data_file = fullfile(dirpath,['TRAC',num2str(var_id),'_n=',num2str(n),'.dat']);
                phi = readOutputFile(data_file,Nx,Nz);    
%                 avgVals(:,:,var_id+1-min(inds)) = avgVals(:,:,var_id+1-min(inds)) + phi;
                avgVals(:,:,var_id+1-min(inds)) = phi;
                

                
                % load in all of the plankton into a matrix
                if plot_type == idx_nums
                    PMAT(:,:,var_id+1-min(inds)) = phi;
                end
            end
    %     end
        avgVals = avgVals/ndt;
        totVals = squeeze(sum(avgVals,3));

        if plot_type == idx_size
            if (bgc_id == idx_nitrate || bgc_id == idx_det)
                pltVals = totVals;
            else
                arevals = sum(avgVals >= tol,3);
                arevals = (arevals >= tol);
                dev = totVals;
                dev(dev <= tol) = 1;
                pltVals = zeros(Nx,Nz);
                if bgc_id == idx_phyto
                    ll = lp;
                elseif bgc_id == idx_zoo
                    ll = lz;
                end
                for ii = 1:Ntracs
                   pltVals = pltVals + avgVals(:,:,ii).*ll(ii)./dev; 
                end
                pltVals = pltVals.*arevals;
            end
        elseif plot_type == idx_nums
            % reference depth
            H_ref = -20; %m


            % Organize data
            for ii = 1:Nx
                [~,ind] = min(abs(H_ref - ZZ_tr(:,ii)));
    %             disp([ZZ_tr(ind,ii), -hb_tr(ii)])
                lxvec(ii) = XX_tr(ii,1);
                if -hb_tr(ii) > H_ref
                    PLTMAT(ii,:) = NaN*ones(1,Ntracs);
                    if plot_roem
                        PLTSSEM(ii,:) = NaN*ones(1,NPs);
                    end
                else
                    PLTMAT(ii,:) = squeeze(PMAT(ii,ind,:));
                    if plot_roem
                        PLTSSEM(ii,:) = squeeze(ssemVals(ii,ind,:));
                    end
                end
            end

            if plot_roem
                pltVals = sum(avgVals,3);
                pltVals = squeeze(pltVals);
                
                
                % count the peaks in the reference run
                addpath ~/Desktop/SSEM/utils/
                [NUMP,sind] = countPeaks(PLTSSEM(end,:),1e-5);
                
                SSEMMAT = zeros(Nx,NP);
                
                for ii = 1:NUMP
                    SSEMMAT(:,ii) = sum( PLTSSEM(:,sind(1,ii):sind(2,ii)), 2);
                end
            end
        else
            % plot the total nutrient concentration
            if (bgc_id == idx_phyto || bgc_id == idx_zoo)
                pltVals = sum(avgVals,3);
            else
                pltVals = avgVals;
            end
            pltVals = squeeze(pltVals);
        end

        nlines = 5;
        if bgc_id == idx_phyto
            minval = 0.05;
            contvec = [0.01 0.025 0.1 0.2 0.5 0.8 1.2];
            cmax = 1.5;
        else
            minval = 0.01;
            contvec = 0.1*[0.01 0.1 0.2 0.5 0.8 1.2];
            cmax = .08;
        end
        dN = (max(max(totVals)) - minval)/nlines;
%         contvec = minval:dN:max(max(totVals));
        
        dp = -75;
        figure(100); clf;

        % load in buoyancy
        if overlay == overlay_b
            idx_buoy = 2;
            ocont = 0:2:22;
        elseif (overlay == overlay_n) 
            idx_buoy = 3;
            ocont = [0 1 2, 4 6 10 10:5:30];
        else
            idx_buoy = 3+NP+NZ+1;
            ocont = [0:100:300 300:200:600 600:500:4500];
        end
        data_file = fullfile(dirpath,['TRAC',num2str(idx_buoy),'_n=',num2str(n),'.dat']);
        buoy_ref = readOutputFile(data_file,Nx,Nz);  
    
        if plot_type == idx_size
            A = subplot(1,3,1);
            pcolor(XX_tr,ZZ_tr,totVals);
            hold on
            [C h] = contour(XX_tr,ZZ_tr,totVals,round(contvec,2),'k'); % plankton values
            plot(xx_psi,-hb_psi,'k')
            hold off
            shading interp
            a = colorbar;
            a.TickLabelInterpreter = 'latex';
            axis([min(min(XX_tr))+100e3 max(max(XX_tr)) dp 0])
            colormap(A,cmap)
            clabel(C,h,'Color','k')
            title('Total Concentration','interpreter','latex')
            ylabel('Depth (m)','interpreter','latex')
    %         xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
            xticks(0:50e3:400e3)
            xticklabels({'400' '350' '300' '250' '200' '150' '100' '50' '0'})
            set(gca,'FontSize',24,'TickLabelInterpreter','latex')


            B = subplot(1,3,3);
            pcolor(XX_tr,ZZ_tr,(log10(pltVals)))
            hold on
            plot(xx_psi,-hb_psi,'k')
            hold off
            shading flat
            a = colorbar;
            a.TickLabelInterpreter = 'latex';
            axis([min(min(XX_tr))+100e3 max(max(XX_tr)) dp 0])
            caxis([min(log10(lp)) 0.6])
            colormap(B,cmocean('rain'));
    %         xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
            xticks(0:50e3:400e3)
            xticklabels({'400' '350' '300' '250' '200' '150' '100' '50' '0'})
            title('Average Size (concentration weighted)','interpreter','latex')
            set(gca,'FontSize',24,'TickLabelInterpreter','latex')


            sumbool = avgVals > tol;
            C = subplot(1,3,2);
            pcolor(XX_tr,ZZ_tr,totVals);
            hold on
            [D h] = contour(XX_tr,ZZ_tr,sum(sumbool,3),0:4,'k');
            clabel(D,h,'Color','k')
            plot(xx_psi,-hb_psi,'k')
            hold off
            shading interp
            a = colorbar;
            a.TickLabelInterpreter = 'latex';
            axis([min(min(XX_tr))+100e3 max(max(XX_tr)) dp 0])
            colormap(C,cmap);
            caxis([0 max(max(totVals))])
            xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
            xticks(0:50e3:400e3)
            xticklabels({'400' '350' '300' '250' '200' '150' '100' '50' '0'})
            title('Number of Classes Present','interpreter','latex')
            set(gca,'FontSize',24,'TickLabelInterpreter','latex')

        elseif plot_type == idx_nums
    %         C = subplot(1,2,1);
            if plot_roem
                A = subplot(2,1,1);
                pcolor(XX_tr,ZZ_tr,pltVals)
                shading interp
                hold on
                if overlay > 0
                    [C h] = contour(XX_tr,ZZ_tr,buoy_ref,ocont,'k');  % buoyancy
                else
                    [C h] = contour(XX_tr,ZZ_tr,totVals,round(contvec,2),'k'); % plankton values
                end
                plot(xx_psi,-hb_psi,'k')
                hold off
                shading interp
                a = colorbar;
                a.TickLabelInterpreter = 'latex';
                axis([min(min(XX_tr))+100e3 max(max(XX_tr)) dp 0])
                colormap(cmap)
                clabel(C,h,'Color','k')
                caxis([0 cmax])
    %             xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
                xticks(0:50e3:400e3)
                ylabel('Depth (m)','interpreter','latex')
                xticklabels({'400' '350' '300' '250' '200' '150' '100' '50' '0'})  
                title(['Concentration (mmol N/m$^3$), $\Delta\ell$ = 0.4, $\tau_{max}$ = 0.1, t = ', num2str(round((n-n0)*dt_s./t1year,2)) ' years'],'interpreter','latex') 
                set(gca,'FontSize',28,'TickLabelInterpreter','latex')

                rnd = 2;
                B = subplot(2,1,2);
                PLTMAT = log10(PLTMAT);
                hold on
                for ii = 1:NP
                   plot(lxvec,(PLTMAT(:,ii)),'LineWidth',3,'Color',cmap(round(ii*255/NP),:))
                   
                end
                
                for ii = 1:NP
%                    plot(lxvec,PLTMAT(:,ii),'LineWidth',3,'Color',cmap(round(ii*255/NP),:))
%                    plot(lxvec,SSEMMAT(:,ii),'LineWidth',3,'Color',cmap(round(ii*255/NP),:))
%                     plot(lxvec,SSEMMAT(:,ii),'--','LineWidth',1,'Color',cmap(round(ii*255/NP),:))
                end
                
                hold off
                if bgc_id == idx_phyto
                    legstr = {};
                    for ii = 1:Ntracs
                        legstr{ii} = ['P',num2str(ii),': log$_{10}$($\ell_p$) = ',num2str(round(log10(lp(ii)),2))];
                    end
                    C = 'P';
                    leg = legend(legstr,'Location','southeast');
                    axis([min(min(XX_tr))+100e3 max(max(XX_tr)) -3 0.5])
%                     axis([min(min(XX_tr))+100e3 max(max(XX_tr)) -4 max(max((PLTMAT)))])
                else
                    legstr = {};
                    for ii = 1:Ntracs
                        legstr{ii} = ['Z',num2str(ii),': log$_{10}$($\ell_z$) = ',num2str(round(lz(ii),2))];
                    end
                    C = 'Z';
                    
                    leg = legend(legstr,'Location','northwest');
                    axis([min(min(XX_tr))+100e3 max(max(XX_tr)) -4 0.5])
%                     axis([min(min(XX_tr))+100e3 max(max(XX_tr)) -4 max(max((PLTMAT)))])
                end
                
                ylabel(['log$_{10}$(',C,') (log$_{10}$(mmol N/m$^3$))'],'interpreter','latex')
                set(gca,'FontSize',28,'TickLabelInterpreter','latex')
                xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
                xticks(0:50e3:400e3)
                xticklabels({'400' '350' '300' '250' '200' '150' '100' '50' '0'})   
                set(leg,'Interpreter','latex','FontSize',21)
%                 figure(101)
%                 plot(lxvec,log10(PLTMAT*lp./sum(PLTMAT,2)))

            else
                if bgc_id == idx_phyto
                    [LX,LP] = ndgrid(lxvec,log10(lp));
                else
                    [LX,LP] = ndgrid(lxvec,log10(lz));
                end
                pcolor(LX,LP,PLTMAT)
                hold on
                [C h] = contour(LX,LP,PLTMAT,contvec,'k');
                hold off
                clabel(C,h,'Color','k')
                shading interp
                a = colorbar;
                a.TickLabelInterpreter = 'latex';
                axis([min(min(XX_tr))+100e3 max(max(XX_tr)) min(min((LP))) max(max((LP)))])
                if bgc_id == idx_phyto
                    colormap(cmocean('speed',20));
                else
                    colormap(cmocean('dense',20));
                end
        %         caxis([0 max(max(max(PLTMAT)),0.35)])
                caxis([0 cmax])
                xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
                xticks(0:50e3:400e3)
                xticklabels({'400' '350' '300' '250' '200' '150' '100' '50' '0'})        
                ylabel('Size log10(ESD)','interpreter','latex')
                title(['Phytoplankton Distribution at ', num2str(H_ref), 'm depth'],'interpreter','latex')
                set(gca,'FontSize',24,'TickLabelInterpreter','latex')
%                 max(max(PLTMAT))
            end



        else
            if (bgc_id == idx_nitrate)
               contvec = [1 2 3 5 10 15 20]; 
            end
            pcolor(XX_tr,ZZ_tr,pltVals)
            shading interp
            hold on
            [C h] = contour(XX_tr,ZZ_tr,totVals,round(contvec,2),'k');
            plot(xx_psi,-hb_psi,'k')
            hold off
            shading interp
            a = colorbar;
            a.TickLabelInterpreter = 'latex';
            axis([min(min(XX_tr))+100e3 max(max(XX_tr)) dp 0])
            colormap(cmap)
            clabel(C,h,'Color','k')
            xlabel('Distance From Coast (km)','FontSize',fs,'interpreter','latex')
            xticks(0:50e3:400e3)
            ylabel('Depth (m)','interpreter','latex')
            xticklabels({'400' '350' '300' '250' '200' '150' '100' '50' '0'})  
            title('Total Phytoplankton Concentration (mmol N/m$^3$)','interpreter','latex') 
            set(gca,'FontSize',24,'TickLabelInterpreter','latex')
        end
        drawnow;

%         FF(n) = getframe;
%     end

%     writerObj = VideoWriter('phyto-conc.mp4');
%     writerObj.FrameRate = 10;
%     open(writerObj);
%     writeVideo(writerObj, FF)
%     close(writerObj);

end

