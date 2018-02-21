%%%
%%% Reads in data from the output of 'Overturning' and makes a movie 
%%% of the velocity field
%%%
%%% local_home_dir specifies the directory in the local system in which
%%% run files are stored. N.B. this function will search within a 
%%% subdirectory called 'run_name' to find the run's output files.
%%%
%%% run_name specifies the name of the run.
%%%
%%% This program is used to analyze the data from MAMEBUS runs.
%%%
function M = analyzeData(local_home_dir,run_name,var_id,velocity)
    %%% Calcuate the velocity field from the residual streamfunction
    %%% Load convenience functions
    addpath ../utils;
    addpath ./redblue

    %%%%%%%%%%%%%%%%%%%%%
    %%%%% VARIABLES %%%%%
    %%%%%%%%%%%%%%%%%%%%% 

    %%% Parameter and data file names
    run_name = strtrim(run_name);
    dirpath = fullfile(local_home_dir,run_name);
    params_file = fullfile(dirpath,[run_name,'_in']);  

    %%% Plotting grid
    [Nx Nx_found] = readparam(params_file,'Nx','%u');
    [Nz Nz_found] = readparam(params_file,'Nz','%u');
    [Lx Lx_found] = readparam(params_file,'Lx','%lf');
    [H H_found] = readparam(params_file,'Lz','%lf');
    [tlength tlength_found] = readparam(params_file,'tlength','%u');
    
    if ((~Nx_found) || (~Nz_found) || (~Lx_found) || (~H_found))
        error('Could not read grid parameters');
    end    
  
    %%% Read grid parameters
    [h_c h_c_found] = readparam(params_file,'h_c','%le');
    [theta_s theta_s_found] = readparam(params_file,'theta_s','%lf');
    [theta_b theta_b_found] = readparam(params_file,'theta_b','%lf');
  
    %%% Read bottom topography
    hb = readDataFile (params_file,dirpath,'topogFile',Nx+2,1,H*ones(Nx+2,1));
    hb_psi = 0.5*(hb(1:end-1)+hb(2:end));  
    hb_tr = hb(2:end-1); %%% Remove "ghost" points
    
    %%% Read wind profile
    tau = readDataFile (params_file,dirpath,'tauFile',Nx,tlength,ones(Nx,Nz));
    tau_max = abs(min(min(tau)));
  
    %%% Parameters related to number of iterations
    dt_s = readparam(params_file,'monitorFrequency','%lf');
    tmax = readparam(params_file,'maxTime','%lf');
  
    %%% For convenience
    t1year = 365*86400; %%% Seconds in one year
    t1day = 86400;
    
    %%% Generate full sigma-coordinate grids
    [XX_tr,ZZ_tr,XX_psi,ZZ_psi,XX_u,ZZ_u,XX_w,ZZ_w] ...
                        = genGrids(Nx,Nz,Lx,h_c,theta_s,theta_b,hb_tr,hb_psi);
                    

    %%% Allocate vectors to hold the data
    sx = Nz;
    dx = zeros(Nx,Nz);
    dz = dx;
    u = dx;
    w = dz;
    
    dx = XX_psi(2:Nx+1,:) - XX_psi(1:Nx,:);
    dz = ZZ_psi(:,2:Nz+1) - ZZ_psi(:,1:Nz);
  
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% PLOTTING LOOP %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%% Max number of time steps is the number of whole time steps that can
    %%% fit into tmax, plus one initial save, plus one final save
    Nt = floor(tmax/dt_s) + 2;
  
    %%% Make a movie of the data - allocate a movie array large enough to
    %%% hold the largest possible number of frames
%     figure(1);
%     clf;
%     axes('FontSize',18);
%     M = moviein(Nt);  
%   
%     %%% Tracks whether we should still read data
%     stillReading = true;
%     counter = 1;
%     n = 0;
  
    % Determines whether or not to make a movie and writes a new file to
    % store the data.
%     if mov_on
%     	vidObj = VideoWriter(mov_name);
%         vidObj.FrameRate = 10;
%         open(vidObj)
%     end
    
     %%% At each time iteration...
%     while (stillReading)
%       
%       %%% Get the time value 
%       t = n*dt_s;
%     
%       %%% Load different streamfunctions      
%       switch (var_id)
%         case 0 %%% Residual streamfunction
%           data_file = fullfile(dirpath,['PSIR_n=',num2str(n),'.dat']);
%         case 1 %%% Mean streamfunction
%           data_file = fullfile(dirpath,['PSIM_n=',num2str(n),'.dat']);
%         case 2 %%% Eddy streamfunction
%           data_file = fullfile(dirpath,['PSIE_n=',num2str(n),'.dat']);
%       end
% 
%       %%% Open the output file for reading    
%       dfid = fopen(data_file,'r');
%       if (dfid == -1)
%         stillReading = false;
%         continue;
%         %%% Ignore any missing data - might just be the end of the
%         %%% computation
%       end
% 
%       %%% Get the phi values on the gridpoints
%       psi = readOutputFile (data_file,Nx+1,Nz+1); 
%       if (size(psi,1)~=Nx+1 || size(psi,2)~=Nz+1)
%         error(['ERROR: Could not find data file: ', data_file]);
%       end    
%       
%       %%% Close data file
%       fclose(dfid);
%       
%       %%% Calculate the zonal and vertial velocities with respect to the
%       %%% stream function.
%       for ii = 1:sx(2)
%         u(:,ii) = (psi(2:end,ii) - psi(1:end-1,ii))./dx(:,ii);
%         w(ii,:) = (psi(ii,2:end) - psi(ii,1:end-1))./dz(ii,:);
%       end
%       
%       %%% Determine which velocity to plot
%       if strcmp(velocity,'zonal')
%         %%% Plot the zonal velocity
%         psi_r_lim = psi;
%         limval = 10e-4;
%         psi_r_lim = min(psi_r_lim,limval);
%         psi_r_lim = max(psi_r_lim,-limval);
% %       [C h] = contourf(XX_psi,ZZ_psi,psi_r_lim,-limval:limval/40:limval,'EdgeColor','k');  
%         figure(1)
%         pcolor(XX_tr,ZZ_tr,u);
%         shading interp;     
%         colormap redblue;
%         h=colorbar;        
%         caxis([-limval limval]);
%         set(h,'FontSize',18);
%         hold on
%         plot(XX_tr,-hb_tr,'k','LineWidth',2)
%         hold off
%       
%         %%% Store the image in the movie buffer  
%         xlabel('x');    
%         ylabel('z','Rotation',0);        
%         axis tight;
%         set(gca,'XTick',0:Lx/5:Lx);
%         set(gca,'YTick',-H:H/5:0);
%         title(strcat(['t=',num2str(round(t/t1day)),' days (',num2str(round(t/t1year)),' yr)']));   
%       else
%         %%% Plot the vertical velocity otherwise
%         psi_r_lim = psi;
%         limval = 0.25;
%         psi_r_lim = min(psi_r_lim,-limval);
%         psi_r_lim = max(psi_r_lim,limval);
% %       [C h] = contourf(XX_psi,ZZ_psi,psi_r_lim,-limval:limval/40:limval,'EdgeColor','k');  
%         figure(1)
%         pcolor(XX_tr,ZZ_tr,w);
%         shading interp;     
%         colormap redblue;
%         h=colorbar;        
%         caxis([-limval limval]);
%         set(h,'FontSize',18);
%         hold on
%         plot(XX_tr,-hb_tr,'k','LineWidth',1)
%         hold off
%       
%         %%% Store the image in the movie buffer  
%         xlabel('x');    
%         ylabel('z','Rotation',0);        
%         axis tight;
%         set(gca,'XTick',0:Lx/5:Lx);
%         set(gca,'YTick',-H:H/5:0);
%         title(strcat(['t=',num2str(round(t/t1day)),' days (',num2str(round(t/t1year)),' yr)']));  
%       end
%     
%       nextframe = getframe(gcf);    
%       M(counter) = nextframe; 
%         
%       counter = counter+1;   
%       n = n + 1;
%       
% %     if mov_on
% %         close(vidObj)
% %     end
%     end
    
    
n = 500;
[H H_found] = readparam(params_file,'Lz','%lf');
t = n*dt_s;

  set(0,'DefaultAxesFontSize',16)
    %%% Find the average and max source depth for the shelf water
    depth_file = fullfile(dirpath,['TRAC',num2str(1),'_n=',num2str(n),'.dat']);
    depth = readOutputFile(depth_file,Nx,Nz);  
    loc = 36;
    numcols = Nx - loc;
    
    col_avg = sum(depth(loc:Nz,:))/(numcols+1);
    fignum = 1;
    figure(fignum)
    contourf(XX_tr,ZZ_tr,depth,15)
    title(['Depth Tracer, t = ' num2str(t/(365*24*60*60)) ' years' ])
    xlabel('Distance (km)')
    ylabel('Depth')
    colorbar
    colormap default;
    caxis([-H 0])
    
    nit_file = fullfile(dirpath,['TRAC',num2str(2),'_n=',num2str(n),'.dat']);
    nitrate = readOutputFile(nit_file,Nx,Nz);

    
    n_colavg = sum(nitrate(loc:Nz,:))/(numcols+1);
    
    surf_d = 0;
    it = 0;
    
    for ii = loc:Nz
        for jj = 1:Nz
            if ZZ_tr(ii,jj) > -50
                it = it + 1;
                surf_d = surf_d + depth(ii,jj);
            end
        end
    end
    
    
    surf_n = nitrate((ZZ_tr > - 50));
    avg_surf = sum(surf_n)/length(surf_n);
    
    fignum = fignum+1;
    figure(fignum)
    pcolor(XX_tr,ZZ_tr,nitrate)
    shading interp
    colorbar
    
    %%% Values to display
    disp('Note: Averages are taken on shelf')
    disp(['Average source depth: ' num2str(sum(col_avg)/Nz)])
    disp(['Average surface source depth: ' num2str(surf_d/it)])
    disp(['Average nitrate concentration: ' num2str(sum(n_colavg)/Nz)])
    disp(['Surface nitrate concentration: ' num2str(avg_surf)])
    disp(['Wind stress maximum strength: ' num2str(tau_max)])
    
    data_file = fullfile(dirpath,['PSIR_n=',num2str(n),'.dat']);
    psi = readOutputFile (data_file,Nx+1,Nz+1); 
    
    buoy_file = fullfile(dirpath,['TRAC',num2str(0),'_n=',num2str(n),'.dat']);
    buoy = readOutputFile(buoy_file,Nx,Nz);  
    
    
    %%% Calculate Uptake on Shelf
    N_shelf = nitrate(loc:Nx,:);
    buoy_shelf = buoy(loc:Nx,:);
    a_temp = 0.6;
    b_temp = 1.066;
    c_temp = 1;
    alpha = 0.025; % s^-1 (W/m^2)^-1
    monod = 0;

    kw = 0.04;
    kc = 0.03;
    efold = 30;
    I0 = 340;
    K = kw + kc.*N_shelf;

    IR = I0*exp(K.*ZZ_tr(loc:Nx,:)/efold);
    t_uptake = a_temp.*b_temp.^(c_temp.*buoy_shelf);
    lk = t_uptake./alpha;
    l_uptake = IR./sqrt(IR.^2 + lk.^2);
              
              
    uptake = (l_uptake.*t_uptake).*N_shelf;
    shelf_uptake = sum(sum(uptake)/(numcols+1))/Nz;
    
    disp(['Average Uptake on the Shelf: ' num2str(shelf_uptake)])
    
    %%% Calculate uptake 100km offshore above 300m
    N_offshelf = zeros(size(nitrate(27:30,:)));
    buoy_offshelf = N_offshelf;
    Zgrid = buoy_offshelf;
    for ii = 1:Nz
        if ZZ_tr(27:30,ii) > -150
            N_offshelf(27:30,ii) = nitrate(27:30,ii);
            buoy_offshelf(27:30,ii) = buoy(27:30,ii);
            Zgrid(27:30,ii) = ZZ_tr(27:30,ii);
        end
    end
%     N_offshelf = nitrate(27:30,:);
%     buoy_offshelf = buoy(27:30,:);
    K = kw + kc.*N_offshelf;

    IR = I0*exp(K.*Zgrid/efold);
    t_uptake = a_temp.*b_temp.^(c_temp.*buoy_offshelf);
    lk = t_uptake./alpha;
    l_uptake = IR./sqrt(IR.^2 + lk.^2);
    
    uptake = (l_uptake.*t_uptake).*N_offshelf;
    offshelf_uptake = sum(sum(uptake)/(numcols+1))/Nz;
    
    disp(['Average Uptake 100km off Shelf above 300m: ' num2str(offshelf_uptake)])
    
    %%% Plotting for uptake
%     N_tot = sum(nitrate,2)/40;
%     K = -ones(size(nitrate))*kw; % + kc.*N_tot;
% 
%     opc = zeros(Nz,1);
%     opc = log(0.01)*efold./K
    
    K = kw + kc*nitrate;
    
    IR = I0*exp(K.*ZZ_tr/efold);
%     
    opc = -115*ones(Nx,1);
%     %%% Find 1 percent light level
%     for ii = 1:Nx
%         for jj = Nz:-1:1
%             if (IR(ii,jj) <= 3.4)
%                 opc(ii) = ZZ_tr(ii,jj+1);
%                 break
%             end
%         end
%     end
    
    t_uptake = a_temp.*b_temp.^(c_temp.*buoy);
    lk = t_uptake./alpha;
    l_uptake = IR./sqrt(IR.^2 + lk.^2);
    
    uptake = (l_uptake.*t_uptake).*nitrate;
    
    
    
    fignum = fignum+1;
    figure(fignum);
    
    pcolor(XX_tr,ZZ_tr,uptake)
    hold on
    plot(XX_tr,opc,'w','LineWidth',2)
    hold off
    shading interp
    colorbar
    colormap jet;
    xlabel('Distance (km)')
    ylabel('Depth (m)')
    title('Uptake (NPP) mmol N/(m^3 d)')
    
    
              
    
    limval = 2;
    
    fignum = fignum+1;
    figure(fignum);
%     subplot(1,2,2)
    colorDepth = 1000;
    colormap(redblue(colorDepth));
    hold on
    pcolor(XX_psi,ZZ_psi,psi)
    shading interp;
    [C,hfigc] = contour(XX_psi,ZZ_psi,psi,(-2:0.5:2));
    set(hfigc, ...
        'LineWidth',1.0, ...
        'Color', [0.2 0.2 0.2]);
    hold off;
    caxis([-limval limval]);
    xlabel('Distance (km)')
    ylabel('Depth (m)')
    title('Residual Streamfunction, $\psi^{\dagger}$','interpreter','latex')
    
%     subplot(1,2,1)
    fignum = fignum +1;
    figure(fignum)
    contourf(XX_tr,ZZ_tr,buoy,0:1:20)
    colorbar
    colormap default
    xlabel('Distance (km)')
    ylabel('Depth (m)')
    title('Buoyancy ($^o$C)','interpreter','latex')
%     colormap;
%     caxis([0 20])


    %%% Calculate EKE
    data_file = fullfile(dirpath,['PSIE_n=',num2str(n),'.dat']);
    psie = readOutputFile (data_file,Nx+1,Nz+1); 
    ue = (psie(2:Nx+1,:) - psie(1:Nx,:))./dx;
    eke = ue.^2*(100)*t1day;
    eke = eke./2;
    
%     fignum = fignum+1;
%     figure(fignum)
%     pcolor(XX_tr,ZZ_tr,eke(:,1:Nz))
%     shading interp
%     colorbar
%     colormap jet;
%     xlabel('Distance (km)')
%     ylabel('Depth (m)')
%     title('EKE m^2/s^2')
    
end