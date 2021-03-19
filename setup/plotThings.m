%
%
% plotThings.m
%
%
%
%

addpath ../utils
addpath ../analysis

Nx = 64;
Nz = Nx;
path = '~/Desktop/MAMEBUS_Data/pg_plts';
addpath(path)
LIN = 'DB_DX_LINEAR_n=0.dat';
CUB = 'DB_DX_CUBIC_n=0.dat';
tracers = true;
final_val = 7300;
tracerval = 1;

gridpath = '~/Desktop/MAMEBUS_Data/resolution-tests';

% [XX_tr,ZZ_tr,XX_psi,ZZ_psi,~,hb_psi,xx_psi] = plotSolution (gridpath,'test64c',tracers,tracerval,0.01);

dbdxl = readOutputFile(LIN,Nx+1,Nz+1);
dbdxc = readOutputFile(CUB,Nx+1,Nz+1);

H = 3000;
Hs = 50;
m1km = 1000; %%% Meters in 1 km    
Lx = 400*m1km; %%% Computational domain width

%%% Initial temperature
Hexp = 150; 
Tmax = 22 - 5*XX_tr/Lx;
Tmin = 4;
buoy_init = Tmin + (Tmax-Tmin).*(exp(ZZ_tr/Hexp+1)-exp(-H/Hexp+1))./(exp(1)-exp(-H/Hexp+1));  

Bx = -(5/Lx).*(exp(ZZ_psi/Hexp+1)-exp(-H/Hexp+1))./(exp(1)-exp(-H/Hexp+1));

figure(101)
A = subplot(1,3,1);
val = abs(max(max(abs(Bx))));
cmap = cmocean('ice');
% Initial buoyancy
bmap = cmocean('thermal');
% [C, h] = contourf(XX_tr,ZZ_tr,Bx,0:2:20);
% clabel(C,h)
pcolor(XX_psi,ZZ_psi,Bx);
shading interp
colormap(A,cmap)
title('Analytical Buoyancy Gradient')
set(gca,'FontSize',22)
xticks([0 1 2 3 4]*1e5)
xticklabels({'400', '300', '200', '100', '0'})
yticks(-3000:200:0)
ylabel('Depth (m)')
axis([min(min(XX_tr)) max(max(XX_tr)) -1000 0])
caxis([-val 0])

% val = abs(max(max(abs(dbdxc))));
cmap = cmocean('ice');
B = subplot(1,3,2);
pcolor(XX_psi,ZZ_psi,dbdxl)
shading interp
caxis([-val 0])
title('Linear Buoyancy Gradient')
colormap(B,cmap)
set(gca,'FontSize',22)
xticks([0 1 2 3 4]*1e5)
xticklabels({'', '300', '200', '100', ''})
xlabel('Distance From the Coast (km)')
yticks(-3000:200:0)
yticklabels({})
axis([min(min(XX_tr)) max(max(XX_tr)) -1000 0])

C = subplot(1,3,3);
pcolor(XX_psi,ZZ_psi,dbdxc)
shading interp
colorbar
title('Cubic Buoyancy Gradient')
caxis([-val 0])
colormap(C,cmap)
set(gca,'FontSize',22)
xticks([0 1 2 3 4]*1e5)
xticklabels({'', '300', '200', '100', ''})
yticks(-3000:600:0)
yticklabels({})
axis([min(min(XX_tr)) max(max(XX_tr)) -1000 0])


figure(102)
D = subplot(1,2,1);
val = max(max(abs(dbdxl(2:end-1,2:end-1) - Bx(2:end-1,2:end-1))));
ermap = cmocean('balance');
pcolor(XX_psi(2:end-1,2:end-1),ZZ_psi(2:end-1,2:end-1),dbdxl(2:end-1,2:end-1) - Bx(2:end-1,2:end-1))
shading interp
caxis([-val 0])
title('Linear Error (Linear - Actual)')
colormap(D,ermap)
set(gca,'FontSize',22)
xticks([0 1 2 3 4]*1e5)
xticklabels({'', '300', '200', '100', ''})
xlabel('Distance From the Coast (km)')
yticks(-3000:200:0)
% yticklabels({})
axis([min(min(XX_tr)) max(max(XX_tr)) -1000 0])
caxis([-val val])
ylabel('Depth (m)')

E = subplot(1,2,2);
% ermap = cmocean('delta');
pcolor(XX_psi(2:end-1,2:end-1),ZZ_psi(2:end-1,2:end-1),dbdxc(2:end-1,2:end-1) - Bx(2:end-1,2:end-1))
shading interp
caxis([-val 0])
title('Cubic Error (Cubic - Actual)')
colormap(E,ermap)
set(gca,'FontSize',22)
xticks([0 1 2 3 4]*1e5)
xticklabels({'', '300', '200', '100', ''})
xlabel('Distance From the Coast (km)')
yticks(-3000:200:0)
yticklabels({})
axis([min(min(XX_tr(2:end-1,2:end-1))) max(max(XX_tr(2:end-1,2:end-1))) -1000 0])
colorbar
caxis([-val val])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% UNCOMMENT FOR PRESSURE GRADIENT TESTS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% tracers = false;
% final_val = 7300;
% varname = 'PSIM';
% if tracers > 0
%     trac_name = [varname,num2str(tracer),'_n=',num2str(final_val),'.dat'];
% else
%     trac_name = [varname,'_n=',num2str(final_val),'.dat'];
% end
% 
% path = '~/Desktop/MAMEBUS_Data/pressure_gradients';
% tracerval = 1;
% years = 10;
% % load in various pressure gradient comparisons
% % [XX_tr,ZZ_tr,XX_psi,ZZ_psi,oneHundreth_cm,~,~] = plotSolution (path,'baro_1e-9',tracers,tracerval,years);
% % [~,~,~,~,moneHundreth_cm,~,~] = plotSolution (path,'baro_m1e-9',tracers,tracerval,years);
% % [~,~,~,~,oneTenth_cm,~,~] = plotSolution (path,'onecm_ssh',tracers,tracerval,years);
% % [~,~,~,~,moneTenth_cm,~,~] = plotSolution (path,'minusonecm_ssh',tracers,tracerval,years);
% % [~,~,~,~,onecm,~,~] = plotSolution (path,'tencm_ssh',tracers,tracerval,years);
% % [~,~,~,~,monecm,~,~] = plotSolution (path,'minus10cm_ssh',tracers,tracerval,years);
% % 
% % [~,~,~,~,ref,hb_psi,xx_psi] = plotSolution (gridpath,'test64c',tracers,tracerval,years);
% 
% % [~,~,~,~,~,hb_psi,xx_psi] = plotSolution (gridpath,'test64c',tracers,tracerval,0.01);
% 
% for ii = 1:10
% figure(ii)
% close
% end
% 
% xplt_psi = [];
% hbplt_psi = [];
% 
% for ii = 1:length(xx_psi) - 1
%     xplt_psi = [xplt_psi, [xx_psi(ii) xx_psi(ii) xx_psi(ii+1) xx_psi(ii+1)] ];
%     hbplt_psi = [hbplt_psi, [-3000 -hb_psi(ii) -hb_psi(ii+1) -3000] ];
% end
% 
% ncol = 2;
% nrow = 2;
% fig = 1;
% pltval = 1;
% 
% if tracers 
%     XX = XX_tr;
%     ZZ = ZZ_tr;
% else
%     XX = XX_psi;
%     ZZ = ZZ_psi;
% end
% 
% mval = 1;
% pltcolor = false;
% ncont = -mval:0.02:mval;
% 
% % close all
% cmap = cmocean('balance');
% 
% figure(100);
% clf;
% % smallest pressure gradients
% subplot(nrow,ncol,fig); fig = fig +1;
% if pltcolor
%     surf(XX,ZZ,moneHundreth_cm)
%     view(2)
% else
%     contourf(XX,ZZ,moneHundreth_cm,ncont)
% end
% shading interp
% caxis([-pltval pltval])
% xticks([0 1 2 3 4]*1e5)
% yticks(-3000:600:0)
% ylabel('Depth (m)')
% xticklabels({})
% set(gca,'FontSize',22)
% colormap(cmap)
% hold on
% fill(xplt_psi,hbplt_psi,[192 192 192]./255,'LineStyle','none');
% hold off
% %
% 
% % subplot(nrow,ncol,fig); fig = fig +1;
% % if pltcolor
% %     pcolor(XX,ZZ,ref)
% % else
% %     contourf(XX,ZZ,ref,ncont)
% % end
% % shading interp
% % xticklabels({})
% % xticks([0 1 2 3 4]*1e5)
% % yticks(-3000:600:0)
% % yticklabels({})
% % caxis([-pltval pltval])
% % set(gca,'FontSize',22)
% % colormap(cmap)
% % hold on
% % fill(xplt_psi,hbplt_psi,[192 192 192]./255,'LineStyle','none');
% % hold off
% %
% 
% subplot(nrow,ncol,fig); fig = fig +1;
% if pltcolor
%     pcolor(XX,ZZ,oneHundreth_cm)
% else
%     contourf(XX,ZZ,oneHundreth_cm,ncont)
% end
% shading interp
% caxis([-pltval pltval])
% colorbar
% xticks([0 1 2 3 4]*1e5)
% xticklabels({})
% yticks(-3000:600:0)
% yticklabels({})
% set(gca,'FontSize',22)
% colormap(cmap)
% hold on
% fill(xplt_psi,hbplt_psi,[192 192 192]./255,'LineStyle','none');
% hold off
% %
% 
% % intermediate pressure gradients
% subplot(nrow,ncol,fig); fig = fig +1;
% if pltcolor
%     pcolor(XX,ZZ,moneTenth_cm)
% else
%     contourf(XX,ZZ,moneTenth_cm,ncont)
% end
% shading interp
% caxis([-pltval pltval])
% xticks([0 1 2 3 4]*1e5)
% yticks(-3000:600:0)
% % yticklabels({})
% ylabel('Depth (m)')
% xticklabels({})
% colormap(cmap)
% set(gca,'FontSize',22)
% hold on
% fill(xplt_psi,hbplt_psi,[192 192 192]./255,'LineStyle','none');
% hold off
% xticks([0 1 2 3 4]*1e5)
% xticklabels({'400', '300', '200', '100', '0'})
% xlabel('Distance From the Coast (km)')
% %
% 
% % subplot(nrow,ncol,fig); fig = fig +1;
% % if pltcolor
% %     pcolor(XX,ZZ,ref)
% % else
% %     contourf(XX,ZZ,ref,ncont)
% % end
% % shading interp
% % caxis([-pltval pltval])
% % xticks([0 1 2 3 4]*1e5)
% % xticklabels({})
% % yticks(-3000:600:0)
% % yticklabels({})
% % set(gca,'FontSize',22)
% % colormap(cmap)
% % hold on
% % fill(xplt_psi,hbplt_psi,[192 192 192]./255,'LineStyle','none');
% % hold off
% % xticks([0 1 2 3 4]*1e5)
% % xticklabels({'400', '300', '200', '100', '0'})
% % xlabel('Distance From the Coast (km)')
% %
% 
% subplot(nrow,ncol,fig); fig = fig +1;
% if pltcolor
%     pcolor(XX,ZZ,oneTenth_cm)
% else
%     contourf(XX,ZZ,oneTenth_cm,ncont)
% end
% shading interp
% caxis([-pltval pltval])
% colorbar
% xticks([0 1 2 3 4]*1e5)
% xticklabels({})
% yticks(-3000:600:0)
% yticklabels({})
% set(gca,'FontSize',22)
% colormap(cmap)
% hold on
% fill(xplt_psi,hbplt_psi,[192 192 192]./255,'LineStyle','none');
% hold off
% xticks([0 1 2 3 4]*1e5)
% xticklabels({'400', '300', '200', '100', '0'})
% xlabel('Distance From the Coast (km)')
% %
% 
% %%% Add in text box
% h1 = annotation('textbox',[0.25 0.925 0.2 0.05],...
%     'string','Northward PGF');
% h2 = annotation('textbox',[0.48 0.925 0.2 0.05],...
%     'string',' ');
% h3 = annotation('textbox',[0.70 0.925 0.2 0.05],...
%     'string','Southward PGF');
% set([h1 h2 h3], 'fitboxtotext','on',...
%     'edgecolor','none','FontSize',24)
% 
% h4 = annotation('textbox',[0 0.72 0.2 0.05],...
%     'string','$$\Delta \eta$$ = .1 cm');
% h5 = annotation('textbox',[0 0.26 0.2 0.05],...
%     'string','$$\Delta \eta$$ = 1 cm');
% set([h4 h5], 'fitboxtotext','on',...
%     'edgecolor','none','FontSize',24,'interpreter','latex')

% largest pressure gradients
% subplot(nrow,ncol,fig); fig = fig +1;
% if pltcolor
%     pcolor(XX,ZZ,monecm)
% else
%     contourf(XX,ZZ,monecm,ncont)
% end
% shading interp
% caxis([-pltval pltval])
% xticks([0 1 2 3 4]*1e5)
% yticks(-3000:600:0)
% % yticklabels({})
% ylabel('Depth (m)')
% xticklabels({'400', '300', '200', '100', '0'})
% set(gca,'FontSize',22)
% colormap(cmap)
% 
% subplot(nrow,ncol,fig); fig = fig +1;
% if pltcolor
%     pcolor(XX,ZZ,ref)
% else
%     contourf(XX,ZZ,ref,ncont)
% end
% shading interp
% caxis([-pltval pltval])
% xticks([0 1 2 3 4]*1e5)
% yticks(-3000:600:0)
% yticklabels({})
% xticklabels({'400', '300', '200', '100', '0'})
% xlabel('Distance From the Coast (km)')
% set(gca,'FontSize',22)
% colormap(cmap)
% 
% subplot(nrow,ncol,fig); fig = fig +1;
% if pltcolor
%     pcolor(XX,ZZ,onecm)
% else
%     contourf(XX,ZZ,onecm,ncont)
% end
% shading interp
% caxis([-pltval pltval])
% yticks(-3000:600:0)
% yticklabels({})
% colorbar
% xticks([0 1 2 3 4]*1e5)
% xticklabels({'400', '300', '200', '100', '0'})
% set(gca,'FontSize',22)
% colormap(cmap)
