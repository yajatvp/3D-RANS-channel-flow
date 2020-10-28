%% FIGURES FOR RANS 3D MIXING LENGTH CODE
global x y z dx dy dz alx1 alx2 alx3 Re damping
%addpath ('D:\Halo data\2020-04-07 Fixed Scan analysis\Functions');
%%
yy=[0 y alx2];xx=[0 x alx1];zz=[0 z alx3];

% Mean streamwise velocity
figure;%set(gcf, 'Position', [1921 -215 1920 1.0048e+03])
set(gcf, 'Position', get(0, 'Screensize'));
plot(yy,[0 mean(mean(w(:,:,:,ntt-1),1),3) 0],'-o');
xlabel('Vertical coordinate');ylabel('Streamwise velocity - U(y)');grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',16)
title(['Mean streamwise velocity profile - U(y), Re = ',num2str(Re),', Damping = ',num2str(damping)]);

%%
% Mean total stress
figure;%set(gcf, 'Position', [1921 -215 1920 1.0048e+03])
set(gcf, 'Position', get(0, 'Screensize'));
plot(y,mean(mean(H_vw(:,:,:),1),3),'-o');
xlabel('Vertical coordinate');ylabel('Total shear stress - \tau_v_w(y)');grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',16);
title(['Mean total shear stress profile - \tau_v_w(y), Re = ',num2str(Re),', Damping = ',num2str(damping)]);

%%
figure;%set(gcf, 'Position', [1921 -215 1920 1.0048e+03])
set(gcf, 'Position', get(0, 'Screensize'));
plot(y,mean(mean(tau_vw(:,:,:),1),3),'-o');
xlabel('Vertical coordinate');ylabel('Mean turbulent shear stress - $\overline{vw}$');grid on;
title(['Mean turbulent shear stress profile - $\overline{vw}$, Re = ',num2str(Re),', Damping = ',num2str(damping)],'Interpreter','Latex');
set(findall(gcf,'-property','FontSize'),'FontSize',16);
%%
figure;%set(gcf, 'Position', [1921 -215 1920 1.0048e+03])
set(gcf, 'Position', get(0, 'Screensize'));
pcolor(z,y,squeeze(mean(p(:,:,:),1)));shading flat;axis normal;
colormap(gca,'hot');c = colorbar;c.Label.String = 'Total pressure';
xlabel('Streamwise coordinate');ylabel('Vertical coordinate');grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',16);
title(['Mean total pressure - P(y,z), Re = ',num2str(Re),', Damping = ',num2str(damping)]);

%%
figure;%set(gcf, 'Position', [1921 -215 1920 1.0048e+03])
set(gcf, 'Position', get(0, 'Screensize'));
contourf(z,y,squeeze(mean(w(:,:,:,ntt-1),1)),80,'edgecolor','none');shading flat;axis normal;
colormap(gca,'jet');c = colorbar;c.Label.String = 'Streamwise velocity';
xlabel('Streamwise coordinate');ylabel('Vertical coordinate');%grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',16);caxis([0 1]);
title(['Mean streamwise velocity - U(y,z), Re = ',num2str(Re),', Damping = ',num2str(damping)]);
%tau_uu=0;tau_uv=0;tau_uw=0;tau_vv=0;tau_vw=0;tau_ww=0;

