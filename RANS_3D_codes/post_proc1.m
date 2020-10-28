%%
u0=mean(mean(w(:,:,:,ntt-1),1),3);
str0=mean(mean(H_vw(:,:,:),1),3);
tstr0=mean(mean(tau_vw(:,:,:),1),3);
%%
u1=mean(mean(w(:,:,:,ntt-1),1),3);
str1=mean(mean(H_vw(:,:,:),1),3);
tstr1=mean(mean(tau_vw(:,:,:),1),3);
yy=[0 y alx2];xx=[0 x alx1];zz=[0 z alx3];
%% Mean streamwise velocity
figure;%set(gcf, 'Position', [1921 -215 1920 1.0048e+03])
set(gcf, 'Position', get(0, 'Screensize'));
plot(yy,[0 u0 0],'-o');hold on;
plot(yy,[0 u1 0],'-o');
xlabel('Vertical coordinate');ylabel('Streamwise velocity - U(y)');grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',16)
title(['Mean streamwise velocity profile - U(y), Re = ',num2str(Re)]);
legend('Without damping','Damping');

%% Mean total stress
figure;%set(gcf, 'Position', [1921 -215 1920 1.0048e+03])
set(gcf, 'Position', get(0, 'Screensize'));
plot(y,str0-tstr0,'-o'); hold on;
plot(y,str1-tstr1,'-o');
xlabel('Vertical coordinate');ylabel('Viscous shear stress - \tau_v_w(y)');grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',16);
title(['Mean viscous shear stress profile - \tau_v_w(y), Re = ',num2str(Re)]);
legend('Without damping','Damping');

%% Mean turbulent stress
figure;%set(gcf, 'Position', [1921 -215 1920 1.0048e+03])
set(gcf, 'Position', get(0, 'Screensize'));
plot(y,tstr0,'-o');hold on;
plot(y,tstr1,'-o');
xlabel('Vertical coordinate');ylabel('Mean turbulent shear stress - $\overline{vw}$','Interpreter','Latex');grid on;
title(['Mean turbulent shear stress profile - $\overline{vw}$, Re = ',num2str(Re)],'Interpreter','Latex');
set(findall(gcf,'-property','FontSize'),'FontSize',16);%%
legend('Without damping','Damping');
