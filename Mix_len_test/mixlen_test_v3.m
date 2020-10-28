clc
clear all
close all

% MECH 6370 - Final project
% 05/08/2020
% Program to test the function mix_len
%% Input
global  dx dy dz n1 n2 n3 x y z alx2_2 Re

load ../q_database_tavg.mat
[n1,n2,n3]=size(w);

Re=4200; % Reynolds number of the validation DNS

x=xs;
y=ys_uniform';
z=zs;
dx=x(2)-x(1); %uniform grid
dy=y(2)-y(1); %uniform grid
dz=y(2)-y(1); %uniform grid

alx2=y(end)+0.5*dy;
alx2_2=alx2/2;


figure()
plot(w(100,:,100)/0.66,y,'-b')
hold on
xlabel('w')
ylabel('y')
set (gca,'FontName','Times','FontSize',16)
grid on
axis([0 1.5 0 2])

%%
% Damping or not - 1 or 0 respectively INPUT
damping=1;

%% Reynolds stresses calculation

ke(:,:,:)=ones(n1,n2,n3);
    
    [tau_uu,tau_uv,tau_uw,tau_vv,tau_vw,tau_ww]=mix_len(u,v,w,ke,damping);
   
    ke(:,:,:)=0.5*(tau_uu(:,:,:)+tau_vv(:,:,:)+tau_ww(:,:,:));

%% Post processing
tau_vw_avg=zeros(n2,1);
for j=1:n2
    for i=1:n1
        for k=1:n3
            tau_vw_avg(j)=tau_vw_avg(j)+tau_vw(i,j,k);
        end
    end
    tau_vw_avg(j)=tau_vw_avg(j)/(n1*n3);
end

figure()
plot(tau_vw_avg,y,'-k');
grid on
hold on
xlabel('{\tau_{vw}}')
ylabel('y')
set(gca,'FontName','Times','FontSize',16)

load ../Re_stress_tavg_DNS.mat


tau_vw_reDNSavg=zeros(n2,1);
for j=1:n2
    for i=1:n1
        for k=1:n3
            tau_vw_reDNSavg(j)=tau_vw_reDNSavg(j)-tau_vws(i,j,k);
        end
    end
    tau_vw_reDNSavg(j)=tau_vw_reDNSavg(j)/(n1*n3);
end

plot(tau_vw_reDNSavg,y,'-r')

% print('Mix_len_damp','-dpng')
% print('Mix_len_nodamp','-dpng')