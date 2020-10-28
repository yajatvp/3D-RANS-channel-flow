clc
close
clear
clearvars -global
format long

% MECH 6370 - Final project 
% Group 3 - Benjamin, Federico, Santosh and Yajat
% 05/08/2020
% 3D - RANS solver using mixing length model; Channel flow case

% VARIABLEs NOMENCLATURE:
% 3 Velocity components: u, v, w ; size=[x,y,z,t] - 4D temporal array
% Pressure: p; size=[x,y,z] - 3D space, temporary w.r.t time
% alx1, alx2, alx3: Domain length in 3 directions
% n1, n2, n3 - No. of nodes in 3 directions
% H_uu - Total shear stress in uu
% Hu - Sum of advective and viscous terms according to equation 7.17 in the
% ...textbook - Ferziger Peric's CFD (2002), 3rd Ed.
% tau_uu - Turbulent shear stress - uu component
% damping - 1 if damping is considered, 0 if not

%% Inputs
global  dx dy dz n1 n2 n3 x y z nt dt alx2_2 Re

%%
% Variable initialization
Re=10^5; % Reynolds number input
tmax=10; % Maximum timestamp for the simulation
tol=10^-3; % Divergence tolerence
utol=10^-5; % Velocity tolerence
CFL=0.6; % CFL initialization
damping=0; % Damping or not

% Grid creation
alx1=1; % spanwise length
alx2=2; % wall normal - ALWAYS = 2, so that half height = 1
alx3=6; % streamwise
alx2_2=alx2/2;
n1=5;
n2=80;
n3=40;

dy=alx2/(n2); y(1:n2)=dy/2:dy:alx2-dy/2;
dx=alx1/(n1); x(1:n1)=dx/2:dx:alx1-dx/2;
dz=alx3/(n3); z(1:n3)=dz/2:dz:alx3-dz/2;

% Velocity, pressure initialization - streamwise velocity=1 and other
% velocity=0 at all points
p=ones(n1,n2,n3); % pressure
w_initial=zeros(n1,n2,n3); 
w_initial(:,:,:)=1;

dt(1)=CFL*dz/max(w_initial(:,:,:,1),[],'all'); % Time step size, CHECK CFL condition
nt=floor(tmax/dt(1));
tim(1)=0; % time array

u=zeros(n1,n2,n3,nt);
v=zeros(n1,n2,n3,nt);
w=zeros(n1,n2,n3,nt);
ke(:,:,:)=ones(n1,n2,n3); % TKE
w(:,:,:,1)=w_initial;

%% Time loop for NS equation discretization

for ntt=1:1:nt-1
    
    % Calculation of viscous stresses using 'VISCOUS' function
    [H_uu,H_uv,H_uw,H_vv,H_vw,H_ww]=viscous(u(:,:,:,ntt),v(:,:,:,ntt),w(:,:,:,ntt),Re);
    
    % Calculation of turbulent stresses using 'MIX_LEN' function
    [tau_uu,tau_uv,tau_uw,tau_vv,tau_vw,tau_ww]=mix_len(u(:,:,:,ntt),v(:,:,:,ntt),w(:,:,:,ntt),ke,damping);
    
    ke=0.5*(tau_uu+tau_vv+tau_ww);
    
    % total stress = viscous + turbulent
    H_uu=H_uu+tau_uu;
    H_uv=H_uv+tau_uv;
    H_uw=H_uw+tau_uw;
    H_vv=H_vv+tau_vv;
    H_vw=H_vw+tau_vw;
    H_ww=H_ww+tau_ww;
        
    % Input the total stress in 'NONLINEAR' function to calculate the H_i
    % components according to eqn. 7.17
    [Hu,Hv,Hw]=nonlinear(u(:,:,:,ntt),v(:,:,:,ntt),w(:,:,:,ntt),H_uu,H_uv,H_uw,H_vv,H_vw,H_ww);

    % Input the H_i components in 'SOLVE_PRESSURE' to get the pressure
    % matrix alongwith the gradients
    [p,px,py,pz]=solve_pressure(Hu,Hv,Hw);
    
    % Calculate the forcing in streamwise direction to maintain the
    % momemtum, using the integration of H_vw on both walls
    Fz=sum(sum(H_vw(:,1,:)))*alx1*alx3 + abs(sum(sum(H_vw(:,n2,:))))*alx1*alx3;

    % N-S equation discretization for next timestep
    u(:,:,:,ntt+1)=u(:,:,:,ntt)+dt(ntt)*(Hu(:,:,:)-px(:,:,:));
    v(:,:,:,ntt+1)=v(:,:,:,ntt)+dt(ntt)*(Hv(:,:,:)-py(:,:,:));
    w(:,:,:,ntt+1)=w(:,:,:,ntt)+dt(ntt)*(Hw(:,:,:)-pz(:,:,:)+(Fz/(n1*n2*n3)));
    
    % Calculation of dt at next timestep using CFL condition
    dt(ntt+1)=CFL*dz/max(w(:,:,:,ntt+1),[],'all');
        
    % Calculation of divergence using 'DIVERG' function
    duvw=diverg(u(:,:,:,ntt+1),v(:,:,:,ntt+1),w(:,:,:,ntt+1));
    
    if max(duvw,[],'all')>tol
        fprintf(['\n Divergence criteria not satisfied, max divergence= ', num2str(max(duvw,[],'all')),'\n']);
        break;
    end
    
    if max((w(:,:,:,ntt+1)-w(:,:,:,ntt)),[],'all')<utol
        fprintf(['\n Convergence satisfied \n']);
        break
    end
    tim(ntt+1)=tim(ntt)+dt(ntt);
    
    ntt
    Fz
    div1=max(duvw,[],'all')
    wmax=max(w(:,:,:,ntt+1),[],'all')

%     du=sum(sum((u(n1,:,:,ntt+1)-u(1,:,:,ntt+1))))*alx2*alx3;
%     dv=sum(sum((v(:,n2,:,ntt+1)-v(:,1,:,ntt+1))))*alx1*alx3;
%     dw=sum(sum((w(:,:,n3,ntt+1)-w(:,:,1,ntt+1))))*alx1*alx2;
%     div2=du+dv+dw
%    fprintf(['\n du=',num2str(du),'; dv=',num2str(dv),'; dw=',num2str(dw),'\n']);
    
end

