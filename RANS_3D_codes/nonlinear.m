function [Hu,Hv,Hw]=nonlinear(u,v,w,tau_uu,tau_uv,tau_uw,tau_vv,tau_vw,tau_ww)
global dx dy dz n1 n2 n3 x y z

%% Hu
[t1,~,~]=gradient1(tau_uu);
[~,t2,~]=gradient1(tau_uv);
[~,~,t3]=gradient1(tau_uw);
% add BCs here
t1(1,:,:)=(tau_uu(2,:,:)-tau_uu(n1,:,:))/(2*dx);
t1(n1,:,:)=(tau_uu(1,:,:)-tau_uu(n1-1,:,:))/(2*dx);
t2(:,1,:)=(tau_uv(:,2,:)-tau_uv(:,1,:))/(dy);
t2(:,n2,:)=(tau_uv(:,n2,:)-tau_uv(:,n2-1,:))/(dy);
t3(:,:,1)=(tau_uw(:,:,2)-tau_uw(:,:,n3))/(2*dz);
t3(:,:,n3)=(tau_uw(:,:,1)-tau_uw(:,:,n3-1))/(2*dz);

t=t1+t2+t3;clear t1 t2 t3
% nonlinear terms
t1(2:n1-1,:,:)=((u(3:n1,:,:).^2-u(1:n1-2,:,:).^2)/(2*dx)); % central, all
t1(n1,:,:)=((u(1,:,:).^2-u(n1-1,:,:).^2)/(2*dx));
t1(1,:,:)=((u(2,:,:).^2-u(n1,:,:).^2)/(2*dx));

t2(:,2:n2-1,:)=(((u(:,3:n2,:).*v(:,3:n2,:))-(u(:,1:n2-2,:).*v(:,1:n2-2,:)))/(2*dy)); % central
t2(:,n2,:)=(0-((u(:,n2,:).*v(:,n2,:)))/(dy*0.5)); % forward
t2(:,1,:)=(((u(:,1,:).*v(:,1,:))-0)/(dy*0.5)); % backward

t3(:,:,2:n3-1)=(((u(:,:,3:n3).*w(:,:,3:n3))-(u(:,:,2:n3-1).*w(:,:,2:n3-1)))/(2*dz)); % central
t3(:,:,1)=(((u(:,:,2).*w(:,:,2))-(u(:,:,n3).*w(:,:,n3)))/(2*dz)); 
t3(:,:,n3)=(((u(:,:,1).*w(:,:,1))-(u(:,:,n3-1).*w(:,:,n3-1)))/(2*dz)); 

Hu(:,:,:)=t-(t1+t2+t3);tu=(t1+t2+t3);
clear t t1 t2 t3
%% Hv
[t1,~,~]=gradient1(tau_uv);
[~,t2,~]=gradient1(tau_vv);
[~,~,t3]=gradient1(tau_vw);
% add BCs here
t1(1,:,:)=(tau_uv(2,:,:)-tau_uv(n1,:,:))/(2*dx);
t1(n1,:,:)=(tau_uv(1,:,:)-tau_uv(n1-1,:,:))/(2*dx);
t2(:,1,:)=(tau_vv(:,2,:)-tau_vv(:,1,:))/(dy);
t2(:,n2,:)=(tau_vv(:,n2,:)-tau_vv(:,n2-1,:))/(dy);
t3(:,:,1)=(tau_vw(:,:,2)-tau_vw(:,:,n3))/(2*dz);
t3(:,:,n3)=(tau_vw(:,:,1)-tau_vw(:,:,n3-1))/(2*dz);

t=t1+t2+t3;clear t1 t2 t3
% nonlinear terms
t1(2:n1-1,:,:)=((u(3:n1,:,:).*v(3:n1,:,:)-u(1:n1-2,:,:).*v(1:n1-2,:,:))/(2*dx)); % Spanwise
t1(n1,:,:)=((u(1,:,:).*v(1,:,:)-u(n1-1,:,:).*v(n1-1,:,:))/(2*dx));
t1(1,:,:)=((u(2,:,:).*v(2,:,:)-u(n1,:,:).*v(n1,:,:))/(2*dx));

t2(:,2:n2-1,:)=(((v(:,3:n2,:).^2)-(v(:,1:n2-2,:).^2))/(2*dy)); % central
t2(:,n2,:)=(0-((v(:,n2,:).^2))/(dy*0.5)); % backward
t2(:,1,:)=(((v(:,1,:).^2)-0)/(dy*0.5)); % forward

t3(:,:,2:n3-1)=(((v(:,:,3:n3).*w(:,:,3:n3))-(v(:,:,2:n3-1).*w(:,:,2:n3-1)))/(2*dz)); % central
t3(:,:,1)=(((v(:,:,2).*w(:,:,2))-(v(:,:,n3).*w(:,:,n3)))/(2*dz)); % forward
t3(:,:,n3)=(((v(:,:,1).*w(:,:,1))-(v(:,:,n3-1).*w(:,:,n3-1)))/(2*dz)); % backward

Hv(:,:,:)=t-(t1+t2+t3);tv=(t1+t2+t3);
clear t t1 t2 t3
%% Hw
[t1,~,~]=gradient1(tau_uw);
[~,t2,~]=gradient1(tau_vw);
[~,~,t3]=gradient1(tau_ww);
% add BCs here
t1(1,:,:)=(tau_uw(2,:,:)-tau_uw(n1,:,:))/(2*dx);
t1(n1,:,:)=(tau_uw(1,:,:)-tau_uw(n1-1,:,:))/(2*dx);
t2(:,1,:)=(tau_vw(:,2,:)-tau_vw(:,1,:))/(dy);
t2(:,n2,:)=(tau_vw(:,n2,:)-tau_vw(:,n2-1,:))/(dy);
t3(:,:,1)=(tau_ww(:,:,2)-tau_ww(:,:,n3))/(2*dz);
t3(:,:,n3)=(tau_ww(:,:,1)-tau_ww(:,:,n3-1))/(2*dz);

t=t1+t2+t3;clear t1 t2 t3
% nonlinear terms
t1(2:n1-1,:,:)=((u(3:n1,:,:).*w(3:n1,:,:)-u(1:n1-2,:,:).*w(1:n1-2,:,:))/(2*dx)); % Spanwise
t1(n1,:,:)=((u(1,:,:).*w(1,:,:)-u(n1-1,:,:).*w(n1-1,:,:))/(2*dx));
t1(1,:,:)=((u(2,:,:).*w(2,:,:)-u(n1,:,:).*w(n1,:,:))/(2*dx));

t2(:,2:n2-1,:)=(((w(:,3:n2,:).*v(:,3:n2,:))-(w(:,1:n2-2,:).*v(:,1:n2-2,:)))/(2*dy)); % central
t2(:,n2,:)=(0-((w(:,n2,:).*v(:,n2,:)))/(0.5*dy)); % backward
t2(:,1,:)=(((w(:,1,:).*v(:,1,:))-0)/(dy*0.5)); % forward

t3(:,:,2:n3-1)=(((w(:,:,3:n3).*w(:,:,3:n3))-(w(:,:,2:n3-1).*w(:,:,2:n3-1)))/(2*dz)); % central
t3(:,:,1)=(((w(:,:,2).*w(:,:,2))-(w(:,:,n3).*w(:,:,n3)))/(2*dz)); 
t3(:,:,n3)=(((w(:,:,1).*w(:,:,1))-(w(:,:,n3-1).*w(:,:,n3-1)))/(2*dz)); 

Hw(:,:,:)=t-(t1+t2+t3);tw=(t1+t2+t3);

