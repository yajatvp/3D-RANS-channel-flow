%% FUNCTIONS
% VISCOUS TERMS CALCULATION
function [tau_uu,tau_uv,tau_uw,tau_vv,tau_vw,tau_ww]=viscous(u,v,w,Re)
global dx dy dz n1 n2 n3

% Note that velocity inputs u,v,w a     re on a centered grid, and must be
% handled accordingly

% 1 - Spanwise, 2 - Normal (Vertical), 3 - Streamwise
% Periodic boundary conditions (1 = n) for directions 1 and 3

% Starting equation for calculating turbulent viscosity

% Calculate fluid shear strain rate tensor - Central difference method
    % Diagonal Components
    S33(:,:,2:n3-1)=((w(:,:,3:n3)-w(:,:,1:n3-2))/(2*dz)); % Streamwise
    S33(:,:,1)=((w(:,:,2)-w(:,:,n3))/(2*dz));
    S33(:,:,n3)=((w(:,:,1)-w(:,:,n3-1))/(2*dz));

    S22(:,2:n2-1,:)=((v(:,3:n2,:)-v(:,1:n2-2,:))/(2*dy)); % Vertical
    S22(:,1,:)=(v(:,1,:)-0)/(0.5*dy); % Add citation preferably - Forward diff
    S22(:,n2,:)=(0-v(:,n2,:))/(0.5*dy); % Backward diff

    S11(2:n1-1,:,:)=((u(3:n1,:,:)-u(1:n1-2,:,:))/(2*dx)); % Spanwise
    S11(n1,:,:)=((u(1,:,:)-u(n1-1,:,:))/(2*dx));
    S11(1,:,:)=((u(2,:,:)-u(n1,:,:))/(2*dx));
    
    % Symmetric components
    t1(2:n1-1,:,:)=((v(3:n1,:,:)-v(1:n1-2,:,:))/(2*dx));t2(:,2:n2-1,:)=(u(:,3:n2,:)-u(:,1:n2-2,:))/(2*dy);
    % Adjust for the edges
    t1(n1,:,:)=((v(1,:,:)-v(n1-1,:,:))/(2*dx)); t2(:,n2,:)=(0-u(:,n2,:))/(0.5*dy);
    t1(1,:,:)=((v(2,:,:)-v(n1,:,:))/(2*dx));    t2(:,1,:)=(u(:,1,:)-0)/(0.5*dy);
    S12(:,:,:)=0.5*(t1+t2);

    clear t1 t2
    % Note that for this component, the boundary condition at the wall
    % (velocity = 0) was taken into account
    t1(:,:,2:n3-1)=((u(:,:,3:n3)-u(:,:,1:n3-2))/(2*dz));t2(2:n1-1,:,:)=((w(3:n1,:,:)-w(1:n1-2,:,:))/(2*dx));
    % Adjust for the edges
    t1(:,:,n3)=((u(:,:,1)-u(:,:,n3-1))/(2*dz));t2(n1,:,:)=((w(1,:,:)-w(n1-1,:,:))/(2*dx));
    t1(:,:,1)=((u(:,:,2)-u(:,:,n3))/(2*dz));t2(1,:,:)=((w(2,:,:)-w(n1,:,:))/(2*dx));
    S13(:,:,:)=0.5*(t1+t2);

    clear t1 t2
    t1(:,:,2:n3-1)=((v(:,:,3:n3)-v(:,:,1:n3-2))/(2*dz));t2(:,2:n2-1,:)=(w(:,3:n2,:)-w(:,1:n2-2,:))/(2*dy);
    % Adjust for the edges
    t1(:,:,n3)=((v(:,:,1)-v(:,:,n3-1))/(2*dz));t2(:,n2,:)=(0-w(:,n2,:))/(0.5*dy);
    t1(:,:,1)=((v(:,:,2)-v(:,:,n3))/(2*dz));t2(:,1,:)=(w(:,1,:)-0)/(0.5*dy);
    S23(:,:,:)=0.5*(t1+t2);
    mu=1/Re;
    % Combine calculated components for turbulent shear stress output
    tau_uu(:,:,:)=2*mu*S11(:,:,:); 
    tau_vv(:,:,:)=2*mu*S22(:,:,:);
    tau_ww(:,:,:)=2*mu*S33(:,:,:);
    tau_uv(:,:,:)=2*mu*S12(:,:,:);
    tau_uw(:,:,:)=2*mu*S13(:,:,:);
    tau_vw(:,:,:)=2*mu*S23(:,:,:);

end
