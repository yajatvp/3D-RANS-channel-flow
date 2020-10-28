%% FUNCTIONS
% MIXING LENGTH

function [tau_uu,tau_uv,tau_uw,tau_vv,tau_vw,tau_ww]=mix_len(u,v,w,ke,damping)
global dx dy dz n1 n2 n3 x y z alx2_2 Re

alx2=alx2_2*2;
% Note that velocity inputs u,v,w a     re on a centered grid, and must be
% handled accordingly

% 1 - Spanwise, 2 - Normal (Vertical), 3 - Streamwise
% Periodic boundary conditions (1 = n) for directions 1 and 3
if damping==0    % Prandtl's Model of mixing length
    l_k=0.217;
    for j=1:floor(n2/2)
        if (y(j)<=l_k*alx2_2)
            mixl(j)=0.41*y(j);

        else
            mixl(j)=0.089*alx2_2;
%             mixl(j)=0.41*ylimit;
        end
    end
    
    for j=(1+floor(n2/2)):n2
        if (abs(alx2-y(j))<=l_k*alx2_2)
            mixl(j)=0.41*(alx2-y(j));
        else
            mixl(j)=0.089*alx2_2;
%             mixl(j)=0.41*ylimit; % Hp: symmetry respect to center line. Otherwise modify ylimit
        end
    end
%     figure()
%     plot(mixl,y)
    
else
    
    wmean(:)=mean(mean(w(:,:,:),3),1);
    tau_wall=0.5*(abs(wmean(1)-0)/(0.5*dy)+abs(wmean(n2)-0)/(0.5*dy)); % mean tau at the wall
    
    for j=1:floor(n2/2)
        yplus1(j)=y(j)*sqrt(Re*tau_wall); % Check multiplication
    end
    
    for j=(1+floor(n2/2)):n2
        yplus2(j-floor(n2/2))=(alx2-y(j))*sqrt(Re*tau_wall);
    end
    
    for j=1:floor(n2/2)
        mixl1(j)=0.41*y(j)*(1-exp(-yplus1(j)/26)); % multiplication check
    end
    
    mixl0=0.089*alx2_2;
    mixl=mixl1;
    mixl(mixl1(:)>mixl0)=mixl0; % 3D mixling length is supposed to be considered - mixl(x,y,z), citation for yplus formulation to be added
    
    for j=1+floor(n2/2):n2
        mixl1q(j-floor(n2/2))=0.41*(alx2-y(j)).*(1-exp(-yplus2(j-floor(n2/2))/26)); % multiplication check
    end
    
    mixl0q=0.089*alx2_2; % same as above, suffix 'q' indicates for upper wall
    mixlq=mixl1q;
    mixlq(mixl1q(:)>mixl0q)=mixl0q;
    
    for j=1:floor(n2/2)
        yplus(j)=yplus1(j);
    end
    for j=1+floor(n2/2):n2
        yplus(j)=yplus2(j-floor(n2/2));
    end
    for j=1:floor(n2/2)
        mixl_f(j)=mixl(j);
    end
    for j=1+floor(n2/2):n2
        mixl_f(j)=mixlq(j-floor(n2/2));
    end
%     
%     figure()
%     plot(mixl_f,y)

end
% Starting equation for calculating turbulent viscosity
if damping==1
    mu_t(:,2:n2-1,:)=[(mixl_f(2:n2-1).^2)].*(abs((w(:,3:n2,:)-w(:,1:n2-2,:))/(2*dy)));
    mu_t(:,n2,:)=(mixl_f(n2)^2)*(abs((0-w(:,n2,:))/(0.5*dy)));
    mu_t(:,1,:)=(mixl_f(1)^2)*(abs((w(:,1,:)-0)/(0.5*dy)));
else
    mu_t(:,2:n2-1,:)=[(mixl(2:n2-1).^2)].*(abs((w(:,3:n2,:)-w(:,1:n2-2,:))/(2*dy)));
    mu_t(:,n2,:)=(mixl(n2).^2)*(abs((0-w(:,n2,:))/(0.5*dy)));
    mu_t(:,1,:)=(mixl(1).^2)*(abs((w(:,1,:)-0)/(0.5*dy)));
end

mu_t_avg=zeros(n2,1);
for j=1:n2
    for k=1:n3
        for i=1:n1
            mu_t_avg(j)=mu_t_avg(j)+mu_t(i,j,k);
        end
    end
    mu_t_avg(j)=mu_t_avg(j)/(n1*n3);
end
% figure()
% plot(y,mu_t_avg,'-k')

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
    
    % Combine calculated components for turbulent shear stress output
    tau_uu(:,:,:)=2*mu_t(:,:,:).*S11(:,:,:)-(2/3)*(ke(:,:,:));
    tau_vv(:,:,:)=2*mu_t(:,:,:).*S22(:,:,:)-(2/3)*(ke(:,:,:));
    tau_ww(:,:,:)=2*mu_t(:,:,:).*S33(:,:,:)-(2/3)*(ke(:,:,:));
    tau_uv(:,:,:)=2*mu_t(:,:,:).*S12(:,:,:);
    tau_uw(:,:,:)=2*mu_t(:,:,:).*S13(:,:,:);
    tau_vw(:,:,:)=2*mu_t(:,:,:).*S23(:,:,:);

end
