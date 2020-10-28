function [duvw]=diverg(u,v,w)
global dx dy dz n1 n2 n3 x y z
[du,~,~]=gradient1(u);
du(1,:,:)=(u(2,:,:)-u(n1,:,:))/(2*dx);
du(n1,:,:)=(u(1,:,:)-u(n1-1,:,:))/(2*dx);

[~,dv,~]=gradient1(v);
dv(:,1,:)=(v(:,2,:)-v(:,1,:))/(dy);
dv(:,n2,:)=(v(:,n2,:)-v(:,n2-1,:))/(dy);

[~,~,dw]=gradient1(w);
dw(:,:,1)=(w(:,:,2)-w(:,:,n3))/(2*dz);
dw(:,:,n3)=(w(:,:,1)-w(:,:,n3-1))/(2*dz);

duvw=du+dv+dw;
end
