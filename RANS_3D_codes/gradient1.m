function [ax,ay,az]=gradient1(a)
global n1 n2 n3 dx dy dz

ax=zeros(n1,n2,n3);
ay=zeros(n1,n2,n3);
az=zeros(n1,n2,n3);

ax(1,:,:)=(a(2,:,:)-a(1,:,:))/dx;
ax(n1,:,:)=(a(n1,:,:)-a(n1-1,:,:))/dx;
ax(2:n1-1,:,:)=(a(3:n1,:,:)-a(1:n1-2,:,:))/(2*dx);

ay(:,1,:)=(a(:,2,:)-a(:,1,:))/dy;
ay(:,n2,:)=(a(:,n2,:)-a(:,n2-1,:))/dy;
ay(:,2:n2-1,:)=(a(:,3:n2,:)-a(:,1:n2-2,:))/(2*dy);

az(:,:,1)=(a(:,:,2)-a(:,:,1))/dz;
az(:,:,n3)=(a(:,:,n3)-a(:,:,n3-1))/dz;
az(:,:,2:n3-1)=(a(:,:,3:n3)-a(:,:,1:n3-2))/(2*dz);

end