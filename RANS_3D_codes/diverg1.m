function duvw=diverg1(u,v,w)
global alx1 alx2 alx3 n1 n2 n3
du=sum(sum((u(n1,:,:)-u(1,:,:))))*alx2*alx3;
dv=sum(sum((v(:,n2,:)-v(:,1,:))))*alx1*alx3;
dw=sum(sum((w(:,:,n3)-w(:,:,1))))*alx1*alx2;

duvw=du+dv+dw;
end