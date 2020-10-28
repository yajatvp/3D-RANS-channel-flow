function [p,px,py,pz]=solve_pressure(Hu,Hv,Hw)
global dx dy dz n1 n2 n3 x y z

%% Create the matrix
% The blocks can be subdivided into smaller components for easier
% combination

% A and Ap blocks of size n1 correspond to each row in the spanwise direction
%    C corresponds to the j+1, j-1
%    C2 corresponds to setting up the equations at j=1, and J=n2 to be
%    forward and backward differencing

% blocks of size n1*n2 correspond to all rows for a given value of k (i.e
%    a plane in the 12 dimension for each value of k)

% The final matrix of size n1*n2*n3 corresponds to all the planes combined
% for every value of k (i.e. the full volume)

%% blocks of size n1

ap = 1;
a = -2*((1/dx)^2+(1/dy)^2+(1/dz)^2);
b = 1/dx^2;
c = 1/dy^2;

% A_plus block
    Ap=diag(ap*ones((n1),1));
    % Add on side diagonals
    %Ap = Ap + diag(b*ones((n1)-1,1),1) + diag(b*ones((n1)-1,1),-1) ;
    % Add boundary condition diagonals
    %Ap = Ap + diag(b*ones(1,1),n1-1) + diag(b*ones(1,1),-n1+1);

% Wall condition block for pressure at top and bottom rows
    C2 =diag(-1*ones((n1),1));

% Wall Condition C block for j+2 and j-2
    C =diag(c*ones((n1),1));
    
% A block (Internal)

    A = diag(a*ones((n1),1));
    % Add on side diagonals
    A = A + diag(b*ones((n1)-1,1),1) + diag(b*ones((n1)-1,1),-1) ;
    % Add boundary condition diagonals
    A = A + diag(b*ones(1,1),n1-1) + diag(b*ones(1,1),-n1+1);

%% blocks of size n1*n2

% D block for k bounds 
d = 1/dz^2;
D = diag(d*ones((n1*n2),1));
% Include pressure corrections at bottom wall
D(1:n1,1:n1) = diag(0*ones((n1),1));
% Include pressure corrections at top wall
D(1+n1*(n2-1):n1*n2,1+n1*(n2-1):n1*n2) = diag(0*ones((n1),1));

% Initialize array for one block of n1*n2 (plane in x,y)
Block_J = zeros(n1*n2);

% assign all x blocks across j to define a set of equations for a plane
for j = 1:n2
    Block_J(1+n1*(j-1):n1*j,1+n1*(j-1):n1*j) = A;
    if j == 1
       % Include Ap block at bounds
       Block_J(1:n1, 1:n1) = Ap;
       % Include C2 block for j+1
       Block_J(1:n1, n1+1:2*n1) = C2;
       % Include C block for j+2
       %Block_J(1:n1, 2*n1+1:3*n1) = C;
    elseif j == n2
       % Include Ap block at bounds
       Block_J(n1*(n2-1)+1:n1*n2, n1*(n2-1)+1:n1*n2) = Ap; 
       % Include C2 block for j-1
       Block_J(n1*(n2-1)+1:n1*n2, n1*(n2-2)+1:n1*(n2-1)) = C2; 
       % Include C block for j-2
       %Block_J(n1*(n2-1)+1:n1*n2, n1*n2-3*(n1)+1:n1*(n2-2)) = C;
    else
        Block_J(1+n1*(j-1):n1*j, 1+n1*j-2*n1:n1*(j-1)) = C;
        Block_J(1+n1*(j-1):n1*j, 1+ n1*j:n1*(1+j)) = C;
    end
end

%% Final combined matrix of size n1*n2*n3

xy = n1*n2;

A_matrix = zeros(n1*n2*n3);
for k = 1:n3
    % Add the Block_J matrix for every value of k (every plane)
    A_matrix(1+xy*(k-1):xy*(k),1+xy*(k-1):xy*(k)) = Block_J;
    
    % Apply uniform pressure outlet
    if k == n3
        A_matrix(1+xy*(k-1):xy*(k),1+xy*(k-1):xy*(k)) = diag(1*ones((n1*n2),1));
    end
end

for  k = 1:n3-1
    % Add the D block for the streamwise portion
    % k + 1
    A_matrix(1+xy*(k-1):xy*(k),1+xy*k:xy*(k+1)) = D;
    % k - 1
    if k ~= n3-1
        A_matrix(1+xy*k:xy*(k+1),1+xy*(k-1):xy*(k)) = D;
    end
end
% Add the D block for periodic BC of k
    % k + 1
    A_matrix(1:xy, n1*n2*(n3-1)+1:n1*n2*n3) = D;
    % k - 1
    %A_matrix(n1*n2*(n3-1)+1:n1*n2*n3, 1:xy) = D;
%% Right size B vector

[Hux,~,~]=gradient1(Hu);
Hux(1,:,:)=(Hu(2,:,:)-Hu(n1,:,:))/(2*dx);
Hux(n1,:,:)=(Hu(1,:,:)-Hu(n1-1,:,:))/(2*dx);

[~,Hvy,~]=gradient1(Hv);
Hvy(:,1,:)=(Hv(:,2,:)-Hv(:,1,:))/(dy);
Hvy(:,n2,:)=(Hv(:,n2,:)-Hv(:,n2-1,:))/(dy);

[~,~,Hwz]=gradient1(Hw);
Hwz(:,:,1)=(Hw(:,:,2)-Hw(:,:,n3))/(2*dz);
Hwz(:,:,n3)=(Hw(:,:,1)-Hw(:,:,n3-1))/(2*dz);

tmp=0;
for k=1:n3
    for j=1:n2
        for i=1:n1
            tmp=tmp+1;
            B(tmp)=Hux(i,j,k)+Hvy(i,j,k)+Hwz(i,j,k);
            if j==1 | j==n2
            B(tmp)=0;
            end
        end
    end
end

B(1+(n1*n2*(n3-1)):end)=1; % Pressure value at outlet plane

% [p11,it]=sor(A_matrix,B');
% 
% p1=linsolve(A_matrix,B');
% 
% [p1,err,it,flag]=sor1(A_matrix,ones(n1*n2*n3),B',1.5);
%%
p1=A_matrix\B';

%%
tmp=0;
for k=1:n3
    for j=1:n2
        for i=1:n1
            tmp=tmp+1;
            p(i,j,k)=p1(tmp);
        end
    end
end

%%
[px,py,pz]=gradient1(p);
px(1,:,:)=(p(2,:,:)-p(n1,:,:))/(2*dx);
px(n1,:,:)=(p(1,:,:)-p(n1-1,:,:))/(2*dx);
py(:,1,:)=(p(:,2,:)-p(:,1,:))/(dy);
py(:,n2,:)=(p(:,n2,:)-p(:,n2-1,:))/(dy);
pz(:,:,1)=(p(:,:,2)-p(:,:,n3))/(2*dz);
pz(:,:,n3)=(p(:,:,1)-p(:,:,n3-1))/(2*dz);
% 
% % %
% figure;%set(gcf, 'Position', get(0, 'Screensize'));
% pcolor(z,y,squeeze(p(10,:,:))); shading flat; axis normal;%caxis([0.8 1]);
% colormap(gca,'jet');c = colorbar;c.Label.String = 'Pressure';
