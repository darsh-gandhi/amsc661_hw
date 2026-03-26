clear; close all;

Nx = 100; % # of pts in x
Ny = 64; % # of pts in y

hx = 2.0*pi/Nx; % step size in x
hy = 2.0/Ny; % step size in y

x = -pi + (0:Nx-1)'*hx; % x mesh
y = (1:Ny)'*hy; % y mesh
y_bc = (0:Ny)'*hy;

N = Nx*Ny; % size of A (NxN)
kidx = @(i,j) (j-1)*Nx + mod(i,Nx) + 1;

%% create A matrix

I = zeros(6*N,1);
J = zeros(6*N,1);
V = zeros(6*N,1);
b = zeros(N,1);
id = 0;

for j = 1:Ny
    for i = 0:Nx-1
        k  = kidx(i, j);
        xi = x(i+1);
        dval = -2/hx^2;

        id=id+1;
        I(id)=k;
        J(id)=kidx(i+1,j);
        V(id)= 1/hx^2;

        id=id+1;
        I(id)=k;
        J(id)=kidx(i-1,j);
        V(id)= 1/hx^2;

        dval = dval - 1/hy^2;
        if j > 1
            id=id+1;
            I(id)=k;
            J(id)=kidx(i,j-1);
            V(id)= 1/hy^2;
        end
        
        dval = dval - 1/hy^2;
        if j < Ny
            id=id+1;
            I(id)=k;
            J(id)=kidx(i,j+1);
            V(id)= 1/hy^2;
        else
            id=id+1;
            I(id)=k;
            J(id)=kidx(i,j-1);
            V(id)= 1/hy^2;
        end

        id=id+1;
        I(id)=k;
        J(id)=k;
        V(id)=dval;

        if -pi/2 <= xi && xi <= pi/2
            b(k) = -cos(xi);
        end
    end
end

A = sparse(I(1:id), J(1:id), V(1:id), N, N);

%% solving system

u = A\b;
u_inner = reshape(u,Nx,Ny)';
U = [zeros(1,Nx);u_inner];

%% plotting
[X,Y] = meshgrid(x,y_bc);

figure;
contourf(X,Y,U,'LineColor','none');
colorbar;
colormap(hot);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
set(gca,'XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',{'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});