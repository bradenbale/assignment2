clear all;
close all;

%%  ELEC 4700 Assignment 2: Finite Difference Method
% ELEC 4700
% Assignment - 1
% Braden Bale (101072763)

%% Question 1(a)
% The matrix form of the Finite Difference Method is used solve for the
% electrostatic potential of a recetangular region with insolating sides.
% Below some varibles for simulation are set up, with F and G being used
% for the matrices.

L = 3;
W = 2;

dx = 0.25;
dy = 0.25;
nx = L/dx; 
ny = W/dy;

V0 = 1;

F = zeros(nx*ny, 1);
G = sparse(nx*ny, nx*ny);

%% 
% Below the Finite Difference Method is used to set up the G and F matrices . 



for i = 1:(nx)
    for j = 1:(ny)
        n = (j + (i-1).*ny);
        if i == 1
            G(n, n) = 1;
            F(n) = V0;
        elseif i == nx
            G(n, n) = 1;
            F(n) = 0;
        elseif j == 1
            G(n, n) = -1/dy^2 + (1/dx^2)*-2;
            G(n, ((j-1) + (i-1).*ny)) = 1/(dx^2);
            G(n, ((j+1) + (i-1).*ny)) = 1/(dx^2);
            G(n, (j + (i).*ny)) = 1/(dy^2);
        elseif j == ny
            G(n, n) = -1/dy^2 + (1/dx^2)*-2;
            G(n, ((j-1) + (i-1).*ny)) = 1/(dx^2);
            G(n, ((j+1) + (i-1).*ny)) = 1/(dx^2);
            G(n, (j + (i-2).*ny)) = 1/(dy^2);
        else
            G(n, n) = (1/dx^2 + 1/dy^2)*-2;
            G(n, ((j-1) + (i-1).*ny)) = 1/(dx^2);
            G(n, ((j+1) + (i-1).*ny)) = 1/(dx^2);
            G(n, (j + (i-2).*ny)) = 1/(dy^2);
            G(n, (j + (i).*ny)) = 1/(dy^2);
        end
    end
end

V = G\F;

V = reshape(V,[],nx)';

%% 
% Finally, in Figure 1 below a 2-D plot of V(x) is shown.


figure(1);
surf(V);
title('2-D Plot of V(x)');
xlabel('x');
ylabel('y');
view(120, 30);

%%

%% Question 1(b)
% When comparing the analytical solution to the solution retrieved from
% using the Finite Difference Method the analytical solution could not be
% found. This code is commented out as it does not correctly complete the
% question.

%% 
% Below the matrices G anf F are set up similarily to Question 1(a).


% for i = 1:(nx)
%     for j = 1:(ny)
%         n = (j + (i-1).*ny);
%         if i == 1
%             G(n, n) = 1;
%             F(n) = V0;
%         elseif i == nx
%             G(n, n) = 1;
%             F(n) = V0;
%         elseif j == 1
%             G(n, n) = 1;
%             F(n) = 0;
%         elseif j == ny
%             G(n, n) = 1;
%             F(n) = 0;
%         else
%             G(n, n) = (1/dx^2 + 1/dy^2)*-2;
%             G(n, ((j-1) + (i-1).*ny)) = 1/(dx^2);
%             G(n, ((j+1) + (i-1).*ny)) = 1/(dx^2);
%             G(n, (j + (i-2).*ny)) = 1/(dy^2);
%             G(n, (j + (i).*ny)) = 1/(dy^2);
%         end
%     end
% end

%% 
% The analytical solution was attempted by building a matrices which was
% the same size as V using xp and xp below and repeating them for the
% length of nx and ny.


% VA = zeros(nx, ny);
% 
% max = 1000;
% 
% a = W;
% b = L;
% 
% xp = linspace(-L, L, ny);
% yp = linspace(0, W, nx);
% 
% 
% 
% for i = 1: nx
%    x(i, :) = xp; 
% end
% 
% for j = 1: ny
%    y(:, j) = yp; 
% end
% 

%% 
% Below the for loop to caluculate the analytical solution is shown.


% for k = 1: max
%     n = k*2 -1;
%     VA = VA + (1/n)*((cosh((n*pi.*x)./a))./(cosh((n.*pi.*b)./a))).*(sin((n.*pi.*y)./a));
% end
% 
% VA = VA.*4.*1./pi;
% 
% 
% V = G\F;
% 
% V = reshape(V,[],nx)';
% 
% figure(1);
% surf(V);
% xlabel('x');
% ylabel('y');
% 
% figure(2);
% surf(VA');
% xlabel('x');
% ylabel('y');

%% Question 2
% From Question 1, conductivity of a surface could be added to the area.
% Two boxes were added to create a bottle-neck to understand conductivity's
% effect on the system. The conductivity in the boxes is set to 10^-2 and
% outside the boxes is set to 1. For this question a conductivity map,
% voltage map, electric field matrix and current density matrix were set up
% and shown.


rho_out_of_box = 1;
rho_in_box = 10^-2;

%% 
% Below the boxees for the conductivity map are set up.


cMap = ones(nx, ny).*rho_out_of_box;
for i = 1: nx
    for j = 1 : ny
        if ((i > ((0.5*nx)- (0.1*nx))) && (i <((0.5*nx)+(0.1*nx))) && (j > ((0.5*ny)+(0.25*ny))))
            cMap(i, j) = rho_in_box;
        elseif ((i > ((0.5*nx)- (0.1*nx))) && (i < ((0.5*nx)+(0.1*nx))) && (j < ((0.5*ny)-(0.25*ny))))
            cMap(i, j) = rho_in_box;
        end
    end
end

G2 = sparse(nx*ny);
B = zeros(1,nx*ny);

%% 
% The conductivity map for the Question is set up using the Finite
% Difference Method.

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G2(n, n) = 1;
            B(n) = 1;
        elseif i == nx
            G2(n, n) = 1;
            B(n) = 0;
        elseif j == 1
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxp+rxm+ryp);
            G2(n,nxp) = rxp;
            G2(n,nxm) = rxm;
            G2(n,nyp) = ryp;
        elseif j ==  ny
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;

            G2(n,n) = -(rxm+rxp+rym);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxm+rxp+rym+ryp);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
            G2(n,nyp) = ryp;
        end

    end
end

%% 
% Then the voltage map is set up below.


V2 = G2\B';

Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        Vmap(i, j) = V2(n);
    end
end

%% 
% Finally, with the voltage map found the electric field can be found.


for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
        elseif i == nx
            Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
        else
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
        elseif j == ny
            Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
        else
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

Jx = cMap.*Ex;
Jy = cMap.*Ey;

%% 
% Figure 2 below shows the conductiviy map.

figure(2)
surf(cMap);
xlabel('x');
ylabel('y');
title('Conductivity Map');

%% 
% Figure 3 below shows the potential map.

figure(3)
surf(Vmap);
xlabel('x');
ylabel('y');
title('Potential Map');

%% 
% Figure 4 below shows the electric field.

figure(4)
quiver(Ex', Ey',  1);
xlabel('x');
ylabel('y');
title('Electric Field');

%% 
% Figure 5 below shows the current density.
% 


figure(5)
quiver(Jx', Jy', 1);
xlabel('x');
ylabel('y');
title('Current Density');

%% 
% From this section onward small changes are made to the code to demonstrate
% different tests to the system.


Jx = 0;
Jy = 0;
Ex = 0;
Ey = 0;
Vmap = 0;
cMap = 0;
V2 = 0;

G2 = sparse(nx*ny);
B = zeros(1,nx*ny);

%% 
% First, a new current density is set. The variable dx and dy are changed
% from 0.25 to 0.10.


dx = 0.10;
dy = 0.10;
nx = L/dx; 
ny = W/dy;

cMap = ones(nx, ny).*rho_out_of_box;
for i = 1: nx
    for j = 1 : ny
        if ((i > ((0.5*nx)- (0.1*nx))) && (i <((0.5*nx)+(0.1*nx))) && (j > ((0.5*ny)+(0.3*ny))))
            cMap(i, j) = rho_in_box;
        elseif ((i > ((0.5*nx)- (0.1*nx))) && (i < ((0.5*nx)+(0.1*nx))) && (j < ((0.5*ny)-(0.3*ny))))
            cMap(i, j) = rho_in_box;
        end
    end
end

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G2(n, n) = 1;
            B(n) = 1;
        elseif i == nx
            G2(n, n) = 1;
            B(n) = 0;
        elseif j == 1
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxp+rxm+ryp);
            G2(n,nxp) = rxp;
            G2(n,nxm) = rxm;
            G2(n,nyp) = ryp;
        elseif j ==  ny
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;

            G2(n,n) = -(rxm+rxp+rym);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxm+rxp+rym+ryp);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
            G2(n,nyp) = ryp;
        end

    end
end



V2 = G2\B';

Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        Vmap(i, j) = V2(n);
    end
end

for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
        elseif i == nx
            Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
        else
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
        elseif j == ny
            Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
        else
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

Jx = cMap.*Ex;
Jy = cMap.*Ey;

%%
% Figure 6 is shown below.

figure(6)
quiver(Jx', Jy', 1);
xlabel('x');
ylabel('y');
title('Current with Mesh Density set to 1/10.');

%% 
% Next, the size of the bottleneck is varied.


Jx = 0;
Jy = 0;
Ex = 0;
Ey = 0;
Vmap = 0;
V2 = 0;

dx = 0.25;
dy = 0.25;
nx = L/dx; 
ny = W/dy;

G2 = sparse(nx*ny);
B = zeros(1,nx*ny);

cMap = ones(nx, ny).*rho_out_of_box;
for i = 1: nx
    for j = 1 : ny
        if ((i > ((0.5*nx)- (0.1*nx))) && (i <((0.5*nx)+(0.1*nx))) && (j > ((0.5*ny)+(0.1*ny))))
            cMap(i, j) = rho_in_box;
        elseif ((i > ((0.5*nx)- (0.1*nx))) && (i < ((0.5*nx)+(0.1*nx))) && (j < ((0.5*ny)-(0.1*ny))))
            cMap(i, j) = rho_in_box;
        end
    end
end

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G2(n, n) = 1;
            B(n) = 1;
        elseif i == nx
            G2(n, n) = 1;
            B(n) = 0;
        elseif j == 1
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxp+rxm+ryp);
            G2(n,nxp) = rxp;
            G2(n,nxm) = rxm;
            G2(n,nyp) = ryp;
        elseif j ==  ny
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;

            G2(n,n) = -(rxm+rxp+rym);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxm+rxp+rym+ryp);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
            G2(n,nyp) = ryp;
        end

    end
end

V2 = G2\B';

Vmap = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        Vmap(i, j) = V2(n);
    end
end

for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
        elseif i == nx
            Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
        else
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
        elseif j == ny
            Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
        else
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

Jx = cMap.*Ex;
Jy = cMap.*Ey;

%%
% Figure 7 is shown below.


figure(7)
quiver(Jx', Jy', 1);
xlabel('x');
ylabel('y');
title('Current with a Larger Bottle-Neck.');

%%

Jx = 0;
Jy = 0;
Ex = 0;
Ey = 0;
Vmap = 0;
cMap = 0;
V2 = 0;

G2 = sparse(nx*ny);
B = zeros(1,nx*ny);

cMap = ones(nx, ny).*rho_out_of_box;
for i = 1: nx
    for j = 1 : ny
        if ((i > ((0.5*nx)- (0.1*nx))) && (i <((0.5*nx)+(0.1*nx))) && (j > ((0.5*ny)+(0.4*ny))))
            cMap(i, j) = rho_in_box;
        elseif ((i > ((0.5*nx)- (0.1*nx))) && (i < ((0.5*nx)+(0.1*nx))) && (j < ((0.5*ny)-(0.4*ny))))
            cMap(i, j) = rho_in_box;
        end
    end
end

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G2(n, n) = 1;
            B(n) = 1;
        elseif i == nx
            G2(n, n) = 1;
            B(n) = 0;
        elseif j == 1
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxp+rxm+ryp);
            G2(n,nxp) = rxp;
            G2(n,nxm) = rxm;
            G2(n,nyp) = ryp;
        elseif j ==  ny
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;

            G2(n,n) = -(rxm+rxp+rym);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxm+rxp+rym+ryp);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
            G2(n,nyp) = ryp;
        end

    end
end

V2 = G2\B';

Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        Vmap(i, j) = V2(n);
    end
end

for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
        elseif i == nx
            Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
        else
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
        elseif j == ny
            Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
        else
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

Jx = cMap.*Ex;
Jy = cMap.*Ey;

%%
% Figure 8 is shown below.


figure(8)
quiver(Jx', Jy', 1);
xlabel('x');
ylabel('y');
title('Current with a Smaller Bottle-Neck.');
%%

% Then, the conductivity in the box was varied.

Jx = 0;
Jy = 0;
Ex = 0;
Ey = 0;
Vmap = 0;
cMap = 0;
V2 = 0;

G2 = sparse(nx*ny);
B = zeros(1,nx*ny);

rho_in_box = 0.1;

cMap = ones(nx, ny).*rho_out_of_box;
for i = 1: nx
    for j = 1 : ny
        if ((i > ((0.5*nx)- (0.1*nx))) && (i <((0.5*nx)+(0.1*nx))) && (j > ((0.5*ny)+(0.25*ny))))
            cMap(i, j) = rho_in_box;
        elseif ((i > ((0.5*nx)- (0.1*nx))) && (i < ((0.5*nx)+(0.1*nx))) && (j < ((0.5*ny)-(0.25*ny))))
            cMap(i, j) = rho_in_box;
        end
    end
end

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G2(n, n) = 1;
            B(n) = 1;
        elseif i == nx
            G2(n, n) = 1;
            B(n) = 0;
        elseif j == 1
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxp+rxm+ryp);
            G2(n,nxp) = rxp;
            G2(n,nxm) = rxm;
            G2(n,nyp) = ryp;
        elseif j ==  ny
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;

            G2(n,n) = -(rxm+rxp+rym);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxm+rxp+rym+ryp);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
            G2(n,nyp) = ryp;
        end

    end
end

V2 = G2\B';

Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        Vmap(i, j) = V2(n);
    end
end

for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
        elseif i == nx
            Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
        else
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
        elseif j == ny
            Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
        else
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
        end
    end
end



Ex = -Ex;
Ey = -Ey;

Jx = cMap.*Ex;
Jy = cMap.*Ey;

%%
% Figure 9 is shown below.


figure(9)
quiver(Jx', Jy', 1);
xlabel('x');
ylabel('y');
title('Current with Conductivity in the Box set to 0.1.');

%%

Jx = 0;
Jy = 0;
Ex = 0;
Ey = 0;
Vmap = 0;
cMap = 0;
V2 = 0;

G2 = sparse(nx*ny);
B = zeros(1,nx*ny);

rho_in_box = 0.5;

cMap = ones(nx, ny).*rho_out_of_box;
for i = 1: nx
    for j = 1 : ny
        if ((i > ((0.5*nx)- (0.1*nx))) && (i <((0.5*nx)+(0.1*nx))) && (j > ((0.5*ny)+(0.25*ny))))
            cMap(i, j) = rho_in_box;
        elseif ((i > ((0.5*nx)- (0.1*nx))) && (i < ((0.5*nx)+(0.1*nx))) && (j < ((0.5*ny)-(0.25*ny))))
            cMap(i, j) = rho_in_box;
        end
    end
end

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G2(n, n) = 1;
            B(n) = 1;
        elseif i == nx
            G2(n, n) = 1;
            B(n) = 0;
        elseif j == 1
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxp+rxm+ryp);
            G2(n,nxp) = rxp;
            G2(n,nxm) = rxm;
            G2(n,nyp) = ryp;
        elseif j ==  ny
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;

            G2(n,n) = -(rxm+rxp+rym);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxm+rxp+rym+ryp);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
            G2(n,nyp) = ryp;
        end

    end
end

V2 = G2\B';

Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        Vmap(i, j) = V2(n);
    end
end

for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
        elseif i == nx
            Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
        else
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
        elseif j == ny
            Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
        else
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

Jx = cMap.*Ex;
Jy = cMap.*Ey;

%%
% Figure 10 is shown below.


figure(10)
quiver(Jx', Jy', 1);
xlabel('x');
ylabel('y');
title('Current with Conductivity in the Box set to 0.5.');

%%

Jx = 0;
Jy = 0;
Ex = 0;
Ey = 0;
Vmap = 0;
cMap = 0;
V2 = 0;

G2 = sparse(nx*ny);
B = zeros(1,nx*ny);

rho_in_box = 0.001;

cMap = ones(nx, ny).*rho_out_of_box;
for i = 1: nx
    for j = 1 : ny
        if ((i > ((0.5*nx)- (0.1*nx))) && (i <((0.5*nx)+(0.1*nx))) && (j > ((0.5*ny)+(0.25*ny))))
            cMap(i, j) = rho_in_box;
        elseif ((i > ((0.5*nx)- (0.1*nx))) && (i < ((0.5*nx)+(0.1*nx))) && (j < ((0.5*ny)-(0.25*ny))))
            cMap(i, j) = rho_in_box;
        end
    end
end

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G2(n, n) = 1;
            B(n) = 1;
        elseif i == nx
            G2(n, n) = 1;
            B(n) = 0;
        elseif j == 1
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxp+rxm+ryp);
            G2(n,nxp) = rxp;
            G2(n,nxm) = rxm;
            G2(n,nyp) = ryp;
        elseif j ==  ny
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;

            G2(n,n) = -(rxm+rxp+rym);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G2(n,n) = -(rxm+rxp+rym+ryp);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
            G2(n,nyp) = ryp;
        end

    end
end

V2 = G2\B';

Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        Vmap(i, j) = V2(n);
    end
end

for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
        elseif i == nx
            Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
        else
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
        elseif j == ny
            Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
        else
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

Jx = cMap.*Ex;
Jy = cMap.*Ey;

%%
% Figure 11 is shown below.


figure(11)
quiver(Jx', Jy', 1);
xlabel('x');
ylabel('y');
title('Current with Conductivity in the Box set to 0.001.');

%%



