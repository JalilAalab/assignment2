%% ELEC 4700 Assignment 2 Jalil (Rohana) Aalab #100995788 Feb/25/2018
%% Part 1 A
% This code is used to provide a visual understanding of the electrostatic
% potential that is present within a given 2D rectangular region. This
% intial section will set up the grid and basic variables prior to
% evaluating the potentials. The first case will assume that the potential
% does not change with respect to the y-coordinates, and therefore will be
% treated as a 1D situation. 

% Inputs
L = 30;
W = 20;

%Grid
x = linspace(0,L);
y = linspace(0,W);
dx = x(2) - x(1);
dy = y(2) - y(1);

% Make matrices
N = L*W;
M = zeros(N,N);
B = zeros(N,1);

% Interior Points
for i = 2:L-1
    for j = 2:W-1
        n = i + (j-1)*L;
        M(n,n) = -4;
        M(n,n-1) = 1;
        M(n,n+1) = 1;
        M(n,n-L) = 1;
        M(n,n+L) = 1;
        B(n,1) = 0;
    end
end

% Left BC point
i = 1;
for j = 1:W
    n = i + (j-1)*L;
    M(n,n) = 1;
     B(n,1) = 1;
end

% Right BC 
i = L;
for j = 1:W
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 0;
end

% Bottom BC
j = 1;
for i = 1:L
    n = i + (j-1)*L;
    M(n,n) = 1;
end

% Top BC
j = W;
for i = 1:L
    n = i + (j-1)*L;
    M(n,n) = 1;
end

% Solve for potential
phi_vec = M\B;

% Convert vector to 2D array
for i = 1 : L
    for j = 1 : W
        n = i + (j-1)*L;
        phi(i,j) = phi_vec(n);
    end
end

% Plot
figure(1);
mesh(phi);
xlabel('x');
ylabel('y');
zlabel('V(x,y)');
title('V(x,y)');   
%%
% The graph shows that the potential (V) is 1 on one end of the region, and
% on the other end it is 0. Between these sides, the potential slowly drops
% as it is distributed from the V = 1 side to the zero potential
% side. This produces the curve that lowers down to a flat line, and this
% represents that the raised potential spreads out to the other side but
% the value of the potential drops as it is moved further away from the
% V = 1 side. 
%% Part 1 B
% This section provides a visual representation of the variance of electric
% potential over a 2D rectangular region. 

% Interior Points
for i = 2:L-1
    for j = 2:W-1
        n = i + (j-1)*L;
        M(n,n) = -4;
        M(n,n-1) = 1;
        M(n,n+1) = 1;
        M(n,n-L) = 1;
        M(n,n+L) = 1;
        B(n,1) = 0;
    end
end

% Left BC point
i = 1;
for j = 1:W
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 1;
end

% Right BC 
i = L;
for j = 1:W
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 1;
end

% Bottom BC
j = 1;
for i = 1:L
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 0;
end

% Top BC
j = W;
for i = 1:L
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 0;
end

% Solve for potential
phi_vec = M\B;

% Convert vector to 2D array
for i = 1 : 30
    for j = 1 : 20
        n = i + (j-1)*L;
        phi(i,j) = phi_vec(n);
    end
end

% Plot
figure(2);
surf(phi);
xlabel('x');
ylabel('y');
zlabel('V(x,y)');
title('V(x,y) Numeric');
%% 
% As opposed to the previous figure, this region has boundary conditions
% such that both ends have a raised potential at V = 1. The areas of the
% region closest to either end are influenced by the raised potential and
% the potential of those areas rises as well. However, the potential drops
% in proportion to the distance away from either end. Therefore, the middle
% of the region is furtherst away from either end and has the lowest
% potential because it is far from the influence of the raised voltage.

% Part 1 b analytical
phi2 = phi;
a = 30;
b = 10;
added = 0;
for i = 1:L
    for j = 1:W
        for n = 1:2:1003
            added = added + ((1/n)*(cosh(n*pi*i/a))*(sin(n*pi*j/a))*(1/(cosh(n*pi*b/a))));
        end
    phi2(i,j) = (4/pi) * added;
    end
end
figure(3);
surf(phi2);
xlabel('x');
ylabel('y');
zlabel('V(x,y)');
title('V(x,y) Analytic');
%%
% This graph is the analytic solution for the electric potential. It is
% derived from using the analytic series solution, which is described above
% by the variable 'added'. It can be seen that this analytic solution is
% not as accurate as the previous numeric solution. This may be due to the
% fact that the variable n is supposed to be all the odd numbers from 1 to
% infinity, however this code only had a maximum n of 1003. Therefore, the
% graph does not properly show the raised voltage of 1 V on both ends of
% the region. For this example, the superior solution is therefore the
% numeric solution. 
%% Part 2 A
% This section provides a visual understanding of the current flow in a
% given 2D rectangular region. The distribution of the conductivity
% (sigma(x,y)), current density(J(x,y)), horizontal electric field (Ex),
% vertical electric field (Ey), and electric potential (V(x,y)) are
% displayed. The region has a 'bottleneck' area which is highly resistive,
% and the effect of this on other variables is explored.

% Interior Points
for i = 2:L-1
    for j = 2:W-1
        n = i + (j-1)*L;
        M(n,n) = -4;
        M(n,n-1) = 1;
        M(n,n+1) = 1;
        M(n,n-L) = 1;
        M(n,n+L) = 1;
        B(n,1) = 0;
    end
end

% Left BC point
i = 1;
for j = 1:W
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 1;
end

% Right BC 
i = L;
for j = 1:W
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 1;
end

% Bottom BC
j = 1;
for i = 1:L
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 0;
end

% Top BC
j = W;
for i = 1:L
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 0;
end

% Solve for potential
phi_vec = M\B;

% Convert vector to 2D array
for i = 1 : 30
    for j = 1 : 20
        n = i + (j-1)*L;
        phi(i,j) = phi_vec(n);
    end
end

% Plot
figure(4);
mesh(phi);
xlabel('X');
ylabel('Y');
zlabel('V(x,y)');
title('V(x,y)');
%% 
% This figure is the same as the previous figure with identical potentials
% and boundary conditions. Other variables and factors are changed in order
% to compare their effects, therefore this figure is a 'control' variable
% to relate the differences of various physical effects. 

% Electric Field 
Vmap = zeros(L,W);
for i = 1:L
    for j = 1:W
        n = j + (i-1)*W;
        Vmap(i,j) = phi_vec(n);
    end
end

for i = 1:L
    for j = 1:W
        if i == 1
            Ex(i,j) = (Vmap(i+1,j) - Vmap(i,j));
        elseif i == L
            Ex(i,j) = (Vmap(i,j) - Vmap(i-1,j));
        else
            Ex(i,j) = (Vmap(i+1,j) - Vmap(i-1,j))*0.5;
        end
        if j == 1
            Ey(i,j) = (Vmap(i,j+1) - Vmap(i,j));
        elseif j == W
            Ey(i,j) = (Vmap(i,j) - Vmap(i,j-1));
        else
            Ey(i,j) = (Vmap(i,j+1) - Vmap(i,j-1))*0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

figure(5);
mesh(Ex);
xlabel('X');
ylabel('Y');
zlabel('Ex');
title('Ex');
%% 
% This is the graph of the horizontal electric field, or the electric field
% that is produced in the 'x' direction. The electric field is the gradient
% of electric potential (Voltage V), and Ex is the gradient in the 'x'
% component.

figure(6);
mesh(Ey);
xlabel('X');
ylabel('Y');
zlabel('Ey');
title('Ey');
%% 
% This is the graph of the vertical electric field, or the electric field
% that is produced in the 'y' direction. The electric field is the gradient
% of electric potential (Voltage V), and Ey is the gradient in the 'y'
% component.

Et = Ex + Ey;
figure(7);
mesh(Et);
xlabel('X');
ylabel('Y');
zlabel('E Total');
title('E Total');
%% 
% This is the graph of the total electric field, which is produced by
% adding components of the electric field in both the 'x' and 'y'
% directions. Therefore, this graph shows the electric field as the
% complete gradient of the 2D electric potential (Voltage V). 

% Sigma
sigma = ones(L,W);
for i = 1:L
    for j = 1:W
        if j <= (W/3) || j >= (W*2/3)
            if i >= (L/3) && i <= (L*2/3)
                sigma(i,j) = 10^-12;
            end
            
        end
    end
end

figure(8);
mesh(sigma);
xlabel('X');
ylabel('Y');
zlabel('Sigma(x,y)');
title('Sigma(x,y)');
%%
% This graph depicts the distribution of the variable sigma over the
% region. The variable sigma is the conductivity of different areas of the
% region, and is the inverse of the resistivity of an area. The graph shows
% the 'bottle-neck' area, where there are 2 'boxes' that contain a
% conductivity that is significantly lower than that of the surrounding
% area. In a physical sense, these 2 boxes represent areas of the region
% that are highly resistive, and which will resist the flow of electric
% current. 

% Current Density
J = sigma .* Et;
figure(9);
mesh(J);
xlabel('X');
ylabel('Y');
zlabel('J(x,y)');
title('J(x,y)');
%%
% The graph shown depicts the distribution of current density (J) over the
% region. From physics, Ohm's Law is J = sigma * E. This formula states
% that the current density is equal to the product of the electric field
% and conductivity. This graph is similar to the previous figure of the
% total electric field (Et), however this graph has 2 'boxes' near the
% middle where there is a near 0 current density. These 2 areas exist
% because the conductivity of those areas is extremely low (meaning that
% the resistance in those areas is extremely high), and this decreases the
% ability of current to flow through these areas. Therefore, the lack of
% current in those areas produces a near 0 current density in those areas.

% Current Flow I = V/R
R = sigma;
for i = 1:L
    for j = 1:W
        if sigma(i,j) == (10^-12)
            R(i,j) = 1 / (10^-12);
        end
    end
end

Current = phi ./ R;
C0 = sum(Current(1,:));
CL = sum(Current(L,:));
c = (C0 + CL) / 2;
figure(10);
mesh(Current);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title(['Avg Current = ', num2str(c), ' Amp']);
%%
% The graph represents the flow of current throughout the region. Another
% form of Ohm's Law is I = V/R, which states that the current is equal to
% the voltage divided by the resistance. This graph is similar the previous
% graph of the distribution of electric potential (V(x,y)), however there
% are 2 'boxes' in the region where the current is 0. These 2 areas are the
% areas with high resistivity. From the formula, a high resistance with
% respect to voltage produces a low (or 0) current. Therefore, these 2
% boxes act as a 'bottle-neck' region which narrows the path of current.
% The average current is stated and is calculated by averaging the current
% at both ends of the region. 

%% Part 2 B
% This section explores the relationship between the mesh size and the
% electric current.

% Inputs
L = 30;
W = 20;

%Grid
x = linspace(0,L);
y = linspace(0,W);
dx = x(2) - x(1);
dy = y(2) - y(1);

% Make matrices
N = L*W;
M = zeros(N,N);
B = zeros(N,1);

%% Part 2 C
% This section investigates the change of the current with respect to the
% change in the width of the bottle-neck region. 

% Interior Points
for i = 2:L-1
    for j = 2:W-1
        n = i + (j-1)*L;
        M(n,n) = -4;
        M(n,n-1) = 1;
        M(n,n+1) = 1;
        M(n,n-L) = 1;
        M(n,n+L) = 1;
        B(n,1) = 0;
    end
end

% Left BC point
i = 1;
for j = 1:W
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 1;
end

% Right BC 
i = L;
for j = 1:W
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 1;
end

% Bottom BC
j = 1;
for i = 1:L
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 0;
end

% Top BC
j = W;
for i = 1:L
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 0;
end

% Solve for potential
phi_vec = M\B;

% Convert vector to 2D array
for i = 1 : 30
    for j = 1 : 20
        n = i + (j-1)*L;
        phi(i,j) = phi_vec(n);
    end
end

% Sigma
sigma1 = ones(L,W);
sigma2 = ones(L,W);
sigma3 = ones(L,W);
sigma4 = ones(L,W);
for i = 1:L
    for j = 1:W
        if j <= (W/4) || j >= (W*3/4)
            if i >= (L/3) && i <= (L*2/3)
                sigma1(i,j) = 10^-12;
            end 
        end
        if j <= (W/3.1) || j >= (W - (W/3.1))
            if i >= (L/3) && i <= (L*2/3)
                sigma2(i,j) = 10^-12;
            end 
        end 
        if j <= (W/2.5) || j >= (W - (W/2.5))
            if i >= (L/3) && i <= (L*2/3)
                sigma3(i,j) = 10^-12;
            end 
        end
        if j <= (W/2.1) || j >= (W - (W/2.1))
            if i >= (L/3) && i <= (L*2/3)
                sigma4(i,j) = 10^-12;
            end 
        end
    end
end

% Current Flow I = V/R
R1 = sigma1;
R2 = sigma2;
R3 = sigma3;
R4 = sigma4;
for i = 1:L
    for j = 1:W
        if sigma1(i,j) == (10^-12)
            R1(i,j) = 1 / (10^-12);
        end
    end
end
for i = 1:L
    for j = 1:W
        if sigma2(i,j) == (10^-12)
            R2(i,j) = 1 / (10^-12);
        end
    end
end
for i = 1:L
    for j = 1:W
        if sigma3(i,j) == (10^-12)
            R3(i,j) = 1 / (10^-12);
        end
    end
end
for i = 1:L
    for j = 1:W
        if sigma4(i,j) == (10^-12)
            R4(i,j) = 1 / (10^-12);
        end
    end
end

Current1 = phi ./ R1;
C01 = sum(Current1(1,:));
CL1 = sum(Current1(L,:));
c1 = (C01 + CL1) / 2;
figure(11);
subplot(2,2,1);
mesh(Current1);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title(['Avg I = ', num2str(c1)]);

Current2 = phi ./ R2;
C02 = sum(Current2(1,:));
CL2 = sum(Current2(L,:));
c2 = (C02 + CL2) / 2;
subplot(2,2,2);
mesh(Current2);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title(['Avg I = ', num2str(c2)]);

Current3 = phi ./ R3;
C03 = sum(Current3(1,:));
CL3 = sum(Current3(L,:));
c3 = (C03 + CL3) / 2;
subplot(2,2,3);
mesh(Current3);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title(['Avg I = ', num2str(c3)]);

Current4 = phi ./ R4;
C04 = sum(Current4(1,:));
CL4 = sum(Current4(L,:));
c4 = (C04 + CL4) / 2;
subplot(2,2,4);
mesh(Current4);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title(['Avg I = ', num2str(c4)]);
%%
% The figure shows 4 different graphs which all show the effect of
% different conductivity regions on the current. The different graphs show
% 'bottle-neck' regions where the path for current to flow through the
% middle region becomes increasingly narrow. However, each graph produces
% the same average current, meaning that both ends of each graph are at the
% same current. This observation would suggest that the narrowing of the
% 'bottle-neck' region has no effect on the current flowing from one end to
% the other. This conclusion can be supported by the fact that the average
% current is calculated from the average of the current at both ends of the
% region, and current is defined as I = V/R. The voltage distribution and
% resistivity of both ends has not changed by the narrowing of the
% 'bottle-neck' in the middle of the region, therefore the currents would
% be identical. 

%% Part 2 D
% This section analyzes the effects of various conductivities (sigma(x,y))
% on the current. 

% Interior Points
for i = 2:L-1
    for j = 2:W-1
        n = i + (j-1)*L;
        M(n,n) = -4;
        M(n,n-1) = 1;
        M(n,n+1) = 1;
        M(n,n-L) = 1;
        M(n,n+L) = 1;
        B(n,1) = 0;
    end
end

% Left BC point
i = 1;
for j = 1:W
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 1;
end

% Right BC 
i = L;
for j = 1:W
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 1;
end

% Bottom BC
j = 1;
for i = 1:L
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 0;
end

% Top BC
j = W;
for i = 1:L
    n = i + (j-1)*L;
    M(n,n) = 1;
    B(n,1) = 0;
end

% Solve for potential
phi_vec = M\B;

% Convert vector to 2D array
for i = 1 : 30
    for j = 1 : 20
        n = i + (j-1)*L;
        phi(i,j) = phi_vec(n);
    end
end

% Sigma
sigma = ones(L,W);
for i = 1:L
    for j = 1:W
        if j <= (W/3) || j >= (W*2/3)
            if i >= (L/3) && i <= (L*2/3)
                sigma(i,j) = 10^-12;
            end
            
        end
    end
end

% Current Flow I = V/R
R1 = sigma;
R2 = sigma;
R3 = sigma;
R4 = sigma;

for i = 1:L
    for j = 1:W
        if sigma(i,j) == (10^-12)
            R1(i,j) = 1 / (10^-3);
            R2(i,j) = 1 / (10^0);
            R3(i,j) = 1 / (10^1);
            R4(i,j) = 1 / (10^3);
        end
    end
end

Current = phi ./ R1;
C0 = sum(Current(1,:));
CL = sum(Current(L,:));
c = (C0 + CL) / 2;
figure(12);
subplot(2,2,1);
mesh(Current);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title(['S = 10^-3 & Avg I = ', num2str(c)]);

Current = phi ./ R2;
C0 = sum(Current(1,:));
CL = sum(Current(L,:));
c = (C0 + CL) / 2;
subplot(2,2,2);
mesh(Current);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title(['S = 10^0 & Avg I = ', num2str(c)]);

Current = phi ./ R3;
C0 = sum(Current(1,:));
CL = sum(Current(L,:));
c = (C0 + CL) / 2;
subplot(2,2,3);
mesh(Current);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title(['S = 10^1 & Avg I = ', num2str(c)]);

Current = phi ./ R4;
C0 = sum(Current(1,:));
CL = sum(Current(L,:));
c = (C0 + CL) / 2;
subplot(2,2,4);
mesh(Current);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title(['S = 10^3 & Avg I = ', num2str(c)]);
%%
% The figure shows 4 different graphs that depict the effect of various
% magnitudes of conductivity on the current. The different graphs show that
% increasing the conductivity of the 2 'boxes' in the 'bottle-neck' region
% proportionally increases the current in those 2 areas. However, all the
% graphs produce the same average current. This observation suggests that
% manipulation of the conductivity of the 'bottle-neck' areas (and
% therefore the resistivity) has no effect on the average current flowing
% from one end to the other. This conclusion can be supported by the fact
% that the average current is calculated as the average of the current on
% both ends of the region. Current is defined as I = V/R, and the change in
% the magnitude of conductivity in the 'bottle-neck' region does not affect
% the voltage distribution or resistivity of either end of the region;
% therefore the currents should be equal. 