clear; close all; clc;

% Code for the paper "On data usage and predictive behavior of data-driven 
% predictive control with 1-norm regularization" by Manuel Klädtke and 
% Moritz Schulze Darup submitted to CDC/L-CSS 2025

% @ 2025 Manuel Klädtke 

% This code uses MPT3 for visualization of Polyhedron and solution of mpQPs
% M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari. Multi-Parametric 
% Toolbox 3.0. In Proc. of the European Control Conference, pages 502–510, 
% Zurich, Switzerland, July 17–19 2013. http://control.ee.ethz.ch/~mpt

% TU colors
tuGreen =[116,181,0]/255;
tuOrange=[255,153,0]/255;
tuGray  =[86,86,86]/255;
%% Figure 1: Visualization of data sets

% Uncomment the following lines and comment line 26 to generate a new
% data set. Otherwise, load the one used in the paper.
% ell = 8;
% D = randn(2, ell);
% save('fig1_data', "ell", "D");

load('fig1_data')

D_pm = [D, -D]; % Mirrored data set
% Compute extreme points via linear programming. Quickhull (convhulln) is
% very slow for large ell. Perhaps due to symmetry of D_pm(?)
D_hull = [];
for i = 1:size(D,2)
    [~, fval] = linprog(ones(1,2*ell), [], [], D_pm, D(:,i), zeros(2*ell,1), inf(2*ell,1));
    if fval >= 1-1e-6 % Small tolerance if points are (close to) colinear
        D_hull = [D_hull, D(:,i)];
    end
end
D_pm_hull = [D_hull, -D_hull];
P = Polyhedron('V', D_pm_hull');

% Plot results
figure
subplot(1,2,1)
P.plot('alpha', 0, 'EdgeAlpha', 0.5)
hold on; box on
plot(D(1,:), D(2,:), 'x', 'Color', tuGreen)
plot(-D(1,:), -D(2,:), 'x', 'Color', tuOrange)
title('(a)')

subplot(1,2,2)
plot(D_hull(1,:), D_hull(2,:), 'x', 'Color', tuGreen)
hold on; box on
P.plot('alpha', 0, 'EdgeAlpha', 0)
plot(-D_hull(1,:), -D_hull(2,:), 'x', 'Color', tuOrange)
title('(b)')

%% Figure 3: Symmetry

% Data-generating dynamics
f = @(x, u) 2*x.^2+2*u.^2-1;

% Plot data to visualize f(x,u) in [-1, 1]x[-1, 1]
x = linspace(-1, 1, 100);
u = linspace(-1, 1, 100);
[XX, UU] = meshgrid(x, u);
FF = f(XX, UU);

% Uncomment the following lines and comment line 74 to generate a new
% data set. Otherwise, load the one used in the paper.
% ell = 20;
% Z = rand(2, ell)*2-1;
% save('fig23_data', "ell", "Z");

load('fig23_data')

D = [Z; f(Z(1,:), Z(2,:))];
D_pm = [D, -D];
D_hull = [];
% Same extreme point calculation via linear programming as in Figure 1
for i = 1:size(D,2)
    [~, fval] = linprog(ones(1,2*ell), [], [], D_pm, D(:,i), zeros(2*ell,1), inf(2*ell,1));
    if fval >= 1-1e-6
        D_hull = [D_hull, D(:,i)];
    end
end
D_pm_hull = [D_hull, -D_hull];
P = Polyhedron('V', D_pm_hull');

% DPC parameters
Q = 1;
lambda = 100;

% Implicit predictor is a PWA function of xi and u, which we compute via MPT3
% Rewrite the mpQP according to following syntax:
% s = [y; a_pm], th = z;
% J(th) = min 0.5*s'*H*s + (pF*th+f)'*s + th'*Y*th + C*th + c
%         s.t.  A*s <= b  + pB*th
%               Ae*s = be + pE*th
%               lb  <= s <= ub

H = blkdiag(Q, zeros(2*ell));
f = [0; lambda*ones(2*ell,1)];
Ae = [[0; 0], D_pm(1:2,:);
        -1, D_pm(3,:)];
be = [0;0;0];
pE = [1, 0;
    0, 1;
    0, 0];
lb = [-inf; zeros(2*ell,1)];
ub = inf(2*ell+1,1);
% Constrain the parameter space to only calculate the visualized domain 
% [-2, 2] x [-2, 2]
A = zeros(4, 1+2*ell);
pB = -[1, 0;
        -1, 0;
        0, 1;
        0, -1];

b = 2*[1; 1; 1; 1];

mpQP = Opt('H', H, 'f', f, 'A', A, 'b', b, 'pB', pB, 'Ae', Ae, 'be', be, 'pE', pE, 'lb', lb, 'ub', ub);
sol_mpQP = mpQP.solve();

% primal optimizer is a PWA function s(th)=s(z). Retrieve only the first
% element y(z) of the optimizer s(z) and add it as a PWA function for
% visualization
for i = 1:length(sol_mpQP.xopt.Set)
    tmp_primal = sol_mpQP.xopt.Set(i).getFunction('primal');
    sol_mpQP.xopt.Set(i).addFunction(AffFunction(tmp_primal.F(1,:),tmp_primal.g(1,:)), 'y')
end

% Importantly, note that we used D_pm instead of the preprocessed data
% D_pm_hull. The interested reader may check that the PWA solution indeed
% does not use any columns of D_pm that are not present in D_pm_hull as
% well.
% Alternatively, one can replace D_pm with D_pm_hull in Ae and check that
% the results indeed match. This also requires to redefine 
% ell = size(D_pm_hull, 2) to match the new (reduced) amount of
% optimization variables.

% Figure 3
figure
% Original data and data generating dynamics
subplot(1,2,1)
surf(XX, UU, FF, 'EdgeColor', 'none', 'FaceColor',tuGreen, 'FaceAlpha',0.3)
hold on; grid on
surf(XX, UU, FF, 'EdgeColor', 'none', 'FaceColor',tuGreen, 'FaceAlpha',0.3)
scatter3(D(1,:), D(2,:), D(3,:), 'Marker', 'x')
view(-26.5, 16)
xlabel('x_0'); ylabel('u'); zlabel('x')
title('(a)')

% Preprocessed data and implicit predictor
subplot(1,2,2)
sol_mpQP.xopt.fplot('y', 'Color', tuOrange, 'alpha', 0.3)
xlim([-1, 1]); ylim([-1, 1])
hold on; grid on
scatter3(D_hull(1,:), D_hull(2,:), D_hull(3,:), 'Marker', 'x')
view(-26.5, 16)
xlabel('x_0'); ylabel('u'); zlabel('x')
title('(b)')

%% Figure 2: Scaling

% Repeat same computation of the implicit predictor as in Figure 3 to cover
% the four cases:
% 1) lambda = 100, Domain [-1, 1]^2
% 2) lambda = 100, Domain [-2, 2]^2
% 3) lambda = 50, Domain [-1, 1]^2
% 4) lambda = 50, Domain [-2, 2]^2

% Case 1) and 2) lambda = 100, Domain [-1, 1]^2 and [-2, 2]^2
% Already computed for Figure 3
sol_mpQP_12 = sol_mpQP;

% Case 3) and 4) lambda = 50, Domain [-1, 1] and [-2, 2]^2
lambda = 50;
f = [0; lambda*ones(2*ell,1)];

mpQP = Opt('H', H, 'f', f, 'A', A, 'b', b, 'pB', pB, 'Ae', Ae, 'be', be, 'pE', pE, 'lb', lb, 'ub', ub);
sol_mpQP_34 = mpQP.solve();

% Figure 2
figure
t = tiledlayout(2, 2);
subplots = [];

% 1) lambda = 100, Domain [-1, 1]^2
subplots = [subplots, nexttile];
plot(sol_mpQP_12.xopt.Set, 'Color', tuGreen, 'Alpha', 0.5)
hold on; box on
xlim([-1, 1]); ylim([-1, 1]) % Zoom in
plot(D_hull(1,:), D_hull(2,:), 'x')
xlim([-1,1]); ylim([-1,1])
title('(a)')

% 3) lambda = 50, Domain [-1, 1]^2
subplots = [subplots, nexttile];
plot(sol_mpQP_34.xopt.Set, 'Color', tuOrange, 'Alpha', 0.5)
hold on; box on
xlim([-1, 1]); ylim([-1, 1]) % Zoom in
plot(D_hull(1,:), D_hull(2,:), 'x')
xlim([-1,1]); ylim([-1,1]) 
title('(b)')

% 2) lambda = 100, Domain [-2, 2]^2
subplots = [subplots, nexttile];
plot(sol_mpQP_12.xopt.Set, 'Color', tuGreen, 'Alpha', 0.5)
hold on; box on
plot(D_hull(1,:), D_hull(2,:), 'x')
xlim(2*[-1,1]); ylim(2*[-1,1]) 

% 4) lambda = 50, Domain [-2, 2]^2
subplots = [subplots, nexttile];
plot(sol_mpQP_34.xopt.Set, 'Color', tuOrange, 'Alpha', 0.5)
hold on; box on
plot(D_hull(1,:), D_hull(2,:), 'x')
xlim(2*[-1,1]); ylim(2*[-1,1]) 

xlabel(t,'x_0')
ylabel(t,'u')
