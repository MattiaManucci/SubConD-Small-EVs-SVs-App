%% ------------------ Fig 1 and 2 related to [1]----------------------------------------------------------
% [1] N. Gugliemi, M. Manucci and E. Mengi,
% [2] Mustafa Kilic, Emre Mengi and E. Alper Yildirim
clc
clearvars
close all
%% Plot settings and Paths inclusions
FS = 15;       % Fontsize
FN = 'times';  % Fontname
LW = 1.6;      % Linewidth
MS = 7.8;      % Markersize
addpath(genpath('./eigopt/'))
addpath(genpath('./Algorithm 2 Func./'))
addpath(genpath('./Plot_Functions/'))
%% The affine decompossition
load('Data_Test_Problems/Data1') %load the matrices
theta = @(x)[exp(x), x];
theta_d = @(x)[exp(x); 1]; %theta_d(x) is the derivative of theta(x)
% Define the parametric domain \mathcal{D}:
bounds.lb = -1; % lower bound
bounds.ub = 3; % upper bound
A{1} = A1; A{2} = A2;
%% Options to run the problem
opts.opt_method = 1; %To use eigopt for subproblems
opts.Rel_Error = 1; %Set 0 if you want to run Algorithm 2 of [1] for H(\mu) or 1 for H_r(\mu). Defaoult is 1.
opts.tol = 1e-8; %\varepsilon in Algorithm 2 of [1]. Default is 1e-4.
opts.RSG_tol = 1e-8; %Relative threshold to include eigenvalues coalescence. Default is 1e-6.
opts.EigOptMaxIt=2000; %Maximum number of function evaluation allowed for EigOpt. Default is 2000.
opts.gamma=-4e3; %Lower bound of the curvature for the target function, see [2], default -4e5.
%% RUN ALGORITHM 2 [1]
[f,Ared,mulist,eiglist,pars] = approx_smallesteig_all(A,theta,theta_d,bounds,opts);
profile off
%% PLOTS OPTIONS
plot_opts.Nmu=100; %Number of parameters for the discrete grid
plot_opts.sparse=issparse(A{1}); %Check if problem is sparse or full
plot_opts.log_space=0; %Set 1 to have a logarithmic distribution in the discrete grid
plot_opts.Hermitian=1; %set 0 if the matrix is not Hermitian
plot_opts.rand = 0; %Set 1 if the test semple has to be randomly chosen
plot_opts.mu=[]; %To include specific points in the plot
%% Computation
[muvec,eigvec] = plot_lambdamin(A,theta,bounds,plot_opts); % Large size problem
[muvec_sub,eigvec_sub] = plot_lambdamin(Ared,theta,bounds,plot_opts); % Reduced size problem
%% Reproduce FIGURES in Example 1 [1]:
% Fig1a
figure
plot(muvec,eigvec,'-b','LineWidth',LW)
hold on
plot(muvec_sub(1:end),eigvec_sub(1:end),'--*r','LineWidth',LW)
xlabel('$\mu$','Interpreter','Latex')
lgd=legend('$\lambda_{{\min}}(\mu)$','$\lambda^{\mathcal{V}}_{{\min}}(\mu)$','Location','best');
set(lgd,'Interpreter','Latex');

set(gca,'Fontname',FN,'Fontsize',FS);
set(gcf, 'Color', 'w');

% Fig1b
figure
semilogy(muvec,abs(eigvec_sub-eigvec)./abs(eigvec_sub),'-b','LineWidth',LW)
xlabel('$\mu$','Interpreter','Latex')
ylabel('$\frac{\lambda^{\mathcal{V}}_{{\min}}(\mu)-\lambda_{{\min}}(\mu)}{|\lambda^{\mathcal{V}}_{{\min}}(\mu)|}$','Interpreter','Latex')

set(gca,'Fontname',FN,'Fontsize',FS);
set(gcf, 'Color', 'w');

% Fig1c
figure
semilogy(1:1:(numel(f)),f,'-ob','LineWidth',LW)
xlabel('$j$','Interpreter','Latex')
ylabel('$H^{(j)}_r(\mu_j)$','Interpreter','Latex')

set(gca,'Fontname',FN,'Fontsize',FS);
set(gcf, 'Color', 'w');

%Fig2 (a and b)
[~,~,~,~,~,~, ~,~] = func_plots(A,theta,theta_d,bounds,opts);