%% ------------------ Fig 3 related to [1]----------------------------------------------------------
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
addpath(genpath('./Algorithm 3 Func./'))
addpath(genpath('./Plot_Functions/'))
profile on
%% The affine decompossition
load('Data_Test_Problems/Data2')
theta = @(x)[x.^2, x];
theta_d = @(x)[2*x;1];
A{1} = A1;
A{2} = A2;
% Define the parametric domain \mathcal{D}:
bounds.lb = -2;
bounds.ub = 4;
%% Options to run the problem
opts.opt_method = 1;
opts.num_init_inter=2;
opts.tol = 1e-8;
opts.RSG_tol = 1e-7;
opts.EigOptMaxIt=2000;
opts.gamma=-4e3;
opts.Rel_Error = 1;
%% RUN ALGORITHM 2 [1]
[f,Ared,mulist,eiglist,pars] = approx_smallesteig_all(A,theta,theta_d,bounds,opts);
profile off
profile viewer
%% PLOT OPTIONS
plot_opts.Nmu=80;
plot_opts.sparse=issparse(A{1});
plot_opts.log_space=0;
plot_opts.rand = 0; %Set 1 if the test semple has to be randomly chosen
plot_opts.mu=mulist;
plot_opts.Hermitian=1;
%% Computation
[muvec,eigvec] = plot_lambdamin(A,theta,bounds,plot_opts); % FP
[muvec_sub,eigvec_sub] = plot_lambdamin(Ared,theta,bounds,plot_opts); % RP
%% Reproduce FIGURES in Example 2 [1]:
%fig3a
figure
[v,ind]=sort(muvec);
plot(muvec(ind),eigvec(ind),'-b','LineWidth',LW)
hold on
plot(muvec_sub(ind),eigvec_sub(ind),'--r','LineWidth',LW)
plot(mulist,eigvec(plot_opts.Nmu+1:end),'*b','LineWidth',LW)

xlabel('$\mu$','Interpreter','Latex')
lgd=legend('$\lambda_{{\min}}(\mu)$','$\lambda^{\mathcal{V}}_{{\min}}(\mu)$', '$\lambda(\mu_j)$','Location','best');
set(lgd,'Interpreter','Latex');

set(gca,'Fontname',FN,'Fontsize',FS);
set(gcf, 'Color', 'w');

%fig3b
figure
semilogy(muvec(1:plot_opts.Nmu),abs(eigvec_sub(1:plot_opts.Nmu)-eigvec(1:plot_opts.Nmu))./abs(eigvec_sub(1:plot_opts.Nmu)),'-b','LineWidth',LW)
xlabel('$\mu$','Interpreter','Latex')
ylabel('$\frac{\lambda^{\mathcal{V}}_{{\min}}(\mu)-\lambda_{{\min}}(\mu)}{|\lambda^{\mathcal{V}}_{{\min}}(\mu)|}$','Interpreter','Latex')

set(gca,'Fontname',FN,'Fontsize',FS);
set(gcf, 'Color', 'w');

