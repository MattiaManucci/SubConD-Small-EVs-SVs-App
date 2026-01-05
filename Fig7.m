%% ------------------ Fig. 7 related to [1]----------------------------------------------------------
% [1] N. Gugliemi, M. Manucci and E. Mengi, 
% [2] Mustafa Kilic, Emre Mengi and E. Alper Yildirim
clc 
clearvars
close all
addpath(genpath('./eigopt/'))
addpath(genpath('./Algorithm 2 Func./'))
addpath(genpath('./Algorithm 3 Func./'))
addpath(genpath('./Plot_Functions/'))
FS = 15;       % Fontsize
FN = 'times';  % Fontname
LW = 1.6;      % Linewidth
MS = 7.8;      % Markersize
%% Define the BS matrix
L=0; S=200; %Boundaries
m=1e3; %Points of the discretization
ds=(S-L)/(m+1);
s=((L+ds):ds:(S-ds))';
%Diffusion term
B(:,3)=(s.^2)/(2*(ds^2)); B(:,1)=(s.^2)./(2*(ds^2)); B(:,2)=-(s.^2)/(ds^2);
D=spdiags(B,[-1,0,1],m,m)'; D=D*0.5*(0.1^2);
%Advection-Reaction term
B(:,3)=-s./(2*ds); B(:,1)=s./(2*ds); B(:,2)=-1;
R=spdiags(B,[-1,0,1],m,m)'; R=1e-2*R;
I=speye(m);
A{1} = D+R;
A{2} = -I;
A{3} = -1i*I;
theta = @(x)[1,x(1),x(2)];
theta_d = @(x)[0, 0; 1, 0; 0,1];
%Define the parametric domain \mathcal{D}:
bounds.lb = [-0.4; 0];
bounds.ub = [0; 0.5];
%% Options to run the problem
opts.opt_method = 1;
opts.tol = 1e-6;
opts.RSG_tol = 1e-7;
opts.Rel_Error = 1;
opts.gamma=-4e3;
Nt=1;
%% RUN ALGORITHM 2 [1]
opts.EigOptMaxIt= 1000;
opts.Nt = Nt;
[ff,curerror,~,Ared,mulist,dimU,dimV, pars] = SSSV_EIGOPT(A,theta,theta_d,bounds,opts);
%% PLOT OPTIONS
plot_opts.sparse=issparse(A{1});
plot_opts.Nmu=[40,40];
plot_opts.log_space=[0,1];
plot_opts.rand = [0,0];
plot_opts.Hermitian=0;
%% Plots
[~,~, eigvec] = plot_lambdamin_MD(A,theta,bounds,plot_opts);
plot_opts.sparse=issparse(Ared{1});
[mu,~,eigvec_sub] = plot_lambdamin_MD(Ared,theta,bounds,plot_opts);
err_Hybrid = (eigvec_sub-eigvec)./(abs(eigvec_sub));
MaxErr=max(err_Hybrid);
err_Hybrid = reshape(err_Hybrid, plot_opts.Nmu(1), plot_opts.Nmu(2));
err_Hybrid(err_Hybrid==0)=eps;
err_Hybrid = err_Hybrid';

%Fig. 7a
figure
Z_log = log10(abs(err_Hybrid));
surf(mu{1}, mu{2}, Z_log)
shading interp; % Smooth colors across surface
view(2); % Top-down view
c = colorbar;
% Customize color bar to show original scale
c.Ticks = log10([1e-15, 1e-14, 1e-13, 1e-12, 1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3]);  % Adjust based on data range
c.TickLabels = {'1e-15', '1e-14', '1e-13', '1e-12', '1e-11','1e-10','1e-9','1e-8','1e-7','1e-6','1e-5','1e-4','1e-3'};
hold on
plot(mulist(1,:),mulist(2,:),'rx')
%% Define the affine decomposition for A(\mu)^*A(\mu)
AA{1} = full(A{1}'*A{1}); AA{2} = I; AA{3} = I; AA{4} = full(-A{1}-A{1}');
AA{5} = full(1i*A{1}-1i*A{1}');
theta_AA = @(x)[1, x(1).^2, x(2).^2, x(1), x(2)];
theta_d_AA = @(x)[0, 0; 2*x(1), 0;  0, 2*x(2); 1, 0; 0, 1];
opts.Nt=40;
[f_2,Ared_2,mulist_2,eiglist_2,pars_2] = subspace_SCMM(AA,theta_AA,theta_d_AA,bounds,opts);
%% Fig. 7b
figure
semilogy(f_2,'b','LineWidth',LW)
xlabel('j')
ylabel('$H^{(j)}_r(\mu)$',Interpreter='latex')

set(gcf, 'Color', 'w');



