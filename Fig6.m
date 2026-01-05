%% ------------------ Fig. 6 related to [1]----------------------------------------------------------
% [1] N. Gugliemi, M. Manucci and E. Mengi, 
% [2] Mustafa Kilic, Emre Mengi and E. Alper Yildirim
clc 
clearvars
close all
addpath(genpath('./eigopt/'))
addpath(genpath('./Algorithm 2 Func./'))
addpath(genpath('./Algorithm 3 Func./'))
addpath(genpath('./Plot_Functions/'))
%% Define the BS affine decomposition
L=0; S=200; %Boundaries
m=2e4; %Points of the discretization
ds=(S-L)/(m+1);
s=((L+ds):ds:(S-ds))';
%Diffusion term
B(:,3)=(s.^2)/(2*(ds^2)); B(:,1)=(s.^2)./(2*(ds^2)); B(:,2)=-(s.^2)/(ds^2);
D=spdiags(B,[-1,0,1],m,m)';
%Advection-Reaction term
B(:,3)=-s./(2*ds); B(:,1)=s./(2*ds); B(:,2)=-1;
R=spdiags(B,[-1,0,1],m,m)';
A{1} = D;
A{2} = R;
theta = @(x)[0.5*x(1)^2, x(2)];
theta_d = @(x)[x(1), 0; 0, 1];
%Define the parametric domain \mathcal{D}:
bounds.lb = [0.05; 1e-3];
bounds.ub = [0.25; 2e-2];
%% Options to run the problem
opts.opt_method = 1;
opts.tol = 1e-4;
opts.RSG_tol = 1e-7;
opts.Rel_Error = 1;
opts.gamma=-4e3;
Nt=1;
%% RUN ALGORITHM 2 [1]
opts.EigOptMaxIt= 1000;
opts.Nt = Nt;
profile on
[ff,curerror,~,Ared,mulist,dimU,dimV, pars] = SSSV_EIGOPT(A,theta,theta_d,bounds,opts);
profile off
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
MaxErr=max(abs(err_Hybrid));
fprintf('Maximum error over the test set is %d\n',MaxErr);
err_Hybrid = reshape(err_Hybrid, plot_opts.Nmu(1), plot_opts.Nmu(2));
err_Hybrid(err_Hybrid==0)=eps;
err_Hybrid = err_Hybrid';
% fig 6
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





