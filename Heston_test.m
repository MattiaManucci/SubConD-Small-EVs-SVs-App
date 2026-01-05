%% ------------------ Heston test problems related to [1]----------------------------------------------------------
% [1] N. Gugliemi, M. Manucci and E. Mengi, 
% [2] Mustafa Kilic, Emre Mengi and E. Alper Yildirim
clc 
clearvars
close all
addpath(genpath('./eigopt/'))
addpath(genpath('./Algorithm 2 Func./'))
addpath(genpath('./Algorithm 3 Func./'))
addpath(genpath('./Plot_Functions/'))
profile off
%% Inizializations
K=100;           % strike
rf=0;            % foregein rate
flag_D=0;
%Space Disc Sizes
m1=100; m2=100; M1 = m1+1; M2 = m2+1; ev = ones(M2,1); E(M1,:) = ev;
ii = 2:M1; jj = 1:m2; mi = max(size(ii)); mj = max(size(jj));  
%% Bulding the OFFLINE matrix 
t1=10;
[Y,Dvm,Dsm,X,Dss,Iv,Ds,Id,Dvv,Dv,Is,Dsmt,Dvmt,Xt,Yt,Dsst,Dvvt,Dst,Dvt,Ivt,u0,S,V,G,s,ds,v] = Heston_Matrix(m1,m2,K,t1);
%% OFFLINE STUCTURES
A0 = kron(Y*Dvm,X*Dsm);
A1_1 = kron(Y,0.5*X*X*Dss);  A1_2=kron(Iv,X*Ds); A1_3=-(0.5)*Id;
A2_1 = kron((0.5)*Y*Dvv,Is); A2_2=kron((Iv)*Dv,Is); A2_4=-kron(Y*Dv,Is); A2_3=-(0.5)*Id;
% Define the parametric domain \mathcal{D}:
sigma_f=4e-1; sigma_in=1.8e-1; r_f=0.2; r_in=0.001; k_f=3; k_in=1.2; eta_f=0.15; eta_in=0.08; rho_f=0.9; rho_in=0.21;
bounds.lb = [sigma_in,r_in,k_in,eta_in,rho_in]; % lower bound
bounds.ub = [sigma_f,r_f,k_f,eta_f,rho_f]; % upper bound
A{1} = A0;
A{2} = A1_1;
A{3} = A1_2+A2_3+A1_3;
A{4} = A2_1;
A{5} = A2_2;
A{6} = A2_4;
theta = @(x)[x(5)*x(1),1,x(2),x(1).^2,x(3)*x(4),x(3)]; 
theta_d = @(x)[x(5), 0, 0, 0,x(1);
                  0, 0, 0, 0, 0;
                  0, 1, 0, 0, 0;
             2*x(1), 0, 0, 0, 0;
                  0, 0,x(4),x(3),0;
                  0, 0, 1, 0, 0];
delete(gcp('nocreate'));
parpool;
%% Options to run the problem
opts.opt_method = 1;
opts.tol = 1e-4;
opts.RSG_tol = 1e-7;
opts.Rel_Error = 1;
opts.gamma=-4e3;
Nt=6;
%% RUN ALGORITHM 2 [1]
opts.EigOptMaxIt= 100;
opts.Nt = Nt;
profile on
[ff,curerror,~,Ared,mulist,dimU,dimV, pars] = SSSV_EIGOPT(A,theta,theta_d,bounds,opts);
profile off
%% PLOT OPTIONS
plot_opts.sparse=issparse(A{1});
plot_opts.Nmu=[7,7,7,7,7];
plot_opts.log_space=[0,0,0,0,0];
plot_opts.rand = [1,1,1,1,1];
plot_opts.Hermitian=0;
%% Computations
[~,~, eigvec] = plot_lambdamin_MD(A,theta,bounds,plot_opts);
plot_opts.sparse=issparse(Ared{1});
[mu,~,eigvec_sub] = plot_lambdamin_MD(Ared,theta,bounds,plot_opts);
err_Hybrid = (eigvec_sub-eigvec)./(abs(eigvec_sub));
[~,~,eigvec_sub_SSCM] = plot_lambdamin_MD(Ared_SSCM,theta,bounds,plot_opts);
err_Hybrid_SSCM = (eigvec_sub_SSCM-eigvec)./(abs(eigvec_sub_SSCM));
nResp = find(err_Hybrid>opts.tol); 
MAXerr_Hybrid3 = max(err_Hybrid);
MEANerr_Hybrid3 = mean(err_Hybrid);
nResp_SSCM = find(err_Hybrid_SSCM>opts.tol); 
MAXerr_Hybrid3_SSCM = max(err_Hybrid_SSCM);
MEANerr_Hybrid3_SSCM = mean(err_Hybrid_SSCM);

