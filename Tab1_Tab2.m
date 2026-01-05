%% ------------------ Test Problems related to [1]----------------------------------------------------------
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
profile on %To monitor the performance
%% The affine decompossition
load('Data_Test_Problems/Data1')
theta = @(x)[exp(x(1)), x(2)];
theta_d = @(x)[exp(x(1)) 0; 0 1];
% Define the parametric domain \mathcal{D}:
bounds.lb = [-10; -10];
bounds.ub = [10; 10];
A{1} = A1;
A{2} = A2;
%% Options to run the problem
opts.opt_method = 1;
opts.tol = 1e-8;
opts.RSG_tol = 1e-7;
opts.Rel_Error = 1;
opts.gamma=-4e3;
%% TAB 1
nt=[10,15,20,25,30,35,40];
nEIGOPT=nt.^2; TOTit=numel(nt);
timeSCM1 = zeros(TOTit,1);
timeEIG1 = zeros(TOTit,1);
nResp1 = zeros(TOTit,1); nResp_SSCM1 = zeros(TOTit,1);
MAXerr_Hybrid1 = zeros(TOTit,1); MAXerr_SSCM1 = zeros(TOTit,1);
MEANerr_Hybrid1 = zeros(TOTit,1); MEANerr_SSCM1 = zeros(TOTit,1);
delete(gcp('nocreate'));
parpool;
%Plot settings
plot_opts.sparse=issparse(A{1});
plot_opts.Nmu=[100,100];
plot_opts.log_space=[0,0];
plot_opts.rand = [1,1];
plot_opts.Hermitian=1;
for i=1:numel(nt)
    opts.Nt = nt(i);
    tic
    [~,~,~,Ared_sub,~,~,~,~,~] = subspace_SCMM(A,theta,theta_d,bounds,opts);
    timeSCM1(i)=toc;
    %% RUN ALGORITHM 2 [1]
    opts.EigOptMaxIt= nEIGOPT(i);
    tic
    [~,Ared,~,~,~] = approx_smallesteig_all(A,theta,theta_d,bounds,opts);
    timeEIG1(i)=toc;
    %% Computations
    [~,~, eigvec] = plot_lambdamin_MD(A,theta,bounds,plot_opts);
    plot_opts.sparse=issparse(Ared{1});
    [~,~,eigvec_sub] = plot_lambdamin_MD(Ared,theta,bounds,plot_opts);
    [~,~,eigvec_sub_SCM] = plot_lambdamin_MD(Ared_sub,theta,bounds,plot_opts);

    err_Hybrid = (eigvec_sub-eigvec)./(abs(eigvec_sub));
    err_SSCM = (eigvec_sub_SCM-eigvec)./(abs(eigvec_sub_SCM));

    MAXerr_Hybrid1(i) = max(err_Hybrid);
    MAXerr_SSCM1(i) = max(err_SSCM);
    MEANerr_Hybrid1(i) = mean(err_Hybrid);
    MEANerr_SSCM1(i) = mean(err_SSCM);

end
%% Tab 2
nt = [50,100,150,200,250,300];
timeSCM2 = zeros(TOTit,1);
timeEIG2 = zeros(TOTit,1);
MAXerr_Hybrid2 = zeros(TOTit,1); MAXerr_SSCM2 = zeros(TOTit,1);
MEANerr_Hybrid2 = zeros(TOTit,1); MEANerr_SSCM2 = zeros(TOTit,1);
for i=1:numel(nt)
    opts.Nt = nt(i);
    tic
    [~,~,~,Ared_sub,~,~,~,~,~] = subspace_SCMM(A,theta,theta_d,bounds,opts);
    timeSCM2(i)=toc;
    %% RUN ALGORITHM 2 [1]
    opts.EigOptMaxIt= 100;
    opts.Nt = opts.Nt;
    tic
    [~,Ared,~,~,~] = approx_smallesteig_all_hybrid(A,theta,theta_d,bounds,opts);
    timeEIG2(i)=toc;
    %% Computations
    [~,~, eigvec] = plot_lambdamin_MD(A,theta,bounds,plot_opts);
    plot_opts.sparse=issparse(Ared{1});
    [~,~,eigvec_sub] = plot_lambdamin_MD(Ared,theta,bounds,plot_opts);
    [~,~,eigvec_sub_SCM] = plot_lambdamin_MD(Ared_sub,theta,bounds,plot_opts);

    err_Hybrid = (eigvec_sub-eigvec)./(abs(eigvec_sub));
    err_SSCM = (eigvec_sub_SCM-eigvec)./(abs(eigvec_sub_SCM));

    MAXerr_Hybrid2(i) = max(err_Hybrid);
    MAXerr_SSCM2(i) = max(err_SSCM);
    MEANerr_Hybrid2(i) = mean(err_Hybrid);
    MEANerr_SSCM2(i) = mean(err_SSCM);

end
