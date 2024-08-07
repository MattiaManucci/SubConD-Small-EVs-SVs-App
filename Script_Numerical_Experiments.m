%% ------------------ Test Problems related to [1]----------------------------------------------------------
% [1] N. Gugliemi, M. Manucci and E. Mengi, 
% [2] Mustafa Kilic, Emre Mengi and E. Alper Yildirim
clc 
clearvars
close all
%% LEGEND OF NUMERICAL EXAMPLES
fprintf('Digit the corresponding integer number to run the test problem:\n\n');
fprintf('1 for full Hermitian matrix with n=100 and 1 parameter\n');
fprintf('2 for full Hermitian matrix with n=100 and 2 parameters\n');
fprintf('3 for full Hermitian matrix with n=2000 and 1 parameter\n');
fprintf('4 for the Thermal Block with n>7000 and 1 parameters\n');
fprintf('5 for Black-Scholes with n=20000 and 2 parameters\n');
fprintf('6 for the pseudospectrum of Black-Schols with n=1000\n');
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
flag=input('Please insert the corresponding number for the problem you want to run\n');
if flag==1
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

end
if flag==2
     %% The affine decompossition
    load('Data_Test_Problems/Data1')
    theta = @(x)[exp(x(1)), x(2)];
    theta_d = @(x)[exp(x(1)) 0; 0 1];
    % Define the parametric domain \mathcal{D}:
    bounds.lb = [-2; -3];
    bounds.ub = [4; 5];
    A{1} = A1;
    A{2} = A2;
    %% Options to run the problem
    opts.opt_method = 1; 
    opts.tol = 1e-8;
    opts.RSG_tol = 1e-7;
    opts.EigOptMaxIt= 900;
    opts.Rel_Error = 1;
    opts.gamma=-4e3;
    %% Comparison with Subspace SCM
    opts.Nt = 40; %Discret grid---> opts.Nt \times opts.Nt
    [f2,curerror,mu_SSCM,Ared_sub,thetalist,mulist,mulist3,eiglist_SSCM,pars_SSCM] = subspace_SCMM(A,theta,theta_d,bounds,opts);
    %% RUN ALGORITHM 2 [1]
    [f,Ared,mulist2,eiglist,pars] = approx_smallesteig_all(A,theta,theta_d,bounds,opts);
    profile off
    %% Reproduce FIGURES in Example 1 [1]:
    %fig3a
    figure
    semilogy(1:1:(numel(f2(1,:))),f2(1,:),'-ob','LineWidth',LW)
    hold on
    %semilogy(1:1:(numel(f2(2,:))),f2(2,:),'-*r','LineWidth',LW)
    semilogy(1:1:(numel(f)),f,'-*c','LineWidth',LW)
    xlabel('$j$','Interpreter','Latex')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');

    
    %fig3b
    figure
    plot(mulist(1,:),mulist(2,:),'bo')
    hold on; % Keep the current plot
    plot(mulist2(1,:),mulist2(2,:),'rx')
    plot(mulist3(1,:),mulist3(2,:),'b+')
    for i = 1:length(mulist(1,:))
        text(mulist(1,i), mulist(2,i), num2str(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
    hold off;
    hold on
    xlabel('$\mu_1$','Interpreter','Latex')
    ylabel('$\mu_2$','Interpreter','Latex')
    lgd=legend('SCM','EigOpt');
    set(lgd,'Interpreter','Latex');

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');

    
end
if flag==3
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
    opts.num_init_inter=1;
    opts.tol = 1e-8;
    opts.RSG_tol = 1e-7;
    opts.EigOptMaxIt=2000;
    opts.gamma=-4e3;
    opts.Rel_Error = 1;
    %% RUN ALGORITHM 2 [1]
    [f,Ared,mulist,eiglist,pars] = approx_smallesteig_all(A,theta,theta_d,bounds,opts);
    profile off
    %% PLOT OPTIONS
    plot_opts.Nmu=80;
    plot_opts.sparse=issparse(A{1});
    plot_opts.log_space=0;
    plot_opts.mu=mulist;
    plot_opts.Hermitian=1;
    %% Computation
    [muvec,eigvec] = plot_lambdamin(A,theta,bounds,plot_opts); % FP
    [muvec_sub,eigvec_sub] = plot_lambdamin(Ared,theta,bounds,plot_opts); % RP
    %% Reproduce FIGURES in Example 2 [1]:
    %fig4a
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
    
    %fig4b
    figure
    semilogy(muvec(1:plot_opts.Nmu),abs(eigvec_sub(1:plot_opts.Nmu)-eigvec(1:plot_opts.Nmu))./abs(eigvec_sub(1:plot_opts.Nmu)),'-b','LineWidth',LW)
    xlabel('$\mu$','Interpreter','Latex')
    ylabel('$\frac{\lambda^{\mathcal{V}}_{{\min}}(\mu)-\lambda_{{\min}}(\mu)}{|\lambda^{\mathcal{V}}_{{\min}}(\mu)|}$','Interpreter','Latex')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');


end
if flag==4
    %% The affine decompossition
    load('Data_Test_Problems/Thermal Block Data')
    theta = @(x)[1, x];
    theta_d = @(x)[0; 1];
    A1=0.2*A1+0.4*A2+0.6*A3+0.8*A4;
    A{1} = A0;
    A{2} = A1;
    % Define the parametric domain \mathcal{D}:
    bounds.lb = 1e-6;
    bounds.ub = 1e2;
    % Define the affine decomposition for A(\mu)^*A(\mu)
    AA{1} = A0'*A0; AA{2} = A1'*A1; AA{3} = A0'*A1+A1'*A0;
    theta_AA = @(x)[1, (x(1))^2,x(1)];
    theta_d_AA = @(x)[0;2*x(1); 1];
    %% Options to run the problem
    % Estimate Condition Number of the Problem
    opts.opt_method = 1; 
    opts.Nt = 1000;
    opts.tol = 1e-2;
    opts.RSG_tol = 1e-7;
    opts.Rel_Error = 1;
    opts.gamma = -4e3; %Lower bound of the curvature for the target function, see [2], default -4e5.
    opts.EigOptMaxIt=2000;
    %% RUN ALGORITHM 3 [1]
    [f,f2,curerror,mu,Ared,thetalist,mulist,eiglist,dimU,dimV,pars] = ...
        approx_smallestsig_all(A,AA,theta,theta_AA,theta_d,theta_d_AA,bounds,opts);
    profile off
    %% PLOT OPTIONS
    plot_opts.Nmu=100;
    plot_opts.sparse=issparse(A{1});
    plot_opts.log_space=1;
    plot_opts.Hermitian=0;
    plot_opts.mu=mulist;
    %% Computations
    [muvec,eigvec] = plot_lambdamin(A,theta,bounds,plot_opts); 
    plot_opts.sparse=issparse(Ared{1});
    [muvec_sub,eigvec_sub] = plot_lambdamin(Ared,theta,bounds,plot_opts); 
    %% Reproduce FIGURES in Thermal Block example [1]:

    %fig6a
    figure
    [v,ind]=sort(muvec);
    semilogx(muvec(ind),eigvec(ind),'-b','LineWidth',LW)
    hold on
    semilogx(muvec_sub(ind),eigvec_sub(ind),'--r','LineWidth',LW)
    semilogx(mulist,eigvec(plot_opts.Nmu+1:end),'*b','LineWidth',LW)

    xlabel('$\mu$','Interpreter','Latex')
    lgd=legend('$\sigma_{{\min}}(\mu)$','$\sigma^{\mathcal{V}}_{{\min}}(\mu)$', '$\sigma(\mu_j)$','Location','best');
    set(lgd,'Interpreter','Latex');

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');
    
    %fig6b
    figure
    loglog(muvec(1:plot_opts.Nmu),(abs(eigvec_sub(1:plot_opts.Nmu)-eigvec(1:plot_opts.Nmu)))./eigvec_sub(1:plot_opts.Nmu),'-b','LineWidth',LW)
    xlabel('$\mu$','Interpreter','Latex')
    ylabel('$\frac{\sigma^{\mathcal{V}}_{{\min}}(\mu)-\sigma_{{\min}}(\mu)}{\sigma^{\mathcal{V}}_{{\min}}(\mu)}$','Interpreter','Latex')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');

    
    %fig6c
    figure
    semilogy(1:1:(numel(f)),f,'-ob','LineWidth',LW)
    xlabel('$j$','Interpreter','Latex')
    ylabel('$S_r^{(j)}(\mu_j)$','Interpreter','Latex')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');


    %fig6d
    cumulativeSum_dimV = cumsum(dimV);
    figure
    plot(1:1:(numel(dimV)),cumulativeSum_dimV,'-ob','LineWidth',LW)
    xlabel('$j$','Interpreter','Latex')
    ylabel('$dim(\mathcal{V}_j)$','Interpreter','Latex')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');

end
if  flag==5
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
    %Define the affine decomposition for A(\mu)^*A(\mu)
    AA{1} = D'*D; AA{2} = R'*R; AA{3} = D'*R+R'*D;
    theta_AA = @(x)[(0.5*x(1)^2)^2, (x(2))^2,0.5*x(2)*x(1)^2];
    theta_d_AA = @(x)[x(1)^3, 0; 0, 2*x(2); x(2)*x(1), 0.5*x(1)^2];
    %% Options to run the problem
    % Estimate Condition Number of the Problem
    opts.Flag_Cond=1;
    opts.opt_method = 1;
    opts.tol = 1e-4;
    opts.RSG_tol = 1e-7;
    opts.Rel_Error = 1;
    opts.gamma = -4e3;
    %% RUN ALGORITHM 3 [1] 
    [f,f2,curerror,~,Ared,thetalist,mulist,eiglist,dimU,dimV,pars] = ...
        approx_smallestsig_all(A,AA,theta,theta_AA,theta_d,theta_d_AA,bounds,opts);
    profile off
    %% PLOT OPTIONS
    plot_opts.sparse=issparse(A{1});
    plot_opts.Nmu=[40,40];
    plot_opts.log_space=[0,0];
    plot_opts.Hermitian=0;
    %% Computations
    [~,~, eigvec] = plot_lambdamin_MD(A,theta,bounds,plot_opts);
    plot_opts.sparse=issparse(Ared{1});
    [mu,~,eigvec_sub] = plot_lambdamin_MD(Ared,theta,bounds,plot_opts);
    %% Reproduce FIGURES in BS example [1]:
    Relative_Error = reshape(abs(eigvec-eigvec_sub)./(eigvec_sub), plot_opts.Nmu(1), plot_opts.Nmu(2));
    Absolute_Error = reshape(abs(eigvec-eigvec_sub), plot_opts.Nmu(1), plot_opts.Nmu(2));
    %fig7a
    figure
    h=gca;
    surf(mu{2},mu{1},Relative_Error)
    xlabel('$r $','Interpreter','Latex')
    ylabel('$\sigma$','Interpreter','Latex')
    set(h,'zscale','log')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');

    %fig7b
    figure
    semilogy(1:1:(numel(f)),f,'-ob','LineWidth',LW)
    xlabel('$j$','Interpreter','Latex')
    ylabel('$S_r^{(j)}(\mu_j)$','Interpreter','Latex')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');


    %Plot absolute error
    figure
    h=gca;
    surf(mu{1},mu{2},Absolute_Error)
    xlabel('$\mu_1$','Interpreter','Latex')
    ylabel('$\mu_2$','Interpreter','Latex')
    set(h,'zscale','log')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');

    
    for j=1:numel(mulist(2,:))
        Int_Val(j)=eiglist{j}(2,1);
    end
   
    figure
    h=gca;
    surf(mu{2},mu{1}, reshape((eigvec), plot_opts.Nmu(1), plot_opts.Nmu(2)),'LineStyle','-','EdgeColor','r','FaceColor','none')
    hold on
    surf(mu{2},mu{1}, reshape((eigvec_sub), plot_opts.Nmu(1), plot_opts.Nmu(2)),'LineStyle','--','EdgeColor','b','FaceColor','none')
    set(h,'zscale','log')
    plot3(mulist(2,:), mulist(1,:), (Int_Val),'kx');
    xlabel('$r$','Interpreter','Latex')
    ylabel('$\sigma$','Interpreter','Latex')
    legend('1','2','3')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');


end
if flag==6
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
    %Define the affine decomposition for A(\mu)^*A(\mu)
    AA{1} = A{1}'*A{1}; AA{2} = I; AA{3} = I; AA{4} = -A{1}-A{1}';
    AA{5} = 1i*A{1}-1i*A{1}';
    theta_AA = @(x)[1, x(1).^2, x(2).^2, x(1), x(2)];
    theta_d_AA = @(x)[0, 0; 2*x(1), 0;  0, 2*x(2); 1, 0; 0, 1];
    %% Options to run the problem
    opts.Flag_Cond=1;
    opts.opt_method = 1;
    opts.Nt = 1000;
    opts.Flag_Cond=1;
    opts.tol = 1e-6;
    opts.RSG_tol = 1e-8;
    opts.Rel_Error = 1;
    opts.SIN = 1; % flag variable for special initialization of the subspaces, see Sec. 7.3 [1]
    opts.EigOptMaxIt = 2000;
    opts.RSG_tol = -4e3;
    %% RUN ALGORITHM 3 [1]
    % Note: go into the code and uncomment some functions oterwhise it
    % takes too much time.
    [f,f2,curerror,mu,Ared,thetalist,mulist,eiglist,dimU,dimV,pars] = ...
        approx_smallestsig_all(A,AA,theta,theta_AA,theta_d,theta_d_AA,bounds,opts);
    profile off
    %% PLOT OPTIONS
    plot_opts.Nmu=[40,40];
    plot_opts.sparse=issparse(A{1});
    plot_opts.log_space=[1,1];
    plot_opts.Hermitian=0;
    %% Computations
    [~,~, eigvec] = plot_lambdamin_MD(A,theta,bounds,plot_opts);
    plot_opts.sparse=issparse(Ared{1});
    [mu,~,eigvec_sub] = plot_lambdamin_MD(Ared,theta,bounds,plot_opts);
    %% Reproduce FIGURES in BS-Pseudospectrum example [1]
    Relative_Error = reshape(abs(eigvec-eigvec_sub)'./(eigvec_sub'), plot_opts.Nmu(1), plot_opts.Nmu(2));
    Absolute_Error = reshape(abs(eigvec-eigvec_sub)', plot_opts.Nmu(1), plot_opts.Nmu(2));
    %fig8a
    figure
    h=gca;
    surf(mu{2},mu{1},Relative_Error)
    xlabel('$i\Im(z)$','Interpreter','Latex')
    ylabel('$\Re(z)$','Interpreter','Latex')
    set(h,'zscale','log')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');
    
    %fig8b
    figure
    semilogy(1:1:numel(f),f,'-ob','LineWidth',LW)
    hold on
    semilogy(numel(f)+1:1:numel(f2),f2(numel(f)+1:1:numel(f2)),'--or','LineWidth',LW)
    xlabel('$j$','Interpreter','Latex')
    ylabel('$S^{(j)}_r(\mu_j)$','Interpreter','Latex')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');


    %Absolute Error
    figure
    h=gca;
    surf(mu{1},mu{2},Absolute_Error)
    xlabel('$\Re(z)$','Interpreter','Latex')
    ylabel('$i\Im(z)$','Interpreter','Latex')
    set(h,'zscale','log')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');

end