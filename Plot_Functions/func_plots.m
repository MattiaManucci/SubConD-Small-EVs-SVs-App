function [ff,curerror,mu,Ared,thetalist, ...
    mulist, eiglist,pars] = func_plots(A,theta,thetap,bounds,options)


FS = 15;       % Fontsize
FN = 'times';  % Fontname
LW = 1.6;      % Linewidth
MS = 7.8;      % Markersize

if isfield(options,'opt_method')

    opt_method = options.opt_method;

else

    opt_method = 1;

end

pars.opt_method = opt_method;


if isfield(options,'eig_method')

    eig_method = options.eig_method;

else

    eig_method = issparse(A{1});

end


if isfield(options,'num_init_inter')
    num_init_inter = options.num_init_inter;
else
    num_init_inter = 4;
end



ff=[];
kappa = length(A);

dim = length(bounds.lb);
% sp = issparse(A1);

pars.gamma = -4e3;
pars.theta = theta;
pars.thetap = thetap;
% pars.gamma = -2*max(eig(A1 + A2));
pars.tol = 10^-6;
pars.minmax = 1;
tol = 1*10^-2;
ne = 1;
curerror = 10000;

opts.maxit=3000;

% Compute the minimum and maximum eigenvalues. If the matrices are sparse
% it uses iterative algorithms, otherwise direct algorithm.
if eig_method

    pars.lambounds = [];

    for j = 1:kappa

        Bmax = eigs(A{j},1,'largestreal',opts);
        if isnan(Bmax)
            Eigen=eig(full(A{j}));
            Bmax=max(Eigen); Bmin=min(Eigen);
        else
            Bmin = eigs(A{j},1,'smallestreal',opts);
            if isnan(Bmin)
                Eigen=eig(full(A{j}));
                Bmin=min(Eigen);
            end
        end

        pars.lambounds = [pars.lambounds; Bmin Bmax];

    end


else

    pars.lambounds = [];

    for j = 1:kappa

        D = eig(A{j});

        pars.lambounds = [pars.lambounds; min(D) max(D)];

    end

end



pars.options = optimoptions('linprog','Display','none');

%mulist = [bounds.lb (bounds.lb+bounds.ub)/2  bounds.ub];

h = bounds.ub - bounds.lb;

for j = 1:num_init_inter
    mulist(:,j) = bounds.lb + 0.2*j.*h;
end

% Just to try
pars.mu=mulist;

P = [];
eiglist = [];
pars.eigvecs=[];
for j = 1:num_init_inter

    mu = mulist(:,j);

    thetanew = theta(mu);
    thetalist(j,:) = thetanew;

    Amu = thetanew(1)*A{1};

    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end
    % Compute first elements of the subspace for the starting selected
    % parameter
    if eig_method

        [V,D] = svds(Amu,ne,'smallestnz',opts);
        if isnan(diag(D))
            [V,D] = svds(Amu,ne,'smallestnz',opts);
            [d,ind] = sort(diag(D));
            D=d(ind(1:ne+1)); D=diag(D); V=V(:,ind(1:ne+1));
        end
        Pext = V(:,1:ne);

        eiglist = [eiglist diag(D(1:ne,1:ne))];
    else
        [V,D] = eig(Amu);
        [~,inds] = sort(diag(D));
        Pext = V(:,inds(1:ne));

        eiglist = [eiglist diag(D(inds(1:ne+1),inds(1:ne+1)))];
    end
    pars.eigvecs = [pars.eigvecs,Pext];
    if j>1
        Pext=Pext-P*(P'*Pext);
    end
    [Pext,~] = qr(Pext,0);
    P=[P,Pext];

end

%[P,~] = qr(P,0);
% pars.eigvecs = P;
pars.premult=pars.eigvecs'*P;

iter = 0;



% Project the problem
for j = 1:kappa
    AP{j} = P'*(A{j}*P);
end

for j = 1:kappa
    PA1=A{j}*P;
    for jj=j:kappa
        PA2=P'*(A{jj})';
        AP2{kappa*(j-1)+jj} = PA2*PA1;
    end
end

for j = 2:kappa
    for jj=1:(j-1)
        AP2{kappa*(j-1)+jj} = AP2{kappa*(jj-1)+j};
    end
end

pars.A = AP;
pars.Afull = A;
pars.ne = ne;
pars.thetalist = thetalist;
pars.eiglist = eiglist;
pars.P = P;



muc=bounds.lb+linspace(0,1,300).*h;
muc=[mulist,muc];
muc=sort(muc);

for i=1:(300+numel(mulist))



    thetanew = theta(muc(:,i));


    Amu = thetanew(1)*A{1};

    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end

    lambda(:,i)=eigs(Amu,2,'smallestreal',opts);
    

    [LBU(i),LB(i),UB(i)] = LB_ULB_UB(muc(i),pars);

end

%% Plots


figure
plot(muc,lambda(1,:),'-b','LineWidth',LW)
hold on
%plot(muc,lambda(2,:),'-c','LineWidth',LW)
plot(muc,UB,'--r','LineWidth',LW)
plot(muc,LB,'--.k','LineWidth',LW)
plot(muc,LBU,'-.c','LineWidth',LW)
plot(mulist,eiglist(1,:),'xr','LineWidth',LW*2)

xlabel('$\mu$','Interpreter','Latex')
lgd=legend('$\lambda_{{\min}}(\mu)$','$\lambda^{\mathcal{V}}_{{UB}}(\mu)$',...
           '$\lambda_{{LB}}(\mu)$ for SCM','$\lambda_{{LB}}(\mu)$ for improved SCM','($\mu_i$, $\lambda_i$)','Location','best');
set(lgd,'Interpreter','Latex');

set(gca,'Fontname',FN,'Fontsize',FS);
set(gcf, 'Color', 'w');

Ared=pars.A;


return
