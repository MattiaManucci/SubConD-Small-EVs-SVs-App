function [ff,curerror,mu,Ared, ...
                mulist,dimU,dimV, pars] = SSSV_EIGOPT(A,theta,thetap,bounds,options)
%% [1] M. Manucci, E. Mengi and N. Guglielmi, arxiv 2024 
if isfield(options,'num_init_inter')
    num_init_inter = options.num_init_inter;
else
    num_init_inter = 1;
end
if isfield(options,'EigOptMaxIt')
    pars.itertol = options.EigOptMaxIt; 
else
    pars.itertol = 2000;
end
if isfield(options,'Flag_Cond')
    pars.flag_cond = options.Flag_Cond; 
else
    pars.flag_cond = 0;
end
if isfield(options,'max_l')
    n = options.max_l; 
else
    n = 30;
end
if isfield(options,'Rel_Error')
    pars.Rel_Error = options.Rel_Error;
else
    pars.Rel_Error = 1;
end
if isfield(options,'Rel_Error')
    pars.Rel_Error = options.Rel_Error;
else
    pars.Rel_Error = 1;
end
if isfield(options,'tol')
    tol = options.tol;
else
    tol = 1e-4;
end
if isfield(options,'RSG_tol')
    options.RSG_tol = options.RSG_tol;
else
    options.RSG_tol = 1e-6;
end
if isfield(options,'gamma')
    pars.gamma = options.gamma;
else
    pars.gamma = -4e5; 
end
if isfield(options,'max_iter')
    max_iter=options.max_iter;
else
    max_iter=1000;
end
            
kappa = length(A);

pars.theta = theta;
pars.thetap = thetap;
pars.tol = tol*0.1;
pars.minmax = 1;

curerror = 10000;
opts.maxit=20000; %Max-Iteration of sparse singular-value solver

sp = issparse(A{1});
p = numel(bounds.lb);

%This is required for the Linear Prgramming of A(mu)^*A(mu)
%% Inizialization of the subspaces
seed=123; rng(seed);
h = bounds.ub - bounds.lb; mulist = zeros(p,num_init_inter);
for j = 1:num_init_inter
    mulist(:,j) = bounds.lb + rand.*h;
end
pars.mu=mulist;

PV = []; PU = []; ne = zeros(max_iter+num_init_inter);

for j = 1:num_init_inter

    ne(j)=1; 
    mu = mulist(:,j);
    thetanew = theta(real(mu));

    Amu = thetanew(1)*A{1};

    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end

    if sp

        [U,D,V] = svds(Amu,ne(j)+1,'smallest',opts);
        while (abs(D(1,1)-D(end,end)))<options.RSG_tol
            ne(j)=ne(j)+1;
            [U,D,V] = svds(Amu,ne(j)+1,'smallest',opts);
        end
        PVext = V(:,(end-ne(j)+1):end);
        PUext = U(:,(end-ne(j)+1):end);
    else
        [U,D,V] = svd(full(Amu));
        [svdAj,inds] = sort(diag(D));
        for i=1:n
            if abs((svdAj(i)- svdAj(i+1)))>options.RSG_tol 
                ne(j)=i;
                break
            end
        end
        PVext = V(:,inds(1:ne(j))); PUext = U(:,inds(1:ne(j)));
    end

    PU = [PU PUext];
    PV = [PV PVext];
    
    pars.eigvecs{j} = PVext;

end
sp=1;
[PV,~] = qr(PV,0); [PU,~] = qr(PU,0);
pars.PU = PU; pars.PV = PV;
PUext = PU;              
PVext = PV;

iter = num_init_inter+1;
PUold=[]; PVold=[];

AP = cell(kappa,1);
AP_V = cell(kappa,1);


for j = 1:kappa
    for jj = 1:kappa
      AP2{kappa*(j-1)+jj} =[];
    end
end
%% MAIN LOOP
ff = zeros(max_iter,1); 
dimU = zeros(max_iter,1);
dimV = zeros(max_iter,1);
while ((curerror > tol)&&(iter<max_iter)) 
    
    ne(iter)=1;

    for j = 1:kappa
        AP_V{j} = [AP_V{j},A{j}*PVext];
    end

    for j = 1:kappa
        if iter==(num_init_inter+1)
            AP_off_diag_1=[];
            AP_off_diag_2=[];
        else
            AP_off_diag_1=(PUext'*A{j})*PVold;
            AP_off_diag_2=PUold'*(A{j}*PVext);
        end
        AP{j} = [AP{j}, AP_off_diag_2;  AP_off_diag_1, PUext'*A{j}*PVext ];
    end
   
   
    pars.AV = AP_V; pars.AU = AP;
    pars.ne = ne;
    pars.PV = PV; pars.PU = PU;
    [curerror,mu,parsout]=eigopt('sigma_error_all',bounds,pars);
    fprintf('EigOpt at iteration %d required %d function evaluations\n',iter,parsout.nfevals);
    fprintf('Current surrogate error is %g \n',curerror);
    ff(iter)=curerror;

    thetanew = theta(mu);
    
    Amu = thetanew(1)*A{1};
    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end
    
    if sp==1
        
        [U,D,V] = svds(Amu,ne(iter)+1,'smallest',opts);
        while (abs(D(1,1)-D(end,end)))<options.RSG_tol
            ne(iter)=ne(iter)+1;
            [U,D,V] = svds(Amu,ne(iter)+1,'smallest',opts);
        end
        PVext = V(:,(end-ne(iter)+1):end);
        PUext = U(:,(end-ne(iter)+1):end);
        
    else
        [U,D,V] = svd(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
            if abs((eigAj(1)- eigAj(i+1)))>options.RSG_tol
                ne(iter)=i;
                break
            end
        end
        PVext = V(:,inds(1:ne(iter)));
        PUext = U(:,inds(1:ne(iter)));

    end
    
    PVext = PVext - PV*(PV'*PVext);
    PVext = PVext - PV*(PV'*PVext);
    PVext = PVext - PV*(PV'*PVext);

    PUext = PUext - PU*(PU'*PUext);
    PUext = PUext - PU*(PU'*PUext);
    PUext = PUext - PU*(PU'*PUext);
    
    [PUext,~] = qr(PUext,0);
    [PVext,~] = qr(PVext,0);
    
    PVold=PV; PUold=PU;
    PU = [PU PUext]; PV = [PV PVext];
    dimU(iter)=size(PUext,2); dimV(iter)=size(PVext,2);

    mulist = [mulist mu];

    pars.mu = [pars.mu, mu];
    iter = iter+1;
end
%% Verify over the grid
% Discreet grid
mu_1d = cell(p,1); Ntrain = options.Nt; flag_log = 0;
for i=1:p
    mu_1d{i}=zeros(p,Ntrain);

    if flag_log==1
        mu_1d{i}=logspace(log10(bounds.lb(i)),log10(bounds.ub(i)),Ntrain);
    else
        nn=Ntrain-3;
        j = 0:1:nn;
        cx = [1,cos(((2*j + 1)/(2*(nn+1)))*pi),-1];
        cx=(cx+1)/2;
        mu_1d{i}=bounds.lb(i)+cx.*(bounds.ub(i)-bounds.lb(i));
    end
end

indices = cell(1, p);
[indices{:}] = ndgrid(1:Ntrain);
allCombinations = cell2mat(cellfun(@(x) x(:), indices, 'UniformOutput', false));
Ntrain = size(allCombinations, 1);

mu_t = zeros(p,Ntrain);
for i = 1:Ntrain
    currentCombination = allCombinations(i, :);
    v=[];
    for j=1:p
        v=[v;mu_1d{j}(currentCombination(j))];
    end
    mu_t(:,i) = v;
end
curerror = tol*10; fff = zeros(Ntrain,1);
while (curerror > tol)
    ne(iter)=1;

    for j = 1:kappa
        AP_V{j} = [AP_V{j},A{j}*PVext];
    end

    for j = 1:kappa
        if iter==(num_init_inter+1)
            AP_off_diag_1=[];
            AP_off_diag_2=[];
        else
            AP_off_diag_1=(PUext'*A{j})*PVold;
            AP_off_diag_2=PUold'*(A{j}*PVext);
        end
        AP{j} = [AP{j}, AP_off_diag_2;  AP_off_diag_1, PUext'*A{j}*PVext ];
    end

    pars.AV = AP_V; pars.AU = AP;
    pars.ne = ne;
    pars.PV = PV; pars.PU = PU;

    parfor ii=1:Ntrain
        [fff(ii),~] = sigma_error_all(mu_t(:,ii),pars);
    end
    [ff(iter),ind]=max(fff); mu=mu_t(:,ind(1)); 
    fprintf('Current surrogate error at iteration %d is \n',iter,ff(iter));
    if ff(iter)<tol
        break
    end
    %% Delate all parameters that verifies the exit condition
    ind_2=find(fff<tol); 
    mu_t(:,ind_2)=[]; Ntrain=Ntrain-numel(ind_2); fff = zeros(Ntrain,1);
    %% ----------------------
    thetanew = theta(mu);
    Amu = thetanew(1)*A{1};
    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end
    
    if sp==1
        
        [U,D,V] = svds(Amu,ne(iter)+1,'smallest',opts);
        while (abs(D(1,1)-D(end,end)))<options.RSG_tol
            ne(iter)=ne(iter)+1;
            [U,D,V] = svds(Amu,ne(iter)+1,'smallest',opts);
        end
        PVext = V(:,(end-ne(iter)+1):end);
        PUext = U(:,(end-ne(iter)+1):end);
        
    else
        [U,D,V] = svd(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
            if abs((eigAj(1)- eigAj(i+1)))>options.RSG_tol
                ne(iter)=i;
                break
            end
        end
        PVext = V(:,inds(1:ne(iter)));
        PUext = U(:,inds(1:ne(iter)));

    end
    
    PVext = PVext - PV*(PV'*PVext);
    PVext = PVext - PV*(PV'*PVext);
    PVext = PVext - PV*(PV'*PVext);

    PUext = PUext - PU*(PU'*PUext);
    PUext = PUext - PU*(PU'*PUext);
    PUext = PUext - PU*(PU'*PUext);
    
    [PUext,~] = qr(PUext,0);
    [PVext,~] = qr(PVext,0);
    
    PVold=PV; PUold=PU;
    PU = [PU PUext]; PV = [PV PVext];
    dimU(iter)=size(PUext,2); dimV(iter)=size(PVext,2);

    mulist = [mulist mu];

    pars.mu = [pars.mu, mu];
    iter = iter+1;

end


Ared = AP_V;

return
