function [ff,Ared,mulist,eiglist,pars] = approx_smallesteig_all_hybrid(A,theta,thetap,bounds,options)
%% [1] M. Manucci, E. Mengi and N. Guglielmi, arxiv 2024 
%% Some useful inizializations
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


ff=[];
kappa = length(A);
n=size(A{1},1);

dim = length(bounds.lb);
sp = issparse(A{1});

pars.theta = theta;
pars.thetap = thetap;

pars.tol = tol*10^-1;
pars.minmax = 1;

curerror = 10000;
opts.maxit=30000;

%% Computing the box-constraints for the LP
if sp==1
    pars.lambounds = [];
    for j = 1:kappa
        Bmax = eigs(A{j},1,'largestreal',opts);
        Bmin = eigs(A{j},1,'smallestreal',opts);
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
%% Loop to compute the starting subspace
seed=123; rng(seed);
h = bounds.ub - bounds.lb;
for j = 1:num_init_inter
    mulist(:,j) = bounds.lb + rand.*h;
end
pars.mu=mulist;
P = [];
eiglist = [];

for j = 1:num_init_inter
    
    ne(j)=1;
    mu = mulist(:,j);
    thetanew = theta(mu);    
    thetalist(j,:) = thetanew;
    Amu = thetanew(1)*A{1};

    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};    
    end

    if sp
        i=1;
        [V,D] = eigs(Amu,10*(ne(j)+1),'smallestreal',opts);
        while (abs(D(1,1)-D(end,end)))<options.RSG_tol
            ne(j)=ne(j)+1;
            [V,D] = eigs(Amu,10*ne(j)+1,'smallestreal',opts);
        end
        Pext = V(:,1:ne(j));
        eiglist = [eiglist diag(D(1:ne(j)+1,1:ne(j)+1))];
       
    else
        [V,D] = eig(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
           if abs((eigAj(i)- eigAj(i+1)))>options.RSG_tol %Check if expression is corrected
               ne(j)=i;
               break
           end
        end
        Pext = V(:,inds(1:ne(j))); eiglist = [eiglist diag(D(inds(1:ne(j)+1),inds(1:ne(j)+1)))];
    end
    P = [P Pext];
    pars.eigvecs{j} = Pext;
end

[P,~] = qr(P,0);
for j=1:numel(pars.eigvecs)
    pars.premult{j}=pars.eigvecs{j}'*P;
end
Pext=P;
iter = num_init_inter+1;
Pold=[];
AP=cell(kappa,1); PA1=cell(kappa,1); PA2{j}=cell(kappa,1);
for j = 1:kappa
    AP{j}=[]; PA1{j}=[]; PA2{j}=[];
end
for j = 1:kappa
    for jj = 1:kappa
        AP2{kappa*(j-1)+jj} = [];
    end
end
%% MAIN LOOP FOR GREEDY SELECTION
while (curerror > tol)    
    ne(iter)=1;
    for j = 1:kappa
        if iter==(num_init_inter+1)
            AP_off_diag=[];
        else
            AP_off_diag=Pold'*(A{j}*Pext);
        end
        AP{j} = [AP{j}, AP_off_diag;  AP_off_diag', Pext'*A{j}*Pext ];
    end
    for j = 1:kappa
        PA1{j}=[PA1{j},A{j}*Pext];
        PA2{j}=[PA2{j};Pext'*A{j}];
    end
    if i>1
        for j=1:kappa
            for jj=j:kappa
                AP2{kappa*(j-1)+jj} = [AP2{kappa*(j-1)+jj},  Pold'*A{jj}*PA1{j}; PA2{jj}*A{j}*Pold' , PA2{jj}*PA1{j} ];
            end
        end
    else
        for j=1:kappa
            for jj=j:kappa
                AP2{kappa*(j-1)+jj} = PA2{jj}*PA1{j};
            end
        end
    end
    
    for j = 2:kappa
        for jj=1:(j-1)
            AP2{kappa*(j-1)+jj} = AP2{kappa*(jj-1)+j}';
        end
    end
    
    pars.A = AP;
    pars.Afull = AP2;
    pars.ne = ne;
    pars.thetalist = thetalist;
    pars.eiglist = eiglist;
    pars.P = P;
    [curerror,mu,parsout]=eigopt('lamin_error_all',bounds,pars); 
    fprintf('EigOpt at iteration %d required %d function evaluations\n',iter,parsout.nfevals);
    fprintf('Current surrogate error is %g \n',curerror);
       
    ff=[ff,curerror];
    pars.tol = curerror*1e-1; %Dynamically adjust the exit tolerance of EigOpt
    thetanew = theta(mu);
    Amu = thetanew(1)*A{1};
    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end
    %% Updating the subspace
    if sp==1
        
        [V,D] = eigs(Amu,10*(ne(iter)+1),'smallestreal',opts);
        while (abs(D(1,1)-D(end,end)))<options.RSG_tol
            ne(iter)=ne(iter)+1;
            [V,D] = eigs(Amu,10*(ne(iter)+1),'smallestreal',opts);
        end
        Pext = V(:,1:ne(iter));
        eiglist = [eiglist diag(D(1:ne(iter)+1,1:ne(iter)+1))];
       
    else
        
        [V,D] = eig(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
           if abs((eigAj(i)- eigAj(i+1)))>options.RSG_tol %Check if expression is corrected
               ne(iter)=i;
               break
           end
        end
        Pext = V(:,inds(1:ne(iter))); eiglist = [eiglist diag(D(inds(1:ne(iter)+1),inds(1:ne(iter)+1)))];
        
    end
    
    pars.eigvecs{iter} =  Pext;
    
    
    Pext = Pext - P*(P'*Pext);
    Pext = Pext - P*(P'*Pext);
    Pext = Pext - P*(P'*Pext);
    
    [Pext,~] = qr(Pext,0);
    
    Pold=P;
    P = [P Pext];

    mulist = [mulist mu];
    thetalist = [thetalist; thetanew];
    pars.mu = [pars.mu, mu];
    for j=1:(numel(pars.eigvecs)-1)
        pars.premult{j}=[pars.premult{j},pars.eigvecs{j}'*Pext];
    end
    pars.premult{numel(pars.eigvecs)}=pars.eigvecs{numel(pars.eigvecs)}'*P;

    iter = iter+1;
end
%% Verify over the grid
% Discreet grid
p = dim;
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
curerror = tol*10; fff = zeros(Ntrain,1); i=1;
while (curerror > tol)
    ne(iter)=1;
    for j = 1:kappa
        if iter==(num_init_inter+1)
            AP_off_diag=[];
        else
            AP_off_diag=Pold'*(A{j}*Pext);
        end
        AP{j} = [AP{j}, AP_off_diag;  AP_off_diag', Pext'*A{j}*Pext ];
    end
    for j = 1:kappa
        PA1{j}=[PA1{j},A{j}*Pext];
        PA2{j}=[PA2{j};Pext'*A{j}];
    end
    if i>1
        for j=1:kappa
            for jj=j:kappa
                AP2{kappa*(j-1)+jj} = [AP2{kappa*(j-1)+jj},  Pold'*A{jj}*PA1{j}; PA2{jj}*A{j}*Pold' , PA2{jj}*PA1{j} ];
            end
        end
    else
        for j=1:kappa
            for jj=j:kappa
                AP2{kappa*(j-1)+jj} = PA2{jj}*PA1{j};
            end
        end
    end
    
    for j = 2:kappa
        for jj=1:(j-1)
            AP2{kappa*(j-1)+jj} = AP2{kappa*(jj-1)+j}';
        end
    end
    
    pars.A = AP;
    pars.Afull = AP2;
    pars.ne = ne;
    pars.thetalist = thetalist;
    pars.eiglist = eiglist;
    pars.P = P;


    parfor ii=1:Ntrain
        [fff(ii),~] = lamin_error_all(mu_t(:,ii),pars);
    end
    [f,ind]=max(fff); mu=mu_t(:,ind); 
    %% Delate all parameters that verifies the exit condition
    ind_2=find(fff<tol); 
    mu_t(:,ind_2)=[]; Ntrain=Ntrain-numel(ind_2); fff = zeros(Ntrain,1);
    %% ----------------------
    curerror=f;
    pars.tol = curerror*1e-1; %Dynamically adjust the exit tolerance of EigOpt
    thetanew = theta(mu);
    Amu = thetanew(1)*A{1};
    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end
    %% Updating the subspace
    if sp==1
        
        [V,D] = eigs(Amu,ne(iter)+1,'smallestreal',opts);
        while (abs(D(1,1)-D(end,end)))<options.RSG_tol
            ne(iter)=ne(iter)+1;
            [V,D] = eigs(Amu,ne(iter)+1,'smallestreal',opts);
        end
        Pext = V(:,1:ne(iter));
        eiglist = [eiglist diag(D(1:ne(iter)+1,1:ne(iter)+1))];
       
    else
        
        [V,D] = eig(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
           if abs((eigAj(i)- eigAj(i+1)))>options.RSG_tol %Check if expression is corrected
               ne(iter)=i;
               break
           end
        end
        Pext = V(:,inds(1:ne(iter))); eiglist = [eiglist diag(D(inds(1:ne(iter)+1),inds(1:ne(iter)+1)))];
        
    end
    
    pars.eigvecs{iter} =  Pext;
    
    
    Pext = Pext - P*(P'*Pext);
    Pext = Pext - P*(P'*Pext);
    Pext = Pext - P*(P'*Pext);
    
    [Pext,~] = qr(Pext,0);
    
    Pold=P;
    P = [P Pext];

    mulist = [mulist mu];
    thetalist = [thetalist; thetanew];
    pars.mu = [pars.mu, mu];
    for j=1:(numel(pars.eigvecs)-1)
        pars.premult{j}=[pars.premult{j},pars.eigvecs{j}'*Pext];
    end
    pars.premult{numel(pars.eigvecs)}=pars.eigvecs{numel(pars.eigvecs)}'*P;

    iter = iter+1;

end


Ared = AP;

return