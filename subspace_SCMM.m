function [ff,curerror,mu,Ared,thetalist, ...
                mulist,mulist2, eiglist,pars] = subspace_SCMM(A,theta,thetap,bounds,options)
% Function for the method in [3]
% [3] P. Sirkovic and D. Kressner, SIMAX 2016


Ntrain=options.Nt;
pars.Rel_Error=options.Rel_Error;
            
if isfield(options,'num_init_inter')
    num_init_inter = options.num_init_inter;
else
    num_init_inter = 1;
end

ff=[]; ff2=[];
kappa = length(A);
n=size(A{1},1);

dim = length(bounds.lb);
sp = issparse(A{1});

pars.gamma = -4e4; 
%pars.gamma = -2*max(eig(A{1} + A{2}));
pars.theta = theta;
pars.thetap = thetap;

pars.tol = 10^-6;
pars.minmax = 1;
tol = options.tol;
pars.Rel_Error=options.Rel_Error;


curerror = 10000;
opts.maxit=3000;

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
seed=123; rng(seed);

h = bounds.ub - bounds.lb;

for j = 1:num_init_inter
    mulist(:,j) = bounds.lb + rand.*h;
end

pars.mu=mulist;

P = [];
eiglist = [];

%% Starting Subspace

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
        
        [V,D] = eigs(Amu,ne+1,'smallestreal',opts);
        while (abs(D(1)-D(2))/abs(D(1)))<options.RSG_tol
            ne(j)=ne(j)+1;
            [V,D] = eigs(Amu,ne+1,'smallestreal',opts);
        end
        Pext = V(:,1:ne);
        eiglist = [eiglist diag(D(1:ne+1,1:ne+1))];
       
    else
        [V,D] = eig(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
           if abs((eigAj(i)- eigAj(i+1)))/abs(eigAj(i))>options.RSG_tol %Check if expression is corrected
               ne(j)=i;
               break
           end
 
        end
        Pext = V(:,inds(1:ne)); eiglist = [eiglist diag(D(inds(1:ne+1),inds(1:ne+1)))];
        
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

for j = 1:kappa
    AP{j}=[];
    PA1{j}=[];
    PA2{j}=[];
end
for j = 1:kappa
    for jj = 1:kappa
      AP2{kappa*(j-1)+jj} =[];
    end
end

% Cheb Points
n=Ntrain-1;
i = 0:1:Ntrain;
cx = cos(((2*i + 1)/(2*(n+1)))*pi);
cx=(cx+1)/2;
cx=linspace(0,1,Ntrain);
mu_t=bounds.lb+cx.*h;
    
if dim>1
    mu_1=kron(mu_t(1,:),ones(1,Ntrain));
    mu_2=kron(ones(1,Ntrain),mu_t(2,:));
    mu_t=[mu_1;mu_2]; Ntrain=Ntrain^2;
end

mulist2=mulist;
%% Main loop
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
    
    % Uncomment to reproduce Figure 3
    % pars.itertol=625;
    % [curerror2,mu2,~]=eigopt('lamin_error_all',bounds,pars); 
    % ff2=[ff2,curerror2];  mulist2 = [mulist2 mu2];

    for ii=1:Ntrain
        [fff(ii),~] = lamin_error_all(mu_t(:,ii),pars);
    end
    [f,ind]=max(fff); mu=mu_t(:,ind);
    curerror=f;

    ff=[ff,curerror];
    display(curerror);

    thetanew = theta(mu);
    
    Amu = thetanew(1)*A{1};
    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end
    
    
    if sp==1
        
        [V,D] = eigs(Amu,ne(iter)+1,'smallestreal',opts);
        while (abs(D(1)-D(2))/abs(D(1)))<options.RSG_tol
            ne(iter)=ne(iter)+1;
            [V,D] = eigs(Amu,ne+1,'smallestreal',opts);
        end
        Pext = V(:,1:ne(iter));
        eiglist = [eiglist diag(D(1:ne(iter)+1,1:ne(iter)+1))];
       
    else
        
        [V,D] = eig(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
           if abs((eigAj(i)- eigAj(i+1)))/abs(eigAj(i))>options.RSG_tol %Check if expression is corrected
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
ff=[ff;ff2]
Ared = AP;

return
