function [ff,ff2,curerror,mu,Ared,thetalist, ...
                mulist, eiglist,dimU,dimV, pars] = approx_smallestsig_all(A,AA,theta,theta_AA,thetap,thetap_AA,bounds,options)


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
if isfield(options,'SIN')
    in_flag=options.SIN;
else
    in_flag=0;
end
if isfield(options,'max_iter')
    max_iter=options.max_iter;
else
    max_iter=1000;
end
            
ff=[];
kappa = length(A);
kappa_2 = length(AA);

pars.theta = theta;
pars.thetap = thetap;
pars.theta_AA = theta_AA;
pars.thetap_AA = thetap_AA;
pars.tol = tol*0.1;
pars.minmax = 1;

curerror = 10000;
opts.maxit=20000; %Max-Iteration of sparse singular-value solver

sp = issparse(A{1});

%This is required for the Linear Prgramming of A(mu)^*A(mu)
if sp==1

    pars.lambounds = [];

    for j = 1:kappa_2
        if j<=kappa
            Bmax = svds(A{j},1,'largest',opts); Bmax=Bmax^2;
            Bmin = svds(A{j},1,'smallest',opts); Bmin=Bmin^2;
        else
            Bmax = eigs(AA{j},1,'largestreal',opts);
            if normest(imag(AA{j}))>eps
                 Bmin = eigs(AA{j},1,'smallestabs',opts);
            else
                 Bmin = eigs(AA{j},1,'smallestreal',opts);
                 %Bmin = eigs(AA{j},1,'smallestabs',opts); %Use this for
                 %test problem 6
            end
        end
        if isnan(Bmin)
            Bmin=0;
        end
        pars.lambounds = [pars.lambounds; Bmin Bmax];
    end

else

    pars.lambounds = [];

    for j = 1:kappa_2
        D = eig(full(AA{j}));
        pars.lambounds = [pars.lambounds; min(D) max(D)];
    end

end
pars.options = optimoptions('linprog','Display','none');
pars.lb = bounds.lb; 
pars.ub = bounds.ub; 
%% Inizialization of the subspaces
seed=123; rng(seed);
h = bounds.ub - bounds.lb;
if in_flag==1

    [EIG]=eigs(A{1},40,'smallestabs');
    EIG=EIG(real(EIG)>=bounds.lb(1));
    EIG=EIG(real(EIG)<=bounds.ub(1)); 
    EIG=EIG(imag(EIG)>=bounds.lb(2));
    EIG=EIG(imag(EIG)<=bounds.ub(2)); 
    mulist=[EIG'; zeros(1,numel(EIG))];
    pars.mu=mulist;
    num_init_inter=numel(EIG);

else

    for j = 1:num_init_inter
        mulist(:,j) = bounds.lb + rand.*h;
    end
    pars.mu=mulist;

end

PV = []; PU = [];
eiglist = [];
%--------------->Uncomment fot the Pseudospectrum test problem
%sp=0;
for j = 1:num_init_inter

    ne(j)=1; 
    mu = mulist(:,j);
    thetanew = theta(real(mu));
    thetalist(j,:) = thetanew;
    thetalist_2(j,:) = theta_AA(mu);

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
        eiglist{j} = [diag(D(1:ne(j)+1,1:ne(j)+1))];
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
        eiglist{j} = flipud(diag(D(inds(1:ne(j)+1),inds(1:ne(j)+1))));
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
dimU(1) = size(PUext,2); dimV(j) = size(PVext,2);

for j=1:numel(pars.eigvecs)
    pars.premult{j}=pars.eigvecs{j}'*PV;
end
iter = num_init_inter+1;
PUold=[]; PVold=[];

for j = 1:kappa
    PA1{j}=[]; PA2{j}=[];
end
for j=1:kappa
     AP{j}=[];
end

for j=1:kappa
    AP_V{j}=[];
end

for j = 1:kappa
    for jj = 1:kappa
      AP2{kappa*(j-1)+jj} =[];
    end
end
%% MAIN LOOP
tol2=tol;
eiglist_2=eiglist;
flag_A3=1;
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
    pars.thetalist = thetalist;
    pars.PV = PV; pars.PU = PU;
    if flag_A3==1
        [curerror,mu,parsout]=eigopt('sigma_error_all',bounds,pars);
        mu=real(mu);
        fprintf('EigOpt at iteration %d required %d function evaluations\n',iter,parsout.nfevals);
        fprintf('Current surrogate error is %g \n',curerror);
        ff=[ff,curerror];
    end
    %% Loop to certify the approximation accuracy for Algorithm 3 of [1] ...
    % (call Algorithm 2 [1] for A^*A with inizialized subspaces provided by Algorithm 3)
    if curerror<tol
        flag_A3=0; ff2=ff;
    end
    if flag_A3==0
        % This could be implemented by updating step by step before entering this loop
        for j = 1:kappa_2
            PA1{j}=AA{j}*PV;
            PA2{j}=PV'*AA{j};
        end

        for j=1:kappa_2
            for jj=j:kappa_2
                AP2{kappa_2*(j-1)+jj} = PA2{jj}*PA1{j};
            end
        end

        for j = 1:kappa_2
            for jj=1:(j-1)
                AP2{kappa_2*(j-1)+jj} = AP2{kappa_2*(jj-1)+j}';
            end
        end

        for j = 1:kappa_2
            AV{j} = PV'*(AA{j}*PV);
        end

        pars.A = AV;
        pars.Afull = AP2;
        pars.APV = PA1;
        pars.ne = ne;
        pars.thetalist = real(thetalist_2);
        for j=1:numel(eiglist)
            pars.eiglist{j} = flipud(eiglist{j}.^2);
        end
        pars.P = PV;
        pars.theta=theta_AA;
        pars.thetap = thetap_AA;
        [curerror2,mu,parsout]=eigopt('lamin_error_all_sig',bounds,pars);
        mu=real(mu);
        fprintf('EigOpt at iteration %d required %d function evaluations\n',iter,parsout.nfevals);
        fprintf('Current surrogate error is %g \n',curerror2);
    
        %% Check the singular-values error
        thetanew = theta_AA(mu);

        Amu = thetanew(1)*AV{1};
        for k = 2:kappa_2
            Amu = Amu + thetanew(k)*AV{k};
        end
        [~,D,~] = eig(Amu);

        if mu==mulist(:,end)
            break
        end
        
        if  pars.Rel_Error
            if curerror2<1
               curerror2=curerror2/2;
               ff2=[ff2,curerror2];
            else
                ff2=[ff2,curerror2];
            end

        else
              curerror2=curerror2/(2*sqrt(D(end,end)));
              ff2=[ff2,curerror2];
        end

        if curerror2<tol2
            break
        else
            curerror=curerror2;
            pars.theta= theta;
            pars.thetap = thetap;
        end
    end
 
    thetanew = theta(mu);
    
    Amu = thetanew(1)*A{1};
    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end
    
    if sp
        
        [U,D,V] = svds(Amu,ne(iter)+1,'smallest',opts);
        while (abs(D(1,1)-D(end,end)))<options.RSG_tol
            ne(iter)=ne(iter)+1;
            [U,D,V] = svds(Amu,ne(iter)+1,'smallest',opts);
        end
        PVext = V(:,(end-ne(iter)+1):end);
        PUext = U(:,(end-ne(iter)+1):end);
        eiglist{iter} = diag(D(1:ne(iter)+1,1:ne(iter)+1));
        
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

        eiglist{iter} =  diag(D(inds(1:ne(iter)+1),inds(1:ne(iter)+1)));
    end

    eiglist_2=eiglist;
    pars.eigvecs{iter} =  PVext;
    
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
    thetalist = [thetalist; thetanew];
    thetalist_2 = [thetalist_2; theta_AA(mu)];

    pars.mu = [pars.mu, mu];
    for j=1:(numel(pars.eigvecs)-1)
        pars.premult{j}=[pars.premult{j},pars.eigvecs{j}'*PVext];
    end
    pars.premult{numel(pars.eigvecs)}=pars.eigvecs{numel(pars.eigvecs)}'*PV;
    iter = iter+1;
end
Ared = AP_V;

return
