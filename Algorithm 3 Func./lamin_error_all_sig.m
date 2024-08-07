 function [f,fd] = lamin_error_all_sig(mu,pars)
% f:  target function of EigOPt evaluated at mu
% fd: derivative of target function of EigOPt evaluated at mu

%Check parameter is real and in the domain range
if abs(imag(mu))>0
    mu=real(mu);
end
for i=1:numel(pars.lb)
   if mu(i)<pars.lb(i)
       mu(i)=pars.lb(i);
   end
   if mu(i)>pars.ub(i)
       mu(i)=pars.ub(i);
   end
end

tol_rho_lambda=eps;

if isfield(pars,'opt_method')
    opt_method = pars.opt_method;    
else    
    opt_method = 1;    
end
% Option for Relative or Absolute Errors
RE=pars.Rel_Error;

ne = pars.ne;
nr=1; %Set r>1 if the problem may have coalescence of eigenvalues  
alength = length(pars.A); 
dim = length(mu); 
thetanew = pars.theta(mu);    

APmu = thetanew(1)*pars.A{1};
for k = 2:alength
    APmu = APmu + thetanew(k)*pars.A{k};    
end

[V,D] = eig(APmu); D=real(D);
[~,inds] = sort(diag(D));
Lu = D(inds(1:nr),inds(1:nr));
%% LB from SirK16---->Algorithm 1 in [1]
options = pars.options;
[~,kappa] = size(pars.eiglist); beta=zeros(kappa,1);
% Computing beta_j(mu)
for j = 1:kappa
    Li = diag(pars.eiglist{j}(1:ne(j)));
    li = pars.eiglist{j}(1);
    lip = pars.eiglist{j}(ne(j)+1);
    beta(j,1) = min(eig((Li-li*eye(ne(j))) - pars.premult{j}*pars.premult{j}'*(Li-lip*eye(ne(j)))));
    b(j,1)=-pars.eiglist{j}(1)-beta(j,1);
end
% Defining LP constrains
M = -pars.thetalist;
%b = -pars.eiglist(1,:)'-beta;
for j = 1:alength
    vecadd = zeros(1,alength);
    vecadd(j) = 1;
    M = [M; vecadd; -vecadd];
    b = [b; pars.lambounds(j,2); -pars.lambounds(j,1)];
end
M = [M;-thetanew]; b = [b;0]; %Only for the sigular value case, i.e. it can't be smaller than zero
%% Linear programming for lower bound
x = linprog(thetanew',M,b,[],[],[],[],options);
eta = thetanew*x;
%% Computing rho(mu)
AAmu=zeros(size(pars.A{1},1));
for j=1:alength
    for jj=1:alength
        AAmu=AAmu+thetanew(j)*thetanew(jj)*pars.Afull{(alength)*(j-1)+jj};
    end
end


temp = real((V(:,inds(1:nr)))'*AAmu*V(:,inds(1:nr)));
rhosq = real(max(eig(temp - Lu*Lu)));

% This control is to avoid numerical instabilities since rho depends from
% A^*AA^*A
if pars.flag_cond==1
    for j=1:numel(pars.mu(1,:))
        if norm(mu-pars.mu(:,j))/norm(mu)<1e-1
            rhosq=0;
        end
    end
end

if abs(rhosq)<tol_rho_lambda
    slb = D(inds(1),inds(1));
else
    slb = min(D(inds(1),inds(1)),eta) - ...
        (2*rhosq)/(abs(D(inds(1),inds(1)) - eta) + sqrt((D(inds(1),inds(1))-eta)^2 + 4*rhosq));
end

if RE
    f = (D(inds(1),inds(1)) - slb)./abs(D(inds(1),inds(1)));
else
    f = (D(inds(1),inds(1)) - slb);
end

if abs(imag(f))>0
    f=real(f);
end

h = 10^-8;

slbh = min(D(inds(1),inds(1)),eta+h) - ...
  (2*rhosq)/(abs(D(inds(1),inds(1)) - eta - h) + sqrt((D(inds(1),inds(1)) - eta - h)^2 + 4*rhosq));
% One should consider that eta constrain depends on mu, this would require
% to solve another LP problem to evaluate the discrete derivative.
% For computational sake we avoid this considering eta depend from (mu)
% only in the theta functions


if opt_method ~= 1
    f = -f;
end

thetanewp = pars.thetap(mu);

for j = 1:dim
    
    APmup = thetanewp(1,j)*pars.A{1};
    for k = 2:alength
        APmup = APmup + thetanewp(k,j)*pars.A{k};    
    end

    if abs(rhosq)<tol_rho_lambda
        fd(j,1) = 0;
    else
        fd(j,1) =  real(V(:,inds(1))'*(APmup*V(:,inds(1)))/norm(V(:,inds(1)))^2) - ((slbh - slb)/h)*thetanewp(:,j)'*x;
        % Compute the derivative of the relative error
        if RE
            fd(j,1) =  (fd(j,1)*(abs(D(inds(1),inds(1))))...
                -((D(inds(1),inds(1))-slb)*sign(D(inds(1),inds(1)))*real(V(:,inds(1))'*(APmup*V(:,inds(1)))/norm(V(:,inds(1)))^2)))....
                /(abs(D(inds(1),inds(1)))^2);
        end
    end

end
end

