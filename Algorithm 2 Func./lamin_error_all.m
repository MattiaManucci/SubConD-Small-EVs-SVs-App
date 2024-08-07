function [f,fd] = lamin_error_all(mu,pars)
% f:  target function of EigOPt evaluated at mu
% fd: derivative of target function of EigOPt evaluated at mu
if isfield(pars,'opt_method')
    opt_method = pars.opt_method;
else
    opt_method = 1;
end

% Option for Relative or Absolute Errors
RE=pars.Rel_Error;

ne = pars.ne;
nr=max(ne);
alength = length(pars.A);
dim = length(mu);
thetanew = pars.theta(mu);
APmu = thetanew(1)*pars.A{1};

for k = 2:alength
    APmu = APmu + thetanew(k)*pars.A{k};
end
[V,D] = eig(APmu);
[~,inds] = sort(diag(D));
%% LB from SirK16---->Algorithm 1 in [1]
Umu = pars.P*V(:,inds(1:nr));
Lu = D(inds(1:nr),inds(1:nr));
options = pars.options;
[~,kappa] = size(pars.eiglist);
%% Computing beta_j(mu)
beta=zeros(kappa,1);
for j = 1:kappa
    Li = diag(pars.eiglist(1:ne(j),j));
    li = pars.eiglist(1,j);
    lip = pars.eiglist(ne(j)+1,j);
    beta(j,1) = min(eig((Li-li*eye(ne(j))) - pars.premult{j}*pars.premult{j}'*(Li-lip*eye(ne(j)))));
end
%% Setting the LP constraints
M = -pars.thetalist;
b = -pars.eiglist(1,:)'-beta;
for j = 1:alength
    vecadd = zeros(1,alength);
    vecadd(j) = 1;
    M = [M; vecadd; -vecadd];
    b = [b; pars.lambounds(j,2); -pars.lambounds(j,1)];
end
%% Linear programming for lower bound
x = linprog(thetanew',M,b,[],[],[],[],options);
eta = thetanew*x;
%% Computation of rho(mu)
AAmu=zeros(size(pars.A{1},1));
for j=1:alength
    for jj=1:alength
        AAmu=AAmu+thetanew(j)*thetanew(jj)*pars.Afull{(alength)*(j-1)+jj};
    end
end

temp = (V(:,inds(1:nr)))'*AAmu*V(:,inds(1:nr));
rhosq = abs(max(eig(temp - Lu*Lu))); %The abs is taken beacuse when rhosq is around machine precision (at iterpolation points for instance) then it could become negative
%% Evaluating f
slb = min(D(inds(1),inds(1)),eta) - ...
    (2*rhosq)/(abs(D(inds(1),inds(1)) - eta) + sqrt((D(inds(1),inds(1))-eta)^2 + 4*rhosq));
if RE
    f = (D(inds(1),inds(1)) - slb)./abs(D(inds(1),inds(1)));
else
    f = (D(inds(1),inds(1)) - slb);
end
if opt_method ~= 1
    f = -f;
end
if abs(imag(f))>0
    display(f)
    fprintf('Complex f should not apper in the Hermitian contest, if imaginary part larger than machine precision then something is going wrong')
end
%% Evaluating df
h = 10^-8;
slbh = min(D(inds(1),inds(1)),eta+h) - ...
    (2*rhosq)/(abs(D(inds(1),inds(1)) - eta - h) + sqrt((D(inds(1),inds(1)) - eta - h)^2 + 4*rhosq));
% One should consider that eta constrain depends on mu, this would require
% to solve another LP problem to evaluate the discrete derivative.
% For computational sake we avoid this considering eta independnt from (mu)
thetanewp = pars.thetap(mu); fd=zeros(dim,1);
for j = 1:dim

    APmup = thetanewp(1,j)*pars.A{1};
    for k = 2:alength
        APmup = APmup + thetanewp(k,j)*pars.A{k};
    end
    % Compute the derivative for the absolute error
    fd(j,1) =  real(V(:,inds(1))'*(APmup*V(:,inds(1)))/norm(V(:,inds(1)))^2) - ((slbh - slb)/h)*thetanewp(:,j)'*x;

    % Compute the derivative for the relative error
    if RE
        fd(j,1) =  (fd(j,1)*(abs(D(inds(1),inds(1))))...
            -((D(inds(1),inds(1))-slb)*sign(D(inds(1),inds(1)))*real(V(:,inds(1))'*(APmup*V(:,inds(1)))/norm(V(:,inds(1)))^2)))....
            /(abs(D(inds(1),inds(1)))^2);
    end

    if opt_method ~= 1
        fd(j,1) = -fd(j,1);
    end

end
end

