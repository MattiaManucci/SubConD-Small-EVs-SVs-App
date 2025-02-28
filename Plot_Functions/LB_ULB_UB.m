function [LBU,LB,UB] = LB_ULB_UB(mu,pars)
%
% Auxiliary Routine; called by eigopt
%
% UPDATED LOWER BOUNDS

if isfield(pars,'opt_method')
    opt_method = pars.opt_method;    
else    
    opt_method = 1;    
end


nr = 1;
ne = pars.ne;
alength = length(pars.A);

dim = length(mu);


thetanew = pars.theta(mu);    

APmu = thetanew(1)*pars.A{1};
for k = 2:alength
    APmu = APmu + thetanew(k)*pars.A{k};    
end


[V,D] = eig(APmu);
[~,inds] = sort(diag(D));
Umu = pars.P*V(:,inds(1:nr));
Lu = D(inds(1:nr),inds(1:nr));
% keyboard

options = pars.options;

[~,kappa] = size(pars.eiglist);
for j = 1:kappa

    Li = diag(pars.eiglist(1:ne,j));
    li = pars.eiglist(1,j);
    lip = pars.eiglist(ne+1,j);
     ViUmu = pars.premult(j,:)*V(:,inds(1:nr));
     beta(j,1) = min(eig((Li-li*eye(ne)) - ViUmu*ViUmu'*(Li-lip*eye(ne))));
end

% [V,D] = eig(APmu);
% [~,inds] = sort(diag(D));
% Umu = pars.P*V(:,inds(1:nr));
% Lu = D(inds(1:nr),inds(1:nr));
% % keyboard
% 
% options = pars.options;
% 
% [~,kappa] = size(pars.eiglist);
% for j = 1:kappa
% 
%     Li = diag(pars.eiglist(1:ne,j));
%     li = pars.eiglist(1,j);
%     lip = pars.eiglist(ne+1,j);
%     ViUmu = pars.eigvecs(:,(j-1)*ne+1:j*ne)'*Umu;
%     beta(j,1) = min(eig((Li-li*eye(ne)) - ViUmu*ViUmu'*(Li-lip*eye(ne))));
% 
% end

% keyboard

M = -pars.thetalist;
b = -pars.eiglist(1,:)'-beta;

for j = 1:alength
    vecadd = zeros(1,alength);
    vecadd(j) = 1;
    M = [M; vecadd; -vecadd];
    b = [b; pars.lambounds(j,2); -pars.lambounds(j,1)];
end

                        
% keyboard
                        
%% cvx_begin quiet
%%   variable x(2)
%%   dual variables y z
%%   minimize( theta * x )
%%   subject to
%%      y : M * x <= b;
%% cvx_end
 
 
% keyboard
% f = min(eig(APmu)) - cvx_optval;

x = linprog(thetanew',M,b,[],[],pars.lambounds(:,1),pars.lambounds(:,2),options);
%x=pinv(M)*b;

eta = thetanew*x;

% AAmu=zeros(size(pars.A{1},1));
% for j=1:alength
%     for jj=1:alength
%         AAmu=AAmu+thetanew(j)*thetanew(jj)*pars.Afull{(alength)*(j-1)+jj};
%     end
% end

Amu = thetanew(1)*pars.Afull{1};
for j = 2:alength
    Amu = Amu + thetanew(j)*pars.Afull{j};
end

temp = Amu*Umu;
rhosq = max(eig(temp'*temp - Lu*Lu));

% temp = (V(:,inds(1:nr)))'*AAmu*V(:,inds(1:nr));
% rhosq = max(eig(temp - Lu*Lu));
LBU = min(D(inds(1),inds(1)),eta) - ...
  (2*rhosq)/(abs(D(inds(1),inds(1)) - eta) + sqrt((D(inds(1),inds(1))-eta)^2 + 4*rhosq));




M = -pars.thetalist;
b = -pars.eiglist(1,:)'-0*beta;

% for j = 1:alength
%     vecadd = zeros(1,alength);
%     vecadd(j) = 1;
%     M = [M; vecadd; -vecadd];
%     b = [b; pars.lambounds(j,2); -pars.lambounds(j,1)];
% end

LB = thetanew*linprog(thetanew',M,b,[],[],pars.lambounds(:,1),pars.lambounds(:,2),options);


UB = D(inds(1),inds(1));

%LBU=max(LB,LBU);

end

