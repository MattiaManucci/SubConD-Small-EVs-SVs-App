function [f,fd] = sigma_error_all(mu,pars)
%% [1] M. Manucci, E. Mengi and N. Guglielmi, arxiv 2024 
%% Output
% f:  target function of EigOPt evaluated at mu
% fd: derivative of target function of EigOPt evaluated at mu
%% --------------------------------------------------------
if isfield(pars,'opt_method')
    opt_method = pars.opt_method;    
else    
    opt_method = 1;    
end

mu=real(mu);
% Option for Relative or Absolute Errors
RE=pars.Rel_Error;

alength = length(pars.AU);
alength_2 = length(pars.AV);

dim = length(mu);

thetanew = pars.theta(mu);   

APUmu = thetanew(1)*pars.AU{1};
APVmu = thetanew(1)*pars.AV{1};

for k = 2:alength
    APUmu = APUmu + thetanew(k)*pars.AU{k};    
end

for k = 2:alength_2
    APVmu = APVmu + thetanew(k)*pars.AV{k}; 
end

[uUB,SUB,vUB]=svds(APVmu,1,'smallest');
SLB=sqrt(vUB'*(APUmu'*APUmu)*vUB);
SLB=real(SLB); 
%% Computation of f
if RE
    f = (SUB-SLB)/SUB;
else
    f = (SUB-SLB);
end
% This can happen numerically when close to machine precision
if f<0
    SUB=SLB; f=0;
end
%% Computation of f'
thetanewJac = pars.thetap(mu); fd=zeros(dim,1);
for j = 1:dim
    
    APVmud = thetanewJac(1,j)*pars.AV{1};
    APUmud = thetanewJac(1,j)*pars.AU{1};
    for k = 2:alength
        APUmud = APUmud + thetanewJac(k,j)*pars.AU{k};   
    end
    for k = 2:alength
        APVmud = APVmud + thetanewJac(k,j)*pars.AV{k};
    end

    dSUB=real(uUB'*APVmud*vUB); dSLB=(1/(2*SLB))*(real(vUB'*(APUmud'*APUmu)*vUB)+real(vUB'*(APUmu'*APUmud)*vUB));

    if RE
        fd(j,1)=((dSUB-dSLB)*SUB-dSUB*(SUB-SLB))./(SUB^2);
    else
        fd(j,1)=dSUB-dSLB;
    end
    
    if opt_method ~= 1
        fd(j,1) = -fd(j,1);
    end
    
end

   
end

