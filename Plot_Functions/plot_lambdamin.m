function [muvec, eigvec] = plot_lambdamin(A,theta,bounds,opts)

Nmu=opts.Nmu;
kappa=numel(A);

if opts.log_space==1
    mu=logspace(log10(bounds.lb),log10(bounds.ub),Nmu);
else
    mu=linspace(bounds.lb,bounds.ub,Nmu);
end

muvec = [];
eigvec = [];

opt.maxit=1000;

for i = 1:Nmu   
       
    thetanew = theta(mu(i));
    Amu = thetanew(1)*A{1};

    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};    
    end
    
    if opts.sparse==1
        if opts.Hermitian==1
            dist = eigs(Amu,1,'smallestreal',opt);
        else
            dist = svds(Amu,1,'smallest',opt);
        end
    else
        if opts.Hermitian==1
            [~,D] = eig(Amu);
            dist=min(diag(D));
        else
            D = svd(Amu);
            dist=min(D);
        end
    end

    muvec = [muvec mu(i)];
    eigvec = [eigvec dist(1)];
   
end

if ~isempty(opts.mu)
    mu=opts.mu; Nmu=size(mu,2);
    
for i = 1:Nmu   
       
    thetanew = theta(mu(i));
    Amu = thetanew(1)*A{1};

    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};    
    end

    if opts.sparse==1
        if opts.Hermitian==1
            dist = eigs(Amu,1,'smallestreal',opt);
        else
            dist = svds(Amu,1,'smallest',opt);
        end
    else
        if opts.Hermitian==1
            [~,D] = eig(Amu);
            dist=min(diag(D));
        else
            D = svd(Amu);
            dist=min(D);
        end
    end

    muvec = [muvec mu(i)];
    eigvec = [eigvec dist(1)];
   
end
end

return;
