function [mu,allCombinations, Eig_Val] = plot_lambdamin_MD(A,theta,bounds,opts)

Nmu=opts.Nmu;
Hermitian=opts.Hermitian;
kappa=numel(A);
p=numel(bounds.ub);
mu=cell(p,1);
opt.maxit=3000;

for i=1:p
    mu{i}=zeros(p,Nmu(i));

    if opts.log_space(i)==1
        mu{i}=logspace(log10(bounds.lb(i)-bounds.lb(i)+1),log10(bounds.ub(i)-bounds.lb(i)+1),Nmu(i));
        mu{i}=mu{i}+(bounds.lb(i))-1;
    else
        mu{i}=linspace(bounds.lb(i),bounds.ub(i),Nmu(i));
    end
end

indices = cell(1, p);
[indices{:}] = ndgrid(1:Nmu(i));
allCombinations = cell2mat(cellfun(@(x) x(:), indices, 'UniformOutput', false));
numCombinations = size(allCombinations, 1);

for i=1:numCombinations
    currentCombination = allCombinations(i, :);
    v=[];
    for j=1:p
        v=[v,mu{j}(currentCombination(j))];
    end

    thetanew = theta(v);
    Amu=thetanew(1)*A{1};
    for k=2:kappa
        Amu=Amu+thetanew(k)*A{k};
    end

    if opts.sparse==1
        if Hermitian
            dist = eigs(Amu,1,'smallestreal',opt);
        else
            dist = svds(Amu,1,'smallest',opt);
        end
    else
        
        if Hermitian
            D = eig(Amu);
            dist=min((D));
        else
            D = svd(Amu);
            dist=min((D));
        end


    end
    Eig_Val(i) = dist(1);
end


end

