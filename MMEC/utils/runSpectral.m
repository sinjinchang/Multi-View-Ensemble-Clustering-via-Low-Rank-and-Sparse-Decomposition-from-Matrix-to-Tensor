function [idx,D,V] = runSpectral(L,K)

[U,S,~] = svd(L,'econ');
S = diag(S);
r = sum(S>1e-4*S(1));
U = U(:,1:r);S = S(1:r);
U = U*diag(sqrt(S));
U = normr(U);
L = (U*U').^4;

D = diag(1./sqrt(sum(L,2)));
L = D*L*D;
[U,~,~] = svd(L);
V = U(:,1:K);
V = D*V;

n = size(V,1);
M = zeros(K,K,20);
for i=1:size(M,3)
    inds = false(n,1);
    while sum(inds)<K
        j = ceil(rand()*n);
        inds(j) = true;
    end
    M(:,:,i) = V(inds,:);
end
idx = kmeans(V,K,'emptyaction','singleton','start',M,'display','off');


end

