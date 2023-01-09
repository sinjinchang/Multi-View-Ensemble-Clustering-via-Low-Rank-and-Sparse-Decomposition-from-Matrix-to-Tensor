function IDX = BasicCluster_RPS(Data,r,K,dist,randKi)
% Only for K-means. Use the Random Parameter Selection Strategy.    

    [n,~] = size(Data);
    IDX = zeros(n,r);
    [n1,~] = size(randKi);
    
    if n1>1
        Ki = randKi; % here randKi is the given Ki   
    elseif randKi==1&&sqrt(n)>K
        Ki = randsample(K:ceil(sqrt(n)),r,true); % here Ki is randomized
    else
        Ki = K*ones(r,1); % here Ki is equal to K
    end
    
    for i=1:r
        IDX(:,i) = litekmeans(Data,Ki(i),'distance',dist, 'replicates',1);
    end
    
end