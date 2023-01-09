function IDX = BasicCluster_RFS2(data,r,K,dist,percent)
% Only for K-means. Use the Random Feature Selection strategy.
    [n,p] = size(data);
    IDX = zeros(n,r);
    num = ceil(p*percent);
    
    for i=1:r
        temp = randperm(p);
        temp_data = data(:,temp(1:num));
        l = (sum(temp_data,2)>0);
        temp_idx = litekmeans(temp_data(l>0,:),K,'distance',dist);
        IDX(:,i) = getIdx(l, temp_idx);
    end
    
end
function result = getIdx(l, temp_idx)
    n  = length(l);
    result = zeros(n,1);
    k=1;
    for i = 1 : n
       if l(i)>0
          result(i) = temp_idx(k); 
          k = k+1;
       end
    end
end
