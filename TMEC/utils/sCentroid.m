function C = sCentroid(idx,n,r,sumKi)
    C = zeros(n,sumKi(r+1));
    for l = 1:n
%         C(l,idx(l,:)+(sumKi(1:r))') =  1;
        for i = 1:r
           if idx(l,i)>0
                C(l,idx(l,i)+sumKi(i)) = 1;
           end
        end
    end
end