function result = getBinaryMatrix(IDX)
    [~,r] = size(IDX);
    result = [];
    for i = 1 : r
        temp = IDX(:,i) + 1;
        Ki = max(temp);
        label =eye(Ki);
        tlabel = label(temp,:);
        result = [result tlabel(:,2:end)];
    end
end
