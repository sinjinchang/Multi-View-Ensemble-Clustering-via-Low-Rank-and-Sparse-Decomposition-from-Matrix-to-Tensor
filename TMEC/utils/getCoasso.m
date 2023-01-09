function reu = getCoasso(IDX, sign)
    if nargin<2
        sign = 0;
    end
    [n,r] = size(IDX);
    if sign==0
        reu = zeros(n,n);
        for i = 1 : r
            temp = IDX(:,i) + 1;
            Ki = max(temp);
            label =eye(Ki);
            tlabel = label(temp,:);
            tlabel = tlabel(:,2:end);
            tlabel(sum(tlabel,2)==0, :) = nan;
            tS = tlabel*tlabel';
            tS(tS==0) = -1;
            tS(isnan(tS)) = 0;
            reu = reu + tS;
        end
        reu = reu ./ r;
    elseif sign==1
        result = [];
        for i = 1 : r
            temp = IDX(:,i) + 1;
            Ki = max(temp);
            label =eye(Ki);
            tlabel = label(temp,:);
            tlabel = tlabel(:,2:end);
            result = [result tlabel];
        end
        tmp = result*result';
        reu = tmp ./ r;
    end
end
