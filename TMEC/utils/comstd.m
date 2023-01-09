function reu = comstd(gnd, K, Z, H, times)
    if nargin<5
        times = 20;
    end
    acc = zeros(1, times);
    nmi = zeros(1, times);
    rn = zeros(1, times);
    
    reuH = zeros(times, 3);
    reuS = zeros(times, 3);
    for i=1:times
        ind = runSpectral(Z, K);
        [ACC1, Rn1, NMI1] = exMeasure(ind,gnd);
        resultS = [ACC1 Rn1 NMI1];
        reuS(i,:) =resultS;
        ind = kmeans(H,K,'replicates',10);
        [ACC2, Rn2, NMI2] = exMeasure(ind,gnd);
        resultK = [ACC2 Rn2 NMI2];
        reuH(i, :) = resultK;
        
        resultM = [resultS; resultK];
        resultM = max(resultM);
        acc(i) = resultM(1);
        rn(i) = resultM(2);
        nmi(i) = resultM(3);
    end
    reu.acc = [mean(acc) std(acc)];
    reu.rn = [mean(rn) std(rn)];
    reu.nmi = [mean(nmi) std(nmi)];
    
    reu.reuS = reuS;
    reu.reuH = reuH;
end

