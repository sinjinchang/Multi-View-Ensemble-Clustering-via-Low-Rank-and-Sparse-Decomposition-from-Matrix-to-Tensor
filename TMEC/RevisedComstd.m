function reu = RevisedComstd(gnd, K, Z, times)
%     if nargin<4
%         times = 20;
%     end
    acc = zeros(1, times);
    ar = zeros(1, times);
    nmi = zeros(1, times);
    F = zeros(1, times);
    P = zeros(1, times);
    R = zeros(1, times);
    
%     reuH = zeros(times, 3);
    reuS = zeros(times, 6);
    for i=1:times
        ind = runSpectral(Z, K);
        [ACC1, AR1, NMI1, F1, P1, R1] = exMeasure(ind,gnd);
        resultS = [ACC1 AR1 NMI1 F1 P1 R1];
        reuS(i,:) =resultS;
%         ind = kmeans(H,K,'replicates',10);
%         [ACC2, Rn2, NMI2] = exMeasure(ind,gnd);
%         resultK = [ACC2 Rn2 NMI2];
%         reuH(i, :) = resultK;
        
        resultM = [resultS];
%         resultM = max(resultM);
        acc(i) = resultM(1);
        ar(i) = resultM(2);
        nmi(i) = resultM(3);
        f(i) = resultM(4);
        p(i) = resultM(5);
        r(i) = resultM(6);
        
    end
    
    
    
    reu.acc = [mean(acc) std(acc)];
    reu.ar = [mean(ar) std(ar)];
    reu.nmi = [mean(nmi) std(nmi)];
    reu.f = [mean(f) std(f)];
    reu.p = [mean(p) std(p)];
    reu.r = [mean(r) std(r)];
    
    reu.reuS = reuS;
%     reu.reuH = reuH;
end

