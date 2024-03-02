clear 
addpath(genpath(cd))
addpath('./BasicPartitions')
dataset='./BasicPartitions/ORL_mtv.mat';
%filename = ['EnsembleClustering_',dataset];
filename = ['ORL_mtv'];
fid1 = fopen(strcat(filename,'.txt'),'a+');  

tmp       = load(dataset);
CellS     = tmp.S';
M         = tmp.M;
K         = tmp.K;
numClust  = K;
num_views = M;
if isfield(tmp, 'truth')
    gnd = tmp.truth;
elseif isfield(tmp, 'truelabel')
    gnd = tmp.truelabel;
else
    gnd = tmp.gnd;
end 

V         = length(CellS);
N         = size(CellS{1},1);
T         = CellS;
T_tensor  = cat(3, CellS{:,:});
t         = T_tensor(:);


veclambda = [0.01:0.01:0.1 0.11:0.01:0.5];
for index = 1:length(veclambda)
    tic;
    lambda = veclambda(index);


    for k=1:V
        Z{k} = zeros(N,N); 
        Y{k} = zeros(N,N);
        E{k} = zeros(N,N); 
    end
    Z_tensor = cat(3, Z{:,:});
    E_tensor = cat(3, E{:,:});

    y = zeros(N*N*V,1);
    dim1 = N;dim2 = N;dim3 = V;
    myNorm = 'tSVD_1';
    sX = [N, N, V];

    tol = 1e-6;
    iter = 0;
    mu = 10e-3; 
    max_mu = 10e10; 
    pho_mu = 2;
    max_iter=50;

    tic;
    while iter < max_iter
        Zpre=Z_tensor;
        Epre=E_tensor;
        fprintf('----processing iter %d--------\n', iter+1);
        %% update Z
        Y_tensor = cat(3, Y{:,:});
        y = Y_tensor(:);
        e = E_tensor(:);

        [z, objV] = wshrinkObj(t - e + 1/mu*y,1/mu,sX,0,3)   ;
        Z_tensor = reshape(z, sX);
        Z{1}=Z_tensor(:,:,1);
        Z{2}=Z_tensor(:,:,2);
        Z{3}=Z_tensor(:,:,3);

        %% update E
        F = [T{1}-Z{1}+Y{1}/mu;T{2}-Z{2}+Y{2}/mu;T{3}-Z{3}+Y{3}/mu];
        [Econcat] = solve_l1l2(F,lambda/mu);
        E{1} = Econcat(1:size(T{1},1),:);
        E{2} = Econcat(size(T{1},1)+1:size(T{1},1)+size(T{2},1),:);
        E{3} = Econcat(size(T{1},1)+size(T{2},1)+1:end,:);
        E_tensor = cat(3, E{:,:});

        for k=1:V
            Y{k} = Y{k} + mu*(T{k}-Z{k}-E{k});
        end

        %% check convergence
        leq = T_tensor-Z_tensor-E_tensor;
        leqm = max(abs(leq(:)));
        difZ = max(abs(Z_tensor(:)-Zpre(:)));
        difE = max(abs(E_tensor(:)-Epre(:)));
        err = max([leqm,difZ,difE]);
        fprintf('iter = %d, mu = %.3f, difZ = %.3f, difE = %.8f,err=%d\n'...
                , iter,mu,difZ,difE,err);
        if err < tol
            break;
        end

        iter = iter + 1;
        mu = min(mu*pho_mu, max_mu);
    end
     time = toc;
    S = zeros(N,N);
    for k=1:num_views
        S = S + Z{k};
    end

    reu = RevisedComstd(gnd, K, S, S, 5);
    fprintf('acc=%.5f, rn=%.5f, nmi=%.5f with time %.5f\n', ...
        reu.acc(1), reu.rn(1), reu.nmi(1), time);
    fprintf(fid1,'lambda = %f: ACC=%f,std=%f, NMI=%f, std=%f, AR=%f, std=%f,time=%f \r\n',lambda,reu.acc(1),reu.acc(2),reu.nmi(1),reu.nmi(2),reu.rn(1),reu.rn(2), time);
end
