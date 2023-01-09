clear
addpath(genpath(cd))
addpath('./data')
dataset='ORL_mtv.mat';
filename = ['EnsembleClusteringTRPCA_',dataset];
fid1 = fopen(strcat(filename,'.txt'),'a+');
rand('seed',42);

%%
tmp       = load(dataset);
CellS     = tmp.S;
M         = tmp.M;
K         = tmp.K;
numClust  = K;
num_views = M;
gnd       = tmp.gnd;
V         = length(CellS);
N         = size(CellS{1},1);
T         = CellS;
T = cat(3, CellS{:,:});
opts.DEBUG    = 1;
opts.eps      = 1e-6;
opts.max_iter = 300;

%%Grid tuning
veclambda = [0.01:0.01:0.3];
for index = 1:length(veclambda)
    tic;
    lambda = veclambda(index);
    
    %%Initialize
    
    tic;
    Z  = RMSC(T, lambda, opts);
       
    time = toc;

    
    reu = RevisedComstd(gnd, K, Z, 5);
    fprintf('acc=%.5f, rn=%.5f, nmi=%.5f with time %.5f\n', ...
        reu.acc(1), reu.ar(1), reu.nmi(1), time);
    fprintf(fid1,'lambda = %f: ACC=%f,std=%f, NMI=%f, std=%f, AR=%f, std=%f, \r\n',lambda,reu.acc(1),reu.acc(2),reu.nmi(1),reu.nmi(2),reu.ar(1),reu.ar(2));
end