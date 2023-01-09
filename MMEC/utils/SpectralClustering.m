function groups = SpectralClustering(W, n, varargin) 
%SPECTRALCLUSTERING Spectral clustering.
%   GROUPS = SpectralClustering(W, N, varargin) clusters data into N groups
%   with affinity W.

% Input Arguments
% W                -- symmetric affinity matrix.
% n                -- number of groups.
% 'Start'          -- initial group for k-means.
%     'sample'(default):
%     the same as the k-means
% 'MaxIter'        -- maximum number of iterations for KMeans 
%     1000(default):
%     positive integer
% 'Replicates'     -- number of replications for KMeans
%     20(default):
%     positive integer
% 'Eig_Solver'     -- eig function of matlab
%     eig(default):
%     eigs

% Algorithm (Ncut)
% Ncut: min_{A_i} sum_i cut(A_i, bar(A)_i) / vol(A_i),
% reformulate as min_{A_i} \sum_i h_i' L h_i 
% s.t. h_ij = 1 / sqrt(Vol(A_j)), if v_i \in A_j
% relax as min_H tr(H' L H) s.t. H' D H = I
% solution is H = D^{-0.5} eig(L_sym).
% in this version, do another row normalization of H.

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

% check data
if ~issymmetric(W)
    error(['In ''' mfilename ''': affinity matrix is not symmetric'])
end
% define defaults options
% Set default 
vararg = {'Start', 'sample', ...
          'MaxIter', 1000, ...
          'Replicates', 20, ...
          'Eig_Solver', 'eig'};
% Overwrite by input
vararg = vararginParser(vararg, varargin);
% Generate variables
for pair = reshape(vararg, 2, []) % pair is {propName;propValue}
   eval([pair{1} '= pair{2};']);
end

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
% The computation is equivalent to:
% - compute the largest eigenvectors of D^{-1} W
% - normalize the rows of the resultant matrix
% - then apply kmeans to the rows.
if strcmpi(Eig_Solver, 'eig')
    [V, D] = eig( cnormalize(full(W), 1)' );
    [~, ord] = sort(real(diag(D)), 'descend');
    kerN = V(:, ord(1:n));
    clear V D;
elseif strcmpi(Eig_Solver, 'eigs')
    [kerN, ~] = eigs( cnormalize(W, 1)', n, 'LR' );
end
kerN = cnormalize_inplace(kerN')';
groups = kmeans(kerN, n, 'Start', Start, ...
                         'MaxIter', MaxIter, ...
                         'Replicates', Replicates, ...
                         'EmptyAction', 'singleton');
end

function vararg = vararginParser(vararg, vararg_in)
%VARARGINPARSER Input parser.

% How to use:
% % Set default 
% vararg = {'firstparameter', 1, 'secondparameter', magic(3)};
% % Overwrite by input
% vararg = vararginParser(vararg, varargin);
% % Generate variables
% for pair = reshape(vararg, 2, []) % pair is {propName;propValue}
%    eval([pair{1} '= pair{2};']);
% end

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

% count arguments
if mod(length(vararg_in), 2) ~= 0
    error('varargin needs propertyName/propertyValue pairs')
end
% overwrite vararg.
optionNames = vararg(1:2:end);
for pair = reshape(vararg_in, 2, []) % pair is {propName;propValue}
    optName = pair{1};
    index = find( strcmpi(optName, optionNames) );
    if ~isempty(index)
        vararg{index * 2} = pair{2};
    else
        error('%s is not a recognized parameter name', optName)
    end
end

end

function [Y, Xnorm] = cnormalize(X, p)
%CNORMALIZE normalizes columns.
% 
% [Y, XNORM] = cnormalize(X, P) normalizes columns of matrix X to unit
% ell_P norm, and returen the norm values to XNORM and data to Y.

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

if ~exist('p', 'var')
    p = 2;
end

if p == Inf
    Xnorm = max(abs(X), [], 1);
else
    Xnorm = sum(abs(X) .^p, 1) .^(1/p);
end
Y = bsxfun(@rdivide, X, Xnorm + eps);
end

function [X, Xnorm] = cnormalize_inplace(X, p)
%CNORMALIZE_INPLACE normalizes columns.
%   This is a inplace version of CNORMALIZE.

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

N = size( X, 2 );
if ~exist('p', 'var')
    p = 2;
end

if nargout > 1, Xnorm = zeros(1, N); end;

for iN = 1:N
    if p == Inf
        cnorm = max(abs(X(:, iN)), [], 1);
    else
        cnorm = sum(abs(X(:, iN)) .^p, 1) .^(1/p);
    end
    X(:, iN) = X(:, iN) / (cnorm + eps);
    
    if nargout > 1, Xnorm(iN) = cnorm; end;
end
end

