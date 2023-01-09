function [Z] = RMSC(S, lambda, opts)
stream = RandStream.getGlobalStream;
reset(stream);

[m,p,n]=size(S);%num of samples
if m~=p
    error('input matrix T must be a square matrix( transitions matrix ).\n');
end
Z=zeros(m,p);
E=randn(m,p,n);
Y=zeros(m,p,n);
B=zeros(m,p);

if isfield(opts,'mu')
    mu=opts.mu;
else
    mu=1e-3;
end
if isfield(opts,'rho')
    rho=opts.rho;
else
    rho=2.0;
end
if isfield(opts,'max_iter')
    max_iter=opts.max_iter;
else
    max_iter=100;
end

step=0;
while(1)
    %tic;
    step=step+1;
    max_inf_norm=-1;
    for i=1:n
        Si=S(:,:,i);
        Ei=E(:,:,i);
        diff=Z+Ei-Si;
        inf_norm=norm(diff,'inf');
        max_inf_norm=max(max_inf_norm,inf_norm);
    end
    if step>1 && max_inf_norm<opts.eps  
        break;
    end
    
    if step > max_iter
         fprintf('reach max iterations %d \n',step);
         break;
    end
 
    %update Z
    for i=1:n
    B=B+(S(:,:,i)-E(:,:,i)-Y(:,:,i)/mu);
    end
    B=B/n;
    C=1/(mu*2*n);
    [U, Sigma, V] = svd(B,'econ');
     Sigma = diag(Sigma);
     svp=length(find(Sigma>C));
      if svp>=1
          Sigma = Sigma(1:svp)-C;
      else
          svp = 1;
          Sigma = 0;
      end
     Z= U(:, 1:svp) *diag(Sigma)*V(:, 1:svp)';

%update Ei
 for i=1:n
          F = S(:,:,i)-Z-Y(:,:,i)/mu;
         [Econcat] = solve_l1l2(F,lambda/mu);
          E(:,:,i) = Econcat;
          Y(:,:,i)=Y(:,:,i)+mu*(Z+E(:,:,i)-S(:,:,i));
 end

%     Z=Z+mu*(P-Q);
    %update mu
    mu=min(rho*mu,1e10);
end

