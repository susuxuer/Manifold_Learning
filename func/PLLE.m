function [O,Res] = PLLE(X,K,d)
warning off;
[D,N] = size(X);
% STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS 
X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
[sorted,index] = sort(distance);
neighborhood = index(2:(1+K),:);
% STEP2: SOLVE FOR RECONSTRUCTION WEIGHTS
if(K>D) 
  tol=1e-3; % regularlizer in case constrained fits are ill conditioned
else
  tol=0;
end
W = zeros(K,N);
for ii=1:N
   z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
   C = z'*z;                                        % local covariance
   C = C + eye(K,K)*tol*trace(C);                   % regularlization (K>D)
   W(:,ii) = C\ones(K,1);                           % solve Cw=1
   W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
end;
% STEP 3: COMPUTE EMBEDDING FROM EIGENVECTS OF COST MATRIX M=(I-W)'(I-W)
M = sparse(1:N,1:N,ones(1,N),N,N,4*K*N); 
for ii=1:N
   w = W(:,ii);
   jj = neighborhood(:,ii);
   M(ii,jj) = M(ii,jj) - w';
   M(jj,ii) = M(jj,ii) - w;
   M(jj,jj) = M(jj,jj) + w*w';
end;

% CALCULATION OF EMBEDDING
options.disp = 0; options.isreal = 1; options.issym = 1; 
[Y,eigenvals] = eigs(M,d+1,0,options);
O = Y(:,1:d)'*sqrt(N);   % bottom evect is [1,1,1,1...] with eval 0
%calculate Z
XX=[];
rr=1;
for ii=1:D-1
    for jj=ii+1:D
     XX(rr,:)=X(ii,:).*X(jj,:);
     rr=rr+1;
    end
end
Z=[ones(1,N);X;X.*X;XX];
Res = O * pinv(Z);