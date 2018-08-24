function [T,Res]= PLE(X, d, k)
%LAPLACIAN_EIGEN Performs non-linear dimensionality reduction using Laplacian Eigenmaps
%
%   [mappedX, mapping] = laplacian_eigen(X, no_dims, k, sigma, eig_impl)
%
% Performs non-linear dimensionality reduction using Laplacian Eigenmaps.
% The data is in matrix X, in which the rows are the observations and the
% columns the dimensions. 
% The reduced data is returned in the matrix mappedX.
sigma = 1.45;
no_dims=d;
[m,N] = size(X');
    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    if ~exist('k', 'var')
        k = 12;
    end
	if ~exist('sigma', 'var')
		sigma = 1;
    end
  
    % Construct neighborhood graph
    disp('Constructing neighborhood graph...');
    if size(X, 1) < 4000
        G = squareform(pdist(X, 'euclidean'));
        % Compute neighbourhood graph
        [tmp, ind] = sort(G); 
        for i=1:size(G, 1)
            G(i, ind((2 + k):end, i)) = 0; 
        end
        G = sparse(double(G));
        G = max(G, G');             % Make sure distance matrix is symmetric
    else
        G = find_nn(X, k);
    end
    G = G .^ 2;
	G = G ./ max(max(G));
    
    % Compute weights (W = G)
    disp('Computing weight matrices...');
    
    % Compute Gaussian kernel (heat kernel-based weights)
    G(G ~= 0) = exp(-G(G ~= 0) / (2 * sigma ^ 2));
        
    % Construct diagonal weight matrix
    D = diag(sum(G, 2));
    
    % Compute Laplacian
    L = D - G;
    L(isnan(L)) = 0; D(isnan(D)) = 0;
	L(isinf(L)) = 0; D(isinf(D)) = 0;
    
    % Construct eigenmaps (solve Ly = labda*Dy)
    disp('Constructing Eigenmaps...');
	tol = 0;
   
        options.disp = 0;
        options.isreal = 1;
        options.issym = 1;
        [mappedX, lambda] = eigs(L, D, no_dims + 1, tol, options);			% only need bottom (no_dims + 1) eigenvectors
     
    % Sort eigenvectors in ascending order
    lambda = diag(lambda);
    [lambda, ind] = sort(lambda, 'ascend');
    lambda = lambda(2:no_dims + 1);
  
    % Final embedding
	mappedX = mappedX(:,ind(2:no_dims + 1));
    T =mappedX';
    XX=[];
rr=1;
newX =X';
for ii=1:m-1
    for jj=ii+1:m
     XX(rr,:)=newX(ii,:).*newX(jj,:);
     rr=rr+1;
    end
end
Z=[ones(1,N);newX;newX.*newX;XX];
Res = T* pinv(Z);
    
    
    
