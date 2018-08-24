% --- ISOMAP Algorithm
% Written by Tenenbaum, de Silva, and Langford (2000). 

function [Y, R] = Isomap_I(D, K, options); 

%%%%% Step 0: Initialization and Parameters %%%%%
N = size(D,1); 
INF =  1000*max(max(D))*N;  %% effectively infinite distance
dims = options.dims; 
comp = 1;  % Assume one component.
Y.coords = cell(length(dims),1); 
R = zeros(1,length(dims)); 
%%%%% Step 1: Construct neighborhood graph %%%%%
[tmp, ind] = sort(D); 
for i=1:N
    D(i,ind((2+K):end,i)) = INF; 
end;
D = min(D,D');    %% Make sure distance matrix is symmetric
%%%%% Step 2: Compute shortest paths %%%%%
for k=1:N
     D = min(D,repmat(D(:,k),[1 N])+repmat(D(k,:),[N 1])); 
end
%%%%% Remove outliers from graph %%%%%
n_connect = sum(~(D==INF));        %% number of points each point connects to
[tmp, firsts] = min(D==INF);       %% first point each point connects to
[comps, I, J] = unique(firsts);    %% represent each connected component once
size_comps = n_connect(comps);     %% size of each connected component
[tmp, comp_order] = sort(size_comps);  %% sort connected components by size
comps = comps(comp_order(end:-1:1));    
size_comps = size_comps(comp_order(end:-1:1)); 
n_comps = length(comps);               %% number of connected components
if (comp>n_comps)                
     comp=1;                              %% default: use largest component
end
Y.index = find(firsts==comps(comp)); 
D = D(Y.index, Y.index); 
N = length(Y.index); 
%%%%% Step 3: Construct low-dimensional embeddings (Classical MDS) %%%%%
opt.disp = 0; 
[vec, val] = eigs(-.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2)), max(dims), 'LR', opt); 
h = real(diag(val)); 
[foo,sorth] = sort(h);  sorth = sorth(end:-1:1); 
val = real(diag(val(sorth,sorth))); 
vec = vec(:,sorth); 
D = reshape(D,N^2,1); 
for di = 1:length(dims)
     if (dims(di)<=N)
         Y.coords{di} = real(vec(:,1:dims(di)).*(ones(N,1)*sqrt(val(1:dims(di)))'))'; 
         r2 = 1-corrcoef(reshape(real(L2_distance(Y.coords{di}, Y.coords{di},0)),N^2,1),D).^2; 
         R(di) = r2(2,1); 
     end
end
clear D; 





 