function [T,Res] = PLTSA(data,d,K,NI)
[m,N] = size(data);  % m是输入样本点的维数. N是输入样本个数
% Step 0:  Neighborhood Index
if nargin<4                   %如果输入参数个数<4
    if length(K)==1           %向量K的长度
        K = repmat(K,[1,N]);  %矩阵K的内容堆叠在[1*N]的矩阵上
    end;
    NI = cell(1,N);            %产生一个1*N的空单元数组
    if m>N
        a = sum(data.*data);  %  .* 用于矩阵中对应元素的相乘,并且竖向相加
        dist2 = sqrt(repmat(a',[1 N]) + repmat(a,[N 1]) - 2*(data'*data)); %计算距离
        for i=1:N
            % Determine ki nearest neighbors of x_j
            [dist_sort,J] = sort(dist2(:,i));   %sort 对dist2 进行排序，dist_sort 返回对每列排序后的结果
            %J返回的是排序后矩阵中每个数在矩阵未排序前每列中的位置
            Ii = J(1:K(i));             %i个近邻点
            NI{i} = Ii;                 %放进单元数组N中
        end;
    else                                %   m<=N
        for i=1:N
            % Determine ki nearest neighbors of x_j
            x = data(:,i); ki = K(i);
            dist2 = sum((data-repmat(x,[1 N])).^2,1);    
            [dist_sort,J] = sort(dist2);  
            Ii = J(1:ki);  
            NI{i} = Ii;
        end;
    end;
else
    K = zeros(1,N);         %1*N的全为0的矩阵
    for i=1:N
        K(i) = length(NI{i});       %length()求得是向量的长度
    end;
end;
% Step 1:  local information
BI = {}; 
%Theta = {};            %!?
for i=1:N 
    % Compute the d largest right singular eigenvectors of the cenantered matrix
    Ii = NI{i}; ki = K(i);
    Xi = data(:,Ii)-repmat(mean(data(:,Ii),2),[1,ki]);  %中心化样本矩阵，得到hatX
    W = Xi'*Xi;                     %协方差矩阵（PSD）
  %(感觉有点多余)  W = (W+W')/2;                   %构造成实对称矩阵
    [Vi,Si] = schur(W);             %S是上三角矩阵（对角线为特征值），V是酉阵
    [s,Ji] = sort(-diag(Si));       %给特征值从大到小来排序
    
    Vi = Vi(:,Ji(1:d));             %找到前d个特征向量
    %----%
    % construct Gi
    Gi = [repmat(1/sqrt(ki),[ki,1]) Vi];  %d+1个列向量
    % compute the local orthogonal projection Bi = I-Gi*Gi' 
    % that has the null space span([e,Theta_i^T]). 
    BI{i} = eye(ki)-Gi*Gi';                 %BI是正交投影矩阵（Gi是论文里面的Qn）
end;
B = speye(N);     %创建单位稀疏矩阵
%K稀疏： 只要出了这K个值之外，其他的值很小，我们就认为向量是稀疏的
for i=1:N
    Ii = NI{i};
    B(Ii,Ii) = B(Ii,Ii)+BI{i};
    B(i,i) = B(i,i)-1;
end;
B = (B+B')/2;
options.disp = 0;                   %
options.isreal = 1;                 %B为实矩阵
options.issym = 1;                  %矩阵B对称
[U,D] = eigs(B,d+2,0,options);          %D为d个最大特征值对角阵，U的列向量为对应的特征向量
lambda = diag(D);                       %将特征值取出来构成一个列向量
[lambda_s,J] = sort(abs(lambda));        %将特征值构成的列向量按升序排列
U = U(:,J);                               %将J放入矩阵U(U wei 全局对齐矩阵)
lambda = lambda(J);                             
T = U(:,2:d+1)';                           %取前(2:d+1）d个特征向量构成全局坐标矩阵
XX=[];
rr=1;
X = data;
for ii=1:m-1
    for jj=ii+1:m
     XX(rr,:)=X(ii,:).*X(jj,:);
     rr=rr+1;
    end
end
Z=[ones(1,N);X;X.*X;XX];
Res = T*pinv(Z);
