% Isomap的函数
% 它是由两个函数组成的，需要相互调用实现Isomap降维方法
% 组成的两个函数为L2_distance.m和Isomap_I。
% Data：输入的高维数据，行数为样本，列数为维数；
% D：输入的降维的维数；
% k：输入的邻近点的个数。
% Y：降维处理后的数据,列数表示维数，行数表示样本数。

function [T, Res] = Isomap_2(X, d, K)
DD = L2_distance(X,X,1);
options.dims = d;
[Y, R] = Isomap_I(DD, K, options);
Y=Y.coords{1}';
T = Y';
[D,N] = size(X);
newX = X;
XX=[];
rr=1;
for ii=1:D-1
    for jj=ii+1:D
     XX(rr,:)=newX(ii,:).*newX(jj,:);
     rr=rr+1;
    end
end
Z=[ones(1,N);newX;newX.*newX;XX];
Res =Y'*pinv(Z);


