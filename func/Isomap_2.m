% Isomap�ĺ���
% ����������������ɵģ���Ҫ�໥����ʵ��Isomap��ά����
% ��ɵ���������ΪL2_distance.m��Isomap_I��
% Data������ĸ�ά���ݣ�����Ϊ����������Ϊά����
% D������Ľ�ά��ά����
% k��������ڽ���ĸ�����
% Y����ά����������,������ʾά����������ʾ��������

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


