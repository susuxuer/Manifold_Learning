function [T,Res] = PLTSA(data,d,K,NI)
[m,N] = size(data);  % m�������������ά��. N��������������
% Step 0:  Neighborhood Index
if nargin<4                   %��������������<4
    if length(K)==1           %����K�ĳ���
        K = repmat(K,[1,N]);  %����K�����ݶѵ���[1*N]�ľ�����
    end;
    NI = cell(1,N);            %����һ��1*N�Ŀյ�Ԫ����
    if m>N
        a = sum(data.*data);  %  .* ���ھ����ж�ӦԪ�ص����,�����������
        dist2 = sqrt(repmat(a',[1 N]) + repmat(a,[N 1]) - 2*(data'*data)); %�������
        for i=1:N
            % Determine ki nearest neighbors of x_j
            [dist_sort,J] = sort(dist2(:,i));   %sort ��dist2 ��������dist_sort ���ض�ÿ�������Ľ��
            %J���ص�������������ÿ�����ھ���δ����ǰÿ���е�λ��
            Ii = J(1:K(i));             %i�����ڵ�
            NI{i} = Ii;                 %�Ž���Ԫ����N��
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
    K = zeros(1,N);         %1*N��ȫΪ0�ľ���
    for i=1:N
        K(i) = length(NI{i});       %length()����������ĳ���
    end;
end;
% Step 1:  local information
BI = {}; 
%Theta = {};            %!?
for i=1:N 
    % Compute the d largest right singular eigenvectors of the cenantered matrix
    Ii = NI{i}; ki = K(i);
    Xi = data(:,Ii)-repmat(mean(data(:,Ii),2),[1,ki]);  %���Ļ��������󣬵õ�hatX
    W = Xi'*Xi;                     %Э�������PSD��
  %(�о��е����)  W = (W+W')/2;                   %�����ʵ�Գƾ���
    [Vi,Si] = schur(W);             %S�������Ǿ��󣨶Խ���Ϊ����ֵ����V������
    [s,Ji] = sort(-diag(Si));       %������ֵ�Ӵ�С������
    
    Vi = Vi(:,Ji(1:d));             %�ҵ�ǰd����������
    %----%
    % construct Gi
    Gi = [repmat(1/sqrt(ki),[ki,1]) Vi];  %d+1��������
    % compute the local orthogonal projection Bi = I-Gi*Gi' 
    % that has the null space span([e,Theta_i^T]). 
    BI{i} = eye(ki)-Gi*Gi';                 %BI������ͶӰ����Gi�����������Qn��
end;
B = speye(N);     %������λϡ�����
%Kϡ�裺 ֻҪ������K��ֵ֮�⣬������ֵ��С�����Ǿ���Ϊ������ϡ���
for i=1:N
    Ii = NI{i};
    B(Ii,Ii) = B(Ii,Ii)+BI{i};
    B(i,i) = B(i,i)-1;
end;
B = (B+B')/2;
options.disp = 0;                   %
options.isreal = 1;                 %BΪʵ����
options.issym = 1;                  %����B�Գ�
[U,D] = eigs(B,d+2,0,options);          %DΪd���������ֵ�Խ���U��������Ϊ��Ӧ����������
lambda = diag(D);                       %������ֵȡ��������һ��������
[lambda_s,J] = sort(abs(lambda));        %������ֵ���ɵ�����������������
U = U(:,J);                               %��J�������U(U wei ȫ�ֶ������)
lambda = lambda(J);                             
T = U(:,2:d+1)';                           %ȡǰ(2:d+1��d��������������ȫ���������
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
