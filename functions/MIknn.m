%% k-nearest neighbor estimation of Mutual Information
% Computes the mutual information between the first column and the
% remaining columns of the input matrix B (whole and single)

function out = MIknn(B,k,metric)

if ~exist('metric','var'), metric='maximum'; end

[N,M]=size(B);

% set subspaces of lower dimension
M_Y = B(:,1);
M_X = B(:,2:end);
for ix = 1:M-1
    eval(['M_X',num2str(ix+1),' = B(:,',num2str(ix+1),');'])
end

% neighbor search in space of higher dimension
atria = nn_prepare(B, metric);
[~, distances] = nn_search(B, atria, (1:N)', k, 0);
dd=distances(:,k);

if ~isempty(M_Y)
    atriaY = nn_prepare(M_Y, metric);
    [count_Y, tmp] = range_search(M_Y, atriaY, (1:N)', dd, 0);
    tmp=tmp(:,2);
    for n=1:length(tmp)
        count_Y(n)=max(k-1,count_Y(n)-sum(tmp{n}==dd(n)));
    end
end

if ~isempty(M_X)
    atriaX = nn_prepare(M_X, metric);
    [count_X, tmp] = range_search(M_X, atriaX, (1:N)', dd, 0);
    tmp=tmp(:,2);
    for n=1:length(tmp)
        count_X(n)=max(k-1,count_X(n)-sum(tmp{n}==dd(n)));
    end
end

for ix = 1:M-1
    M_X_tmp = eval(['M_X',num2str(ix+1)]);
    M_YX_tmp = eval(['[M_Y M_X',num2str(ix+1),']']);
    
    if ~isempty(M_X_tmp)
         atria_X = nn_prepare(M_X_tmp, metric);
         [count_X_tmp, tmp] = range_search(M_X_tmp, atria_X, (1:N)', dd, 0);
         tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
         for n=1:length(tmp)-1
             count_X_tmp(n)=max(k-1,count_X_tmp(n)-sum(tmp{n}==dd(n)));
         end
     else
         count_X_tmp=(Ni(iy)-1)*ones(Ni(iy),1);
    end 
    
    if ~isempty(M_YX_tmp)
         atria_YX = nn_prepare(M_YX_tmp, metric);
         [count_YX_tmp, tmp] = range_search(M_YX_tmp, atria_YX, (1:N)', dd, 0);
         tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
         for n=1:length(tmp)-1
             count_YX_tmp(n)=max(k-1,count_YX_tmp(n)-sum(tmp{n}==dd(n)));
         end
     else
         count_YX_tmp=(Ni(iy)-1)*ones(Ni(iy),1);
    end 
    
    IYX_tmp = psi(N) + (1/N)*(sum(psi(count_YX_tmp+1) - psi(count_X_tmp+1) -psi(count_Y+1)));
    
    eval(['out.I_YX',num2str(ix),'=IYX_tmp;']);
    
end


IYX = psi(k) + psi(N) - (1/N)*(sum(psi(count_X+1) + psi(count_Y+1))) ; 
out.I_YX = IYX;
    
end

    
