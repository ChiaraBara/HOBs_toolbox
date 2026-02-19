%% Compute prediction error time series by using a local linear estimation strategy

%%% INPUT
% B: data matrix (in column)
% jv: index column of the target variable
% iv: index column of the source variable(s)
% k: number of neighbors used to locally describe dynamics; if vectorial, the search of the optimal k is performed
% opt_figure: if equal to 1, the localizzation of the optimal value of k is depicted
%%% OUTPUT
% eUn: time series collecting the residuals of the prediction for each point of the multivariate series
% k_opt: (optimal) value of neighbors used to locally describe dynamics

function [eUn,k_opt] = nonlin_loclin_pred(B,jv,iv,k,opt_figure)

    if nargin < 5, opt_figure=0; end

    N = size(B,1);
    X = B(:,iv); Y = B(:,jv);
    eUn_tmp = nan(N,length(k));
    for ik = 1:length(k)
        parfor i = 1:N
            X_tmp = X; Y_tmp = Y;
            X_tmp(i,:) = []; Y_tmp(i,:) = [];
            idx_tmp = knnsearch(X_tmp,X(i,:),'K',k(ik),'Distance','chebychev');  
            Xn_tmp = X_tmp(idx_tmp,:); Yn_tmp = Y_tmp(idx_tmp,:);
            Bn_tmp=  [Yn_tmp Xn_tmp]; 
            out=lrp_LinRegStatic(Bn_tmp,1:size(Yn_tmp,2),size(Yn_tmp,2)+1:size(Yn_tmp,2)+size(Xn_tmp,2));
            eAm=out.eA;
            eYn_tmp = X(i,:)*eAm;
            eUn_tmp(i,ik) = Y(i,:)-eYn_tmp;
        end
        var_eUn_tmp(ik) = var(eUn_tmp(:,ik));
    end
    [~,opt_idx] = min(var_eUn_tmp);
    k_opt = k(opt_idx);
    if opt_figure == 1
        figure; scatter(k,var_eUn_tmp,'ok');
        hold on; plot(k(opt_idx),var_eUn_tmp(opt_idx),'or');
    end
    eUn = eUn_tmp(:,opt_idx);
end