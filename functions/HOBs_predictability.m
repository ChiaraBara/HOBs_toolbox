%% Compute WMS via predictability measures

%%% INPUT
% B: data matrix (in column)
% jv: index column of the target variable
% iv: two element matrix with index column of the source variables
% kopt (optional): structure indicating a fixed value of k used to find neighbors on which perfom regression; if not present, the search of the optimal value of k is performed
%%% OUTPUT
% Delta_MP: WMS metric computed in terms of predictability measures 
% k_opt: (optimal) values of k used to obtain the three metrics of mutual predictability


function [Delta_MP, k_opt] = HOBs_predictability(B,jv,iv,kopt)

    B_tmp = B(:,[jv iv]);

    [N,M]=size(B_tmp);
    for m=1:M
        B_tmp(:,m)=B_tmp(:,m)-mean(B_tmp(:,m));
    end

    if nargin < 4
        veck_X1X2 = 2:N;
        veck_X1 = veck_X1X2; veck_X2 = veck_X1X2;
    else
        veck_X1X2 = kopt.kopt_X1X2;
        veck_X1 = kopt.kopt_X1;
        veck_X2 = kopt.kopt_X2;
    end
   
    Y = B_tmp(:,1);

    [eUn_X1X2, kopt_X1X2] = nonlin_loclin_pred(B_tmp,1,[2 3],veck_X1X2); 
    [eUn_X1, kopt_X1] = nonlin_loclin_pred(B_tmp,1,2,veck_X1); 
    [eUn_X2, kopt_X2] = nonlin_loclin_pred(B_tmp,1,3,veck_X2); 

    var_Y = var(Y);
    var_eUn_X1X2 = var(eUn_X1X2); MP_X1X2 = var_Y - var_eUn_X1X2;
    var_eUn_X1 = var(eUn_X1); MP_X1 = var_Y - var_eUn_X1;
    var_eUn_X2 = var(eUn_X2); MP_X2 = var_Y - var_eUn_X2;

    Delta_MP = MP_X1X2 - MP_X1 - MP_X2;
    k_opt.kopt_X1X2 = kopt_X1X2;
    k_opt.kopt_X1 = kopt_X1;
    k_opt.kopt_X2 = kopt_X2;

end