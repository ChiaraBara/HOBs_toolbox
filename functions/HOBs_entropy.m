%% Compute WMS via entropy measures

%%% INPUT
% B: data matrix (in column)
% jv: index column of the target variable
% iv: two element matrix with index column of the source variables
% k: number of neighbors
%%% OUTPUT
% Delta_MI = WMS metric computed in terms of entropy measures 

function Delta_MI = HOBs_entropy(B,jv,iv,k)

    B_tmp = B(:,[jv iv]);

    out = MIknn(zscore(B_tmp),k);
    MI_X1X2 = out.I_YX;
    MI_X1 = out.I_YX1;
    MI_X2 = out.I_YX2;

    Delta_MI = MI_X1X2 - MI_X1 - MI_X2;

end