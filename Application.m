clear all; close all; clc;

load([pwd,'\physio_data_1subj.mat']);
addpath([pwd,'\functions\'])

%%%%%%% parameters
pfilter=0.94;   %filter

sel_series = label_data;
first_triplet = [3 1 2];
second_triplet = [2 4 5];
triplet_label = {'\{MAP; SAP, DAP\}', '\{DAP; PEP, LVET\}'};
jv = 1; % target index in the triplets
iv = [2 3]; % sources index in the triplets
k = 5;

numsurr = 100;
alpha = 0.05;
minshift = 20;

for itriplet = 1:2

    if itriplet == 1
        triplet = first_triplet;
    elseif itriplet == 2
        triplet = second_triplet;
    end

    data = data_1subj(:,triplet);
    
    [N,M]=size(data);

    for m=1:M
        B(:,m)=AR_filter(data(:,m),1,pfilter); % filtered series
        B(:,m)=B(:,m)-mean(B(:,m));
        B(:,m) = filloutliers(B(:,m),'spline','quartile');
    end
    
    [Delta_MP(1,triplet), kopt] = HOBs_predictability(B,jv,iv);
    Delta_MI(1,triplet) = HOBs_entropy(B,jv,iv,k);

    for isurr = 1:numsurr

        maxshift=N-minshift;
        lagshift=fix(rand(1,2)*(maxshift-minshift+1)+minshift);% shift casuale tra minshift e maxshift 
        X1_surr=circshift(B(:,iv(1)),lagshift(1));%is-esimo shift dei valori del candidato scelto  
        X2_surr=circshift(B(:,iv(2)),lagshift(2));
        Y = B(:,jv);

        B_surr = [Y, X1_surr, X2_surr];

        Delta_MP_surr(isurr,triplet) = HOBs_predictability(B_surr,1,[2 3], kopt);
        Delta_MI_surr(isurr,triplet) = HOBs_entropy(B_surr,1,[2 3],k);

    end

end

figure;
for itriplet = 1:2

    subplot(2,2,itriplet);
    scatter(1,Delta_MP(1,itriplet),'xr');
    hold on;
    scatter(2*ones(numsurr,1),Delta_MP_surr(:,itriplet),'ob');
    th=prctile(Delta_MP_surr(:,itriplet),100*(1-alpha));
    yline(th,'--k');
    xlim([0.5 2.5]); xticks([]); 
    title(triplet_label(itriplet));
    ylabel('\Delta_{MP} [a.u.]');

    subplot(2,2,2+itriplet);
    scatter(1,Delta_MI(1,itriplet),'xr');
    hold on;
    scatter(2*ones(numsurr,1),Delta_MI_surr(:,itriplet),'ob');
    th=prctile(Delta_MI_surr(:,itriplet),100*(1-alpha));
    yline(th,'--k');
    xlim([0.5 2.5]); xticks([]); 
    ylabel('\Delta_{MI} [nats]');

end