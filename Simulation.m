clear all; close all; clc;

load([pwd,'\data_sim.mat']);
addpath([pwd,'\functions\'])

%%%%%%% parameters

jv = 1; % target index in the triplets
iv = [2 3]; % sources index in the triplets
k = 5;

numsurr = 100;
alpha = 0.05;
minshift = 50;

for isim = 1:size(data_sim,2)

    B_tmp = data_sim{1,isim};
    B(:,jv) = B_tmp(2:end,jv);
    B(:,iv) = B_tmp(1:end-1,iv);

    N = size(B,1);
    
    [Delta_MP(1,isim), kopt] = HOBs_predictability(B,jv,iv);
    Delta_MI(1,isim) = HOBs_entropy(B,jv,iv,k);

    for isurr = 1:numsurr

        maxshift=N-minshift;
        lagshift=fix(rand(1,2)*(maxshift-minshift+1)+minshift);% shift casuale tra minshift e maxshift 
        X1_surr=circshift(B(:,iv(1)),lagshift(1));%is-esimo shift dei valori del candidato scelto  
        X2_surr=circshift(B(:,iv(2)),lagshift(2));
        Y = B(:,jv);

        B_surr = [Y, X1_surr, X2_surr];

        Delta_MP_surr(isurr,isim) = HOBs_predictability(B_surr,1,[2 3], kopt);
        Delta_MI_surr(isurr,isim) = HOBs_entropy(B_surr,1,[2 3],k);

    end

end

figure;
for isim = 1:size(data_sim,2)
    subplot(2,4,isim);
    scatter(1,Delta_MP(1,isim),'xr');
    hold on;
    scatter(2*ones(numsurr,1),Delta_MP_surr(:,isim),'ob');
    th=prctile(Delta_MP_surr(:,isim),100*(1-alpha));
    yline(th,'--k');
    xlim([0.5 2.5]); xticks([]); 
    title(['Sim ',num2str(isim)]);
    ylabel('\Delta_{MP} [a.u.]');

    subplot(2,4,4+isim);
    scatter(1,Delta_MI(1,isim),'xr');
    hold on;
    scatter(2*ones(numsurr,1),Delta_MI_surr(:,isim),'ob');
    th=prctile(Delta_MI_surr(:,isim),100*(1-alpha));
    yline(th,'--k');
    xlim([0.5 2.5]); xticks([]); 
    ylabel('\Delta_{MI} [nats]');
end