% Check the non-employment of GKOS2021 models
% This script is only used to double-check certain aspects
% Not used for anything in the actual paper

useModel=6;
nSigmaz=4;

load(['./SavedOutput/BasicDiscretization/DiscretizedEarnings',num2str(useModel),'nSigmaz',num2str(nSigmaz),'.mat'])

% Compare to Figure 12-E (pg 2327; 25 of 37)

% The process: pi_zupsilon_J
simoptions=struct();
mcmomentsoptions.n_z=n_zupsilon;
[mean,variance,autocorrelation,statdist]=MarkovChainMoments_FHorz(zupsilon_grid_J,pi_zupsilon_J,jequaloneDistzupsilon,simoptions,mcmomentsoptions);

% To 
simoptions.nsims=10^5;
simoptions.n_z=n_zupsilon;
zupsilon_panel=MarkovChainSimulation_FHorz(zupsilon_grid_J,pi_zupsilon_J,jequaloneDistzupsilon,simoptions);
% For each simulation, count years employed (upsilon=0)
yearsemp=zeros(simoptions.nsims,1);
parfor ii=1:simoptions.nsims
    temp=zupsilon_panel(ii,:,2); % 2 is the non-employment shock
    for tt=1:38
        if temp(tt)==0 % if employed
            yearsemp(ii)=yearsemp(ii)+1;
        end
    end
end

[cnt_unique_yearsemp, unique_yearsemp] = hist(yearsemp,unique(yearsemp));

yearsemp_cdf=cumsum(cnt_unique_yearsemp)/sum(cnt_unique_yearsemp);

plot(unique_yearsemp,yearsemp_cdf)

[unique_yearsemp';yearsemp_cdf]



%% Check out the upsilon transition probabilities
% This it to try and figure out whay the mixture probabilities for model 4
% appear as negative for some grid values

Params.mixprobs_eta=@(agej,zlag) (-0.474+1.961*(agej/10)-3.183*zlag-0.187*(agej/10)*zlag)/100; % GKOS2021 denote this p_z

mixprobs_eta=nan(n_z,41);
for ii=1:n_z
    for jj=2:41
        zlag=z_grid_J(ii,jj-1);
        agej=Params.agej(jj);
        mixprobs_eta(ii,jj)=Params.mixprobs_eta(agej,zlag);
    end
end







