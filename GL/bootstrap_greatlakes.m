%fix PCAs to the original per-fold PCA data
%for each fold, resample with replacement from training data, and separately resample with replacement from testing data
%run BBS CV on each of the 8 models (PE, PEw/G, PEw/enrich, PEw/grades, PEw/school, PEw/parent, PEw/PF10, PEw/all5) to get mean correlation
%calculate the PE-X for each of the 7 pairs of models and store that so we get a distribution of differences for each model
%repeat 10000 times
addpath /home/mangstad/repos/Misc_utils

%need: folds, Aa, phenotypes, nuisance
NumComps = [400];

rng(SEED);
Exp = '/home/mangstad/work/ser';

load([Exp '/Data/data.mat']);

nFold = numel(unique(fs));
n = size(ns,1);
nSamples = 10;

% Tests = {
%     'PE';
%     'PE w/ G';
%     'PE w/ enrich';
%     'PE w/ grades';
%     'PE w/ family';
%     'PE w/ school';
%     'PE w/ PF10';
%     'PE w/ all';
%     };

%tmp = mc_bbs(fms,aps(:,1),[ns],fs,NumComps,'NoComponents',1,'Scores',Aas);

outcorr = zeros(nSamples,7);
outcomps = zeros(nSamples,8);

for iSample = 1:nSamples
    %resample within fold
    if (SEED==0)
        Aas_rs = Aas;
        samp = [1:n]';
    else
        samp = zeros(n,1);
        for iFold = 1:nFold
            idx = fs==iFold;
            pos = find(idx);
            nf = sum(idx);
            s = randsample(nf,nf,1);
            samp(idx) = pos(s);
        end
        Aas_rs = [];
        for iFold = 1:nFold
            Aas_rs{iFold} = Aas{iFold}(samp,:);
        end
    end
tic
    tmpPE = mc_bbs(ones(n,1),aps(samp,1),[ns(samp,:)],fs,NumComps,'NoComponents',1,'Scores',Aas_rs);
toc
    tmpPE_G = mc_bbs(ones(n,1),aps(samp,1),[ns(samp,:) aps(samp,2)],fs,NumComps,'NoComponents',1,'Scores',Aas_rs);
    tmpPE_enrich = mc_bbs(ones(n,1),aps(samp,1),[ns(samp,:) aps(samp,3)],fs,NumComps,'NoComponents',1,'Scores',Aas_rs);
    tmpPE_grades = mc_bbs(ones(n,1),aps(samp,1),[ns(samp,:) aps(samp,4)],fs,NumComps,'NoComponents',1,'Scores',Aas_rs);
    tmpPE_family = mc_bbs(ones(n,1),aps(samp,1),[ns(samp,:) aps(samp,5)],fs,NumComps,'NoComponents',1,'Scores',Aas_rs);
    tmpPE_school = mc_bbs(ones(n,1),aps(samp,1),[ns(samp,:) aps(samp,6)],fs,NumComps,'NoComponents',1,'Scores',Aas_rs);
%    tmpPE_pf10 = mc_bbs(ones(n,1),aps(samp,1),[ns(samp,:) aps(samp,7)],fs,NumComps,'NoComponents',1,'Scores',Aas_rs);
    tmpPE_all6 = mc_bbs(ones(n,1),aps(samp,1),[ns(samp,:) aps(samp,2:6)],fs,NumComps,'NoComponents',1,'Scores',Aas_rs);
    outcorr(iSample,:) = [tmpPE.mean_corr tmpPE_G.mean_corr tmpPE_enrich.mean_corr tmpPE_grades.mean_corr tmpPE_family.mean_corr tmpPE_school.mean_corr tmpPE_all6.mean_corr];
%    outcomps(iSample,:) = [mean(tmpPE.bestcomps) mean(tmpPE_G.bestcomps) mean(tmpPE_enrich.bestcomps) mean(tmpPE_grades.bestcomps) mean(tmpPE_family.bestcomps) mean(tmpPE_school.bestcomps) mean(tmpPE_pf10.bestcomps) mean(tmpPE_all5.bestcomps)];
    
end

filename = [Exp '/' sprintf('Results/corr_%05d.txt',SEED)];
save(filename,'outcorr','-ASCII');

%filename = [Exp '/' sprintf('Results/comps_%05d.txt',SEED)];
%save(filename,'outcomps','-ASCII');
