NumComps = [50:50:1000];
%mainpheno = [dat.ses_fac dat.EdYearsAverage dat.Income2Needs dat.reshist_addr1_adi_wsum];
results = mc_bbs(featuremat,mainpheno(:,1),nuisance,folds,10,'LOSOPheno',0,'NoComponents',1);
results_ses = mc_bbs(featuremat,mainpheno(:,1),nuisance,folds,NumComps,'LOSOPheno',0,'NoComponents',1,'Scores',results.Aa);
results_pe = mc_bbs(featuremat,mainpheno(:,2),[nuisance mainpheno(:,[3 4])],folds,NumComps,'LOSOPheno',0,'NoComponents',1,'Scores',results.Aa);
results_inc = mc_bbs(featuremat,mainpheno(:,3),[nuisance mainpheno(:,[2 4])] ,folds,NumComps,'LOSOPheno',0,'NoComponents',1,'Scores',results.Aa);
results_adi = mc_bbs(featuremat,mainpheno(:,4),[nuisance mainpheno(:,[2 3])],folds,NumComps,'LOSOPheno',0,'NoComponents',1,'Scores',results.Aa);

clear Aa;
for i = 1:19
    Aa{i} = results.Aa{i}(:,1:2000);
end
save('Results/perfold_Aa.mat','Aa','-v7.3');

results_ses = rmfield(results_ses,'Aa');
results_pe = rmfield(results_pe,'Aa');
results_inc = rmfield(results_inc,'Aa');
results_adi = rmfield(results_adi,'Aa');

save('Results/bbs_cv_models.mat','results_ses','results_pe','results_inc','results_adi');

load('Results/perfold_Aa.mat');
load('Results/bbs_cv_models.mat');

Model = {'SES','PE','INC','ADI'}';
MeanComps = [mean(results_ses.bestcomps) ...
mean(results_pe.bestcomps) ...
mean(results_inc.bestcomps) ...
mean(results_adi.bestcomps)]'
MaxComps = [max(results_ses.bestcomps) ...
max(results_pe.bestcomps) ...
max(results_inc.bestcomps) ...
max(results_adi.bestcomps)]';
SDComps = [std(results_ses.bestcomps) ...
std(results_pe.bestcomps) ...
std(results_inc.bestcomps) ...
std(results_adi.bestcomps)]';
BBSCV_r = [results_ses.mean_corr ...
results_pe.mean_corr ...
results_inc.mean_corr ...
results_adi.mean_corr]';

out = table(Model,MeanComps,MaxComps,SDComps,BBSCV_r);
out


%BBS permutations
bestcomps = floor(mean(results_ses_nog.bestcomps));

tic;
results_ses_nog_perm = mc_bbs_perm(featuremat,mainpheno(:,1),nuisance,folds,bestcomps,10000,'Scores',Aa,'Perms',pset);
toc;
fold_p_ses = (sum(bsxfun(@ge,squeeze(results_ses_perm.fold_corr_perm),results_ses_perm.fold_corr'-(100*eps)),2))./10000;

bestcomps_pe = floor(mean(results_pe.bestcomps));
bestcomps_inc = floor(mean(results_inc.bestcomps));
bestcomps_adi = floor(mean(results_adi.bestcomps));

results_inc_perm = mc_bbs_perm(featuremat,mainpheno(:,3),[nuisance mainpheno(:,[2 4])],folds,bestcomps_inc,10000,'Scores',Aa,'Perms',pset(:,2:1001));
fold_p_inc = (sum(bsxfun(@ge,squeeze(results_inc_perm.fold_corr_perm),results_inc_perm.fold_corr'),2))./10000;

results_pe_perm = mc_bbs_perm(featuremat,mainpheno(:,2),[nuisance mainpheno(:,[3 4])],folds,bestcomps_pe,10000,'Scores',Aa,'Perms',pset(:,2:1001));
fold_p_pe = (sum(bsxfun(@ge,squeeze(results_pe_perm.fold_corr_perm),results_pe_perm.fold_corr'),2))./10000;

results_adi_perm = mc_bbs_perm(featuremat,mainpheno(:,4),[nuisance mainpheno(:,[2 3])],folds,bestcomps_adi,10000,'Scores',Aa,'Perms',pset(:,2:1001));
fold_p_adi = (sum(bsxfun(@ge,squeeze(results_adi_perm.fold_corr_perm),results_adi_perm.fold_corr'),2))./10000;
