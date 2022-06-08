
%check all models in the same subset of subjects
allpheno = [dat.EdYearsAverage dat.G_lavaan dat.fes_p_ss_int_cult_sum_pr dat.avg_grades dat.crpbi_y_ss_parent dat.srpf_y_ss_ses];
good = ~any(isnan(allpheno),2) & ~any(isnan(nuisance),2);

aps = allpheno(good,:);
ns = nuisance(good,:);
fs = folds(good); %checked to confirm still same 19 sites
[c,~,~,l] = crosstab(fs);

Aas = {};
for i = 1:19
    Aas{i} = Aa{i}(good,:);
end

save('GL/Data/data.mat','aps','ns','fs','Aas','-v7.3');

Tests = {
    'PE';
    'PE w/ G';
    'PE w/ enrich';
    'PE w/ grades';
    'PE w/ family';
    'PE w/ school';
    'PE w/ all';
    };

NumComps = [50:50:1000];
results_pes = mc_bbs(fms,aps(:,1),[ns],fs,NumComps,'NoComponents',1,'Scores',Aas);
results_pes_g = mc_bbs(fms,aps(:,1),[ns aps(:,2)],fs,NumComps,'NoComponents',1,'Scores',Aas);
results_pes_enrich = mc_bbs(fms,aps(:,1),[ns aps(:,3)],fs,NumComps,'NoComponents',1,'Scores',Aas);
results_pes_grades = mc_bbs(fms,aps(:,1),[ns aps(:,4)],fs,NumComps,'NoComponents',1,'Scores',Aas);
results_pes_family = mc_bbs(fms,aps(:,1),[ns aps(:,5)],fs,NumComps,'NoComponents',1,'Scores',Aas);
results_pes_school = mc_bbs(fms,aps(:,1),[ns aps(:,6)],fs,NumComps,'NoComponents',1,'Scores',Aas);
results_pes_pf10 = mc_bbs(fms,aps(:,1),[ns aps(:,7)],fs,NumComps,'NoComponents',1,'Scores',Aas);
results_pes_all = mc_bbs(fms,aps(:,1),[ns aps(:,2:7)],fs,NumComps,'NoComponents',1,'Scores',Aas);

mean_corr = zeros(8,1);
mean_comp = zeros(8,1);
std_comp = zeros(8,1);
mean_corr(1) = results_pes.mean_corr;
mean_corr(2) = results_pes_g.mean_corr;
mean_corr(3) = results_pes_enrich.mean_corr;
mean_corr(4) = results_pes_grades.mean_corr;
mean_corr(5) = results_pes_family.mean_corr;
mean_corr(6) = results_pes_school.mean_corr;
mean_corr(7) = results_pes_pf10.mean_corr;
mean_corr(8) = results_pes_all.mean_corr;

mean_comp(1) = mean(results_pes.bestcomps);
mean_comp(2) = mean(results_pes_g.bestcomps);
mean_comp(3) = mean(results_pes_enrich.bestcomps);
mean_comp(4) = mean(results_pes_grades.bestcomps);
mean_comp(5) = mean(results_pes_family.bestcomps);
mean_comp(6) = mean(results_pes_school.bestcomps);
mean_comp(7) = mean(results_pes_pf10.bestcomps);
mean_comp(8) = mean(results_pes_all.bestcomps);

std_comp(1) = std(results_pes.bestcomps);
std_comp(2) = std(results_pes_g.bestcomps);
std_comp(3) = std(results_pes_enrich.bestcomps);
std_comp(4) = std(results_pes_grades.bestcomps);
std_comp(5) = std(results_pes_family.bestcomps);
std_comp(6) = std(results_pes_school.bestcomps);
std_comp(7) = std(results_pes_pf10.bestcomps);
std_comp(8) = std(results_pes_all.bestcomps);

out = table(Tests,mean_corr,mean_comp,std_comp);
out

cv_corr = load('GL/Results/cv_corr_00000.txt');
cv_comps = load('GL/Results/cv_comps_perfold_00000.txt');
mean_corr = cv_corr';
mean_comp = mean(cv_comps)';
std_comp = std(cv_comps)';
out = table(Tests,mean_corr,mean_comp,std_comp);
out
%run bootstrapping on HPC
%load results

bs_corr = load('GL/Results/allcorr.txt');

obs = bs_corr(1,:);
bs_corr = bs_corr(11:end,:);

pvals = [sum(bsxfun(@gt,bsxfun(@minus,bs_corr,bs_corr(:,1)),0))/10000]';
pvals(1) = NaN;
pvals(pvals==0) = 0.0001;
mean_corr = obs';
diff = mean_corr - mean_corr(1);
out = table(Tests([1:6 8]),mean_corr,diff,pvals);
out

shuf_order = [1,3,5,6,4,2,7];
out(shuf_order,:)
