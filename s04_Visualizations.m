[sitesize,results_ses.fold_corr']
figure;
scatter(sitesize,results_ses.fold_corr);

[~,si] = sort(results_ses.fold_corr,'descend');
for i = 1:4
    figure;
    scatter(results_ses.pheno_residualized(folds==si(i)),results_ses.pheno_predict(folds==si(i)));
end

subjectkey = dat.subjectkey;
site_id_l = dat.site_id_l;
actual_ses = results_ses.pheno_residualized;
predicted_ses = results_ses.pheno_predict;
out = table(subjectkey,site_id_l,actual_ses,predicted_ses);
writetable(out,'Results/scatter_data.csv');


figure;
hist(results_pe_perm.mean_corr_perm(2:end),30)
hold on;
plot([results_pe_perm.mean_corr results_pe_perm.mean_corr],[0 1000],'r-');


