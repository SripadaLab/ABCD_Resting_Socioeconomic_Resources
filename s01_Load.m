%%
%utility funcions
addpath /nfs/corenfs/psych-freewill-data/Users/mangstad/repos/Misc_utils

Exp = '/nfs/corenfs/psych-freewill-data/Data/ABCD/CIFTI/Scripts/ABCD_Resting_Socioeconomic_Resources/';

%%
%Setup files
DataFile = [Exp '/Data/ABCD_rest.csv'];
CorrTemplate = [Exp '/Data/Gordon_Sub_Cere/[Subject].txt'];
CorrTemplate = '/nfs/corenfs/psych-freewill-data/Data/ABCD/CIFTI/Rest/Gordon_Sub_Cere/[Subject].txt';
NetsFile = [Exp '/Data/gordon_sub_cere_parcels.csv'];

%%
%setup
dat = readtable(DataFile);

N = size(dat,1);
nROI = 418;
P = (nROI*(nROI-1))/2;

netsfile = readtable(NetsFile);
netsfile = netsfile(1:nROI,:);
nets = netsfile.NetworkNumber;

%%
%load connectomes
featuremat = zeros(N,P);

parfor iSubject = 1:N
    Subject = dat.Subject{iSubject};
    fprintf(1,'%d\n',iSubject);
    file = strrep(CorrTemplate,'[Subject]',Subject);
    try
        tmp = load(file);
        tmp = tmp(1:418,1:418);
        tmp = mc_flatten_upper_triangle(tmp);
        featuremat(iSubject,:) = mc_FisherZ(tmp);
    catch
    end
end

%%
%setup folds
u = unique(dat.abcd_site_num);
nFold = numel(u);
fold_site_conversion = [[1:nFold]' u];
folds = zeros(size(dat.abcd_site_num));
sitesize = zeros(nFold,1);
for iFold = 1:nFold
    folds(dat.abcd_site_num==u(iFold)) = iFold;
    sitesize(iFold) = sum(folds==iFold);
end

mainpheno = [dat.ses_fac dat.EdYearsAverage dat.Income2Needs dat.reshist_addr1_adi_wsum];

%%
%get nuisance variables
u = unique(dat.RaceEthnicity);
re = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.RaceEthnicity,u{i});
    re(idx,i) = 1;
    if (strcmp(u{i},'NaN'))
        re(idx,i) = NaN;
    end
end
s = sum(re);
i = isnan(s);
re(isnan(re(:,i)),:) = NaN;
%re(:,i) = [];
[~,i] = max(nansum(re));
re(:,5) = [];

u = unique(dat.Sex);
sex = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.Sex,u{i});
    sex(idx,i) = 1;
    if (strcmp(u{i},'NaN'))
        sex(idx,i) = NaN;
    end
end
s = sum(sex);
i = isnan(s);
sex(isnan(sex(:,i)),:) = NaN;
%sex(:,i) = [];
sex(:,2) = [];

%normal and expanded nuisance
nuisance = [orth_poly(dat.Age,2) sex orth_poly(dat.fd,2) re];

%generate exchangibility blocks and permutation orderings
addpath /nfs/corenfs/psych-freewill-data/Data/ABCD/sojoshi/Scripts/repos/ABCD_winkler/share
addpath /nfs/corenfs/psych-freewill-data/Data/ABCD/sojoshi/Scripts/repos/PALM

File = '/nfs/corenfs/psych-freewill-data/Data/ABCD/CIFTI/Scripts/rest_classify/Data/ABCD_zygosity_fornow.csv';

subs = dat.subjectkey;
EB = abcd2blocks(File,'',[100 10],subs,1);
tmp = [EB(:,1) -folds EB(:,2:end)];
pset = palm_quickperms([],tmp,10000,1,0,[],[]);

save('Results/pset.mat','pset');
