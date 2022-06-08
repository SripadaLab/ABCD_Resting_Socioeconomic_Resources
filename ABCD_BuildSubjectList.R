rm(list=ls())
source("setup.R")
library(dplyr)

mergecols = c("subjectkey","eventname")

#load family data
fam = read.abcd(file.path(ABCDDataDir,"acspsw03.txt"))
fam = fam[,c(mergecols,"rel_family_id")]

#load site data
site= read.abcd(file.path(ABCDDataDir,"abcd_lt01.txt"))
site = site[,c(mergecols,"site_id_l")]

#load meanFD data and DOF totals
motion = read.csv(file.path(DataDir,"meanfd.csv"))
motion$dofavg = motion$dof/motion$runs
motion$censoredavg = motion$censored/motion$runs
motion$subjectkey = gsub("NDAR","NDAR_",motion$subjectkey)
motion$eventname = "baseline_year_1_arm_1"

#generate nuisance covariates from raw ABCD files
pdem = read.abcd(file.path(ABCDDataDir,"pdem02.txt"))
cols = c(mergecols,"interview_age","sex","demo_sex_v2",
         "demo_comb_income_v2","demo_prnt_marital_v2",
         "demo_prnt_ed_v2","demo_prtnr_ed_v2","demo_ethn_v2",
         paste0("demo_race_a_p___",c(10:25,77,99)))
pdem = pdem[,cols]
pdem$re_white = as.numeric(pdem$demo_race_a_p___10==1)
pdem$re_black = as.numeric(pdem$demo_race_a_p___11==1)
pdem$re_hisp = as.numeric(pdem$demo_ethn_v2==1)
pdem$re_hisp[is.na(pdem$re_hisp)] = 0

pdem$re_asian = as.numeric(rowSums(pdem[,c(paste0("demo_race_a_p___",c(14:24)))])>0)
pdem$re_other = as.numeric(rowSums(pdem[,c(paste0("demo_race_a_p___",c(12:13,25)))])>0)

temp = with(pdem,re_asian + 2*re_black + 4*re_hisp + 8*re_other + 16*re_white)
temp[temp==1] = "Asian"
temp[temp==2] = "Black"
temp[temp==4] = "Hispanic"
temp[temp==8] = "Other"
temp[temp==16] = "White"
b=1:31
d = c(1:31)[bitwAnd(b,4)>0]
temp[temp %in% d] = "Hispanic"
temp[!(temp %in% c("Asian","Black","Hispanic","Other","White"))] = "Other"
pdem$RaceEthnicity = as.factor(temp)

pdem$Sex = as.factor(pdem$sex)

pdem$Age = pdem$interview_age/12

pdem = pdem[,c(mergecols,"Age","Sex","RaceEthnicity")]

ses = read.csv(file.path(DataDir,"ABCD_ses.csv"))

#grab other variables for grades/enrichment/etc
#g
#pf10

#recode fes eventname so that we can use the 2 year data during baseline
fes = read.abcd(file.path(ABCDDataDir,"abcd_sscep01.txt"))
fes = fes[,c(mergecols,"fes_p_ss_int_cult_sum_pr")]
fes = fes[fes$eventname=="2_year_follow_up_y_arm_1",]
fes$eventname = "baseline_year_1_arm_1"

sag1 = read.abcd(file.path(ABCDDataDir,"abcd_ysaag01.txt"))
sag1 = sag1[,c(mergecols,"sag_grades_last_yr")]
sag1 = sag1[sag1$eventname=="2_year_follow_up_y_arm_1",]
sag1$eventname = "baseline_year_1_arm_1"

sag2 = read.abcd(file.path(ABCDDataDir,"abcd_saag01.txt"))
sag2 = sag2[,c(mergecols,"sag_grade_type")]
sag2 = sag2[sag2$eventname=="2_year_follow_up_y_arm_1",]
sag2$eventname = "baseline_year_1_arm_1"
sag = merge(sag1,sag2)
sag$sag_grade_type = na_if(sag$sag_grade_type,-1)
sag$sag_grade_type = na_if(sag$sag_grade_type,777)
sag$sag_grades_last_yr = na_if(sag$sag_grades_last_yr,-1)
sag$sag_grades_last_yr = na_if(sag$sag_grades_last_yr,777)

grades_dict = c((100+97) / 2,(96+93) / 2,(92+90) / 2,(89+87) / 2,(86+83) / 2,(82+80) / 2,(79+77) / 2,(76+73) / 2,(72+70) / 2,(69+67) / 2,(66+65) / 2,65)
sag$g1_midpoint = grades_dict[sag$sag_grade_type]
sag$g2_midpoint = grades_dict[sag$sag_grades_last_yr]
sag$avg_grades = rowMeans(sag[,c("g1_midpoint","g2_midpoint")],na.rm=T)

sscey = read.abcd(file.path(ABCDDataDir,"abcd_sscey01.txt"))
sscey = sscey[,c(mergecols,"crpbi_y_ss_parent","srpf_y_ss_ses")]

g = read.csv(file.path(DataDir,"ABCD_lavaan_gfactor_loso.csv"))
g$eventname = "baseline_year_1_arm_1"
g = g[,c(mergecols,"G")]
names(g)[3] = "G_lavaan"

p = read.csv(file.path(DataDir,"ABCD_lavaan_pfactor_loso.csv"))
p$eventname = "baseline_year_1_arm_1"
p = p[,c(mergecols,"PF10_lavaan")]

#subjective financial and safety measures
pdem02 = read.abcd(file.path(ABCDDataDir,"pdem02.txt"))
BNU.cols = paste0("demo_fam_exp",1:7,"_v2")
pdem02[,BNU.cols]=apply(pdem02[,BNU.cols], 2, function(x) gsub("777", NA, x) )
pdem02[,BNU.cols]=apply(pdem02[,BNU.cols], 2, function(x) as.numeric(x))
pdem02$basic_needs_sum=rowSums(pdem02[,BNU.cols])
pdem02 = pdem02[,c(mergecols,"basic_needs_sum")]

pnsc = read.abcd(file.path(ABCDDataDir,"abcd_pnsc01.txt"))
crime.cols = paste0("neighborhood",1:3,"r_p")
#pnsc[,crime.cols]=apply(pnsc[,crime.cols], 2, function(x) gsub("777", NA, x) )
pnsc[,crime.cols]=apply(pnsc[,crime.cols], 2, function(x) as.numeric(x))
pnsc$sub_crime_sum = rowSums(pnsc[,crime.cols])
pnsc = pnsc[,c(mergecols,"sub_crime_sum")]

#merge everything together
data = multi.merge(fam,site,motion,pdem,ses,fes,sag,sscey,pdem02,pnsc,g,p,by=mergecols)

data = data[data$eventname=="baseline_year_1_arm_1",]

tr = 0.8
data$GoodTime = (data$TRs - data$censored)*tr/60
data$Include.rest = data$GoodTime >= 8 & data$runs>=2

#exclude for any data problems, like NaN ROIs etc
nansubs = read.csv(file.path(DataDir,"nan_subs.csv"))
data$Include.data = !(data$subjectkey %in% nansubs$subjectkey)

data$Include = data$Include.rest & data$Include.data

sum(data$Include,na.rm=T)

data = data[data$Include==T & !is.na(data$Include),]
data$Subject = gsub("NDAR_","NDAR",data$subjectkey)

sum(data$Include)

#exclude for missing ses variables
good = !is.na(data$Income2Needs) & !is.na(data$EdYearsAverage) & !is.na(data$reshist_addr1_adi_wsum) & !is.na(data$ses_fac)

data = data[good,]

#check for families that cross site and drop
t = table(data$rel_family_id,data$site_id_l)
t = t>0
sum(rowSums(t)>1)
fams = rownames(t)[rowSums(t)>1]
data = data[!(data$rel_family_id %in% fams),]

data$abcd_site = as.character(data$site_id_l)
table(data$abcd_site)
data$abcd_site_num = as.numeric(gsub('site','',data$abcd_site))

t = table(data$abcd_site)
sites = names(t)[t>=75]
data = data[data$abcd_site %in% sites,]

dim(data)

sum(rowSums(is.na(data[,c("EdYearsAverage","G_lavaan","Age","Sex","RaceEthnicity","fd","fes_p_ss_int_cult_sum_pr","avg_grades","srpf_y_ss_ses","crpbi_y_ss_parent","basic_needs_sum","sub_crime_sum")]))==0)

write.csv(data,file.path(DataDir,"ABCD_rest.csv"),row.names=FALSE,na="NaN")
