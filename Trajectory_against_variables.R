library(gplots)
library(RColorBrewer)

peardist <- function(x,distance){
      as.dist(1-cor(t(x)))
}

has_subvec = function(x,s){
  xL = length(x); sL = length(s)
  if(xL < sL) return(FALSE)
  any(sapply(1:(xL-sL+1),function(i){
    isTRUE(all.equal(x[i:(i+sL-1)],s))
  }))
}
#Models--------
#mdrop0<-read.delim('Data_for_trajectory_tests_3.txt',sep='\t')
infection
mdrop<-read.delim('Data_for_trajectory_tests_8.txt',sep='\t')
mdrop<-read.delim('Data_for_trajectory_tests_8b.txt',sep='\t')
mdrop$Traj<-mdrop$trajec
mdrop$T1<-as.numeric(mdrop$Traj)
mdrop$T1[mdrop$Traj!='1']<-0
mdrop$T2<-as.numeric(mdrop$Traj)
mdrop$T2[mdrop$Traj!='2a']<-0
mdrop$T2[mdrop$Traj=='2a']<-1
mdrop$T3<-as.numeric(mdrop$Traj)
mdrop$T3[mdrop$Traj!='2c']<-0
mdrop$T3[mdrop$Traj=='2c']<-1
mdrop$T4<-as.numeric(mdrop$Traj)
mdrop$T4[mdrop$Traj!='2b']<-0
mdrop$T4[mdrop$Traj=='2b']<-1
mdrop$T5<-as.numeric(mdrop$Traj)
mdrop$T5[mdrop$Traj!='3']<-0
mdrop$T5[mdrop$Traj=='3']<-1

# mtest<-mdrop[,c('birthmode_RecievingIAP','sibs','familly_NbSiblingsFt',
#                 'gestage','env_HabitationType','env_FurOrFeathers','inf_AgeSolidsIntroduced',
#                 "m_probio","birth_season",'exclusiveBF_6months','exclusiveBF_12months',
#                 #"m_FattyAcid","m_FolicAcid",
#                 "m_AgeDelivery",'antibiotic3m',#'antibiotic6m',
#                 'probiotic3m','probiotic6m',
#                 'inf_Prob_Lactobacillus3m','inf_Prob_Lactobacillus6m',
#                 'inf_Prob_Bifidobacterium3m','inf_Prob_Bifidobacterium6m',
#                 'T1','T2','T3')]

mtest<-mdrop[,c('birthmode_RecievingIAP','sibs','familly_NbSiblingsFt',
                'gestage','env_HabitationType','env_FurOrFeathers',
                'inf_AgeSolidsIntroduced',#'inf_AgeFormulaIntroduced',
                "m_probio","birth_season",'exclusiveBF_6months','exclusiveBF_12months',
                #"m_FattyAcid","m_FolicAcid",
                "m_AgeDelivery",'antibiotic3m',#'antibiotic6m',
                'probiotic3m','probiotic6m',
                'inf_Prob_Lactobacillus3m','inf_Prob_Lactobacillus6m',
                'inf_Prob_Bifidobacterium3m','inf_Prob_Bifidobacterium6m',
                'T1','T2','T3','T4','T5')]


mtest$inf_AgeSolidsIntroduced[mtest$inf_AgeSolidsIntroduced>6]<-6

mtest[,c(3,4,8,12:15)]<-scale(mtest[,c(3,4,8,12:15)])
#mtest[,c(3,4,9,13:16)]<-scale(mtest[,c(3,4,9,13:16)])
mtest$birthmode_RecievingIAP[mtest$birthmode_RecievingIAP=='Vaginal_NoAB']<-'a_Vaginal_NoAB'
mtest<-na.omit(mtest)


#Do the tests again but for all the different birthmodes. 
#Try the partition tree again with the new trajectories and limited variables
#Look at the probiotics in the first 3 months like I did for AB. Also try the Prob sepertaly

# T1mod<- step(glm(T1~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
#                     exclusiveBF_6months+exclusiveBF_12months+
#                     #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
#                     #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
#                     gestage+env_HabitationType+env_FurOrFeathers+
#                     probiotic3m+#probiotic6m+
#                     m_probio+#birth_season+
#                     m_AgeDelivery+
#                     inf_AgeSolidsIntroduced+ inf_AgeFormulaIntroduced+
#                     antibiotic3m+#antibiotic6m+
#                     inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
#                     inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
#                     #feeding_bottle_ever+feeding_formula_ever
#                     ,
#                   data=mtest,family="binomial"))
T1mod_ind_pvalue<-list()
T1mod_ind_est<-list()
T2mod_ind_pvalue<-list()
T2mod_ind_est<-list()
T3mod_ind_pvalue<-list()
T3mod_ind_est<-list()
T4mod_ind_pvalue<-list()
T4mod_ind_est<-list()
T5mod_ind_pvalue<-list()
T5mod_ind_est<-list()
mod_ind_names<-list()
for (ii in 1:19) {
  mod_ind_names[[ii]]<-rownames(summary(glm(as.formula(paste("T1 ~", names(mtest)[ii])),data=mtest,family="binomial"))$coef)[-1]
  T1mod_ind_pvalue[[ii]]<-as.numeric(summary(glm(as.formula(paste("T1 ~", names(mtest)[ii])),data=mtest,family="binomial"))$coef[-1,4])
  T1mod_ind_est[[ii]]<-as.numeric(summary(glm(as.formula(paste("T1 ~", names(mtest)[ii])),data=mtest,family="binomial"))$coef[-1,1])
  
  T2mod_ind_pvalue[[ii]]<-(summary(glm(as.formula(paste("T2 ~", names(mtest)[ii])),data=mtest,family="binomial"))$coef[-1,4])
  T2mod_ind_est[[ii]]   <-as.numeric(summary(glm(as.formula(paste("T2 ~", names(mtest)[ii])),data=mtest,family="binomial"))$coef[-1,1])
  
  T3mod_ind_pvalue[[ii]]<-as.numeric(summary(glm(as.formula(paste("T3 ~", names(mtest)[ii])),data=mtest,family="binomial"))$coef[-1,4])
  T3mod_ind_est[[ii]]   <-as.numeric(summary(glm(as.formula(paste("T3 ~", names(mtest)[ii])),data=mtest,family="binomial"))$coef[-1,1])
  
  T4mod_ind_pvalue[[ii]]<-as.numeric(summary(glm(as.formula(paste("T4 ~", names(mtest)[ii])),data=mtest,family="binomial"))$coef[-1,4])
  T4mod_ind_est[[ii]]   <-as.numeric(summary(glm(as.formula(paste("T4 ~", names(mtest)[ii])),data=mtest,family="binomial"))$coef[-1,1])
  
  T5mod_ind_pvalue[[ii]]<-as.numeric(summary(glm(as.formula(paste("T5 ~", names(mtest)[ii])),data=mtest,family="binomial"))$coef[-1,4])
  T5mod_ind_est[[ii]]   <-as.numeric(summary(glm(as.formula(paste("T5 ~", names(mtest)[ii])),data=mtest,family="binomial"))$coef[-1,1])
  }


IndMods_est<-as.data.frame(cbind(unlist(mod_ind_names),
               unlist(T1mod_ind_est),
               unlist(T2mod_ind_est),
               unlist(T3mod_ind_est),
               unlist(T4mod_ind_est),
               unlist(T5mod_ind_est)
               ))
IndMods_pval<-as.data.frame(cbind(unlist(mod_ind_names),
               as.numeric(unlist(T1mod_ind_pvalue)),
               as.numeric(unlist(T2mod_ind_pvalue)),
               as.numeric(unlist(T3mod_ind_pvalue)),
               as.numeric(unlist(T4mod_ind_pvalue)),
               as.numeric(unlist(T5mod_ind_pvalue))
               ))

names(IndMods_est)<-c('ModelResp',
                  'T1 est',
                  'T2 est',
                  'T3 est',
                  'T4 est',
                  'T5 est')
names(IndMods_pval)<-c('ModelResp',
                  'T1 pvalue','T2 pvalue','T3 pvalue',
                  'T4 pvalue','T5 pvalue')
rownames(IndMods_est)<-IndMods_est$ModelResp
rownames(IndMods_pval)<-IndMods_est$ModelResp
IndMods_est<-IndMods_est[,-which(names(IndMods_est) %in% c('ModelResp'))]
IndMods_pval<-IndMods_pval[,-which(names(IndMods_pval) %in% c('ModelResp'))]

for (ii in 1:nrow(IndMods_pval)) {
  for (jj in 1:ncol(IndMods_pval)) {
    IndMods_est[ii,jj]<-as.numeric(IndMods_est[ii,jj])
    if (as.numeric(IndMods_pval[ii,jj])>0.05) {
      IndMods_est[ii,jj]<-NA
    }
  }
}


IndMods_est <- mutate_all(IndMods_est, function(x) as.numeric(as.character(x)))
IndMods_est<-IndMods_est[rowSums(is.na(IndMods_est)) == 0, ]
df<-IndMods_est
df<-df[rowSums(is.na(df)) != ncol(df), ]

df %>% drop_na() 

T1mod<- step(glm(T1~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
                    exclusiveBF_6months+exclusiveBF_12months+
                    #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
                    #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
                    gestage+env_HabitationType+env_FurOrFeathers+
                    probiotic3m+#probiotic6m+
                    #m_probio+#birth_season+
                    m_AgeDelivery+
                    inf_AgeSolidsIntroduced+#inf_AgeFormulaIntroduced+
                    antibiotic3m+#antibiotic6m+
                    inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
                    inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
                    #feeding_bottle_ever+feeding_formula_ever
                    ,
                  data=mtest,family="binomial"))


T2mod<- step(glm(T2~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
                    exclusiveBF_6months+
                    exclusiveBF_12months+
                    #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
                    #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
                    gestage+env_HabitationType+env_FurOrFeathers+
                    probiotic3m+#probiotic6m+
                    m_probio+#birth_season+
                    m_AgeDelivery+
                    inf_AgeSolidsIntroduced+#inf_AgeFormulaIntroduced+
                    antibiotic3m+#antibiotic6m+
                    inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
                    inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
                    #feeding_bottle_ever+feeding_formula_ever
                    #feeding_bottle_ever+feeding_formula_ever
                    ,
                  data=mtest,family="binomial"))

T3mod<- step(glm(T3~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
                    exclusiveBF_6months+
                    exclusiveBF_12months+
                    #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
                    #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
                    gestage+env_HabitationType+env_FurOrFeathers+
                    probiotic3m+#probiotic6m+
                    m_probio+#birth_season+
                    m_AgeDelivery+
                    inf_AgeSolidsIntroduced+#inf_AgeFormulaIntroduced+
                    antibiotic3m+#antibiotic6m+
                    inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
                    inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
                    #feeding_bottle_ever+feeding_formula_ever
                    #feeding_bottle_ever+feeding_formula_ever
                    ,
                  data=mtest,family="binomial"))

T4mod<- step(glm(T4~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
                    #exclusiveBF_6months+
                    exclusiveBF_12months+
                    #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
                    #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
                    gestage+env_HabitationType+env_FurOrFeathers+
                    probiotic3m+#probiotic6m+
                    m_probio+#birth_season+
                    m_AgeDelivery+
                    inf_AgeSolidsIntroduced+#inf_AgeFormulaIntroduced+
                    antibiotic3m+#antibiotic6m+
                    inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
                    inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
                    #feeding_bottle_ever+feeding_formula_ever
                    #feeding_bottle_ever+feeding_formula_ever
                    ,
                  data=mtest,family="binomial"))

T5mod<- step(glm(T5~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
                    #exclusiveBF_6months+
                    exclusiveBF_12months+
                    #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
                    #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
                    gestage+env_HabitationType+env_FurOrFeathers+
                    probiotic3m+#probiotic6m+
                    m_probio+#birth_season+
                    m_AgeDelivery+
                    inf_AgeSolidsIntroduced+#inf_AgeFormulaIntroduced+
                    antibiotic3m+#antibiotic6m+
                    inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
                    inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
                    #feeding_bottle_ever+feeding_formula_ever
                    #feeding_bottle_ever+feeding_formula_ever
                    ,
                  data=mtest,family="binomial"))

#old models.-----
# T3mod <- step(glm(T3~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
#                     exclusiveBF_6months+exclusiveBF_12months+
#                     #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
#                     #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
#                     gestage+env_HabitationType+env_FurOrFeathers+
#                     probiotic3m+#probiotic6m+
#                     m_probio+#birth_season+
#                     m_AgeDelivery+inf_AgeSolidsIntroduced+
#                     antibiotic3m+#antibiotic6m+
#                     inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
#                     inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
#                     #feeding_bottle_ever+feeding_formula_ever
#                     ,
#                   data=mtest,family="binomial"))
# 
# T2mod_to1 <- step(glm(T2~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
#                     exclusiveBF_6months+exclusiveBF_12months+
#                     #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
#                     #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
#                     gestage+env_HabitationType+env_FurOrFeathers+
#                     probiotic3m+#probiotic6m+
#                     m_probio+#birth_season+
#                     m_AgeDelivery+inf_AgeSolidsIntroduced+
#                     antibiotic3m+#antibiotic6m+
#                     inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
#                     inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
#                     #feeding_bottle_ever+feeding_formula_ever
#                     ,
#                   data=mtest[mtest$T1==1|mtest$T2==1,],family="binomial"))
# 
# T3mod_to1 <- step(glm(T3~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
#                     exclusiveBF_6months+#exclusiveBF_12months+
#                     #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
#                     #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
#                     gestage+env_HabitationType+env_FurOrFeathers+
#                     probiotic3m+#probiotic6m+
#                     m_probio+#birth_season+
#                     m_AgeDelivery+inf_AgeSolidsIntroduced+
#                     #antibiotic3m+#antibiotic6m+
#                     inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
#                     inf_Prob_Bifidobacterium3m#+inf_Prob_Bifidobacterium6m
#                     #feeding_bottle_ever+feeding_formula_ever
#                     ,
#                   data=mtest[mtest$T1==1|mtest$T3==1,],family="binomial"))
# 
# T2mod_to3 <- step(glm(T2~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
#                     exclusiveBF_6months+exclusiveBF_12months+
#                     #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
#                     #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
#                     gestage+env_HabitationType+env_FurOrFeathers+
#                     probiotic3m+#probiotic6m+
#                     m_probio+#birth_season+
#                     m_AgeDelivery+inf_AgeSolidsIntroduced+
#                     antibiotic3m+#antibiotic6m+
#                     inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
#                     inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
#                     #feeding_bottle_ever+feeding_formula_ever
#                     ,
#                   data=mtest[mtest$T2==1|mtest$T3==1,],family="binomial"))


# T1mod_CsecVab<-step(glm(T1~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
#                     exclusiveBF_6months+exclusiveBF_12months+
#                     #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
#                     #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
#                     gestage+env_HabitationType+env_FurOrFeathers+
#                     probiotic3m+#probiotic6m+
#                     m_probio+#birth_season+
#                     m_AgeDelivery+inf_AgeSolidsIntroduced+
#                     antibiotic3m+#antibiotic6m+
#                     inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
#                     inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
#                     #feeding_bottle_ever+feeding_formula_ever
#                     ,
#                   data=mtest[mtest$birthmode_RecievingIAP=='Vaginal_AB'| mtest$birthmode_RecievingIAP=='C-section',],family="binomial"))
# 
# T2mod_CsecVab<-step(glm(T2~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
#                     exclusiveBF_6months+#exclusiveBF_12months+
#                     #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
#                     #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
#                     gestage+env_HabitationType+env_FurOrFeathers+
#                     probiotic3m+#probiotic6m+
#                     m_probio+#birth_season+
#                     m_AgeDelivery+inf_AgeSolidsIntroduced+
#                     antibiotic3m+#antibiotic6m+
#                     inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
#                     inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
#                     #feeding_bottle_ever+feeding_formula_ever
#                     ,
#                   data=mtest[mtest$birthmode_RecievingIAP=='Vaginal_AB'| mtest$birthmode_RecievingIAP=='C-section',],family="binomial")
#                   )
# 
# 
# T3mod_CsecVab<-step(glm(T3~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
#                     exclusiveBF_6months+exclusiveBF_12months+
#                     #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
#                     #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
#                     gestage+env_HabitationType+env_FurOrFeathers+
#                     probiotic3m+#probiotic6m+
#                     m_probio+#birth_season+
#                     m_AgeDelivery+inf_AgeSolidsIntroduced+
#                     antibiotic3m+#antibiotic6m+
#                     inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
#                     inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
#                     #feeding_bottle_ever+feeding_formula_ever
#                     ,
#                   data=mtest[mtest$birthmode_RecievingIAP=='Vaginal_AB'| mtest$birthmode_RecievingIAP=='C-section',],family="binomial"))
# 
# T2mod_to1_CsecVab<-step(glm(T2~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
#                     exclusiveBF_6months+exclusiveBF_12months+
#                     #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
#                     #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
#                     gestage+env_HabitationType+env_FurOrFeathers+
#                     probiotic3m+#probiotic6m+
#                     m_probio+#birth_season+
#                     m_AgeDelivery+inf_AgeSolidsIntroduced+
#                     antibiotic3m+#antibiotic6m+
#                     inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
#                     inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
#                     #feeding_bottle_ever+feeding_formula_ever
#                     ,
#                   data=mtest[(mtest$T2==1|mtest$T1==1) & (mtest$birthmode_RecievingIAP=='Vaginal_AB'| mtest$birthmode_RecievingIAP=='C-section'),],family="binomial"))
# 
# T3mod_to1_CsecVab<-step(glm(T3~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
#                     exclusiveBF_6months+exclusiveBF_12months+
#                     #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
#                     #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
#                     gestage+env_HabitationType+env_FurOrFeathers+
#                     probiotic3m+#probiotic6m+
#                     m_probio+#birth_season+
#                     m_AgeDelivery+inf_AgeSolidsIntroduced+
#                     #antibiotic3m+#antibiotic6m+
#                     inf_Prob_Lactobacillus3m+#inf_Prob_Lactobacillus6m+
#                     inf_Prob_Bifidobacterium3m#+inf_Prob_Bifidobacterium6m
#                     #feeding_bottle_ever+feeding_formula_ever
#                     ,
#                   data=mtest[(mtest$T3==1|mtest$T1==1) & (mtest$birthmode_RecievingIAP=='Vaginal_AB'| mtest$birthmode_RecievingIAP=='C-section'),],family="binomial"))
# 
# 
# T2mod_to3_CsecVab  <- step(glm(T2~birthmode_RecievingIAP+sibs+familly_NbSiblingsFt+
#                     exclusiveBF_6months+exclusiveBF_12months+
#                     #inf_AgeFormulaIntroduced+inf_FirstAntibiotic_age+
#                     #m_alc+m_edu+p_edu+m_SmokingPrior+m_FattyAcid+m_FolicAcid+
#                     gestage+env_HabitationType+env_FurOrFeathers+
#                     probiotic3m+#probiotic6m+
#                     m_probio+#birth_season+
#                     m_AgeDelivery+inf_AgeSolidsIntroduced+
#                     antibiotic3m+#antibiotic6m+
#                     inf_Prob_Lactobacillus3m+inf_Prob_Lactobacillus6m+
#                     inf_Prob_Bifidobacterium3m+inf_Prob_Bifidobacterium6m
#                     #feeding_bottle_ever+feeding_formula_ever
#                     ,
#                   data=mtest[(mtest$T3==1|mtest$T2==1) & (mtest$birthmode_RecievingIAP=='Vaginal_AB'|mtest$birthmode_RecievingIAP=='C-section'),],family="binomial"))
# 
# 


#Cseec and Vab
# t1CsecVab<-as.data.frame(summary(T1mod_CsecVab)$coefficients)
# names(t1CsecVab)[1]<-'T1CsecVABest'
# names(t1CsecVab)[4]<-'T1CsecVABp'
# t1CsecVab<-t1CsecVab[2:nrow(t1CsecVab),c(1,4)]
# t1CsecVab$variables<-rownames(t1CsecVab)
# 
# t2CsecVab<-as.data.frame(summary(T2mod_CsecVab)$coefficients)
# names(t2CsecVab)[1]<-'T2CsecVABest'
# names(t2CsecVab)[4]<-'T2CsecVABp'
# t2CsecVab<-t2CsecVab[2:nrow(t2CsecVab),c(1,4)]
# t2CsecVab$variables<-rownames(t2CsecVab)
# 
# t3CsecVab<-as.data.frame(summary(T3mod_CsecVab)$coefficients)
# names(t3CsecVab)[1]<-'T3CsecVABest'
# names(t3CsecVab)[4]<-'T3CsecVABp'
# t3CsecVab<-t3CsecVab[2:nrow(t3CsecVab),c(1,4)]
# t3CsecVab$variables<-rownames(t3CsecVab)
# 
# t2_1CsecVab<-as.data.frame(summary(T2mod_to1_CsecVab)$coefficients)
# names(t2_1CsecVab)[1]<-'T2_1CsecVABest'
# names(t2_1CsecVab)[4]<-'T2_1CsecVABp'
# t2_1CsecVab<-t2_1CsecVab[2:nrow(t2_1CsecVab),c(1,4)]
# t2_1CsecVab$variables<-rownames(t2_1CsecVab)
# 
# t3_1CsecVab<-as.data.frame(summary(T3mod_to1_CsecVab)$coefficients)
# names(t3_1CsecVab)[1]<-'T3_1CsecVABest'
# names(t3_1CsecVab)[4]<-'T3_1CsecVABp'
# t3_1CsecVab<-t3_1CsecVab[2:nrow(t3_1CsecVab),c(1,4)]
# t3_1CsecVab$variables<-rownames(t3_1CsecVab)
# 
# t2_3CsecVab<-as.data.frame(summary(T2mod_to3_CsecVab)$coefficients)
# names(t2_3CsecVab)[1]<-'T2_3CsecVABest'
# names(t2_3CsecVab)[4]<-'T2_3CsecVABp'
# t2_3CsecVab<-t2_3CsecVab[2:nrow(t2_3CsecVab),c(1,4)]
# t2_3CsecVab$variables<-rownames(t2_3CsecVab)
# 
#PLotting----
#All
t1c<-as.data.frame(summary(T1mod)$coefficients)
names(t1c)[1]<-'T1est'
names(t1c)[4]<-'T1p'
t1c<-t1c[2:nrow(t1c),c(1,4)]
t1c$variables<-rownames(t1c)
t2c<-as.data.frame(summary(T2mod)$coefficients)
names(t2c)[1]<-'T2est'
names(t2c)[4]<-'T2p'
t2c<-t2c[2:nrow(t2c),c(1,4)]
t2c$variables<-rownames(t2c)
t3c<-as.data.frame(summary(T3mod)$coefficients)
names(t3c)[1]<-'T3est'
names(t3c)[4]<-'T3p'
t3c<-t3c[2:nrow(t3c),c(1,4)]
t3c$variables<-rownames(t3c)
t4c<-as.data.frame(summary(T4mod)$coefficients)
names(t4c)[1]<-'T4est'
names(t4c)[4]<-'T4p'
t4c<-t4c[2:nrow(t4c),c(1,4)]
t4c$variables<-rownames(t4c)
t5c<-as.data.frame(summary(T5mod)$coefficients)
names(t5c)[1]<-'T5est'
names(t5c)[4]<-'T5p'
t5c<-t5c[2:nrow(t5c),c(1,4)]
t5c$variables<-rownames(t5c)

df_list <- list(t1c, t2c, t3c, t4c, t5c)

#merge all data frames in list
dfm<-Reduce(function(x, y) merge(x, y, all=TRUE,by='variables'), df_list)
labs<-dfm[, -grep("est", colnames(dfm))]
dfm<-dfm[, -grep("p", colnames(dfm))]
dfm[is.na(dfm)]<-0



rownames(dfm)<-dfm$variables
rownames(labs)<-labs$variables
dfm<-dfm[ , -which(names(dfm) %in% c("variables"))]
labs<-labs[ , -which(names(labs) %in% c("variables"))]
rownames(dfm)[rownames(dfm)=="antibiotic3m"]<-'Antibiotics in first 3 months'
rownames(dfm)[rownames(dfm)=="antibiotic6m"]<-'Antibiotics in first 6 months'
rownames(dfm)[rownames(dfm)=="birthmode_RecievingIAPC-section"]<-"C-section birth"
rownames(dfm)[rownames(dfm)=="birthmode_RecievingIAPVaginal_AB"]<-"Vaginal birth (intrapartum antibiotics)"
rownames(dfm)[rownames(dfm)=="sibs"]<-"Number of siblings"
rownames(dfm)[rownames(dfm)=="feeding_bottle_everyes"]<-"Bottle feeding"
rownames(dfm)[rownames(dfm)=="feeding_formula_everyes"]<-"Formula feeding"
rownames(dfm)[rownames(dfm)=="inf_AgeSolidsIntroduced"]<-"Age solids introduced"
rownames(dfm)[rownames(dfm)=="inf_Prob_Bifidobacterium6m"]<-'Probiotic (Bifidobacterium) in first 6 months'
rownames(dfm)[rownames(dfm)=="inf_Prob_Lactobacillus3m"]<-'Probiotic (Lactobacillus) in first 3 months'
rownames(dfm)[rownames(dfm)=="inf_Prob_Lactobacillus6m"]<-'Probiotic (Lactobacillus) in first 6 months'
rownames(dfm)[rownames(dfm)=="m_AgeDelivery"]<-'Mothers age at delivery'
rownames(dfm)[rownames(dfm)=="m_probio"]<-'Mother used probiotics'
rownames(dfm)[rownames(dfm)=="probiotic6m"]<-'Probiotic (all) use in first 6 months'
rownames(dfm)[rownames(dfm)=="probiotic3m"]<-'Probiotic (all) use in first 3 months'
rownames(dfm)[rownames(dfm)=="sibsyes"]<-'Siblings present'
rownames(dfm)[rownames(dfm)=="env_FurOrFeathersTRUE"]<-'Pets present'
rownames(dfm)[rownames(dfm)=="exclusiveBF_12monthsyes"]<-'No formula in first 12 months'
rownames(dfm)[rownames(dfm)=="exclusiveBF_6monthsyes"]<-'No formula in first 6 months'
rownames(dfm)[rownames(dfm)=="inf_Prob_Bifidobacterium3m"]<-'Probiotic (Bifidobacterium) in first 3 months'
rownames(dfm)[rownames(dfm)=="gestage"]<-'Gestational age'
rownames(dfm)[rownames(dfm)=="familly_NbSiblingsFt"]<-'Number of siblings'

dfm[dfm>0]<-dfm[dfm>0]+0.1
dfm[dfm<0]<-dfm[dfm<0]-0.1

for (ii in seq(nrow(dfm))) {
  for (jj in seq(ncol(dfm))) {
    if (dfm[ii,jj]>0){dfm[ii,jj]<-dfm[dfm>0,]+0.25}
    if (dfm[ii,jj]<0){dfm[ii,jj]<-dfm[ii,jj]-0.25}
  }
}


colnames(dfm)<-c('T1','T2','T3','T4','T5')
l0<-labs
labs[l0<=0.1]<-'.'
labs[l0<0.05]<-'*'
labs[l0<0.01]<-'**'
labs[l0<0.001]<-'***'
labs[l0>=0.05 |is.na(l0)]<-''
labs<-as.matrix(labs)

pdf('Trajectory_heatmap_with5Trajectories_v2.pdf',width=12,height = 6)
#quartz(width=8,height=6);
par(oma=c(5,1,1,15))
 heatmap.2(as.matrix(dfm),dendrogram='none',
           Rowv=FALSE, Colv=FALSE,trace='none',
           col=rev(brewer.pal(n=11,'RdBu')),
           #lhei=c(2,6), lwid=c(2,2.5),
           lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2, 6), lwid=c(1, 2.5, 1),
           density.info = 'none',
           cellnote=labs, notecol="black", notecex=2,
           key.xtickfun = function() {
            breaks = pretty(parent.frame()$breaks)
            breaks = breaks[c(1:length(breaks))]
            list(at = parent.frame()$scale01(breaks),
                 labels = breaks)
          })
dev.off()


#Percentage heatmap

mdrop<-read.delim('Data_for_trajectory_tests_8b.txt',sep='\t')

mtest<-mdrop[,c('birthmode_RecievingIAP','sibs','familly_NbSiblingsFt',
                'env_HabitationType','env_FurOrFeathers','inf_AgeSolidsIntroduced','inf_AgeFormulaIntroduced',
                "m_probio",'exclusiveBF_6months','exclusiveBF_12months',
                "m_AgeDelivery",'antibiotic3m',#'antibiotic6m',
                'probiotic3m','probiotic6m',
                'inf_Prob_Lactobacillus3m','inf_Prob_Lactobacillus6m',
                'inf_Prob_Bifidobacterium3m','inf_Prob_Bifidobacterium6m',
                'trajec')]


mtest$inf_AgeSolidsIntroduced[mtest$inf_AgeSolidsIntroduced>6]<-6

mtest$Csection<-0
mtest$Csection[mtest$birthmode_RecievingIAP=='C-section']<-1
mtest$VABX<-0
mtest$VABX[mtest$birthmode_RecievingIAP=='Vaginal_AB']<-1
mtest$Vag<-0
mtest$Vag[mtest$birthmode_RecievingIAP=='Vaginal_NoAB']<-1

mtest$Apartment<-0
mtest$Apartment[mtest$env_HabitationType=='Appartment']<-1
mtest$Apartment<-0
mtest$Apartment[mtest$env_HabitationType=='Appartment']<-1
mtest$House<-0
mtest$House[mtest$env_HabitationType=='House']<-1
mtest$Townhouse<-0
mtest$Townhouse[mtest$env_HabitationType=='Townhouse / Semi-detached house']<-1

mtest$inf_AgeSolidsIntroduced[mtest$inf_AgeSolidsIntroduced==6]<-'6'
mtest$inf_AgeSolidsIntroduced[mtest$inf_AgeSolidsIntroduced<6]<-'5'
mtest$inf_AgeFormulaIntroduced[mtest$inf_AgeFormulaIntroduced==999]<-NA
mtest$trajec[mtest$trajec=='1']<-1
mtest$trajec[mtest$trajec=='2a']<-2
mtest$trajec[mtest$trajec=='2b']<-4
mtest$trajec[mtest$trajec=='3']<-5
mtest$trajec[mtest$trajec=='2c']<-3

mtest$T1<-0
mtest$T1[mtest$trajec==1]<-1
mtest$T1[is.na(mtest$trajec)]<-NA
chisqresT1 <- data.frame(var=names(mtest)[c(-19,-26,-27,-28,-29,-30)],T1=rep(NA,24))
for(i in chisqresT1$var){
  j=1
    if (i!='inf_AgeFormulaIntroduced') {
      chisqresT1[chisqresT1$var==i,paste("T",j,sep="")] <- chisq.test(table(mtest[,i],mtest$T1)[,c(1,2)])$p.value
    }
    if (i=='inf_AgeFormulaIntroduced') {
      chisqresT1[chisqresT1$var==i,paste("T",j,sep="")] <- chisq.test(tapply(mtest$inf_AgeFormulaIntroduced,mtest$T1,median,na.rm=T)[c(j,2)])$p.value
    }
  }

rownames(chisqresT1)<-chisqresT1$var
chisqresT1<-chisqresT1[-which(rownames(chisqresT1) %in% c('m_AgeDelivery','probiotic3m','inf_Prob_Lactobacillus3m','inf_Prob_Bifidobacterium3m')),]
#chisqresT1$T1[chisqresT1$T1>0.05]<-NA

mtest$T2<-0
mtest$T2[mtest$trajec==2]<-1
mtest$T2[is.na(mtest$trajec)]<-NA
chisqresT2 <- data.frame(var=names(mtest)[c(-19,-26,-27,-28,-29,-30)],T2=rep(NA,24))
for(i in chisqresT2$var){
  j=1
    if (i!='inf_AgeFormulaIntroduced') {
      chisqresT2[chisqresT2$var==i,2] <- chisq.test(table(mtest[,i],mtest$T2)[,c(1,2)])$p.value
    }
    if (i=='inf_AgeFormulaIntroduced') {
      chisqresT2[chisqresT2$var==i,2] <- chisq.test(tapply(mtest$inf_AgeFormulaIntroduced,mtest$T2,median,na.rm=T)[c(j,2)])$p.value
    }
  }

rownames(chisqresT2)<-chisqresT2$var
chisqresT2<-chisqresT2[-which(rownames(chisqresT2) %in% c('m_AgeDelivery','probiotic3m','inf_Prob_Lactobacillus3m','inf_Prob_Bifidobacterium3m')),]
#chisqresT2$T2[chisqresT2$T2>0.05]<-NA

mtest$T3<-0
mtest$T3[mtest$trajec==3]<-1
mtest$T3[is.na(mtest$trajec)]<-NA
chisqresT3 <- data.frame(var=names(mtest)[c(-19,-26,-27,-28,-29,-30)],T3=rep(NA,24))
for(i in chisqresT3$var){
  j=1
    if (i!='inf_AgeFormulaIntroduced') {
      chisqresT3[chisqresT3$var==i,2] <- chisq.test(table(mtest[,i],mtest$T3)[,c(1,2)])$p.value
    }
    if (i=='inf_AgeFormulaIntroduced') {
      chisqresT3[chisqresT3$var==i,2] <- chisq.test(tapply(mtest$inf_AgeFormulaIntroduced,mtest$T3,median,na.rm=T)[c(j,2)])$p.value
    }
  }

rownames(chisqresT3)<-chisqresT3$var
chisqresT3<-chisqresT3[-which(rownames(chisqresT3) %in% c('m_AgeDelivery','probiotic3m','inf_Prob_Lactobacillus3m','inf_Prob_Bifidobacterium3m')),]
#chisqresT3$T3[chisqresT3$T3>0.05]<-NA


mtest$T4<-0
mtest$T4[mtest$trajec==4]<-1
mtest$T4[is.na(mtest$trajec)]<-NA
chisqresT4 <- data.frame(var=names(mtest)[c(-19,-26,-27,-28,-29,-30)],T4=rep(NA,24))
for(i in chisqresT4$var){
  j=1
    if (i!='inf_AgeFormulaIntroduced') {
      chisqresT4[chisqresT4$var==i,2] <- chisq.test(table(mtest[,i],mtest$T4)[,c(1,2)])$p.value
    }
    if (i=='inf_AgeFormulaIntroduced') {
      chisqresT4[chisqresT4$var==i,2] <- chisq.test(tapply(mtest$inf_AgeFormulaIntroduced,mtest$T4,median,na.rm=T)[c(j,2)])$p.value
    }
  }

rownames(chisqresT4)<-chisqresT4$var
chisqresT4<-chisqresT4[-which(rownames(chisqresT4) %in% c('m_AgeDelivery','probiotic3m','inf_Prob_Lactobacillus3m','inf_Prob_Bifidobacterium3m')),]
#chisqresT4$T4[chisqresT4$T4>0.05]<-NA

mtest$T5<-0
mtest$T5[mtest$trajec==5]<-1
mtest$T5[is.na(mtest$trajec)]<-NA
chisqresT5_ <- data.frame(var=names(mtest)[c(-19,-26,-27,-28,-29,-30)],T5=rep(NA,24))
for(i in chisqresT5_$var){
  j=1
    if (i!='inf_AgeFormulaIntroduced') {
      chisqresT5_[chisqresT5_$var==i,2] <- chisq.test(table(mtest[,i],mtest$T5)[,c(1,2)])$p.value
    }
    if (i=='inf_AgeFormulaIntroduced') {
      chisqresT5_[chisqresT5_$var==i,2] <- chisq.test(tapply(mtest$inf_AgeFormulaIntroduced,mtest$T5,median,na.rm=T)[c(j,2)])$p.value
    }
  }

rownames(chisqresT5_)<-chisqresT5_$var
chisqresT5_<-chisqresT5_[-which(rownames(chisqresT5_) %in% c('m_AgeDelivery','probiotic3m','inf_Prob_Lactobacillus3m','inf_Prob_Bifidobacterium3m')),]
#chisqresT5_$T5[chisqresT5_$T5>0.05]<-NA




chisqres <- data.frame(var=names(mtest)[c(-19,-26,-27,-28,-29,-30)],T2=rep(NA,24),T3=rep(NA,24),T4=rep(NA,24),T5=rep(NA,24))
for(i in chisqres$var){
  for(j in c(2:5)){
    if (i!='inf_AgeFormulaIntroduced') {
      chisqres[chisqres$var==i,paste("T",j,sep="")] <- chisq.test(table(mtest[,i],mtest$trajec)[,c(1,j)])$p.value
    }
    if (i=='inf_AgeFormulaIntroduced') {
      chisqres[chisqres$var==i,paste("T",j,sep="")] <- chisq.test(tapply(mtest$inf_AgeFormulaIntroduced,mtest$trajec,median,na.rm=T)[c(j,5)])$p.value
    }
  }
}
rownames(chisqres)<-chisqres$var
chisqres<-chisqres[ ,-which(names(chisqres) %in% c('var'))]
chisqres$T1<-NA
chisqres<-chisqres[,c(5,1,2,3,4)]
chisqres<-chisqres[-which(rownames(chisqres) %in% c('probiotic3m','inf_Prob_Lactobacillus3m','inf_Prob_Bifidobacterium3m')),]

chisqresT5 <- data.frame(var=names(mtest)[c(-19,-26,-27,-28,-29,-30)],T2=rep(NA,24),T3=rep(NA,24),T4=rep(NA,24),T1=rep(NA,24))
chisqresT5['inf_AgeSolidsIntroduced',var]<-NA
for(i in chisqresT5$var){
  for(j in c(1:4)){
    if (i!='inf_AgeFormulaIntroduced') {
      chisqresT5[chisqresT5$var==i,paste("T",j,sep="")] <- chisq.test(table(mtest[,i],mtest$trajec)[,c(j,5)])$p.value
    }
    if (i=='inf_AgeFormulaIntroduced') {
      chisqresT5[chisqresT5$var==i,paste("T",j,sep="")] <- chisq.test(tapply(mtest$inf_AgeFormulaIntroduced,mtest$trajec,median,na.rm=T)[c(j,5)])$p.value
    }
  }
}
rownames(chisqresT5)<-chisqresT5$var
chisqresT5<-chisqresT5[ ,-which(names(chisqresT5) %in% c('var'))]
chisqresT5<-chisqresT5[,c(4,1:3)]
chisqresT5$T5<-NA

#comgining all the different test results
chi<-cbind(chisqresT1,chisqresT2$T2,chisqresT3$T3,chisqresT4$T4,chisqresT5_$T5)
chi<-chi[ ,-which(names(chi) %in% c('var'))]
names(chi)<-c('T1','T2','T3','T4','T5')
for (ii in rownames(chi)) {
  for (jj in 1:5) {
    chi[ii,jj]<-as.numeric(min(chi[ii,jj],chisqres[ii,jj],chisqresT5[ii,jj],na.rm = T))
  }
}






perper<-mtest

perper$trajec[mtest$trajec=='1']<-'T1'
perper$trajec[mtest$trajec=='2']<-'T2'
perper$trajec[mtest$trajec=='3']<-'T3'
perper$trajec[mtest$trajec=='4']<-'T4'
perper$trajec[mtest$trajec=='5']<-'T5'

perper$Csection<-0
perper$Csection[perper$birthmode_RecievingIAP=='C-section']<-1
perper$VABX<-0
perper$VABX[perper$birthmode_RecievingIAP=='Vaginal_AB']<-1
perper$Vag<-0
perper$Vag[perper$birthmode_RecievingIAP=='Vaginal_NoAB']<-1

perper$sibs[perper$sibs=='yes']<-1
perper$sibs[perper$sibs=='no']<-0
perper$Apartment<-0
perper$Apartment[perper$env_HabitationType=='Appartment']<-1
perper$Apartment<-0
perper$Apartment[perper$env_HabitationType=='Appartment']<-1
perper$House<-0
perper$House[perper$env_HabitationType=='House']<-1
perper$Townhouse<-0
perper$Townhouse[perper$env_HabitationType=='Townhouse / Semi-detached house']<-1


perper$inf_AgeFormulaIntroduced[perper$inf_AgeFormulaIntroduced==999]<-NA
perper$m_probio[perper$m_probio>1]<-2


pp<-data.frame(matrix(nrow = 5, ncol = 21))
for (ii in c(2,5,6,8:10,12:18,20:25)) {
  pp[,ii]<- prop.table(table(perper$trajec,perper[,ii]),margin = 1)[,2]
  names(pp)[ii]<-names(perper)[ii]
  if (ii==6) {
    pp[,ii]<- prop.table(table(perper$trajec,perper[,ii]),margin = 1)[,1]
    names(pp)[ii]<-names(perper)[ii]
  }
}
rownames(pp)<-c('T1','T2','T3','T4','T5')
pp$X7<-tapply(perper$inf_AgeFormulaIntroduced,perper$trajec,median,na.rm=T)
names(pp)[names(pp)=='X7']<-'Median age of formula introduced'
pp<-t(pp)
pp<-pp[-which(rownames(pp) %in% c('probiotic3m','inf_Prob_Lactobacillus3m','inf_Prob_Bifidobacterium3m')),]

#P value for chisqres and chisqrest5
#chisqres
# pp_pvalue<-pp
# pp_pvalue[is.numeric(pp_pvalue)]<-NA
# 
# for (ii in rownames(chisqres)) {
#   for (jj in c('T2','T3','T4','T5'))
#   if (ii %in% rownames(pp_pvalue)) {
#      pp_pvalue[ii,jj]<-chisqres[ii,jj]
#   }
# }
# 
# #chisqres 5
# pp_pvalueT5<-pp
# pp_pvalueT5[!is.na(pp_pvalueT5)]<-NA
# 
# for (ii in rownames(chisqresT5)) {
#   for (jj in c('T1','T2','T3','T4'))
#   if (ii %in% rownames((pp_pvalueT5))) {
#      pp_pvalueT5[ii,jj]<-chisqresT5[ii,jj]
#   }
# }

pp_pvalue<-chi
#Renaming------

rownames(pp)[rownames(pp) == 'sibs']             <- 'Siblings'
rownames(pp)[rownames(pp) == 'env_FurOrFeathers']            <- 'Household has pets'
rownames(pp)[rownames(pp) == 'inf_AgeSolidsIntroduced']      <- 'Solid foods introduced before 6 mo'
rownames(pp)[rownames(pp) == 'inf_AgeFormulaIntroduced']          <- 'Formula was introduced at 3 wk'
rownames(pp)[rownames(pp) == 'm_probio']                    <- 'Mother took probiotic before birth'
rownames(pp)[rownames(pp) == 'exclusiveBF_6months']        <- 'No formula before 6 mo'
rownames(pp)[rownames(pp) == 'exclusiveBF_12months']            <- 'No formula before 12 mo'
rownames(pp)[rownames(pp) == 'antibiotic3m']            <- 'Antiobiotic recieved before 3 mo'
rownames(pp)[rownames(pp) == 'probiotic3m'] <- 'Probiotic (any) recieved before 3 mo'
rownames(pp)[rownames(pp) == 'probiotic6m'] <- 'Probiotic (any) recieved before 6 mo'
rownames(pp)[rownames(pp) == 'inf_Prob_Lactobacillus3m'] <- 'Probiotic (Lactobacillus) recieved before 3 mo'
rownames(pp)[rownames(pp) == 'inf_Prob_Lactobacillus6m'] <- 'Probiotic (Lactobacillus) recieved before 6 mo'
rownames(pp)[rownames(pp) == 'inf_Prob_Bifidobacterium3m'] <- 'Probiotic (Bifidobacterium) recieved before 3 mo'
rownames(pp)[rownames(pp) == 'inf_Prob_Bifidobacterium6m'] <- 'Probiotic (Bifidobacterium) recieved before 6 mo'
rownames(pp)[rownames(pp) == 'Csection']           <- 'C-section delivery'
rownames(pp)[rownames(pp) == 'VABX']            <- 'Vaginal delivery intrapartum antibiotics (VABX)'
rownames(pp)[rownames(pp) == 'Vag']          <- 'Vaginal delivery'
rownames(pp)[rownames(pp) == 'Apartment']          <- 'Residence: apartment'
rownames(pp)[rownames(pp) == 'House']          <- 'Residence: detached house'
rownames(pp)[rownames(pp) == 'Townhouse']          <- 'Residence: townhouse/semi-detached house'

rownames(pp_pvalue)[rownames(pp_pvalue) == 'sibs']             <- 'Siblings'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'env_FurOrFeathers']            <- 'Household has pets'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'inf_AgeSolidsIntroduced']      <- 'Solid foods introduced before 6 mo'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'inf_AgeFormulaIntroduced']          <- 'Formula was introduced at 3 wk'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'm_probio']                    <- 'Mother took probiotic before birth'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'exclusiveBF_6months']        <- 'No formula before 6 mo'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'exclusiveBF_12months']            <- 'No formula before 12 mo'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'antibiotic3m']            <- 'Antiobiotic recieved before 3 mo'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'probiotic3m'] <- 'Probiotic (any) recieved before 3 mo'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'probiotic6m'] <- 'Probiotic (any) recieved before 6 mo'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'inf_Prob_Lactobacillus3m'] <- 'Probiotic (Lactobacillus) recieved before 3 mo'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'inf_Prob_Lactobacillus6m'] <- 'Probiotic (Lactobacillus) recieved before 6 mo'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'inf_Prob_Bifidobacterium3m'] <- 'Probiotic (Bifidobacterium) recieved before 3 mo'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'inf_Prob_Bifidobacterium6m'] <- 'Probiotic (Bifidobacterium) recieved before 6 mo'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'Csection']           <- 'C-section delivery'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'VABX']            <- 'Vaginal delivery intrapartum antibiotics (VABX)'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'Vag']          <- 'Vaginal delivery'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'Apartment']          <- 'Residence: apartment'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'House']          <- 'Residence: detached house'
rownames(pp_pvalue)[rownames(pp_pvalue) == 'Townhouse']          <- 'Residence: townhouse/semi-detached house'


#PLotting----

# pp_pvalueT50<-pp_pvalueT5
# l0 <- pp_pvalueT5
# #labs<-l0
# pp_pvalueT5[l0 <= 0.1] <- '.'
# pp_pvalueT5[l0 < 0.05] <- '+'
# pp_pvalueT5[l0 < 0.01] <- '++'
# pp_pvalueT5[l0 < 0.001] <- '+++'
# pp_pvalueT5[l0 >= 0.05 | is.na(l0)] <- ''
# pp_pvalueT5 <- as.matrix(pp_pvalueT5)
# labsT5<-pp_pvalueT5
# labsT5<-labsT5[-which(rownames(labsT5) %in% c('X1','X3','X4','X11','X19')),]

#combining multiple sign values
# labsC<-labs
# for (ii in 1:nrow(labs)) {
#   for (jj in 1:ncol(labsT5)) {
#     labsC[ii,jj]<-paste(labs[ii,jj],labsT5[ii,jj],sep='\n')
#   }
# }

pp<-pp[-which(rownames(pp) %in% c('X1','X3','X4','X11','X19',
                                  'birthmode_RecievingIAP',
                                   'familly_NbSiblingsFt',
                                   'env_HabitationType',
                                   'Solid foods introduced before 6 mo',
                                   "Mother took probiotic before birth" ,
                                   "Probiotic (any) recieved before 6 mo",
                                   "Probiotic (Bifidobacterium) recieved before 6 mo",
                                   "Vaginal delivery",
                                   "Residence: detached house",
                                   "Residence: townhouse/semi-detached house",
                                  "Median age of formula introduced"  
                                  )),]
pp_pvalue<-pp_pvalue[-which(rownames(pp_pvalue) %in% c('birthmode_RecievingIAP',
                                             'familly_NbSiblingsFt',
                                             'env_HabitationType',
                                             'Solid foods introduced before 6 mo',
                                             "Mother took probiotic before birth" ,
                                             "Probiotic (any) recieved before 6 mo",
                                             "Probiotic (Bifidobacterium) recieved before 6 mo",
                                             "Vaginal delivery",
                                             "Residence: detached house",
                                             "Residence: townhouse/semi-detached house" ,
                                             "Median age of formula introduced",  
                                             "Formula was introduced at 3 wk"
                                             )),]


pp_pvalue0<-pp_pvalue
l0 <- pp_pvalue
#labs<-l0
pp_pvalue[l0 <= 0.1] <- '.'
pp_pvalue[l0 < 0.05] <- '*'
pp_pvalue[l0 < 0.01] <- '**'
pp_pvalue[l0 < 0.001] <- '***'
pp_pvalue[l0 >= 0.05 | is.na(l0)] <- ''
pp_pvalue <- as.matrix(pp_pvalue)
labs<-pp_pvalue
#labs<-labs[-which(rownames(labs) %in% c('X1','X3','X4','X11','X19')),]


# write.table(l0,'trajectory_personal_chisq_pvalues.txt')
# write.table(pp,'trajectory_personal_chisq_estimates.txt')
# pp<-read.delim('trajectory_personal_chisq_estimates.txt',sep = ' ')
# l0<-pp_pvalue<-read.delim('trajectory_personal_chisq_pvalues.txt',sep = ' ');#pp_pvalue<-l0
# 

pdf('TrajectoryPersonal_percentages_T_switched_pvalueT5_3.pdf',width=18,height = 16)
  par(oma=c(1,1,1,50))
  heatmap.2(as.matrix(pp),scale='row',dendrogram='none',
            Colv=FALSE,trace='none',distfun=peardist,
            col=rev(brewer.pal(n=9,'RdBu')),density.info = 'none',
            cexRow=4,cexCol=4,#srtRow = 315,
            lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2, 16), lwid=c(1, 2, 1),
            cellnote=labs, notecol="black", notecex=8,
            key.xtickfun = function() {
               breaks = pretty(parent.frame()$breaks)
               breaks = breaks[c(1:length(breaks))]
               list(at = parent.frame()$scale01(breaks),labels = breaks)
            }
  )
dev.off()   


