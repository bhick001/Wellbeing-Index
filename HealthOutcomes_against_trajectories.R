library(gplots)
peardist <- function(x, distance) {
  as.dist(1 - cor(t(log(x + min(
    x[x > 0]
  ) / 10))))
}

infection_0to3m <-
  read.delim('infectiondata/0to3Mo_aggregated.csv', sep = ',')
colnames(infection_0to3m)[2:13] <-
  paste(colnames(infection_0to3m)[2:13], "0to3m", sep = "_")
infection_0to6m <-
  read.delim('infectiondata/0to6Mo_aggregated.csv', sep = ',')
colnames(infection_0to6m)[2:13] <-
  paste(colnames(infection_0to6m)[2:13], "0to6m", sep = "_")
infection_6to12m <-
  read.delim('infectiondata/6to12Mo_aggregated.csv', sep = ',')
colnames(infection_6to12m)[2:13] <-
  paste(colnames(infection_6to12m)[2:13], "6to12m", sep = "_")
infection_12to24m <-
  read.delim('infectiondata/12to24Mo_aggregated.csv', sep = ',')
colnames(infection_12to24m)[2:13] <-
  paste(colnames(infection_12to24m)[2:13], "12to24m", sep = "_")

infection <-
  merge(infection_0to3m,
        infection_0to6m,
        by = 'Family_ID',
        all = T)
infection <- merge(infection, infection_6to12m, by = 'Family_ID', all =
                     T)
infection <- merge(infection, infection_12to24m, by = 'Family_ID', all =
                     T)
colnames(infection)[1] <- 'id'
id <- infection$id
infection <- infection[, -c(1)]
infection <- Filter(function(x)
  sd(x, na.rm = T) != 0, infection)
infection[infection > 1] <- 1
infection$id <- id


GI3w <-
  read.delim('Average_CryingGI_19.11.23.csv', sep = ',')


persondata <- read.delim('helmi_persondata5.txt') #extra cases in T1

m <- read.delim('meta.txt')
m <- m[m$use == 1, ]
m$healthyref <-
  rowSums(m[, c(grep(names(m), pattern = "_0to2yr|_2yr$|_5yr", value = T)[c(7, 8, 10, 12, 14, 15, 19, 26)])], na.rm =
            T)
m$healthyref[abs(m$BMI.for.age) > 2 &
               !is.na(m$BMI.for.age) |
               abs(m$height.for.age) > 2 &
               !is.na(m$height.for.age)] <-
  m$healthyref[abs(m$BMI.for.age) > 2 &
                 !is.na(m$BMI.for.age) |
                 abs(m$height.for.age) > 2 & !is.na(m$height.for.age)] + 1
m$healthyref[abs(m$BMI.for.age_5yr) > 2 &
               !is.na(m$BMI.for.age_5yr) |
               abs(m$height.for.age_5yr) > 2 &
               !is.na(m$height.for.age_5yr)] <-
  m$healthyref[abs(m$BMI.for.age_5yr) > 2 &
                 !is.na(m$BMI.for.age_5yr) |
                 abs(m$height.for.age_5yr) > 2 & !is.na(m$height.for.age_5yr)] + 1

mnew <-
  read.delim(
    "meta.txt"
  )
muse <- mnew[mnew$use == 1, ]
mdrop <- read.delim('Data_for_trajectory_tests.txt', sep = '\t')
smoking <- muse[, c('id', 'm_SmokingPrior', 'p_SmokingDuring')]
smoking <-
  unique(smoking[, c('id', 'm_SmokingPrior', 'p_SmokingDuring')])

house <- muse[, c('id', 'env_HabitationType')]
house <- unique(house[, c('id', 'env_HabitationType')])

persondata$m_SmokingPrior <- NA
persondata$p_SmokingDuring <- NA
persondata$env_HabitationType <- NA
for (fid in unique(m$id)) {
  persondata$m_SmokingPrior[persondata$id == fid] <-
    smoking$m_SmokingPrior[smoking$id == fid]
  persondata$p_SmokingDuring[persondata$id == fid] <-
    smoking$p_SmokingDuring[smoking$id == fid]
  persondata$HealthReference[persondata$id == fid] <-
    unique(m$healthyref[m$Family_ID == fid])
  persondata$env_HabitationType[persondata$id == fid] <-
    house$env_HabitationType[house$id == fid]
}


persondata <- merge(persondata, infection, by = 'id', all.x = T)

#-------

response <-
  c("excessive_crying",
    grep(names(muse), pattern = "_0to2yr|_2yr$|_5yr", value = T)[-c(3, 4, 8, 27:31)])
response_infect <-
  c(grep(names(persondata), pattern = "0to3m|0to6m|0to6m|12to24m", value =
           T)[-c(3, 4, 8, 27:31)])
response <- c(response, response_infect)

#Against All -----
Traj3_clustScores_healthtests <- data.frame(
  var = response,cl1 = NA,cl2a = NA,cl2b = NA,cl2c = NA,cl3 = NA
  )
Traj3_clustScores_healthtestsE <- Traj3_clustScores_healthtests
for (j in na.omit(unique(as.character(persondata$trajec)))) {
  persondata$traj <-as.character(persondata$trajec)
  persondata$traj[persondata$trajec != j] <- "0"
  for (i in unique(Traj3_clustScores_healthtests$var)) {
    reference<-persondata[persondata$HealthReference==0,]
    gr <- persondata$HealthReference
    gr[gr>0&persondata[,i]==1]<-1
    gr[gr>0&persondata[,i]==0|gr>0&is.na(persondata[,i])]<-NA
    persondata$gr<-gr
    mod <-glm(as.formula(paste('gr', "~m_BMI+m_HadAllergyDisease+f_HadAllergyDisease+m_edu+p_edu+
          inf_Sex+siblings+env_FurOrFeathers+household_size+
          env_HabitationType+
          p_SmokingDuring + m_SmokingPrior+
          m_weightGain+m_DiagnosisPregnancyType_GestationalDiabetes
         as.factor(traj)"
        ) ), data = persondata, family = "binomial")
    if ("as.factor(traj)1" %in% rownames(summary(mod)$coef))
      Traj3_clustScores_healthtests[Traj3_clustScores_healthtests$var == i, "cl1"] <-
        summary(mod)$coef["as.factor(traj)1", 4]
    if ("as.factor(traj)2a" %in% rownames(summary(mod)$coef))
      Traj3_clustScores_healthtests[Traj3_clustScores_healthtests$var == i, "cl2a"] <-
        summary(mod)$coef["as.factor(traj)2a", 4]
    if ("as.factor(traj)2b" %in% rownames(summary(mod)$coef))
      Traj3_clustScores_healthtests[Traj3_clustScores_healthtests$var == i, "cl2b"] <-
        summary(mod)$coef["as.factor(traj)2b", 4]
    if ("as.factor(traj)2c" %in% rownames(summary(mod)$coef))
      Traj3_clustScores_healthtests[Traj3_clustScores_healthtests$var == i, "cl2c"] <-
        summary(mod)$coef["as.factor(traj)2c", 4]
    if ("as.factor(traj)3" %in% rownames(summary(mod)$coef))
      Traj3_clustScores_healthtests[Traj3_clustScores_healthtests$var == i, "cl3"] <-
        summary(mod)$coef["as.factor(traj)3", 4]
    
    if ("as.factor(traj)1" %in% rownames(summary(mod)$coef))
      if (summary(mod)$coef["as.factor(traj)1", 4] < 0.05)
        Traj3_clustScores_healthtestsE[Traj3_clustScores_healthtestsE$var == i, "cl1"] <-
          summary(mod)$coef["as.factor(traj)1", 1]
      if ("as.factor(traj)2a" %in% rownames(summary(mod)$coef))
        if (summary(mod)$coef["as.factor(traj)2a", 4] < 0.05)
          Traj3_clustScores_healthtestsE[Traj3_clustScores_healthtestsE$var == i, "cl2a"] <-
            summary(mod)$coef["as.factor(traj)2a", 1]
        if ("as.factor(traj)2b" %in% rownames(summary(mod)$coef))
          if (summary(mod)$coef["as.factor(traj)2b", 4] < 0.05)
            Traj3_clustScores_healthtestsE[Traj3_clustScores_healthtestsE$var == i, "cl2b"] <-
              summary(mod)$coef["as.factor(traj)2b", 1]
          if ("as.factor(traj)2c" %in% rownames(summary(mod)$coef))
            if (summary(mod)$coef["as.factor(traj)2c", 4] < 0.05)
              Traj3_clustScores_healthtestsE[Traj3_clustScores_healthtestsE$var == i, "cl2c"] <-
                summary(mod)$coef["as.factor(traj)2c", 1]
            if ("as.factor(traj)3" %in% rownames(summary(mod)$coef))
              if (summary(mod)$coef["as.factor(traj)3", 4] < 0.05)
                Traj3_clustScores_healthtestsE[Traj3_clustScores_healthtestsE$var == i, "cl3"] <-
                  summary(mod)$coef["as.factor(traj)3", 1]
  }
}
#Traj3_clustScores_healthtestsE
rownames(Traj3_clustScores_healthtestsE) <-Traj3_clustScores_healthtestsE$var
Traj3_clustScores_healthtestsE<-Traj3_clustScores_healthtestsE[,c(1:3,5,4,6)]
colnames(Traj3_clustScores_healthtestsE) <-c('var', 'T1', 'T2', 'T3', 'T4', 'T5')


Traj3_clustScores_healthtestsE_pvalue <-Traj3_clustScores_healthtests
rownames(Traj3_clustScores_healthtestsE_pvalue) <-Traj3_clustScores_healthtestsE_pvalue$var
Traj3_clustScores_healthtestsE_pvalue<-Traj3_clustScores_healthtestsE_pvalue[,c(1:3,5,4,6)]
colnames(Traj3_clustScores_healthtestsE_pvalue) <-c('var', 'T1', 'T2', 'T3', 'T4', 'T5')
#Traj3_clustScores_healthtestsE_pvalue[Traj3_clustScores_healthtestsE_pvalue>0.05]<-NA

df_list <- list(Traj3_clustScores_healthtestsE)
dfm <- Reduce(function(x, y) merge(x, y, all = TRUE, by = 'var'), df_list)
rownames(dfm) <- dfm$var
dfm <- dfm[,-which(names(dfm) %in% c("var", 'T1 vs T1', 'T1 vs T3', "T3 vs T3"))]

df_list <- list(Traj3_clustScores_healthtestsE_pvalue)
labs <- Reduce(function(x, y) merge(x, y, all = TRUE, by = 'var'), df_list)
rownames(labs) <- labs$var
labs <-labs[,-which(names(labs) %in% c("var", 'T1 vs T1', 'T1 vs T3', "T3 vs T3"))]
glm_pvalue<-labs
labs <- labs[rowSums(is.na(dfm)) < ncol(dfm), ]
dfm <- dfm[rowSums(is.na(dfm)) < ncol(dfm), ]


l0 <- labs
labs[l0 <= 0.1] <- '.'
labs[l0 < 0.05] <- '*'
labs[l0 < 0.01] <- '**'
labs[l0 < 0.001] <- '***'
labs[l0 >= 0.05 | is.na(l0)] <- ''
labs <- as.matrix(labs)


#-------------------------
#Chi Squared testing-----
perper<-persondata[,c(2:30,42)]

#Chi square test
perper$T1<-0
perper$T1[perper$trajec=='T1']<-1
perper$T1[is.na(perper$trajec)]<-NA
chisqresT1 <- data.frame(var=names(perper)[c(-30:-50)],T1=rep(NA,29))
for(i in chisqresT1$var){
  gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    
    perper$gr
    
    perper$gr[perper$gr==1]<-i

    j=1
    if (dim(t(table(perper$gr,perper$T1)))[2]>1) {
      chisqresT1[chisqresT1$var==i,paste("T",j,sep="")] <- chisq.test(table(perper$gr,perper$T1)[,c(1,2)])$p.value
    }
  }

rownames(chisqresT1)<-chisqresT1$var

perper$T2<-0
perper$T2[perper$trajec=='T2']<-1
perper$T2[is.na(perper$trajec)]<-NA
chisqresT2 <- data.frame(var=names(perper)[c(-30:-100)],T2=rep(NA,29))
for(i in chisqresT2$var){
  j=1
  gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    if (dim(t(table(perper$gr,perper$T2)))[2]>1) {
      chisqresT2[chisqresT2$var==i,2] <- chisq.test(table(perper$gr,perper$T2)[,c(1,2)])$p.value
    }
    
  }

rownames(chisqresT2)<-chisqresT2$var

perper$T3<-0
perper$T3[perper$trajec=='T3']<-1
perper$T3[is.na(perper$trajec)]<-NA
chisqresT3 <- data.frame(var=names(perper)[c(-30:-100)],T3=rep(NA,29))
for(i in chisqresT3$var){
  j=1
   gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    if (dim(t(table(perper$gr,perper$T3)))[2]>1) {
      chisqresT3[chisqresT3$var==i,2] <- chisq.test(table(perper$gr,perper$T3)[,c(1,2)])$p.value
    }
      }
rownames(chisqresT3)<-chisqresT3$var

perper$T4<-0
perper$T4[perper$trajec=='T4']<-1
perper$T4[is.na(perper$trajec)]<-NA
chisqresT4 <- data.frame(var=names(perper)[c(-30:-100)],T4=rep(NA,29))
for(i in chisqresT4$var){
  gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    j=1
    if (dim(t(table(perper$gr,perper$T4)))[2]>1) {
      chisqresT4[chisqresT4$var==i,2] <- chisq.test(table(perper$gr,perper$T4)[,c(1,2)])$p.value
    }
    
  }
rownames(chisqresT4)<-chisqresT4$var

perper$T5<-0
perper$T5[perper$trajec=='T5']<-1
perper$T5[is.na(perper$trajec)]<-NA
chisqresT5_ <- data.frame(var=names(perper)[c(-30:-100)],T5=rep(NA,29))
for(i in chisqresT5_$var){
  j=1
  gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    if (dim(t(table(perper$gr,perper$T5)))[2]>1) {
      chisqresT5_[chisqresT5_$var==i,2] <- chisq.test(table(perper$gr,perper$T5)[,c(1,2)])$p.value
    }
    
  }
rownames(chisqresT5_)<-chisqresT5_$var

chisqres <- data.frame(var=names(perper)[c(-30:-100)],T2=rep(NA,29),T3=rep(NA,29),T4=rep(NA,29),T5=rep(NA,29))
for(i in chisqres$var){
  for(j in c(2:5)){
    gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    if (dim(t(table(perper$gr,perper[,paste("T",j,sep="")])))[2]>1) {
      chisqres[chisqres$var==i,paste("T",j,sep="")] <- chisq.test(table(perper$gr,perper$trajec)[,c(1,j)])$p.value
    }
    
  }
}
rownames(chisqres)<-chisqres$var
chisqres<-chisqres[ ,-which(names(chisqres) %in% c('var'))]
chisqres$T1<-NA
chisqres<-chisqres[,c(5,1,2,3,4)]


chisqresT5 <- data.frame(var=names(perper)[c(-30:-100)],T2=rep(NA,29),T3=rep(NA,29),T4=rep(NA,29),T1=rep(NA,29))
for(i in chisqresT5$var){
  for(j in c(1:4)){
    gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    if (i!='inf_AgeFormulaIntroduced') {
      chisqresT5[chisqresT5$var==i,paste("T",j,sep="")] <- chisq.test(table(perper$gr,perper$trajec)[,c(j,5)])$p.value
    }
    
  }
}
rownames(chisqresT5)<-chisqresT5$var
chisqresT5$T5<-NA
chisqresT5<-chisqresT5[ ,-which(names(chisqresT5) %in% c('var'))]
chisqresT5<-chisqresT5[,c(4,1:3,5)]


#comgining all the different test results
chi<-cbind(chisqresT1,chisqresT2$T2,chisqresT3$T3,chisqresT4$T4,chisqresT5_$T5)
chi<-chi[ ,-which(names(chi) %in% c('var'))]
names(chi)<-c('T1','T2','T3','T4','T5')
for (ii in rownames(chi)) {
  for (jj in 1:5) {
    chi[ii,jj]<-min(chi[ii,jj],chisqres[ii,jj],chisqresT5[ii,jj],glm_pvalue[ii,jj],na.rm = T)
  }
}
pp_pvalue<-chi


#percentage tables health reference
perper<-perper0
pp<-data.frame(matrix(nrow = 5, ncol = 42))
for (ii in 1:42) {
  gr <- persondata$HealthReference
  gr[gr>0 & perper[,ii]==1]<-1
  gr[gr>0 & perper[,ii]==0 | gr>0 &is.na(perper[,ii])]<-NA
  perper$gr<-gr
  if (length(unique(perper[!is.na(perper[,ii]),ii]))>1) {
      pp[,ii]<- prop.table(table(perper$trajec,perper$gr),margin = 1)[,2]
  }
  names(pp)[ii]<-names(perper)[ii]
}
rownames(pp)<-c('T1','T2','T3','T4','T5')
pp<-t(pp)
pp_infections<-pp


#-----------------------------------------
#-----------------------------------------
#-----------------------------------------
#INfection data----
#Chi square test

perper<-persondata[,c(53:94,42)]
perper$trajec[perper$trajec=='1']<-'T1'
perper$trajec[perper$trajec=='2a']<-'T2'
perper$trajec[perper$trajec=='2b']<-'T4'
perper$trajec[perper$trajec=='2c']<-'T3'
perper$trajec[perper$trajec=='3']<-'T5'
perper0<-perper

perper$T1<-0
perper$T1[perper$trajec=='T1']<-1
perper$T1[is.na(perper$trajec)]<-NA
chisqresT1 <- data.frame(var=names(perper)[c(-43:-100)],T1=rep(NA,42))
for(i in chisqresT1$var){
  gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    j=1
    if (dim(t(table(perper$gr,perper$T1)))[2]>1) {
      chisqresT1[chisqresT1$var==i,paste("T",j,sep="")] <- chisq.test(table(perper$gr,perper$T1)[,c(1,2)])$p.value
    }
    
  }
rownames(chisqresT1)<-chisqresT1$var

perper$T2<-0
perper$T2[perper$trajec=='T2']<-1
perper$T2[is.na(perper$trajec)]<-NA
chisqresT2 <- data.frame(var=names(perper)[c(-43:-100)],T2=rep(NA,42))
for(i in chisqresT2$var){
  j=1
  gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    if (dim(t(table(perper$gr,perper$T2)))[2]>1) {
      chisqresT2[chisqresT2$var==i,2] <- chisq.test(table(perper$gr,perper$T2)[,c(1,2)])$p.value
    }
    
  }

rownames(chisqresT2)<-chisqresT2$var

perper$T3<-0
perper$T3[perper$trajec=='T3']<-1
perper$T3[is.na(perper$trajec)]<-NA
chisqresT3 <- data.frame(var=names(perper)[c(-43:-100)],T3=rep(NA,42))
for(i in chisqresT3$var){
  j=1
   gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    if (dim(t(table(perper$gr,perper$T3)))[2]>1) {
      chisqresT3[chisqresT3$var==i,2] <- chisq.test(table(perper$gr,perper$T3)[,c(1,2)])$p.value
    }
    
  }
rownames(chisqresT3)<-chisqresT3$var


perper$T4<-0
perper$T4[perper$trajec=='T4']<-1
perper$T4[is.na(perper$trajec)]<-NA
chisqresT4 <- data.frame(var=names(perper)[c(-43:-100)],T4=rep(NA,42))
for(i in chisqresT4$var){
  gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    j=1
    if (dim(t(table(perper$gr,perper$T4)))[2]>1) {
      chisqresT4[chisqresT4$var==i,2] <- chisq.test(table(perper$gr,perper$T4)[,c(1,2)])$p.value
    }
    
  }
rownames(chisqresT4)<-chisqresT4$var

perper$T5<-0
perper$T5[perper$trajec=='T5']<-1
perper$T5[is.na(perper$trajec)]<-NA
chisqresT5_ <- data.frame(var=names(perper)[c(-43:-100)],T5=rep(NA,42))
for(i in chisqresT5_$var){
  j=1
  gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    if (dim(t(table(perper$gr,perper$T5)))[2]>1) {
      chisqresT5_[chisqresT5_$var==i,2] <- chisq.test(table(perper$gr,perper$T5)[,c(1,2)])$p.value
    }
    
  }

rownames(chisqresT5_)<-chisqresT5_$var


chisqres <- data.frame(var=names(perper)[c(-43:-100)],T2=rep(NA,42),T3=rep(NA,42),T4=rep(NA,42),T5=rep(NA,42))
for(i in chisqres$var){
  for(j in c(2:5)){
    gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    if (dim(t(table(perper$gr,perper[,paste("T",j,sep="")])))[2]>1) {
      chisqres[chisqres$var==i,paste("T",j,sep="")] <- chisq.test(table(perper$gr,perper$trajec)[,c(1,j)])$p.value
    }
    
  }
}
rownames(chisqres)<-chisqres$var
chisqres<-chisqres[ ,-which(names(chisqres) %in% c('var'))]
chisqres$T1<-NA
chisqres<-chisqres[,c(5,1,2,3,4)]


chisqresT5 <- data.frame(var=names(perper)[c(-43:-100)],T2=rep(NA,42),T3=rep(NA,42),T4=rep(NA,42),T1=rep(NA,42))
for(i in chisqresT5$var){
  for(j in c(1:4)){
    gr <- persondata$HealthReference
    gr[gr>0 & persondata[,i]==1]<-1
    gr[gr>0 & persondata[,i]==0 | gr>0 &is.na(persondata[,i])]<-NA
    perper$gr<-gr
    if (i!='inf_AgeFormulaIntroduced') {
      chisqresT5[chisqresT5$var==i,paste("T",j,sep="")] <- chisq.test(table(perper$gr,perper$trajec)[,c(j,5)])$p.value
    }
    
  }
}
rownames(chisqresT5)<-chisqresT5$var
chisqresT5$T5<-NA
chisqresT5<-chisqresT5[ ,-which(names(chisqresT5) %in% c('var'))]
chisqresT5<-chisqresT5[,c(4,1:3,5)]


#comgining all the different test results
chi<-cbind(chisqresT1,chisqresT2$T2,chisqresT3$T3,chisqresT4$T4,chisqresT5_$T5)
chi<-chi[ ,-which(names(chi) %in% c('var'))]
names(chi)<-c('T1','T2','T3','T4','T5')
for (ii in rownames(chi)) {
  for (jj in 1:5) {
    chi[ii,jj]<-min(chi[ii,jj],chisqres[ii,jj],chisqresT5[ii,jj],glm_pvalue[ii,jj],na.rm = T)
  }
}
pp_pvalue_infections<-chi



#percentage tables health reference
perper<-perper0
pp<-data.frame(matrix(nrow = 5, ncol = 42))
for (ii in 1:42) {
  gr <- persondata$HealthReference
  gr[gr>0 & perper[,ii]==1]<-1
  gr[gr>0 & perper[,ii]==0 | gr>0 &is.na(perper[,ii])]<-NA
  perper$gr<-gr
  if (length(unique(perper[!is.na(perper[,ii]),ii]))>1) {
      pp[,ii]<- prop.table(table(perper$trajec,perper$gr),margin = 1)[,2]
  }
  names(pp)[ii]<-names(perper)[ii]
}
rownames(pp)<-c('T1','T2','T3','T4','T5')
pp<-t(pp)
pp_infections<-pp



#------------
#Joining the two

pp_pvalue<-rbind(pp_pvalue,pp_pvalue_infections)
pp<-rbind(pp_health,pp_infections)

pp_pvalue<-pp_pvalue[-which(rownames(pp_pvalue) %in% c('MilkAllergyDiagnosis_0to2yr','AllergyDiagnosis_0to2yr',
                                'FoodAllergyDiagnosis_0to2yr','Obesity_5yr',
                                'BodyFat_5yr','Underweight_5yr')),]

pp<-pp[-which(rownames(pp) %in% c('MilkAllergyDiagnosis_0to2yr','AllergyDiagnosis_0to2yr',
                                'FoodAllergyDiagnosis_0to2yr','Obesity_5yr',
                                'BodyFat_5yr','Underweight_5yr')),]



#plotting----
pp_pvalue0<-pp_pvalue
l0 <- pp_pvalue
pp_pvalue[l0 <= 0.1] <- '.'
pp_pvalue[l0 < 0.05] <- '*'
pp_pvalue[l0 < 0.01] <- '**'
pp_pvalue[l0 < 0.001] <- '***'
pp_pvalue[l0 >= 0.05 | is.na(l0)] <- ''
pp_pvalue <- as.matrix(pp_pvalue)

pp2<-pp
pp_pvalue2<-pp_pvalue

pp2 <- pp2[!apply(pp_pvalue2, 1, function(x) all(x=="")),] 
pp_pvalue2 <- pp_pvalue2[!apply(pp_pvalue2, 1, function(x) all(x=="")),] 
 

pdf('TrajectoryHealthOutcomes.pdf',width=18,height = 16)
par(oma=c(5,1,1,55))
heatmap.2(as.matrix(pp3),dendrogram='none',scale='row',
          Colv=FALSE,Rowv=FALSE,trace='none',distfun=peardist,
          col=rev(brewer.pal(n=8,'RdBu')),density.info = 'none',
          cexRow=4,cexCol=4,#key=FALSE,
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2, 20), lwid=c(1, 3.5, 1),
          cellnote=pp_pvalue3, notecol="black", notecex=8,
          key.xtickfun = function() {
             breaks = pretty(parent.frame()$breaks)
             breaks = breaks[c(1:length(breaks))]
             list(at = parent.frame()$scale01(breaks),labels = breaks)
          }
)
dev.off()




