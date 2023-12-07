#wellbeing index based on health-associated bacteria----------------------------------------

library(gplots)
library(RColorBrewer)
library(ggplot2)

peardist <- function(x,distance){
  as.dist(1-cor(t(log(x + min(x[x>0])))))
}

m <- read.delim("meta.txt")
traj <- read.delim("trajectories.txt")
for(i in unique(traj$Family_ID)) m$RDP_Trajectory_final[m$Family_ID==i] <- unique(traj$TrajectoryRDP_final[traj$Family_ID==i])

m$healthyref <- rowSums(m[,c(grep(names(m),pattern="_0to2yr|_2yr$|_5yr",value=T)[c(7,8,10,12,14,15,19,26)])],na.rm=T)
m$healthyref[abs(m$BMI.for.age)>2&!is.na(m$BMI.for.age)|abs(m$height.for.age)>2&!is.na(m$height.for.age)] <- m$healthyref[abs(m$BMI.for.age)>2&!is.na(m$BMI.for.age)|abs(m$height.for.age)>2&!is.na(m$height.for.age)] +1
m$healthyref[abs(m$BMI.for.age_5yr)>2&!is.na(m$BMI.for.age_5yr)|abs(m$height.for.age_5yr)>2&!is.na(m$height.for.age_5yr)] <- m$healthyref[abs(m$BMI.for.age_5yr)>2&!is.na(m$BMI.for.age_5yr)|abs(m$height.for.age_5yr)>2&!is.na(m$height.for.age_5yr)] +1
table(m$healthyref,m$age,m$Trajectory2)



g <- read.delim("genus.txt")


mgen <- cbind(m,g/rowSums(g))#read.delim("final_helmi_metagenus.txt")
mgenlog <- cbind(m,log((g+1)/rowSums(g)))#read.delim("final_helmi_metagenus.txt")
mgenlog <- mgenlog[mgenlog$use==1,]

top20_6<-names(g[,colSums(g[m$age_weeks<30,]>0)>0.3*nrow(g[m$age_weeks<30,])])
top20_24<-names(g[,colSums(g[m$age_weeks>30,]>0)>0.3*nrow(g[m$age_weeks>30,])])
top20<-unique(c(top20_6,top20_24))

topcor <- cor(mgen[,top20])
topcor[topcor==1] <- 0
apply(topcor,FUN=max,MARGIN = 2)[order(apply(topcor,FUN=max,MARGIN = 2))]

quartz();
gplots::heatmap.2(topcor[names(apply(topcor,FUN=max,MARGIN = 2)[apply(topcor,FUN=max,MARGIN = 2)>0.45]),
                          names(apply(topcor,FUN=max,MARGIN = 2)[apply(topcor,FUN=max,MARGIN = 2)>0.45])],
                   col = brewer.pal(n=9,name="BrBG"),trace=NULL,margins = c(15,15),
          cexRow = 0.5,cexCol = 0.5)


#final model:

top20_6<-names(g[,colSums(g[m$age_weeks<30,]>0)>0.3*nrow(g[m$age_weeks<30,])])
top20_24<-names(g[,colSums(g[m$age_weeks>30,]>0)>0.3*nrow(g[m$age_weeks>30,])])
top20_24<-names(sort(colMeans(g[,top20_24]/rowSums(g[,top20_24])),decreasing = T))[c(1:65)]
#top20_24<- top20_24[-45]

top20<-unique(c(top20_6,top20_24))

gdata2<-mgenlog
gdata2<-gdata2[,names(gdata2) %in% top20]
gdata2$Trajectory<-NA
gdata2$Trajectory[mgenlog$Trajectory2==1 & mgenlog$healthyref==0]<-1
gdata2$Trajectory[mgenlog$Trajectory2!=1 & mgenlog$healthyref==0]<-0.5
gdata2$Trajectory[mgenlog$healthyref>0]<-0
gdata2$age_weeks<-mgenlog$age_weeks

g11.glm.healthyref_24mo_d2 <- step(glm(formula = Trajectory ~ (age_weeks+I(age_weeks^2))*.,
                                       data = gdata2[gdata2$age_weeks>30,c(top20_24,"age_weeks","Trajectory")],
                                       family =binomial),k=3)


top20_1<-names(g[,colSums(g[m$age_weeks<40,]>0)>0.3*nrow(g[m$age_weeks<40,])])
top20_2<-names(g[,colSums(g[m$age_weeks>20,]>0)>0.3*nrow(g[m$age_weeks>20,])])
top20 <- unique(c(top20_1,top20_2))

gdata2<-mgenlog
gdata2<-gdata2[,names(gdata2) %in% top20]
gdata2$Trajectory<-NA
gdata2$Trajectory[mgenlog$Trajectory2==1 & mgenlog$healthyref==0]<-1
gdata2$Trajectory[mgenlog$Trajectory2!=1 & mgenlog$healthyref==0]<-0.5
gdata2$Trajectory[mgenlog$healthyref>0]<-0
gdata2$age_weeks<-mgenlog$age_weeks
gdata2$Family_ID <-mgenlog$Family_ID
  
g11.glm.healthyref_6mo_d <- step(glm(formula = Trajectory ~ (age_weeks+I(age_weeks^2))*.,
                                     data = gdata2[gdata2$age_weeks<40,c(top20_1,"age_weeks","Trajectory")],
                                     family =binomial))







tdata$index <-NA
tdata$index[tdata$age_weeks<30] <- predict(g11.glm.healthyref_6mo_d,newdata=tdata[tdata$age_weeks<30,],type="response")
tdata$index[tdata$age_weeks>30] <- predict(g11.glm.healthyref_24mo_d2,newdata=tdata[tdata$age_weeks>30,],type="response")


for(j in unique(tdata$age)){
  tdata[tdata$age==j,"index"] <- (tdata[tdata$age==j,"index"]-min(tdata[tdata$age==j,"index"],na.rm=T))/(max(tdata[tdata$age==j,"index"],na.rm=T)-min(tdata[tdata$age==j,"index"],na.rm=T))
}




tdata$health_group <- NA
tdata$health_group[tdata$excessive_crying==1]<-"Crying" 
tdata$health_group[tdata$healthyref==0]<-"Control" 
tdata$health_group[tdata$PoorHealth_0to2yr==1]<-"Poor health" 
tdata$health_group[abs(tdata$BMI.for.age)>2|abs(tdata$height.for.age)>2|abs(tdata$BMI.for.age_5yr)>2|abs(tdata$height.for.age_5yr)>2]<-"Growth" 
tdata$health_group[tdata$AllergyDiagnosis_2yr==1|tdata$AtopySymptoms_2yr==1|tdata$AsthmaDiagnosis_2yr==1|tdata$AtopyDiagnosis_5yr==1|
                     tdata$AsthmaDiagnosis_5yr==1|tdata$RhinitisDiagnosis_5yr==1]<-"Allergic" 


tdata$group <-as.character(tdata[,"health_group"]) 
tdata$group[tdata$group=="Crying"]<-NA
tdata$group[tdata$group=="Poor health"]<-NA
tdata$group <-as.factor(tdata$group)  

p <- ggplot(tdata[!is.na(tdata$group)&tdata$group!="NA"&tdata$group!="3",], 
            aes(y = index, x = age,color=group,fill=group)) +
  stat_smooth(method = "loess", formula = y ~x, se = T) + 
  geom_point(position="jitter",alpha=0.2)+
  theme_bw() + 
  xlab("Age (weeks)") + 
  ylab("Microbiota Wellbeing Index") +
  scale_color_manual(values=c("#E41A1C", "#377EB8" ,"green")[1:length(levels(tdata$group))],
                     aesthetics = c("colour","fill"),guide="none")+
  scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

quartz();plot(p)

response <- c("excessive_crying",grep(names(m),pattern="_0to2yr|_2yr$|_5yr",value=T)[-c(3,4,8,27:31)])

p <- list()

for(i in response){
  tdata$group <-as.character(tdata[,i]) 
  tdata$group[is.na(tdata[,i])]<-"NA"
  tdata$group[tdata$healthyref>0&tdata$group==0]<-"3"
  tdata$group <-as.factor(tdata$group)  
  
  p[[i]] <- ggplot(tdata[!is.na(tdata$group)&tdata$group!="NA"&tdata$group!="3",], 
                   aes(y = index, x = age,color=group,fill=group)) +
    stat_smooth(method = "loess", formula = y ~x, se = T) + 
    geom_point(position="jitter",alpha=0.2)+
    #geom_boxplot()+
    theme_bw() + 
    xlab("Age (weeks)") + 
    ylab("Microbiota Wellbeing Index") +
    scale_color_manual(values=c("#E41A1C", "#377EB8" ,"green")[1:length(levels(tdata$group))],
                       aesthetics = c("colour","fill"),guide="none")+
    scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
    #theme(legend.position="none")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  pdf(paste(i,"indexplots_new_updated.pdf"),width=4,height = 3); plot(p[[i]]);dev.off()
}


tdata$BM <- tdata$birthmode_RecievingIAP
tdata$BM[tdata$BM=="Vaginal_NoAB"]<-"V"
tdata$BM[tdata$BM=="Vaginal_AB"]<-"V-ABX"
tdata$BM[tdata$BM=="C-section"]<-"CS"
tdata$BM[tdata$BM=="CS"|tdata$BM=="V-ABX"]<-"CS/V-ABX"

tdata$var <- paste(tdata$BM,tdata$siblings)

tdata$sibs <- "no"
tdata$sibs[tdata$siblings=="(0,4]"]<-"yes"

tdata$diet3 <- tdata$diet2
tdata$diet3[tdata$diet2=="BF_F"|tdata$diet2=="F_S"]<-"F"
tdata$diet3[tdata$diet2=="BF_S"]<-"BF"

for(i in c("diet3","sibs")){
  tdata$group <-as.character(tdata[,i]) 
  tdata$group <-as.factor(tdata$group)  
  #tdata<-tdata[!is.na(tdata$group)&tdata$group!="NA",]
  p[[i]] <- ggplot(tdata[!is.na(tdata$group),], 
                   aes(y = index, x = age,color=birthmode,fill=birthmode, linetype=factor(group))) + 
    stat_smooth(method = "loess", formula = y ~ x, se = T) + 
    theme_bw() + 
    xlab("Age (weeks)") + 
    ylab("Index") +
    #ylim(c(0.1,0.3))+
    scale_color_manual(values=c("royalblue","red","orange"),
                       aesthetics = c("colour","fill"),guide="none")+
    scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  pdf(paste(i,"indexplots_new.pdf"),width=5,height = 4); plot(p[[i]]);dev.off()
}


p[[i]] <- ggplot(tdata[!is.na(tdata$birthmode)&tdata$diet3=="BF",], 
                 aes(y = index, x = age,color=birthmode,fill=birthmode)) + 
  stat_smooth(method = "loess", formula = y ~ x, se = T) + 
  theme_bw() + 
  xlab("Age (weeks)") + 
  ylab("Index") +
  coord_cartesian(ylim = c(0,0.7))+
  scale_color_manual(values=c("royalblue","red","cyan"),
                     aesthetics = c("colour","fill"),guide="none")+
  scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
quartz(height = 3,width = 3);plot(p[[i]])


p[[i]] <- ggplot(tdata[!is.na(tdata$birthmode)&tdata$diet3!="BF",], 
                 aes(y = index, x = age,color=birthmode,fill=birthmode)) + 
  stat_smooth(method = "loess", formula = y ~ x, se = T) + 
  theme_bw() + 
  xlab("Age (weeks)") + 
  ylab("Index") +
  coord_cartesian(ylim = c(0,0.5))+
  scale_color_manual(values=c("royalblue","red","cyan"),
                     aesthetics = c("colour","fill"),guide="none")+
  scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
quartz(height = 3,width = 3);plot(p[[i]])







p[[i]] <- ggplot(tdata[!is.na(tdata$birthmode)&tdata$sibs=="yes"&tdata$diet3=="BF",], 
                 aes(y = index, x = age,color=birthmode,fill=birthmode)) + 
  stat_smooth(method = "loess", formula = y ~ x, se = T) + 
  theme_bw() + 
  xlab("Age (weeks)") + 
  ylab("Index") +
  coord_cartesian(ylim = c(0,0.5))+
  scale_color_manual(values=c("royalblue","red","cyan"),
                     aesthetics = c("colour","fill"),guide="none")+
  scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
quartz(height = 3,width = 3);plot(p[[i]])

p[[i]] <- ggplot(tdata[!is.na(tdata$birthmode)&tdata$sibs=="no"&tdata$diet3=="BF",], 
                 aes(y = index, x = age,color=birthmode,fill=birthmode)) + 
  stat_smooth(method = "loess", formula = y ~ x, se = T) + 
  theme_bw() + 
  xlab("Age (weeks)") + 
  ylab("Index") +
  coord_cartesian(ylim = c(0,0.5))+
  scale_color_manual(values=c("royalblue","red","cyan"),
                     aesthetics = c("colour","fill"),guide="none")+
  scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
quartz(height = 3,width = 3);plot(p[[i]])

p[[i]] <- ggplot(tdata[!is.na(tdata$birthmode)&tdata$sibs=="yes"&tdata$diet3=="F",], 
                 aes(y = index, x = age,color=birthmode,fill=birthmode)) + 
  stat_smooth(method = "loess", formula = y ~ x, se = T) + 
  theme_bw() + 
  xlab("Age (weeks)") + 
  ylab("Index") +
  coord_cartesian(ylim = c(0,0.5))+
  scale_color_manual(values=c("royalblue","red","cyan"),
                     aesthetics = c("colour","fill"),guide="none")+
  scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
quartz(height = 3,width = 3);plot(p[[i]])

p[[i]] <- ggplot(tdata[!is.na(tdata$birthmode)&tdata$sibs=="no"&tdata$diet3=="F",], 
                 aes(y = index, x = age,color=birthmode,fill=birthmode) ) + 
  stat_smooth(method = "loess", formula = y ~ x, se = T) + 
  theme_bw() + 
  xlab("Age (weeks)") + 
  ylab("Index") +
  coord_cartesian(ylim = c(0,0.5))+
  scale_color_manual(values=c("royalblue","red","cyan"),
                     aesthetics = c("colour","fill"),guide="none")+
  scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
quartz(height = 3,width = 3);plot(p[[i]])


p[[i]] <- ggplot(tdata[!is.na(tdata$BM)&tdata$sibs=="no"&tdata$diet3=="F",], 
                 aes(y = index, x = age,color=BM,fill=BM) ) + 
  stat_smooth(method = "loess", formula = y ~ x, se = T) + 
  theme_bw() + 
  xlab("Age (weeks)") + 
  ylab("Index") +
  coord_cartesian(ylim = c(0,0.5))+
  scale_color_manual(values=c("royalblue","red","cyan"),
                     aesthetics = c("colour","fill"),guide="none")+
  scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
quartz(height = 3,width = 3);plot(p[[i]])

p[[i]] <- ggplot(tdata[!is.na(tdata$BM)&tdata$sibs=="no"&tdata$diet3=="BF",], 
                 aes(y = index, x = age,color=BM,fill=BM) ) + 
  stat_smooth(method = "loess", formula = y ~ x, se = T) + 
  theme_bw() + 
  xlab("Age (weeks)") + 
  ylab("Index") +
  coord_cartesian(ylim = c(0,0.5))+
  scale_color_manual(values=c("royalblue","red","cyan"),
                     aesthetics = c("colour","fill"),guide="none")+
  scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
quartz(height = 3,width = 3);plot(p[[i]])

p[[i]] <- ggplot(tdata[!is.na(tdata$BM)&tdata$sibs=="yes"&tdata$diet3=="F",], 
                 aes(y = index, x = age,color=BM,fill=BM) ) + 
  stat_smooth(method = "loess", formula = y ~ x, se = T) + 
  theme_bw() + 
  xlab("Age (weeks)") + 
  ylab("Index") +
  coord_cartesian(ylim = c(0,0.5))+
  scale_color_manual(values=c("royalblue","red","cyan"),
                     aesthetics = c("colour","fill"),guide="none")+
  scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
quartz(height = 3,width = 3);plot(p[[i]])

p[[i]] <- ggplot(tdata[!is.na(tdata$BM)&tdata$sibs=="yes"&tdata$diet3=="BF",], 
                 aes(y = index, x = age,color=BM,fill=BM) ) + 
  stat_smooth(method = "loess", formula = y ~ x, se = T) + 
  theme_bw() + 
  xlab("Age (weeks)") + 
  ylab("Index") +
  coord_cartesian(ylim = c(0,0.5))+
  scale_color_manual(values=c("royalblue","red","cyan"),
                     aesthetics = c("colour","fill"),guide="none")+
  scale_x_sqrt(breaks=(c(3,6,13,26,39,52,78,104)),labels=c(3,6,13,26,39,52,78,104))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
quartz(height = 3,width = 3);plot(p[[i]])


#health against index------------------------------------------------------------------------------------
indexmodsadj2 <- matrix(ncol=length(response),nrow=8)
colnames(indexmodsadj2) <- response
rownames(indexmodsadj2) <- unique(tdata$age)
indexmodsadj2p <- indexmodsadj2


resp_select <- response[c(1,3:5,10:13,15,17:23,25:27)]

for(i in response[c(1:5,10:13,15:27)]){
  for(j in unique(tdata$age)){
    mod <- step(glm(as.formula(paste(i," ~ m_edu+m_BMI+m_weightGain+m_AgeDelivery+m_DiagnosisPregnancyType_GestationalDiabetes+
                                          m_FattyAcid+m_FolicAcid+m_SmokingPrior+m_HadAllergyDisease+f_HadAllergyDisease+
                                          p_edu+p_SmokingDuring+env_FurOrFeathers+inf_Sex+index")),
                    family="binomial",
                    data=na.omit(tdata[tdata$age==j&tdata$healthyref==0|
                                         tdata$age==j&tdata$healthyref>0&tdata[,i]==1,
                                       c(i,"index","m_edu","m_BMI","m_weightGain","m_AgeDelivery",
                                         "m_DiagnosisPregnancyType_GestationalDiabetes",
                                         "m_FattyAcid","m_FolicAcid","m_SmokingPrior","m_HadAllergyDisease","f_HadAllergyDisease",
                                         "p_edu","p_SmokingDuring","env_FurOrFeathers","inf_Sex")])))
    
    if("index"%in%rownames(summary(mod)$coef)){
      indexmodsadj2[as.character(j),i] <- summary(mod)$coef["index",1]
      indexmodsadj2p[as.character(j),i] <- summary(mod)$coef["index",4]
    }
    
  }
}


indexmodsadj2[is.na(indexmodsadj2)]<-0

indexmodsadj2_sym <- indexmodsadj2p
indexmodsadj2_sym[indexmodsadj2p<0.05]<-"*"
indexmodsadj2_sym[indexmodsadj2p<0.01]<-"**"
indexmodsadj2_sym[indexmodsadj2p<0.001]<-"***"
indexmodsadj2_sym[indexmodsadj2p>0.05|is.na(indexmodsadj2p)]<-""

indexmodsadj2_minimumAdjusted<-indexmodsadj2
indexmodsadj2_minimumAdjusted[indexmodsadj2_minimumAdjusted==min(indexmodsadj2_minimumAdjusted)]<-
  min(indexmodsadj2_minimumAdjusted[indexmodsadj2_minimumAdjusted>min(indexmodsadj2_minimumAdjusted)])
quartz();heatmap.2(as.matrix(t(indexmodsadj2[c("3","6","13","26","39","52","78","104"),resp_select])),
                   col = colorRampPalette(colors=c("#377EB8","white","firebrick")),
                   dendrogram = "none",trace = "none",Colv = "none",Rowv="none",
                   cellnote = as.matrix(t(indexmodsadj2_sym[c("3","6","13","26","39","52","78","104"),resp_select])),
                   notecol = "black",margins=c(10,10),
                   main='indexmodsadj2')



colnames(indexmodsadj2_minimumAdjusted)

pdf('IndexMOdel_Adjusted2.pdf',width=18,height=16)
par(oma=c(4,1,1,30))
heatmap.2(as.matrix(t(indexmodsadj2_minimumAdjusted[c("3","6","13","26","39","52","78","104"),resp_select])),
          col = colorRampPalette(colors=c("#377EB8","white","firebrick")),
          dendrogram = "none",trace = "none",Colv = "none",Rowv="none",
          cellnote = as.matrix(t(indexmodsadj2_sym[c("3","6","13","26","39","52","78","104"),resp_select])),
          notecol = "black",#margins=c(10,10)
          cexRow = 4,cexCol = 4,notecex = 4,
          lmat=rbind(c(5,4,2),c(6,1,3)),lhei=c(2,16),lwid=c(1,2,1)
          )
dev.off() 




indexmods <- matrix(ncol=length(response),nrow=8)
colnames(indexmods) <- response
rownames(indexmods) <- unique(tdata$age)
indexmodsp <- indexmods

for(i in response){
  for(j in unique(tdata$age)){
    mod <- (glm(as.formula(paste(i," ~ index")),
                family="binomial",
                data=na.omit(tdata[tdata$age==j&tdata$healthyref==0|
                                     tdata$age==j&tdata$healthyref>0&tdata[,i]==1,
                                   c(i,"index")])))
    
    if("index"%in%rownames(summary(mod)$coef)){
      indexmods[as.character(j),i] <- summary(mod)$coef["index",1]
      indexmodsp[as.character(j),i] <- summary(mod)$coef["index",4]
    }
    
  }
}

indexmods_sym <- indexmodsp
indexmods_sym[indexmodsp<0.05]<-"*"
indexmods_sym[indexmodsp<0.01]<-"**"
indexmods_sym[indexmodsp<0.001]<-"***"
indexmods_sym[indexmodsp>0.05|is.na(indexmodsp)]<-""

quartz();heatmap.2(as.matrix(t(indexmods[c("3","6","13","26","39","52","78","104"),resp_select])),
                   col = colorRampPalette(colors=c("#377EB8","white","firebrick")),
                   dendrogram = "none",trace = "none",Colv = "none",Rowv='none',
                   cellnote = as.matrix(t(indexmods_sym[c("3","6","13","26","39","52","78","104"),resp_select])),
                   notecol = "black",margins=c(10,10),
                   main='indexmods')


#mediation of birth mode effects--------------------------------------------------------------------------
summary(lm(index~birthmode,data=tdata[tdata$age==3,]))

bmmods <-list()
for(i in response){
  bmmods[[i]] <- glm(as.formula(paste(i," ~ birthmode")),family="binomial",
            data=tdata[tdata$age==3&tdata$healthyref==0|
                         tdata$age==3&tdata$healthyref>0&tdata[,i]==1,])
}


bmmodsadj <-list()
for(i in response){
  bmmodsadj[[i]] <- step(glm(as.formula(paste(i," ~ m_edu+m_BMI+m_weightGain+m_AgeDelivery+m_DiagnosisPregnancyType_GestationalDiabetes+
                                          m_FattyAcid+m_FolicAcid+m_SmokingPrior+m_HadAllergyDisease+f_HadAllergyDisease+
                                          p_edu+p_SmokingDuring+env_FurOrFeathers+inf_Sex+household_size+birthmode")),family="binomial",
                             data=na.omit(tdata[tdata$age==3&tdata$healthyref==0|
                                                  tdata$age==3&tdata$healthyref>0&tdata[,i]==1,
                                                c(i,"index","m_edu","m_BMI","m_weightGain","m_AgeDelivery","m_DiagnosisPregnancyType_GestationalDiabetes",
                                                  "m_FattyAcid","m_FolicAcid","m_SmokingPrior","m_HadAllergyDisease","f_HadAllergyDisease",
                                                  "p_edu","p_SmokingDuring","env_FurOrFeathers","inf_Sex","household_size","birthmode")])))
  
}

bmmodsadjindex <-list()
for(i in response){
  bmmodsadjindex[[i]] <- step(glm(as.formula(paste(i," ~ index+ m_edu+m_BMI+m_weightGain+m_AgeDelivery+m_DiagnosisPregnancyType_GestationalDiabetes+
                                          m_FattyAcid+m_FolicAcid+m_SmokingPrior+m_HadAllergyDisease+f_HadAllergyDisease+
                                          p_edu+p_SmokingDuring+env_FurOrFeathers+inf_Sex+household_size+birthmode")),family="binomial",
                                   data=na.omit(tdata[tdata$age==3&tdata$healthyref==0|
                                                        tdata$age==3&tdata$healthyref>0&tdata[,i]==1,
                                                      c(i,"index","m_edu","m_BMI","m_weightGain","m_AgeDelivery","m_DiagnosisPregnancyType_GestationalDiabetes",
                                                        "m_FattyAcid","m_FolicAcid","m_SmokingPrior","m_HadAllergyDisease","f_HadAllergyDisease",
                                                        "p_edu","p_SmokingDuring","env_FurOrFeathers","inf_Sex","household_size","birthmode")])))
  
}


#associations with background variables----------------------------------------
indextable <- data.frame(variable=c("index","m_edu","m_BMI","m_weightGain","m_AgeDelivery","m_DiagnosisPregnancyType_GestationalDiabetesTRUE",
                                    "m_FattyAcidTRUE","m_FolicAcidTRUE","m_SmokingPriorTRUE","m_HadAllergyDiseaseYes","f_HadAllergyDiseaseYes",
                                    "p_edu","p_SmokingDuringTRUE","env_FurOrFeathersTRUE","inf_SexGirl","household_size","birthmodeC-section","birthmodeVaginal_AB"))
rownames(indextable)<-indextable$variable
for(i in response){
  indextable[intersect(rownames(summary(bmmodsadjindex2[[i]])$coef),indextable$variable),i] <- summary(bmmodsadjindex[[i]])$coef[intersect(rownames(summary(bmmodsadjindex[[i]])$coef),indextable$variable),1]
}
indextable <- scale(indextable[,-1],center=F)
indextable[is.na(indextable)]<-0

indextable_p <- data.frame(variable=c("index","m_edu","m_BMI","m_weightGain","m_AgeDelivery","m_DiagnosisPregnancyType_GestationalDiabetesTRUE",
                                    "m_FattyAcidTRUE","m_FolicAcidTRUE","m_SmokingPriorTRUE","m_HadAllergyDiseaseYes","f_HadAllergyDiseaseYes",
                                    "p_edu","p_SmokingDuringTRUE","env_FurOrFeathersTRUE","inf_SexGirl","household_size","birthmodeC-section","birthmodeVaginal_AB"))
rownames(indextable_p)<-indextable_p$variable
for(i in response){
  indextable_p[intersect(rownames(summary(bmmodsadjindex[[i]])$coef),indextable_p$variable),i] <- summary(bmmodsadjindex2[[i]])$coef[intersect(rownames(summary(bmmodsadjindex[[i]])$coef),indextable_p$variable),4]
}

indextable_sym <- indextable_p
indextable_sym[indextable_p<0.05]<-"*"
indextable_sym[indextable_p<0.01]<-"**"
indextable_sym[indextable_p<0.001]<-"***"
indextable_sym[indextable_p>0.05|is.na(indextable_p)]<-""

quartz();heatmap.2(as.matrix(indextable),col = colorRampPalette(colors=c("#377EB8","white","firebrick")),
                 dendrogram = "none",trace = "none",cellnote = as.matrix(indextable_sym[,-1]))


tdata$par_edu <- apply(tdata[,c("m_edu","p_edu")],FUN = max,MARGIN = 1)
tdata$m_suppl <- tdata$m_FattyAcid
tdata$m_suppl[tdata$m_FolicAcid=="TRUE"]<-"TRUE"
tdata$hh_smoke <- tdata$p_SmokingDuring
tdata$hh_smoke[tdata$m_SmokingPrior=="TRUE"]<-"TRUE"
tdata$par_allergies <- tdata$m_HadAllergyDisease
tdata$par_allergies[tdata$f_HadAllergyDisease=="Yes"]<-"TRUE"

tdata$m_weigth <- cut(tdata$m_BMI,breaks=c(0,21,25,30,100),labels = c("under","a","over","obese"))

indexmod <-list()
for(i in as.factor(unique(tdata$age))){
  indexmod[[i]] <- step(lm(as.formula(paste("index ~ diet2+m_edu+m_weigth+#m_weightGain+m_AgeDelivery+m_DiagnosisPregnancyType_GestationalDiabetes+
                                          m_suppl+m_SmokingPrior+m_HadAllergyDisease+f_HadAllergyDisease+
                                          p_edu+p_SmokingDuring+env_FurOrFeathers+inf_Sex+household_size+birthmode+Antibiotic+inf_Prob_Bifidobacterium+inf_Prob_Lactobacillus")),
                                   data=na.omit(tdata[tdata$age==i,
                                                          c("index","diet2","m_edu","m_weigth",
                                                            "m_suppl","m_SmokingPrior","m_HadAllergyDisease","f_HadAllergyDisease",
                                                            "p_edu","p_SmokingDuring","env_FurOrFeathers","inf_Sex","household_size","birthmode","Antibiotic","inf_Prob_Bifidobacterium","inf_Prob_Lactobacillus")])))
  
}

indextable <- data.frame(variable=c("diet2BF_F","diet2BF_S","diet2F_S","m_edu","m_weigthunder","m_weigthover","m_weigthobese",#"m_weightGain","m_AgeDelivery","m_DiagnosisPregnancyType_GestationalDiabetesTRUE",
                                    "m_supplTRUE","m_SmokingPriorTRUE","m_HadAllergyDiseaseYes","f_HadAllergyDiseaseYes",
                                    "p_edu","p_SmokingDuringTRUE","env_FurOrFeathersTRUE","inf_SexGirl","household_size","birthmodeC-section","birthmodeVaginal_AB",
                                    "Antibiotic","inf_Prob_BifidobacteriumTRUE","inf_Prob_LactobacillusTRUE"))
rownames(indextable)<-indextable$variable
for(i in as.factor(unique(tdata$age))){
  indextable[intersect(rownames(summary(indexmod[[i]])$coef),indextable$variable),i] <- summary(indexmod[[i]])$coef[intersect(rownames(summary(indexmod[[i]])$coef),indextable$variable),1]
}
indextable <- scale(indextable[,-1],center=F)
indextable[is.na(indextable)]<-0

indextable_p <- data.frame(variable=c("diet2BF_F","diet2BF_S","diet2F_S","m_edu","m_weigthunder","m_weigthover","m_weigthobese",#"m_weightGain","m_AgeDelivery","m_DiagnosisPregnancyType_GestationalDiabetesTRUE",
                                    "m_supplTRUE","m_SmokingPriorTRUE","m_HadAllergyDiseaseYes","f_HadAllergyDiseaseYes",
                                    "p_edu","p_SmokingDuringTRUE","env_FurOrFeathersTRUE","inf_SexGirl","household_size","birthmodeC-section","birthmodeVaginal_AB",
                                    "Antibiotic","inf_Prob_BifidobacteriumTRUE","inf_Prob_LactobacillusTRUE"))
rownames(indextable_p)<-indextable_p$variable
for(i in as.factor(unique(tdata$age))){
  indextable_p[intersect(rownames(summary(indexmod[[i]])$coef),indextable_p$variable),i] <- summary(indexmod[[i]])$coef[intersect(rownames(summary(indexmod[[i]])$coef),indextable_p$variable),4]
}
indextable_sym <- indextable_p
indextable_sym[indextable_p<0.05]<-"*"
indextable_sym[indextable_p<0.01]<-"**"
indextable_sym[indextable_p<0.001]<-"***"
indextable_sym[indextable_p>0.05|is.na(indextable_p)]<-""

indextable[abs(rowMeans(indextable[,-1]))>0,]

quartz();
pdf("index_background_variables.pdf")
heatmap.2(as.matrix(indextable[abs(rowMeans(indextable[,-1]))>0,c("3","6","13","26","39","52","78","104")]),
                   col = colorRampPalette(colors=c("#377EB8","white","firebrick")),
                   Rowv=F,Colv = F,notecol = "black",
                    dendrogram = "none",trace = "none",margins = c(10,20),
                   cellnote = as.matrix(indextable_sym[abs(rowMeans(indextable[,-1]))>0,c("3","6","13","26","39","52","78","104")]))
dev.off()

#END HERE--------------------------------------------------------

#predictions based on index----------------------------------
library(pROC)
ylrd_ramp <-colorRampPalette(colors=c("white","firebrick"))

predmod <- list()
predroc <- list()
indexp <-list()
aucs <- list()
est <- list()



h="index"

tdata$index <-NA
tdata$index[tdata$age_weeks<30] <- predict(update(g11.glm.healthyref_6mo_d2,.~.,data=tdata[rownames(tdata[tdata$Family_ID%in%trainingD&tdata$age_weeks<30,]),]),newdata=tdata[tdata$age_weeks<30,],type="response")
tdata$index[tdata$age_weeks>30] <- predict(update(g11.glm.healthyref_24mo_d2,.~.,data=tdata[rownames(tdata[tdata$Family_ID%in%trainingD&tdata$age_weeks>30,]),]),newdata=tdata[tdata$age_weeks>30,],type="response")

tdata$index <-NA
tdata$index[tdata$age_weeks<30] <- predict(update(g11.glm.healthyref_6mo_d,.~.,data=tdata[rownames(tdata[tdata$Family_ID%in%trainingD&tdata$age_weeks<40,]),]),newdata=tdata[tdata$age_weeks<30,],type="response")
tdata$index[tdata$age_weeks>30] <- predict(update(g11.glm.healthyref_24mo_d,.~.,data=tdata[rownames(tdata[tdata$Family_ID%in%trainingD&tdata$age_weeks>20,]),]),newdata=tdata[tdata$age_weeks>30,],type="response")


for(h in 1:10){
indexp[[h]] <- data.frame(BMIgroup1=rep(NA,8))
rownames(indexp[[h]]) <- c("3","6","13","26","39","52","78","104")
aucs[[h]] <- data.frame(BMIgroup1=rep(NA,8))
rownames(aucs[[h]]) <- c("3","6","13","26","39","52","78","104")
est[[h]] <- data.frame(BMIgroup1=rep(NA,8))
rownames(est[[h]]) <- c("3","6","13","26","39","52","78","104")

trainingD<-sample(unique(tdata$Family_ID),length(unique(tdata$Family_ID))*0.75)#sample(rownames(tdata),nrow(tdata)*0.75)
testD <- setdiff(unique(tdata$Family_ID),trainingD)
  
rownames(tdata[tdata$Family_ID%in%trainingD&tdata$age_weeks<30,])

tdata$index <-NA
tdata$index[tdata$age_weeks<30] <- predict(update(g11.glm.healthyref_6mo_d,.~.,data=tdata[rownames(tdata[tdata$Family_ID%in%trainingD&tdata$age_weeks<30,]),]),newdata=tdata[tdata$age_weeks<30,],type="response")
tdata$index[tdata$age_weeks>30] <- predict(update(g11.glm.healthyref_24mo_d,.~.,data=tdata[rownames(tdata[tdata$Family_ID%in%trainingD&tdata$age_weeks>30,]),]),newdata=tdata[tdata$age_weeks>30,],type="response")


for(j in unique(tdata$Family_ID)){
  for(i in unique(tdata$age)){
    tdata$index[tdata$age==i&tdata$Family_ID==j] <- tdata$index[tdata$age==i&tdata$Family_ID==j]+(0-median(tdata$index[tdata$age==i&tdata$healthUse==1]))
  }
} 

for(j in unique(tdata$age)){
  tdata[tdata$age==j,"index"] <- (tdata[tdata$age==j,"index"]-min(tdata[tdata$age==j,"index"],na.rm=T))/(max(tdata[tdata$age==j,"index"],na.rm=T)-min(tdata[tdata$age==j,"index"],na.rm=T))
}

for(j in response){
  tmp <- tdata[rownames(tdata[tdata$Family_ID%in%testD,]),]
  tmp$group <-tmp[,j]
  tmp$group[tmp$healthyref>0&tmp$group==0]<-3
  if(nrow(tmp[tmp$group==1&!is.na(tmp$group),])>30){
  for(i in rownames(indexp[[h]])){
    indexp[[h]][i,j] <- summary(lm(index~group,data=tmp[tmp$age==i&tmp$group!="3",]))$coef["group",4]
    predmod[[paste(h,i,j)]] <- glm(as.numeric(group)~index,family="binomial",
                                   data=tmp[tmp$age==i&tmp$group!="3"&tmp$group!="NA",])
    predroc[[paste(h,i,j)]] <- roc(response=as.numeric(tmp$group[!is.na(tmp$index)&
                                                                   !is.na(tmp$group)&tmp$group!="NA"&
                                                                   tmp$group!="3"&tmp$age==i]), 
                                   predictor = predict(predmod[[paste(h,i,j)]],type="response"))
    aucs[[h]][i,j] <-predroc[[paste(h,i,j)]]$auc
    est[[h]][i,j]  <- summary(predmod[[paste(h,i,j)]])$coefficients[2,1]
  }  
  }
}
}

auc_data <- array(unlist(aucs) , c(8,29,3) )
auc_mean <- apply(auc_data, 1:2 , mean,na.rm=T )
colnames(auc_mean) <- response#colnames(aucs[[1]])[-1]
rownames(auc_mean) <- rownames(aucs[[1]])

indexp_data <- array(unlist(indexp) , c(8,29,3) )
indexp_mean <- apply(indexp_data, 1:2 , mean,na.rm=T )
colnames(indexp_mean) <- response#colnames(indexp[[1]])[-1]
rownames(indexp_mean) <- rownames(indexp[[1]])


indexsym <- indexp_mean
indexsym[indexp_mean<0.1]<-"."
indexsym[indexp_mean<0.05]<-"*"
indexsym[indexp_mean<0.01]<-"**"
indexsym[indexp_mean<0.001]<-"***"
indexsym[indexp_mean>0.1]<-""  

quartz();
heatmap.2(as.matrix(t(auc_mean)),Colv = NULL,dendrogram = "row",trace = "none",
          col=c(ylrd_ramp(50)),margins = c(5,15),density.info = "none",
          key.title ="AUC",key.xlab = "AUC",key.ylab = NULL,keysize = 1,
          cellnote = t(indexsym),notecol="black")




auc_mean <- matrix(nrow = 8,ncol=18)
for(i in 1:8){
  for(j in response){
    if(j%in%colnames(aucs[[h]])){
  auc_mean[i,j] <- unlist(lapply(aucs[[h]][i,j],mean,na.rm=T))  
  }
  } 
  }

indexp_mean <- matrix(nrow = 8,ncol=18)
for(i in 1:8) indexp_mean[i,] <- unlist(lapply(indexp[[h]][i,],median,na.rm=T))
colnames(indexp_mean) <- names(lapply(indexp[[h]][i,],mean,na.rm=T))
rownames(indexp_mean) <- rownames(indexp[[h]])
indexp_mean <- indexp_mean[,-1]
 

