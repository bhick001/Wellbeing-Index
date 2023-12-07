library(plyr)

########################################
#Calculating Stability between time points
m<-read.delim('meta.txt')
m$cluster_use<- as.numeric(ordered(m$pearson_genus_kmeans))

meta<-m
#meta<-read.delim('final_helmi_meta_forAnalysis.txt')
f<-read.delim('genus.txt')
f<-f[meta$use==1,]
reltaxa <- (f+1)/meta$ReadCount
meta<-meta[meta$use==1,]

time <- meta$age

for (i in names(table(meta[, 'Family_ID'])[table(meta[, 'Family_ID']) > 1])) {
    time[meta[,'Family_ID'] == i] <- as.numeric(ordered(time[meta[,'Family_ID'] == i]))
    for (j in time[meta[,'Family_ID'] == i][order(time[meta[,'Family_ID'] == i])][-1]) {
        meta$IntraindividualPearson[meta[,'Family_ID'] == i &  time == j] <- cor((log(t(reltaxa[meta[,'Family_ID'] ==  i & time == j, ]))), (log(t(reltaxa[meta[,'Family_ID'] ==  i & time == (j - 1), ]))))

    }
}

#stability boxplots for time and trajectory groups
quartz();par(mfrow=c(1,5))
boxplot(meta$IntraindividualPearson[meta$Trajectory_use=='1'] ~ meta$age_weeks[meta$Trajectory_use=='1'],
        main=paste('Group 1 N =',length(meta$IntraindividualPearson[meta$Trajectory_use=='1'])),
        notch=T,col='red',xlab='',ylab='')
boxplot(meta$IntraindividualPearson[meta$Trajectory_use=='2a'] ~ meta$age_weeks[meta$Trajectory_use=='2a'],
        main=paste('Group 2 N =',length(meta$IntraindividualPearson[meta$Trajectory_use=='2a'])),
        notch=T,col='green',xlab='',ylab='')
boxplot(meta$IntraindividualPearson[meta$Trajectory_use=='2c'] ~ meta$age_weeks[meta$Trajectory_use=='2c'],
        main=paste('Group 3 N =',length(meta$IntraindividualPearson[meta$Trajectory_use=='3'])),
        notch=T,col='blue',xlab='',ylab='')
boxplot(meta$IntraindividualPearson[meta$Trajectory_use=='2b'] ~ meta$age_weeks[meta$Trajectory_use=='2b'],
        main=paste('Group 4 N =',length(meta$IntraindividualPearson[meta$Trajectory_use=='2b'])),
        notch=T,col='cyan',xlab='',ylab='')
boxplot(meta$IntraindividualPearson[meta$Trajectory_use=='3'] ~ meta$age_weeks[meta$Trajectory_use=='3'],
        main=paste('Group 5 N =',length(meta$IntraindividualPearson[meta$Trajectory_use=='3'])),
        notch=T,col='orange3',xlab='',ylab='')




#distance from eaach other within each time point
f3<-f[meta$age_weeks==3,]
f3d <- peardist(f3) 
f6<-f[meta$age_weeks==6,]
f6d <- peardist(f6) 
f12<-f[meta$age_weeks==12,]
f12d <- peardist(f12) 
f24<-f[meta$age_weeks==24,]
f24d <- peardist(f24) 
f36<-f[meta$age_weeks==36,]
f36d <- peardist(f36) 
f48<-f[meta$age_weeks==48,]
f48d <- peardist(f48) 
f72<-f[meta$age_weeks==72,]
f72d <- peardist(f72) 
f96<-f[meta$age_weeks==96,]
f96d <- peardist(f96) 

df<-data.frame(X=c(f3d,f6d,f12d,f24d,f36d,f48d,f72d,f96d),
               Age=rep(c(3,6,12,24,36,48,72,96),
                       times=c(length(f3d),length(f6d),length(f12d),
                               length(f24d),length(f36d),length(f48d),
                               length(f72d),length(f96d)
                              )
                       )
               )
quartz();boxplot(X~Age,data=df)
meandist<-data.frame(matrix(ncol=1, nrow=8))
meandist$age_weeks<-c(sort(unique(meta$age_weeks)))
meandist$meanAge<-c(mean(f3d),mean(f6d),mean(f12d),mean(f24d),
               mean(f36d),mean(f48d),mean(f72d),mean(f96d))
meandist$sdAge<-c(sd(f3d),sd(f6d),sd(f12d),sd(f24d),
               sd(f36d),sd(f48d),sd(f72d),sd(f96d))
quartz();boxplot(meandist$meanAge~meandist$age_weeks)



#Static to compare other Ts at each age to T1
meta$T1<-meta$Trajectory_use
meta$T2<-meta$Trajectory_use
meta$T3<-meta$Trajectory_use
meta$T4<-meta$Trajectory_use
meta$T5<-meta$Trajectory_use

meta$T2[meta$T2=='1']<-0
meta$T2[meta$T2=='2a']<-1 
meta$T2[meta$Trajectory_use=='2b'|meta$Trajectory_use=='2c'|meta$Trajectory_use=='3']<-NA
meta$T2<-as.numeric(meta$T2)
meta$T3[meta$T3=='1']<-0
meta$T3[meta$T3=='2c']<-1 
meta$T3[meta$Trajectory_use=='2b'|meta$Trajectory_use=='2a'|meta$Trajectory_use=='3']<-NA
meta$T3<-as.numeric(meta$T3)
meta$T4[meta$T4=='1']<-0
meta$T4[meta$T4=='2b']<-1 
meta$T4[meta$Trajectory_use=='2a'|meta$Trajectory_use=='2c'|meta$Trajectory_use=='3']<-NA
meta$T4<-as.numeric(meta$T4)

meta$T5[meta$T5=='1']<-0
meta$T5[meta$T5=='3']<-1 
meta$T5[meta$Trajectory_use=='2b'|meta$Trajectory_use=='2c'|meta$Trajectory_use=='2a']<-NA
meta$T5<-as.numeric(meta$T5)


m6<-meta[meta$age==6,]
m13<-meta[meta$age==13,]
m26<-meta[meta$age==26,]
m39<-meta[meta$age==39,]
m52<-meta[meta$age==52,]
m78<-meta[meta$age==78,]
m104<-meta[meta$age==104,]

m2<-data.frame(age=sort(unique(meta$age))[-1],
               pvalue=rep(NA,length(sort(unique(meta$age))[-1])),
               estimate=rep(NA,length(sort(unique(meta$age))[-1]))
               )
for (ii in sort(unique(meta$age))[-1]) {
  print(ii)
  mm<-meta[meta$age==ii,]
  m<-lm(IntraindividualPearson~T2,data=mm)
  s<-summary(m)
  m2$pvalue[m5$age==ii]<-s$coefficients[2,4]
  m2$estimate[m5$age==ii]<-s$coefficients[2,1]
}
m2$pvalue[m2$pvalue>0.05]<-NA

m3<-data.frame(age=sort(unique(meta$age))[-1],
               pvalue=rep(NA,length(sort(unique(meta$age))[-1])),
               estimate=rep(NA,length(sort(unique(meta$age))[-1]))
               )
for (ii in sort(unique(meta$age))[-1]) {
  print(ii)
  mm<-meta[meta$age==ii,]
  m<-lm(IntraindividualPearson~T3,data=mm)
  s<-summary(m)
  m3$pvalue[m5$age==ii]<-s$coefficients[2,4]
  m3$estimate[m5$age==ii]<-s$coefficients[2,1]
}
m3$pvalue[m3$pvalue>0.05]<-NA


m4<-data.frame(age=sort(unique(meta$age))[-1],
               pvalue=rep(NA,length(sort(unique(meta$age))[-1])),
               estimate=rep(NA,length(sort(unique(meta$age))[-1]))
               )
for (ii in sort(unique(meta$age))[-1]) {
  mm<-meta[meta$age==ii,]
  m<-lm(IntraindividualPearson~T4,data=mm)
  s<-summary(m)
  m4$pvalue[m5$age==ii]<-s$coefficients[2,4]
  m4$estimate[m5$age==ii]<-s$coefficients[2,1]
}
m4$pvalue[m4$pvalue>0.05]<-NA


m5<-data.frame(age=sort(unique(meta$age))[-1],
               pvalue=rep(NA,length(sort(unique(meta$age))[-1])),
               estimate=rep(NA,length(sort(unique(meta$age))[-1]))
               )
for (ii in sort(unique(meta$age))[-1]) {
  mm<-meta[meta$age==ii,]
  m<-lm(IntraindividualPearson~T5,data=mm)
  s<-summary(m)
  m5$pvalue[m5$age==ii]<-s$coefficients[2,4]
  m5$estimate[m5$age==ii]<-s$coefficients[2,1]
}
m5$pvalue[m5$pvalue>0.05]<-NA


df<-meta[,c('age','IntraindividualPearson',"Trajectory_use")]
df$Trajectory_use[meta$Trajectory_use=='3']<-5
df$Trajectory_use[meta$Trajectory_use=='2a']<-2
df$Trajectory_use[meta$Trajectory_use=='2b']<-4
df$Trajectory_use[meta$Trajectory_use=='2c']<-3
names(df)<-c('Age','IntraindividualDistance','Trajectory')
df<-df[!is.na(df$Trajectory),]


pdf('Stability_by_trajectory.pdf',width=8,height = 4)
ggplot(df[df$Age>3,], aes(x=as.factor(Age), y=IntraindividualDistance, fill=Trajectory)) + 
    geom_boxplot(outlier.colour = NULL)+
  theme(
      text = element_text(size=20),
      plot.title = element_text(size=20),
      axis.line = element_line(size=1.5,colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.text=element_text(size=16)
    )+
  coord_cartesian(ylim=c(0, 1))+
   #labs(size="Age (weeks)", shape='', colour="")+
  ylab("Intraindividual Distance")+ xlab("Age (weeks)")

dev.off()