#GMBT-----
library(scales)
library('gbmt')
library('plyr')
library(lessR)
library(gplots)
library(RColorBrewer)
setwd('/Users/kkorpela/Documents/brandon/analysis/wellbeing_index/helmi/helmiRDP/HelmiRDPCombined')

peardist <- function(x,distance){
as.dist(1-cor(t(x)))
}
`%notin%` <- Negate(`%in%`)

f<-read.delim('family.txt',sep='\t')
familyData<-read.delim('meta.txt') 
familyData$sample<-rownames(f)

familyData_keep_firsttimepoints <- list()
familyData_drop<-list()
for (fid in unique(familyData$Family_ID)) {
  dd <- familyData[familyData$Family_ID == fid, ]
  ddage<-sort(dd$age_weeks)
  if (any(3==ddage |6==ddage) & any(12==ddage |24==ddage)) {
    familyData_keep_firsttimepoints[fid] <- fid
  } else {familyData_drop[fid]<-fid }
}

familyData_keep_FirstLasttimepoints <- list()
familyData_drop<-list()
for (fid in unique(familyData$Family_ID)) {
  dd <- familyData[familyData$Family_ID == fid, ]
  ddage<-sort(dd$age_weeks)
  if (any(3==ddage |6==ddage) & any(12==ddage |24==ddage) & any(36==ddage |48==ddage |72==ddage |96==ddage)) {
    familyData_keep_FirstLasttimepoints[fid] <- fid
  } else {familyData_drop[fid]<-fid }
}


familyData$pearson_genus_c<-familyData$pearson_genus_kmeans
f<-f[familyData$use==1,]
reltaxa <- (f)/rowSums(f)

familyData<-familyData[familyData$use==1,]

familyData$pearson_genus_kmeans_score1_01<-familyData$pearson_genus_kmeans_score1
familyData$pearson_genus_kmeans_score1_01[familyData$pearson_genus_kmeans_score1_01>0]<-1
familyData$pearson_genus_kmeans_score2_01<-familyData$pearson_genus_kmeans_score2
familyData$pearson_genus_kmeans_score2_01[familyData$pearson_genus_kmeans_score2_01>0]<-1
familyData$pearson_genus_kmeans_score3_01<-familyData$pearson_genus_kmeans_score3
familyData$pearson_genus_kmeans_score3_01[familyData$pearson_genus_kmeans_score3_01>0]<-1
familyData$pearson_genus_kmeans_score4_01<-familyData$pearson_genus_kmeans_score4
familyData$pearson_genus_kmeans_score4_01[familyData$pearson_genus_kmeans_score4_01>0]<-1

familyData$pearson_genus_kmeans_score1_plus<-familyData$pearson_genus_kmeans_score1
familyData$pearson_genus_kmeans_score1_plus[familyData$pearson_genus_kmeans_score1_plus>0]<-familyData$pearson_genus_kmeans_score1_plus[familyData$pearson_genus_kmeans_score1_plus>0]+0.05
familyData$pearson_genus_kmeans_score2_plus<-familyData$pearson_genus_kmeans_score2
familyData$pearson_genus_kmeans_score2_plus[familyData$pearson_genus_kmeans_score2_plus>0]<-familyData$pearson_genus_kmeans_score2_plus[familyData$pearson_genus_kmeans_score2_plus>0]+0.05
familyData$pearson_genus_kmeans_score3_plus<-familyData$pearson_genus_kmeans_score3
familyData$pearson_genus_kmeans_score3_plus[familyData$pearson_genus_kmeans_score3_plus>0]<-familyData$pearson_genus_kmeans_score3_plus[familyData$pearson_genus_kmeans_score3_plus>0]+0.05
familyData$pearson_genus_kmeans_score4_plus<-familyData$pearson_genus_kmeans_score4
familyData$pearson_genus_kmeans_score4_plus[familyData$pearson_genus_kmeans_score4_plus>0]<-familyData$pearson_genus_kmeans_score4_plus[familyData$pearson_genus_kmeans_score4_plus>0]+0.05

familyData$pearson_genus_kmeans_score1_plus2<-familyData$pearson_genus_kmeans_score1
familyData$pearson_genus_kmeans_score1_plus2[familyData$pearson_genus_kmeans_score1_plus2>0]<-familyData$pearson_genus_kmeans_score1_plus2[familyData$pearson_genus_kmeans_score1_plus2>0]+0.5
familyData$pearson_genus_kmeans_score2_plus2<-familyData$pearson_genus_kmeans_score2
familyData$pearson_genus_kmeans_score2_plus2[familyData$pearson_genus_kmeans_score2_plus2>0]<-familyData$pearson_genus_kmeans_score2_plus2[familyData$pearson_genus_kmeans_score2_plus2>0]+0.5
familyData$pearson_genus_kmeans_score3_plus2<-familyData$pearson_genus_kmeans_score3
familyData$pearson_genus_kmeans_score3_plus2[familyData$pearson_genus_kmeans_score3_plus2>0]<-familyData$pearson_genus_kmeans_score3_plus2[familyData$pearson_genus_kmeans_score3_plus2>0]+0.5
familyData$pearson_genus_kmeans_score4_plus2<-familyData$pearson_genus_kmeans_score4
familyData$pearson_genus_kmeans_score4_plus2[familyData$pearson_genus_kmeans_score4_plus2>0]<-familyData$pearson_genus_kmeans_score4_plus2[familyData$pearson_genus_kmeans_score4_plus2>0]+0.5

familyData_reduced<-familyData[familyData$Family_ID %in% names(familyData_keep),]
familyData_reduced<-familyData[familyData$Family_ID %in% names(familyData_keep_firsttimepoints),]
familyData_reduced<-familyData[familyData$Family_ID %in% names(familyData_keep_FirstLasttimepoints),]

familyData$clustDist1[familyData$pearson_genus_c!=1]<-0
familyData$clustDist2[familyData$pearson_genus_c!=2]<-0
familyData$clustDist3[familyData$pearson_genus_c!=3]<-0
familyData$clustDist4[familyData$pearson_genus_c!=4]<-0

#Trajectories on Cluster goodness/distance
varNames<-c("pearson_genus_kmeans_score1" , "pearson_genus_kmeans_score2","pearson_genus_kmeans_score3","pearson_genus_kmeans_score4"         )
varNames<-c("pearson_genus_kmeans_score1_01" , "pearson_genus_kmeans_score2_01","pearson_genus_kmeans_score3_01","pearson_genus_kmeans_score4_01"         )
varNames<-c("pearson_genus_kmeans_score1_plus" , "pearson_genus_kmeans_score2_plus","pearson_genus_kmeans_score3_plus","pearson_genus_kmeans_score4_plus"         )
varNames<-c("pearson_genus_kmeans_score1_plus2" , "pearson_genus_kmeans_score2_plus2","pearson_genus_kmeans_score3_plus2","pearson_genus_kmeans_score4_plus2"         )

#clusters scores
m_famAll2d3_genus_clust_plus2<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=2, data=familyData_reduced, scaling=0)
m_famAll3d3_genus_clust_plus2<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=3, data=familyData_reduced, scaling=0)
m_famAll4d3_genus_clust_plus2<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=4, data=familyData_reduced, scaling=0)
m_famAll5d3_genus_clust_plus2<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=5, data=familyData_reduced, scaling=0)
m_famAll6d3_genus_clust_plus2<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=6, data=familyData_reduced, scaling=0)
m_famAll7d3_genus_clust_plus2<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=7, data=familyData_reduced, scaling=0)
m_famAll8d3_genus_clust_plus2<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=8, data=familyData_reduced, scaling=0)
m_famAll9d3_genus_clust_plus2<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=9, data=familyData_reduced, scaling=0)
m_famAll10d3_genus_clust_plus2<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=10, data=familyData_reduced, scaling=0)

gbtm_scores<-as.data.frame(rbind(m_famAll2d3_genus_clust_plus2$ic,m_famAll3d3_genus_clust_plus2$ic,m_famAll4d3_genus_clust_plus2$ic,m_famAll5d3_genus_clust_plus2$ic,
           m_famAll6d3_genus_clust_plus2$ic,m_famAll7d3_genus_clust_plus2$ic,m_famAll8d3_genus_clust_plus2$ic,m_famAll9d3_genus_clust_plus2$ic,m_famAll10d3_genus_clust_plus2$ic
                     ))
rownames(gbtm_scores)<-c('g2','g3','g4','g5','g6','g7','g8','g9','g10')
gbtm_scores$group<-c('g2','g3','g4','g5','g6','g7','g8','g9','g10')

colnames(gbtm_scores)<-c('aic','bic','caic',
                         'ssbic','hqic')
quartz();par(mfrow=c(2,2))
plot(x=seq(from=2,to=10),y=gbtm_scores$caic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=10),labels = T)
plot(x=seq(from=2,to=10),y=gbtm_scores$bic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=10),labels = T)
plot(x=seq(from=2,to=10),y=gbtm_scores$hqic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=10),labels = T)
plot(x=seq(from=2,to=10),y=gbtm_scores$ssbic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=10),labels = T)

fitted3d3<-as.data.frame(m_famAll5d3_genus_clust_plus2$assign)
fitted3d3$Family_ID<-rownames(fitted3d3)
familyDataN<-join(familyData,fitted3d3,by='Family_ID')


fitted3d3<-as.data.frame(m_famAll5d3_genus_clust_plus$assign)
fitted3d3$Family_ID<-rownames(fitted3d3)
familyDataN<-join(familyData_reduced,fitted3d3,by='Family_ID')
names(familyDataN)[228]<-'Traj5_clustScores_plus1'
#familyDataN$Traj5_clustScores_plus1[familyDataN$Traj5_clustScores_plus2==3]<-1


fitted3d3<-as.data.frame(m_famAll7d3_genus_clust_plus2$assign)
fitted3d3[fitted3d3==6]<-5
fitted3d3[fitted3d3==7]<-6
fitted3d3$Family_ID<-rownames(fitted3d3)
familyDataN<-join(familyDataN,fitted3d3,by='Family_ID')
names(familyDataN)[217]<-'Traj7_clustScores_plus2_2timepoints'
names(familyDataN)[216]<-'Traj5_clustScores_plus2_3timepoints'

fitted3d3<-as.data.frame(m_famAll5d3_genus_clust_plus2$assign)
fitted3d3$Family_ID<-rownames(fitted3d3)
familyDataN<-join(familyDataN,fitted3d3,by='Family_ID')
names(familyDataN)[230]<-'Traj5_clustScores_plus2'
familyDataN$Traj5_clustScores_plus2[familyDataN$Traj5_clustScores_plus2==3]<-1

#clusters scores 
m_famAll2d3_genus_clust<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=2, data=familyData_reduced, scaling=0)
m_famAll3d3_genus_clust<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=3, data=familyData_reduced, scaling=0)
m_famAll4d3_genus_clust<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=4, data=familyData_reduced, scaling=0)
m_famAll5d3_genus_clust<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=5, data=familyData_reduced, scaling=0)
m_famAll6d3_genus_clust<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=6, data=familyData_reduced, scaling=0)
m_famAll7d3_genus_clust<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=7, data=familyData_reduced, scaling=0)
m_famAll8d3_genus_clust<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=8, data=familyData_reduced, scaling=0)
m_famAll9d3_genus_clust<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=9, data=familyData_reduced, scaling=0)
m_famAll10d3_genus_clust<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=10, data=familyData_reduced, scaling=0)

gbtm_scores<-as.data.frame(rbind(m_famAll2d3_genus_clust$ic,m_famAll3d3_genus_clust$ic,m_famAll4d3_genus_clust$ic,m_famAll5d3_genus_clust$ic,
           m_famAll6d3_genus_clust$ic,m_famAll7d3_genus_clust$ic,m_famAll8d3_genus_clust$ic,m_famAll9d3_genus_clust$ic,m_famAll10d3_genus_clust$ic
                     ))
rownames(gbtm_scores)<-c('g2','g3','g4','g5','g6','g7','g8','g9','g10')
gbtm_scores$group<-c('g2','g3','g4','g5','g6','g7','g8','g9','g10')

colnames(gbtm_scores)<-c('aic','bic','caic',
                         'ssbic','hqic')
quartz();par(mfrow=c(2,2))
plot(x=seq(from=2,to=10),y=gbtm_scores$caic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=10),labels = T)
plot(x=seq(from=2,to=10),y=gbtm_scores$bic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=10),labels = T)
plot(x=seq(from=2,to=10),y=gbtm_scores$hqic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=10),labels = T)
plot(x=seq(from=2,to=10),y=gbtm_scores$ssbic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=10),labels = T)

fitted3d3<-as.data.frame(m_famAll4d3_genus_clust$assign)
fitted3d3$Family_ID<-rownames(fitted3d3)
familyDataN<-join(familyDataN,fitted3d3,by='Family_ID')
names(familyDataN)[231]<-'Traj4_clustScores'

#clusters 0 or 1
m_famAll2d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=2, data=familyData_reduced, scaling=0)
m_famAll3d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=3, data=familyData_reduced, scaling=0)
m_famAll4d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=4, data=familyData_reduced, scaling=0)
m_famAll5d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=5, data=familyData_reduced, scaling=0)
m_famAll6d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=6, data=familyData_reduced, scaling=0)
m_famAll7d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=7, data=familyData_reduced, scaling=0)
m_famAll8d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=8, data=familyData_reduced, scaling=0)
m_famAll9d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=9, data=familyData_reduced, scaling=0)
m_famAll10d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=10, data=familyData_reduced, scaling=0)
m_famAll11d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=11, data=familyData_reduced, scaling=0)
m_famAll12d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=12, data=familyData_reduced, scaling=0)
m_famAll13d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=13, data=familyData_reduced, scaling=0)
m_famAll14d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=14, data=familyData_reduced, scaling=0)
m_famAll15d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=15, data=familyData_reduced, scaling=0)
m_famAll16d3_genus_clust01<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=16, data=familyData_reduced, scaling=0)


gbtm_scores<-as.data.frame(rbind(m_famAll3d3_genus_clust01$ic,m_famAll4d3_genus_clust01$ic,m_famAll5d3_genus_clust01$ic,
           m_famAll6d3_genus_clust01$ic,m_famAll7d3_genus_clust01$ic,m_famAll8d3_genus_clust01$ic,m_famAll9d3_genus_clust01$ic,m_famAll10d3_genus_clust01$ic,
           m_famAll11d3_genus_clust01$ic,m_famAll12d3_genus_clust01$ic,m_famAll13d3_genus_clust01$ic,m_famAll14d3_genus_clust01$ic,
           m_famAll16d3_genus_clust01$ic
          ))
rownames(gbtm_scores)<-c('g3','g4','g5','g6','g7','g8','g9','g10','g11','g12','g13','g14','g16')
gbtm_scores$group<-c('g3','g4','g5','g6','g7','g8','g9','g10','g11','g12','g13','g14','g16')

colnames(gbtm_scores)<-c('aic','bic','caic',
                         'ssbic','hqic')
quartz();par(mfrow=c(2,2))
plot(x=seq(from=3,to=15),y=gbtm_scores$caic, type="l",col='red')
axis(side = 1, at = seq(from=3,to=15),labels = T)
plot(x=seq(from=3,to=15),y=gbtm_scores$bic, type="l",col='red')
axis(side = 1, at = seq(from=3,to=15),labels = T)
plot(x=seq(from=3,to=15),y=gbtm_scores$hqic, type="l",col='red')
axis(side = 1, at = seq(from=3,to=15),labels = T)
plot(x=seq(from=3,to=15),y=gbtm_scores$ssbic, type="l",col='red')
axis(side = 1, at = seq(from=3,to=15),labels = T)


fitted3d3<-as.data.frame(m_famAll10d3_genus_clust01$assign)
fitted3d3$Family_ID<-rownames(fitted3d3)
familyDataN<-join(familyDataN,fitted3d3,by='Family_ID')
names(familyDataN)[232]<-'Traj10_clust01'

read.table(familyDataN,'final_helmi_meta 10 with Trajectories.txt',sep='\t')

write.table(familyDataN,'final_helmi_meta 10v2 with Trajectories.txt',sep='\t')
read.csv()
# Trajectories for Health related bacteria
spcoaearHealth<-as.data.frame(spcoaPearHealth$sites[,  c("MDS1", "MDS2","MDS3", "MDS4","MDS5", "MDS6")])
names(spcoaearHealth)<-c("MDS1_pearGenusHealth", "MDS2_pearGenusHealth","MDS3_pearGenusHealth", "MDS4_pearGenusHealth","MDS5_pearGenusHealth", "MDS6_pearGenusHealth")
varNames<-c("MDS1_pearGenusHealth", "MDS2_pearGenusHealth","MDS3_pearGenusHealth", "MDS4_pearGenusHealth","MDS5_pearGenusHealth", "MDS6_pearGenusHealth")
familyData2 <- cbind(familyData,spcoaearHealth)

m_famAll2d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=2, data=familyData2, scaling=0)
m_famAll3d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=3, data=familyData2, scaling=0)
m_famAll4d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=4, data=familyData2, scaling=0)
m_famAll5d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=5, data=familyData2, scaling=0)
m_famAll6d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=6, data=familyData2, scaling=0)
m_famAll7d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=7, data=familyData2, scaling=0)
m_famAll8d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=8, data=familyData2, scaling=0)
m_famAll9d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=9, data=familyData2, scaling=0)
m_famAll10d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=10, data=familyData2, scaling=0)
m_famAll11d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=11, data=familyData2, scaling=0)
m_famAll12d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=12, data=familyData2, scaling=0)
m_famAll13d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=13, data=familyData2, scaling=0)
m_famAll14d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=14, data=familyData2, scaling=0)
m_famAll15d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=15, data=familyData2, scaling=0)
m_famAll16d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=16, data=familyData2, scaling=0)
m_famAll17d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=17, data=familyData2, scaling=0)
m_famAll18d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=18, data=familyData2, scaling=0)
m_famAll19d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=19, data=familyData2, scaling=0)
m_famAll20d3_genus_bacteriaHealth<-gbmt(x.names=varNames, unit="Family_ID", time="age_weeks", d=3, ng=20, data=familyData2, scaling=0)

gbtm_scores<-as.data.frame(rbind(m_famAll2d3_genus_bacteriaHealth$ic,m_famAll3d3_genus_bacteriaHealth$ic,m_famAll4d3_genus_bacteriaHealth$ic,m_famAll5d3_genus_bacteriaHealth$ic,
           m_famAll6d3_genus_bacteriaHealth$ic,m_famAll7d3_genus_bacteriaHealth$ic,m_famAll8d3_genus_bacteriaHealth$ic,m_famAll9d3_genus_bacteriaHealth$ic,m_famAll10d3_genus_bacteriaHealth$ic,
           m_famAll11d3_genus_bacteriaHealth$ic,m_famAll12d3_genus_bacteriaHealth$ic,m_famAll13d3_genus_bacteriaHealth$ic,m_famAll14d3_genus_bacteriaHealth$ic,m_famAll15d3_genus_bacteriaHealth$ic,
           m_famAll16d3_genus_bacteriaHealth$ic,m_famAll17d3_genus_bacteriaHealth$ic,m_famAll18d3_genus_bacteriaHealth$ic,m_famAll19d3_genus_bacteriaHealth$ic,m_famAll20d3_genus_bacteriaHealth$ic
           ))
rownames(gbtm_scores)<-c('g2','g3','g4','g5','g6','g7','g8','g9','g10',
                         'g11','g12','g13','g14','g15','g16','g17','g18','g19','g20')
gbtm_scores$group<-c('g2','g3','g4','g5','g6','g7','g8','g9','g10',
                         'g11','g12','g13','g14','g15','g16','g17','g18','g19','g20')

colnames(gbtm_scores)<-c('aic','bic','caic',
                         'ssbic','hqic')
quartz();par(mfrow=c(2,2))
plot(x=seq(from=2,to=20),y=gbtm_scores$caic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=20),labels = T)
plot(x=seq(from=2,to=20),y=gbtm_scores$bic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=20),labels = T)
plot(x=seq(from=2,to=20),y=gbtm_scores$hqic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=20),labels = T)
plot(x=seq(from=2,to=20),y=gbtm_scores$ssbic, type="l",col='red')
axis(side = 1, at = seq(from=2,to=20),labels = T)

fitted6d3_healthpeargenus<-as.data.frame(m_famAll6d3_genus_bacteriaHealth$assign)
write.table(fitted6d3_healthpeargenus,'Trajectories_genusHealthBacteria_6groups.txt',sep='\t')

fitted6d3<-as.data.frame(m_famAll6d3_genus_bacteriaHealth$assign)
fitted6d3$Family_ID<-rownames(fitted6d3)
familyDataN<-join(familyData,fitted6d3,by='Family_ID')
names(familyDataN)[119]<-'group'


