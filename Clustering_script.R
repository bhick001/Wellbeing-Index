library(RColorBrewer)
library(cluster)
library(vegan)
library(scales)

cols = brewer.pal(9, "GnBu")
pal = colorRampPalette(cols)

m$order = findInterval(m$age, sort(m$age))

#parents data
md <- read.delim("meta.txt")
gd <- read.delim("organised_combined_genus_table.txt")
g<-read.delim('genus.txt')
f<-read.delim('family.txt')
mparents <- md[md$age_weeks>100&md$ReadCount>50000&!is.na(md$Sample_type),]
gparents <- gd[md$age_weeks>100&md$ReadCount>50000&!is.na(md$Sample_type),]

gall <- rbind(g,gparents)

#plot pcoa with clusters
palette((c("#4DAF4A","#377EB8","#FC8D62","#E41A1C","gray")))

quartz();par(mfrow=c(2,1))
barplot(prop.table(table(m$pam_cl[m$healthUse==1],m$age_weeks[m$healthUse==1]),margin=2),col=1:10,main = "PEARPAM 5- Reference Children",)
barplot(prop.table(table(m$pam_cl[m$birthmode_RecievingIAP=='C-section'],m$age_weeks[m$birthmode_RecievingIAP=='C-section']),margin=2),col=1:10,main = "PEARPAM 5- C-section Children",)


#clustering data -----------------------------------------------------------
pearson_genus <- peardist(log((gall+1)/rowSums(gall))) 
pcoa_pearson_genus <- vegan::capscale(pearson_genus ~ 1)
spcoa<-summary(pcoa_pearson_genus)
pcm_all<-data.frame(spcoa$sites[,c('MDS1','MDS2')],age=c(m$age,rep(10000,nrow(f))))
pcm_all$kmeans <- kmeans(x=pearson_genus ,centers = 5)$cluster


