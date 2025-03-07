#BIOL/CONS 314
#Lecture XX: caclulating phylogenetic community structure
#Jonathan Davies 06/03/2025

#install and load the libraries
library(ape)
library(picante)

#Example datasest in the Picante r package for illustrating calculations of
# phylogenetic community structure (MPD and MNTD)

#laod the data
data(phylocom)

#extract the tree
tree.phylocom<-phylocom$phylo

#exract the community data matrix
sp<-phylocom$sample

#select the first site (row)
site1<-sp[1,]
#gete the names of the species NOT in the site (0 abundance)
sp2<-site1[site1==0]

#pune tht species for the tree NOT in the site
tree<-drop.tip(tree.phylocom, names(sp2))

#calculate MPD and MNTD
plot(tree)
mat<-cophenetic(tree)
mat
dist.mat<-as.dist(mat)
dist.mat
my.mpd<-mean(dist.mat)

mat[mat==0]<-NA
mat.min<-apply(mat, 2, min, na.rm = T)
my.mntd<-mean(mat.min)


##############################################################################
##############################################################################
#Now lets try for the real forst plot data!
##############################################################################

#read in the community matrix
x<-read.table("forest.data.txt", header = T)

#read in the tree
tree<-read.tree("forest.tre")
plot(tree)

#extract data for site in column 4 (site LJ) 
y<-(x[,4])
names(y)<-x[,1]

#get list of taxa NOT in the site (0 abundance)
z<-y[y==0]

#prune the tree
y.tree<-drop.tip(tree, names(z))
plot(y.tree)

#calculate MPD and MNTD
mat<-cophenetic(y.tree)
dist.mat<-as.dist(mat)
my.mpd<-mean(dist.mat)

mat[mat==0]<-NA
mat.min<-apply(mat, 2, min, na.rm = T)
my.mntd<-mean(mat.min)


##############################################################################
##############################################################################
#Now lets calculate the Standard Effect Sizes (SES), also referred to as z-scores, again starting the phylocom data  
##############################################################################

library(picante)

data(phylocom)
comm.tree<-phylocom$phylo
comm.data<-phylocom$sample

#calculate SES values for MPD and MNTD
test.mpd<-ses.mpd(comm.data, cophenetic(comm.tree),null.model="taxa.labels", iterations = 1000)
test.mpd
test.mntd<-ses.mntd(comm.data, cophenetic(comm.tree),null.model="taxa.labels", iterations = 1000)
test.mntd
#generate a box plot for comparison
boxplot(cbind(test.mpd$mpd.obs.z,test.mntd$mntd.obs.z), names=c("MPD", "MNTD"))

#calculate SES values for MPD and MNTD weigting by abundance
test.abund.mpd<-ses.mpd(comm.data, cophenetic(comm.tree),null.model="taxa.labels", abundance.weighted=TRUE, iterations = 1000)
test.abund.mpd
test.abund.mntd<-ses.mntd(comm.data, cophenetic(comm.tree),null.model="taxa.labels", abundance.weighted=TRUE, iterations = 1000)
test.abund.mntd
#generate a box plot for comparison
boxplot(cbind(test.abund.mpd$mpd.obs.z,test.abund.mntd$mntd.obs.z), names=c("weighted MPD", "weighted MNTD"))

#More boxplots:
boxplot(cbind(test.mpd$mpd.obs.z,test.abund.mpd$mpd.obs.z), names=c("MPD", "weighted MPD"))
boxplot(cbind(test.mntd$mntd.obs.z,test.abund.mntd$mntd.obs.z), names=c("MNTD", "weighted MNTD"))


##############################################################################
##############################################################################
#Your task is to calculate SES values for the Tiawanese forest plots
##############################################################################

library(picante)
tree<-read.tree("forest.tre")

x<-read.table("forest.data.txt", header = T, row.names = 1)
#don't forget to transpose your matrix
x.mat<-t(x)

#The rest should be easy!

###################################################
###################################################
###################################################
###################################################

#calculate SES values for MPD and MNTD
test.mpd<-ses.mpd(x.mat, cophenetic(tree),null.model="taxa.labels", iterations = 1000)
test.mpd
test.mntd<-ses.mntd(x.mat, cophenetic(tree),null.model="taxa.labels", iterations = 1000)
test.mntd
#generate a box plot for comparison
boxplot(cbind(test.mpd$mpd.obs.z,test.mntd$mntd.obs.z), names=c("MPD", "MNTD"))

#calculate SES values for MPD and MNTD weigting by abundance
test.abund.mpd<-ses.mpd(x.mat, cophenetic(tree),null.model="taxa.labels", abundance.weighted=TRUE, iterations = 1000)
test.abund.mpd
test.abund.mntd<-ses.mntd(x.mat, cophenetic(tree),null.model="taxa.labels", abundance.weighted=TRUE, iterations = 1000)
test.abund.mntd
#generate a box plot for comparison
boxplot(cbind(test.abund.mpd$mpd.obs.z,test.abund.mntd$mntd.obs.z), names=c("weighted MPD", "weighted MNTD"))

#More boxplots:
boxplot(cbind(test.mpd$mpd.obs.z,test.abund.mpd$mpd.obs.z), names=c("MPD", "weighted MPD"))
boxplot(cbind(test.mntd$mntd.obs.z,test.abund.mntd$mntd.obs.z), names=c("MNTD", "weighted MNTD"))

