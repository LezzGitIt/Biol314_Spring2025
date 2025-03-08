#BIOL/CONS 314
#Lecture XX: Caclulating phylogenetic community structure
# Aaron Skinner 03/08/2025

# Load libraries ----------------------------------------------------------
library(ape)
library(picante)

# Phylogenetic community structure ----------------------------------------
# >Example -----------------------------------------------------------------
#Example datasest in the Picante r package for illustrating calculations of phylogenetic community structure (MPD and MNTD)

#laod the data
data(phylocom)

#extract the tree
tree.phylocom<-phylocom$phylo

#exract the community data matrix
sp<-phylocom$sample

#select the first site (row)
site1<-sp[1,]
#get the names of the species NOT in the site (0 abundance)
sp2<-site1[site1==0]

#prune the species for the tree NOT in the site
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


# >Bird data--------------------------------------------------------------------
# Bring in community data -- currently a dataframe
Meta_birds <- read.xlsx("Derived/Meta_314.xlsx", sheetIndex = 1) %>% tibble()
head(Meta_birds)

# Turn to matrix and transpose
Birds_mat <- Meta_birds %>% column_to_rownames("Id_muestreo") %>% 
  as.matrix()
Birds_t <- Birds_mat %>% t()

## Phylogenetic tree
# Set path to pull in inputs from primary R project
path <- "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/PhD/Analysis/Colombia-SCR-Rd3/Derived/"

## Bring in single phylogenetic tree (BirdTree taxonomy) with just the species from my project 
Birds_tree <- read.tree(paste0(path, "Single_tree.tre"))

# Visualize phylogeny
plot(Birds_tree)

# Extract data for site in column 4 (`C-MB-M-C_03`)
y <- Birds_t[,4]
names(y)<- rownames(Birds_t)

# Get list of taxa NOT in the site (0 abundance)
z<-y[y==0]

# Prune the tree
y.tree<-drop.tip(Meta_tree, names(z))
plot(y.tree)

## Calculate MPD and MNTD
mat<-cophenetic(y.tree)
dist.mat<-as.dist(mat)
my.mpd<-mean(dist.mat)

mat[mat==0]<-NA
mat.min<-apply(mat, 2, min, na.rm = T)
my.mntd<-mean(mat.min)


# Standard Effect Sizes (SES) ---------------------------------------------
#Now lets calculate the Standard Effect Sizes (SES), also referred to as z-scores

# >Example ---------------------------------------------------------------- 
# Start with the phylocom data  

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


# >Bird data --------------------------------------------------------------

#calculate SES values for MPD and MNTD
test.mpd<-ses.mpd(Birds_mat, cophenetic(Birds_tree),null.model="taxa.labels", iterations = 1000)
test.mpd
test.mntd<-ses.mntd(Birds_mat, cophenetic(Birds_tree),null.model="taxa.labels", iterations = 1000)
test.mntd
#generate a box plot for comparison
boxplot(cbind(test.mpd$mpd.obs.z,test.mntd$mntd.obs.z), names=c("MPD", "MNTD"))

#calculate SES values for MPD and MNTD weigting by abundance
test.abund.mpd<-ses.mpd(Birds_mat, cophenetic(Birds_tree),null.model="taxa.labels", abundance.weighted=TRUE, iterations = 1000)
test.abund.mpd
test.abund.mntd<-ses.mntd(Birds_mat, cophenetic(Birds_tree),null.model="taxa.labels", abundance.weighted=TRUE, iterations = 1000)
test.abund.mntd
#generate a box plot for comparison
boxplot(cbind(test.abund.mpd$mpd.obs.z,test.abund.mntd$mntd.obs.z), names=c("weighted MPD", "weighted MNTD"))

#More boxplots:
boxplot(cbind(test.mpd$mpd.obs.z,test.abund.mpd$mpd.obs.z), names=c("MPD", "weighted MPD"))
boxplot(cbind(test.mntd$mntd.obs.z,test.abund.mntd$mntd.obs.z), names=c("MNTD", "weighted MNTD"))