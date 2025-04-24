#BIOL/CONS 314
## Caclulating phylogenetic community structure
# Aaron Skinner 03/08/2025

# Load libraries ----------------------------------------------------------
library(ape)
library(picante)
library(tidyverse)
library(xlsx)
library(janitor)


# To do -------------------------------------------------------------------
 # Compare overall (and controlling for species richness) PD (millions of years) captured in each habitat type. Control for species richness? 

# Background -------------------------------------------------------------------

# Phylogenetic under-dispersion, or clustering (less scattered than random), is usually caused by environmental filtering 
# Phylogenetic over-dispersion (more scattered than random), is usually caused by competition  

## MPD vs MNTD 
# What are the processes causing phylogenetic structure? And where have these processes occured (ie tree tips or ancestral)?
# Mean phylogenetic distance (MPD) is likely better at detecting the process of filtering, & is able to detect processes happening deeper in the tree. For example, an ancestral trait of 'drought resistance' (whatever that looks like) that is shared among a clade is a 'deeper process'
# MNTD likely better at detecting (places more emphasis on) the process of competition, as closely related species that are in competition often speciated recently

# Load data--------------------------------------------------------------------
# Bring in community data -- currently a dataframe
Meta_birds <- read.xlsx("Derived/Meta_314.xlsx", sheetIndex = 1) %>% tibble()

# Turn to matrix and transpose
Birds_mat <- Meta_birds %>% column_to_rownames("Id_muestreo") %>% 
  as.matrix()
Birds_t <- Birds_mat %>% t()

# Bring in covariates, one row per site
co_vars <- read.xlsx("Derived/Site_covs314.xlsx", sheetIndex = 1)

# Functional traits
Traits <- read.xlsx("Derived/Morphology_ft.xlsx", sheetIndex = 1)

## Phylogenetic tree
# Bring in single phylogenetic tree (BirdTree taxonomy) with just the species from my project 
Birds_tree <- read.tree("Derived/Meta_tree314.tre")

# Standard Effect Sizes (SES) ---------------------------------------------
# Now lets calculate the Standard Effect Sizes (SES). Essentially, we are comparing the observed phylogenetic structure to communities generated at random

# Calculate SES values for MPD and MNTD - column of interest: mpd.obs.z
mpd<-ses.mpd(Birds_mat, cophenetic(Birds_tree),null.model="sample.pool", iterations = 100)
# Box plot
boxplot(mpd$mpd.obs.z, names=c("MPD"))

## Weighting by abundance
abund.mpd<-ses.mpd(Birds_mat, cophenetic(Birds_tree),null.model="sample.pool", abundance.weighted=TRUE, iterations = 10)
boxplot(cbind(mpd$mpd.obs.z, abund.mpd$mpd.obs.z), names=c("MPD", "weighted MPD"))

## Join MPD metric to site covariates 
# Create function to subset, join, and rename
merge_mpd_covs <- function(df){
  df %>% rownames_to_column("Id_muestreo") %>% 
    select(Id_muestreo, mpd.obs.z) %>% 
    left_join(co_vars) %>% 
    rename(MPD = mpd.obs.z) %>% 
    tibble()
}
co_vars_all <- merge_mpd_covs(mpd) 

## Compare MPD by habitat types
# Riparian forest is clustered, suggesting there may be filtering for riparian specialists 
boxplot(MPD~Habitat_cons, data = co_vars_all)
aov.test<-aov(MPD~Habitat_cons, data = co_vars_all)
TukeyHSD(aov.test)

# Forest only -------------------------------------------------------------
Traits %>% tabyl(Habitat)
Forest_spp <- Traits %>% 
  filter(Habitat %in% c("Forest", "Woodland", "Riverine")) %>% 
  pull(Species_bt)
TF <- colnames(Birds_mat) %in% Forest_spp

Forest <- Birds_mat[,TF]
dim(Forest) #177 species retained

mpd<-ses.mpd(Forest, cophenetic(Birds_tree),null.model="sample.pool", iterations = 1000)
abund.mpd<-ses.mpd(Forest, cophenetic(Birds_tree),null.model="sample.pool", abundance.weighted=TRUE, iterations = 100)
boxplot(cbind(mpd$mpd.obs.z, abund.mpd$mpd.obs.z), names=c("MPD", "weighted MPD"))

co_vars_forest <- merge_mpd_covs(mpd)

boxplot(MPD~Habitat_cons, data = co_vars_forest, main = "Forest birds", xlab = "Habitat")
aov.test<-aov(MPD~Habitat_cons, data = co_vars_forest)
#Tukey's honest significance test
TukeyHSD(aov.test) # plot()

# >Abiotic factors ---------------------------------------------------------
plot(co_vars_forest$MPD~ co_vars_forest$Elev)
model<-lm(co_vars_forest$MPD~co_vars_forest$Elev)
summary(model) # Negative

plot(co_vars_forest$MPD~co_vars_forest$Avg_temp)
model<-lm(co_vars_forest$MPD~co_vars_forest$Avg_temp)
summary(model) # Positive

## Way oo correlated to include in a single model -- coefficients can't be trusted
cor(co_vars_forest$Avg_temp, co_vars_forest$Elev)
model<-lm(MPD~Avg_temp+Elev, data = co_vars_forest)
summary(model) # Elevation has switched signs 
car::vif(model)

# Forest v Riparian -------------------------------------------------------
ID_Bosque<-subset(co_vars, Habitat_cons == "Bosque")[,1]
ID_Bosque_r<-subset(co_vars, Habitat_cons == "Bosque ripario")[,1]
ID_Bosque_comb<-c(ID_Bosque,ID_Bosque_r) 
Birds_Bosque<-Birds_mat[rownames(Birds_mat) %in% ID_Bosque_comb,]
dim(Birds_Bosque) #115 sites

Bosque.mpd<-ses.mpd(Birds_Bosque, cophenetic(Birds_tree),null.model="sample.pool", abundance.weighted=F, iterations = 100)

sites<-rownames(Bosque.mpd)
habitat<-NULL
for (i in 1:length(sites)){
for (j in 1:length(co_vars$Habitat_cons)){
if (sites[i]==co_vars$Id_muestreo[j]){
habitat[i]<-co_vars$Habitat_cons[j]
break
}
}
}

boxplot(Bosque.mpd$mpd.obs.z~habitat,  main = "Bosque")
t.test(Bosque.mpd$mpd.obs.z~habitat)

# Pasture only ------------------------------------------------------------
Pasture_spp <- Traits %>% 
  filter(Habitat %in% c("Grassland","Human modified","Shrubland","Wetland")) %>%
  pull(Species_bt)
TF <- colnames(Birds_mat) %in% Pasture_spp
Pasture <- Birds_mat[,TF]
dim(Pasture) #60 species retained

Pastizales.mpd<-ses.mpd(Pasture, cophenetic(Birds_tree),null.model="sample.pool", abundance.weighted=F, iterations = 100)

co_vars_pasture <- merge_mpd_covs(Pastizales.mpd)

boxplot(MPD~Habitat_cons, data = co_vars_pasture )
aov.test<-aov(MPD~Habitat_cons, data = co_vars_pasture )
#Tukey's honest significance test
TukeyHSD(aov.test) # plot()