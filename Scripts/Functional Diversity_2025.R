#BIOL/CONS 314
#Lecture 19: Functional Diversity
#Jonathan Davies 26/03/2025

# Background --------------------------------------------------------------

# The Functional richness metric will increase with species richness , so to account for that it's important to create random communities to compare to your observed functional richness 

# Libraries ---------------------------------------------------------------
library(fundiversity)
library(vegan)
library(venneuler)
ggplot2::theme_set(cowplot::theme_cowplot())


# Bring in data -----------------------------------------------------------
# Functional traits - Jonathan turned into a csv
bird_traits<-read.csv("Data/Morphology_ft.csv", row.names = 1)

# Bring in spp x site df & turn to matrix
Meta_birds <- read.xlsx("Derived/Meta_314.xlsx", sheetIndex = 1) 
Birds_mat <- Meta_birds %>% column_to_rownames("Id_muestreo") %>% 
  as.matrix()

# bird.com = habitat x species matrix with 1s and 0s
#bird.com<-read.csv("Data/bird_comm.txt")

# fd_fric() function produces NA at sites w/ fewer species than traits
Traits_red <- bird_traits %>% select(where(is.numeric)) %>% 
  select(Beak.Length_Culmen, Tarsus.Length, Wing.Length, Tarsus.Length, Tail.Length, Mass)

Fric_tbl <- fd_fric(traits = Traits_red, sp_com = Birds_mat) %>% 
  rename_with(~ c("Id_muestreo", "F_ric")) %>% 
  tibble()

Fric_tbl %>% filter(is.na(F_ric)) #43

# Subset
birds<-bird_traits[1:5,1:4]
head(birds)

# Convert to distance matrix
my.dist<-dist(birds)
my.dist
sum(my.dist)

# Turn distance matrix into a dendrogram 
cl<-hclust(my.dist)
plot(cl)
treeheight(cl)

# Functional richness -----------------------------------------------------
## For each row (site) calculate the cumulative phylogenetic diversity
# The function(...) collects all columns for each row into a list.

#Functional Richness (FRic) represents the total amount of functional space filed by a community in a dataset
fd_fric(birds)

birds<-bird_traits[,1:6]

FD<-fd_fric(birds, bird.com)
FD

FD$log.FD<-log(FD$FRic)
FD

# Standardization rescales to divide the functional richness in each habitat by the TOTAL functional richness across all habitats 
#Because the convex hull volume depends on the number and the units of the traits used, it is difficult to compare across datasets, that is why it has been suggested to standardize its value by the total
#volume comprising all species in the dataset
FD.stnd<-fd_fric(birds, bird.com, stand = TRUE)
FD.stnd

#scl.birds<-scale(birds)
#FD.scl<-fd_fric(scl.birds, bird.com)

plot(log(FD$FRic), log(FD.stnd$FRic))

## Compare overlap between habitat types 
#Sometimes you?re interested in the shared functional volumes between pairs of sites
#lets look at functional overlap between the first two sites
#scl.birds<-scale(birds)


# Functional overlap ------------------------------------------------------

## Grass human
#fd_fric_intersect(birds, bird.com[2:3,], stand = TRUE)#Forest and Grassland
Grass_Human<-fd_fric_intersect(birds, bird.com[c(3,4),], stand = TRUE)#Grassland and Human modified
plot(venneuler(c(Human_modified = Grass_Human[3,3],
                 Grassland = Grass_Human[2,3],
                 "Human_modified&Grassland" = Grass_Human[1,3]
)))

# CAUSES ERROR
For_Human<-fd_fric_intersect(birds, bird.com[c(2,4),], stand = TRUE)#Forest and Human modified

#Jac<-Grass_Human[1,3]/(Grass_Human[2,3]+Grass_Human[3, 3]-Grass_Human[1,3])
#Jac
#Jac<-For_Human[1,3]/(For_Human[2,3]+For_Human[3, 3]-For_Human[1,3])
#Jac

#Now calculate functional richness overlap across all pairs of sites
#fd_fric_intersect(birds, bird.com)
                 
plot(venneuler(c(Human_modified = For_Human[3,3],
                 Forest = For_Human[2,3],
"Human_modified&Forest" = For_Human[1,3]
                 )))



# CALC PER SITE -----------------------------------------------------------
# CODE PASTED FROM HISTORY FOR PD , BUT COULD BE ADAPTED TO CALCULATE FUNCTIONAL RICHNESS FOR EACH SITE
Spp_site <- pmap(Forest, function(..., row = list(...)) {
  row_vec <- unlist(row)
  names(row_vec)[row_vec > 0]
})

Pd_site <- map_dbl(Spp_site, \(Spp_vec){
  Tree_site <- keep.tip(Birds_tree, Spp_vec)
  sum(Tree_site$edge.length)
})

Pd_my <- tibble(Id_muestreo = names(Pd_site), Pd_my = Pd_site, Spp_rich) %>%
  left_join(co_vars)
                 

# Non-continuous traits? --------------------------------------------------

#Do not panic. You can still compute the above-mentioned functional diversity indices.
#However, as all indices need continuous descriptors for all considered species, 
#you need to transform the non-continuous trait data into a continuous form. 
#The general idea is to obtain from the trait table a table of quantitative descriptions 
#by defining specific dissimilarity and projecting species dissimilarities onto quantitative 
#space using Principal Coordinates Analysis (PCoA). The framework is fully described in Maire et al. (2015).

#To compute dissimilarity with non-continuous traits you can user Gower?s distance (Gower 1971) or its 
#following adaptations (Pavoine et al. 2009; Podani 1999). You can use the following functions:
#cluster::daisy(), FD::gowdis(), ade4::dist.ktab(), or vegan::vegdist().
#Then you can project these dissimilarities with Principal Coordinates using ape::pcoa() for example.
#You can then select the first dimensions that explains the most variance and use theses as the input ?traits?
#to compute functional diversity indices.
#https://cran.r-project.org/web/packages/fundiversity/vignettes/fundiversity.html