## BIOL 314 -- Phylogeny & Functional Traits birds Meta ## 

# Contents: 
# 1) Load in phylogeny
# 2) Load in traits from Avonet 
# 3) Subset both to just the species we observed in Meta


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(xlsx)
library(phytools)
library(ggtree)
library(conflicted)
conflicts_prefer(purrr::map)

# Bring in data
Meta_birds <- read.xlsx("Derived/Meta_314.xlsx", sheetIndex = 1)
spp_obs <- names(Meta_birds)[-1]
# Traits from Avonet
path <- "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/PhD/Analysis/Colombia-SCR-Rd3/Derived/"
Traits <- read.xlsx(paste0(path, "Excels/Avo_traits_final.xlsx"), 
                    sheetIndex = "Traits")
load("Rdata/Biol314_inputs_02.26.25.Rdata")

# Phylogenetic tree -------------------------------------------------------
## Bring in single phylogenetic tree (BirdTree taxonomy) with just the species from my project 
Meta_tree <- read.tree(paste0(path, "Single_tree.tre"))

# Prune tree --------------------------------------------------------------
## Prune(?) the tree down to just the species we have 
phylo.obs <- ape::keep.tip(Meta_tree, spp_obs, trim.internal = TRUE)

# Export tree -------------------------------------------------------------

write.tree(phylo.obs, "Derived/Meta_tree314.tre")

# Taxonomy ----------------------------------------------------------------
## Generally I have been using the Ayerbe (2018) taxonomy , but the phylogeny was downloaded from Bird Tree
## Bring in data on orders or families 
Spp_join_bt <- Tax_df3 %>% as_tibble() %>%
  distinct(Species_ayerbe, Species_bt, order_gbif, family_gbif) 

Tax_join_meta <- Spp_join_bt %>% 
  select(-Species_ayerbe) %>% 
  distinct() %>% 
  mutate(
    Species_bt = str_replace(Species_bt, " ", "_"), 
    order_gbif = if_else(Species_bt == "Ortalis_guttata", "Galliformes", order_gbif),
    family_gbif = if_else(Species_bt == "Ortalis_guttata", "Galliformes", family_gbif)
    )

# Functional traits -------------------------------------------------------
Morphology_ft <- Spp_join_bt %>% select(-c(order_gbif, family_gbif)) %>% 
  left_join(Traits) %>% 
  mutate(Species_bt = str_replace(Species_bt, " ", "_")) %>% 
  right_join(tibble(spp_obs), by = join_by("Species_bt" == "spp_obs")) %>% 
  slice_head(n = 1, by = Species_bt) %>% 
  distinct(Species_bt, pick(c(7:17)))

# Export functional traits ------------------------------------------------
Morphology_ft %>% as.data.frame() %>%
  write.xlsx("Derived/Morphology_ft.xlsx", row.names = FALSE, showNA = FALSE)
