#BIOL/CONS 314
## Caclulating phylogenetic community structure
# Aaron Skinner 03/08/2025

# Load libraries ----------------------------------------------------------
library(ape)
library(picante)
library(tidyverse)
library(xlsx)
library(janitor)
library(ggpubr)
ggplot2::theme_set(cowplot::theme_cowplot())

# To do -------------------------------------------------------------------

# Background --------------------------------------------------------------

## The key assumptions underlying all measures of phylogenetic diversity are that (functional) traits structuring communities are evolutionary conserved and match to a model of evolution that is well approximated by Brownian motion (i.e. a gradual mode of evolution)

## Under & Overdispersion 
# Phylogenetic under-dispersion, or clustering (less scattered than random), is usually caused by environmental filtering 
# Phylogenetic over-dispersion (more scattered than random), is usually caused by competition  

## MPD vs MNTD 
# What are the processes causing phylogenetic structure? And where have these processes occured (ie tree tips or ancestral)?
# Mean phylogenetic distance (MPD) is likely better at detecting the process of filtering, & is able to detect processes happening deeper in the tree. For example, an ancestral trait of 'drought resistance' (whatever that looks like) that is shared among a clade is a 'deeper process'
# MNTD likely better at detecting (places more emphasis on) the process of competition, as closely related species that are in competition often speciated recently
# NOTE:: B/c these are mean measures the species richness shouldn't matter 

## What about the relationship between FD & PD? 
# Depends on several things, including the taxonomic scale in question. Within closely related species, PD & FD are fairly highly correlated, but within disparate groups convergent evolution may wash out the signal between PD & FD. I.e., Very evolutionarily distinct groups may have similar FD 

# Load data--------------------------------------------------------------------
# Bring in community data -- currently a dataframe
Meta_birds <- read.xlsx("Derived/Meta_314.xlsx", sheetIndex = 1) %>% tibble()

# Turn to matrix and transpose
Birds_mat <- Meta_birds %>% column_to_rownames("Id_muestreo") %>% 
  as.matrix()
Birds_t <- Birds_mat %>% t()

# Bring in covariates, one row per site
co_vars <- read.xlsx("Derived/Site_covs314.xlsx", sheetIndex = 1) %>%
  rename(Habitat = Habitat_cons) %>%
  mutate(Habitat = factor(Habitat, levels = c("Bosque ripario", "Bosque", "Ssp", "Pastizales")))

# Functional traits
Traits <- read.xlsx("Derived/Morphology_ft.xlsx", sheetIndex = 1) 
Traits_red <- Traits %>% 
  column_to_rownames("Species_bt") %>% #select(where(is.numeric)) %>% 
  select(Beak.Length_Culmen, Tarsus.Length, Wing.Length, Tarsus.Length, Tail.Length, Mass) %>% 
  mutate(across(everything(.), as.numeric)) 

## Phylogenetic tree
# Bring in single phylogenetic tree (BirdTree taxonomy) with just the species from my project 
Birds_tree <- read.tree("Derived/Meta_tree314.tre")

# All birds ---------------------------------------------------------------
# >Standard Effect Sizes (SES) ---------------------------------------------
# Use the ses.##() functions to generate both scaled & unscaled values of PD & MPD
# Standard Effect Sizes (SES) compare the observed phylogenetic structure to communities generated at random

## Phylogenetic diversity (PD)
Pd_ses_tbl <- ses.pd(Birds_mat, Birds_tree, null.model="sample.pool", iterations = 10)

## Key function -- Join PD metric- (pd or mpd) to site covariates 
# Create function to subset, join, and rename
merge_pd_covs <- function(df, var){
  var_int <- str_split_i({{ var }}, "\\.", i = 1) # Variable of interest
  col_name_obs <- paste0(var_int, ".obs") # Observed pd value
  col_name_z <- paste0(var_int, ".z") # z- standardized observed pd value
  
  df %>% rownames_to_column("Id_muestreo") %>% 
    select(Id_muestreo, ntaxa, col_name_obs, {{ var }}) %>% 
    left_join(co_vars) %>% 
    rename(!!col_name_z := {{ var }}) %>% 
    tibble()
}

Pd_tbl <- merge_pd_covs(Pd_ses_tbl, var = "pd.obs.z") 

## Mean Phylogenetic Diversity (MPD)
# Calculate SES values for MPD and MNTD - column of interest: mpd.obs.z
Mpd_ses_tbl <- ses.mpd(samp = Birds_mat, dis = cophenetic(Birds_tree), null.model="sample.pool", iterations = 100)
# Box plot
boxplot(Mpd_ses_tbl$mpd.obs.z, names=c("MPD"))

## Weighting by abundance 
# NOTE:: Things don't change much by including abundance
abund.mpd<-ses.mpd(Birds_mat, cophenetic(Birds_tree),null.model="sample.pool", abundance.weighted=TRUE, iterations = 10)
boxplot(cbind(Mpd_ses_tbl$mpd.obs.z, abund.mpd$mpd.obs.z), names=c("MPD", "weighted MPD"))

# Join PD & MPD
Mpd_tbl <- merge_pd_covs(Mpd_ses_tbl, var = "mpd.obs.z")
Pd_mpd <- Pd_tbl %>% left_join(Mpd_tbl) %>% 
  relocate(ntaxa:pd.z, .after = Habitat) %>% 
  rename(spp_rich = ntaxa)

# >Functional richness ----------------------------------------------------
# fd_fric() function produces NA at sites w/ fewer species than traits
Fric_tbl <- fd_fric(traits = Traits_red, sp_com = Birds_mat) %>% 
  rename_with(~ c("Id_muestreo", "F_ric")) %>% 
  tibble()

# Join
Fd_pd_mpd <- Pd_mpd %>% left_join(Fric_tbl)

# >Visualize --------------------------------------------------------------
# There is a non-linear (saturating) positive relationship between species richness and phylogenetic diversity 
Fd_pd_mpd %>% ggplot(aes(x = spp_rich, y = pd.obs)) + 
  geom_point(aes(color = Habitat)) + 
  geom_smooth()

# Custom function pivot_scale
pivot_scale <- function(df, var1, var2 = NULL, var3 = NULL, scale = FALSE){
  if(scale){
    df  <- df  %>%
      mutate(across(c({{ var1 }}, {{ var2 }}, {{ var3 }}), ~ scale(.)))
  }
  df %>% pivot_longer(cols = c({{ var1 }}, {{ var2 }}, {{ var3 }}), 
                           names_to = "Metric", values_to = "Value")
}

# Prep dataframes for plotting via pivot
df1 <- pivot_scale(df = Pd_mpd, var1 = mpd.z, var2 = pd.z, scale = FALSE)
df2 <- pivot_scale(df = Fd_pd_mpd, var1 = spp_rich, 
                   var2 = pd.obs, var3 = F_ric, 
                   scale = TRUE)

# Custom function to plot PD
plot_pd <- function(df, hline = FALSE, label = NULL, title = NULL){
  dodge <- position_dodge(width = 0.75)
  p <- ggplot(df, aes(x = Habitat, y = Value, fill = Metric, 
                      group = interaction(Habitat, Metric))) 
  if (hline) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  }
  p <- p + geom_boxplot(outlier.shape = NA, position = dodge) +
    geom_point(alpha = .15, position = position_jitterdodge(jitter.width = .3))
  if(!is.null(title)) {
    p <- p + labs(y = NULL, title = title) 
  }
  if(!is.null(label)) {
    p <- p + scale_fill_discrete(labels = label)
  }
  return(p)
}

## Plot
# a) Important differences in overall amount of phylogenetic diversity captured, but this is primarily driven by differences in species richness
# NOTE: Could log both SR & PD , and try to have 2 y axis lables, one on both side 
# b) Riparian forest is clustered, suggesting there may be filtering for riparian specialists 
Mpd_pd_p <- plot_pd(df1, hline = TRUE, label = c("MPD", "PD"), title = "Z-standardized")
Pd_sr_p <- plot_pd(df2, label = c("Functional \nrichness", "PD", "Species \nrichness"), title = "Meta - All species")
ggarrange(Pd_sr_p, Mpd_pd_p)

#NOTE:: Several NAs
df1 %>% filter(Metric %in% c("mpd.z", "pd.z") & is.na(Value)) 
Fric_tbl %>% filter(is.na(F_ric)) #43 NAs

# >Models -----------------------------------------------------------------
## Compare MPD by habitat types
# Run an ANOVA & then a Tukey HSD test
aov_tukey <- function(data, x, y){
  form <- reformulate(x, response = y)
  aov.test <- aov(form, data = data)   
  TukeyHSD(aov.test)                   
}
aov_tukey(data = Pd_mpd, x = "Habitat", y = "mpd.z")
aov_tukey(data = Pd_mpd, x = "Habitat", y = "pd.z")
aov_tukey(data = Pd_mpd, x = "Habitat", y = "pd.obs")

# Forest only -------------------------------------------------------------
Traits %>% tabyl(Habitat)
Forest_spp <- Traits %>% 
  # Riverine includes kingfishers, hoazin, drab water tyrant, swallow
  filter(Habitat %in% c("Forest", "Woodland")) %>% #"Riverine"
  pull(Species_bt)
TF <- colnames(Birds_mat) %in% Forest_spp

Forest_mat <- Birds_mat[,TF]
dim(Forest_mat) #177 species retained

# >Standard Effect Sizes (SES) ---------------------------------------------
## Phylogenetic diversity (PD)
Forest_ses_pd <- ses.pd(Forest_mat, Birds_tree, null.model="sample.pool", iterations = 10)

Forest_pd <- merge_pd_covs(Forest_ses_pd, var = "pd.obs.z") 

## Mean Phylogenetic Diversity (MPD)
Forest_ses_mpd <- ses.mpd(samp = Forest_mat, dis = cophenetic(Birds_tree), null.model="sample.pool", iterations = 100)
# Box plot
boxplot(Forest_ses_mpd$mpd.obs.z, names=c("MPD"))

## Weighting by abundance 
# NOTE:: Things don't change much by including abundance
abund.mpd<-ses.mpd(Forest_mat, cophenetic(Birds_tree),null.model="sample.pool", abundance.weighted=TRUE, iterations = 10)
boxplot(cbind(Mpd_ses_tbl$mpd.obs.z, abund.mpd$mpd.obs.z), names=c("MPD", "weighted MPD"))

# Join PD & MPD
Forest_mpd <- merge_pd_covs(Forest_ses_mpd, var = "mpd.obs.z")
Forest_pd_mpd <- Forest_pd %>% left_join(Forest_mpd) %>% 
  relocate(ntaxa:pd.z, .after = Habitat) %>% 
  rename(spp_rich = ntaxa)

# >Functional richness ----------------------------------------------------
# fd_fric() function produces NA at sites w/ fewer species than traits
Forest_fric <- fd_fric(traits = Traits_red, sp_com = Forest_mat) %>% 
  rename_with(~ c("Id_muestreo", "F_ric")) %>% 
  tibble()

Fric_tbl %>% filter(is.na(F_ric)) #43 NAs

# Join
Forest_fd_pd_mpd <- Forest_pd_mpd %>% left_join(Forest_fric)

# >Visualize --------------------------------------------------------------
# There is a non-linear (saturating) positive relationship between species richness and phylogenetic diversity 
Forest_fd_pd_mpd %>% ggplot(aes(x = spp_rich, y = pd.obs)) + 
  geom_point(aes(color = Habitat)) + 
  geom_smooth()

# Prep dataframes for plotting via pivot
df1 <- pivot_scale(df = Forest_fd_pd_mpd, var1 = mpd.z, var2 = pd.z)
df2 <- pivot_scale(Forest_fd_pd_mpd, var1 = spp_rich, 
                   var2 = pd.obs, var3 = F_ric, 
                   scale = TRUE)

## Plot
# a) Important differences in overall amount of phylogenetic diversity captured, but this is primarily driven by differences in species richness
# NOTE: Could log both SR & PD , and try to have 2 y axis lables, one on both side 
# b) Riparian forest is clustered, suggesting there may be filtering for riparian specialists 
Forest_mpd_pd_p <- plot_pd(df1, hline = TRUE, label = c("MPD", "PD"), title = "Z-standardized")
Forest_pd_sr_p <- plot_pd(df2, label = c("Functional \nrichness", "PD", "Species \nrichness"), title = "Meta - Forest species")
ggarrange(Forest_pd_sr_p, Forest_mpd_pd_p)
df2 %>% filter(Metric == "F_ric") %>% select(Value) #%>% 
  #filter(is.na(.))

#NOTE:: Several NAs
df1 %>% filter(Metric %in% c("mpd.z", "pd.z") & is.na(Value)) # 58 NAs
Forest_fric %>% filter(is.na(F_ric)) #82 NAs

# >Models -----------------------------------------------------------------
## Compare MPD by habitat types
# Run an ANOVA & then a Tukey HSD test
aov_tukey(data = Forest_pd_mpd, x = "Habitat", y = "mpd.z")
aov_tukey(data = Forest_pd_mpd, x = "Habitat", y = "pd.z")
aov_tukey(data = Forest_pd_mpd, x = "Habitat", y = "pd.obs")
aov_tukey(data = Forest_pd_mpd, x = "Habitat", y = "spp_rich")


# Visualize forest & all species ------------------------------------------
# Would be nice to have a single legend with 3 colors , & PD having same color across all plots 
Pd_sr_plots <- ggarrange(Pd_sr_p, Forest_pd_sr_p, common.legend = TRUE, nrow = 1)
Mpd_pd_plots <- ggarrange(Mpd_pd_p, Forest_mpd_pd_p, common.legend = TRUE, nrow = 1)
ggarrange(Pd_sr_plots, Mpd_pd_plots, nrow = 2)
ggsave('PD_plot_meta.png', bg = "white")

# EXTRAS ------------------------------------------------------------------
# Could examine Forest v Riparian

# >Abiotic factors ---------------------------------------------------------
# Can do with just forest species or all species 
plot(co_vars_forest$MPD~ co_vars_forest$Elev)
model<-lm(co_vars_forest$MPD~co_vars_forest$Elev)
summary(model) # Negative

plot(co_vars_forest$MPD~co_vars_forest$Avg_temp)
model<-lm(co_vars_forest$MPD~co_vars_forest$Avg_temp)
summary(model) # Positive

## Way too correlated to include in a single model -- coefficients can't be trusted
cor(co_vars_forest$Avg_temp, co_vars_forest$Elev)
model<-lm(MPD~Avg_temp+Elev, data = co_vars_forest)
summary(model) # Elevation has switched signs 
car::vif(model)