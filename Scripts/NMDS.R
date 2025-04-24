# NMDS to understand community composition

#Non-metric multidimensional scaling (NMDS) is a tool that can be used to convey multidimensional species abundance/occurrence data in 2 or 3 dimensions. 
#NMDS in R: https://www.youtube.com/watch?v=OMrtxobDhrM&ab_channel=EcologicalApplicationsinR
#Vegan package in R lecture: https://www.youtube.com/watch?v=sDYPX-CutOY


# Libraries ---------------------------------------------------------------
library(devtools)
library(vegan)
library(ggvegan)
library(ggordiplots)
library(tidyverse)
library(xlsx)
library(plotly)
ggplot2::theme_set(cowplot::theme_cowplot())


# To do -------------------------------------------------------------------

# Explore possibility of plotting at the level of the genus or family for easier visualization


# Mite data ---------------------------------------------------------------
?mite
data(mite)
str(mite) #Counts of 35 species (columns) of mites at 70 sites (rows) 
data(mite.env)
data(mite.pcnm) #principal components of neighboring matrices. The first few (3 in the video) vectors correspond to broad spatial scale autocorrelation=
data(mite.xy)

#Redundancy analysis. Would need to rewatch this to understand exactly what's going on.
rda(mite.env[,1:2]) #If you just include a single set of vars it runs a normal PCA
#When you add spp. data as response variable this runs an RDA. 
rda <- rda(mite ~ SubsDens + WatrCont, data = mite.env) 
#When interpreting output "Inertia' is equivalent to Variance. 'Constrained' is the variance that is explained and 'unconstrained' is residuals 
rda
rda2 <- rda(mite ~ . + as.matrix(mite.pcnm[1:3]), data = mite.env) #The '.' means include all of the variables in the data frame (i.e., mite.env data frame). 
#Difficult to interpret (remove some vars for better visualization), but shows how the sites are located in space and the arrows show predictor vars. I think the RDA axes must somehow represent the spp. communities? 
autoplot(rda2, arrows = T) 


#NMDS Key steps (see Quant methods note): 
#Step 1) calculates the ‘real-world’ distance matrix (based on specified distance metric and all axes (which are species); there is only 1 right answer here). Step 2) tries to reduce the number of axes to the number of dimensions specified (usually 2-3). This then produces the ‘ordination space’ matrix, which describes the euclidian distance between sites in ordination space. Step 3) iteratively repeats step 2 trying to maximize the rank-order correlation between the ‘real-world’ distance and ‘ordination space’ distance & reduce the amount of stress.
?metaMDS

# Load bird data--------------------------------------------------------------
# Bring in community data -- currently a dataframe
Meta_birds <- read.xlsx("Derived/Meta_314.xlsx", sheetIndex = 1) %>% tibble()

# Turn to matrix and transpose
# Remove sites where no species were observed
Birds_mat <- Meta_birds %>% arrange(Id_muestreo) %>% 
  column_to_rownames("Id_muestreo") %>% 
  filter(rowSums(.) != 0) %>%
  as.matrix()

# Bring in covariates, one row per site
Site_covs <- read.xlsx("Derived/Site_covs314.xlsx", sheetIndex = 1)

# Functional traits
Traits <- read.xlsx("Derived/Morphology_ft.xlsx", sheetIndex = 1)

#NMDS to reduce the 35 axes (one for each spp.) into just 2 (or 3) dimensions). Based on bray curtis distances, where a score of 0 means that 2 sites are identical in composition and abundance, and 1 means two sites are entirely different (pure turnover). This can help visualize differences in alpha & beta diversity!
#Note metaMDS does rotate the axes so that axis 1 contains the greatest variance (similar to PCA)

Site_covs2 <- Site_covs %>% arrange(Id_muestreo) %>% 
  filter(Id_muestreo %in% rownames(Birds_mat))


# Run NMDS ----------------------------------------------------------------
# Hellinger transformation minimizes impact of super abundant or super rare data
nmds <- metaMDS(Birds_mat, k = 3, try = 20, trymax = 20, autotransform = T)

nmds #Counts are transformed by default. Stress value shows how much the data had to be contorted to display it in only k = 2 (or 3) dimensions. metaMDS is iteratively changing the points to try to reduce that stress. Rule of thumb is stress < 0.2 means you can interpret the axes, stress < 0.1 is ideal. If stress is too high can increase k = 3 (3 dimensions), or can increase trymax (try too?). 
autoplot(nmds, geom = "text", legend = "none") #Red numbers are sites, blue words are bird taxa. The sites w/ positive scores on NMDS1 have higher abundances of the species w/ positive scores for NMDS1 
col <- c("red", "blue", "green", "orange")
shp <- c(18, 20, 15, 21)
#Will want to do this by LC type
#Base R option (no thanks!)
plot(nmds$points, col = col[Site_covs2$Habitat_cons]) #pch = shp[mite.env$Shrub]
#ordispider(nmds, groups = mite.env$Shrub, label = T)

#GGplot option
gg_ordiplot(nmds, groups = Site_covs2$Habitat_cons, kind = "sd", spiders = T) #hull = T

#Create a plot similar to Jo's Fig 19c
#I think key to really customizing these plots is to set plot = F and using the internal data frame

#See Karp 2019 paper Fig 2 for ideas on how you might organize this plot. Can color code sites by continuous gradients (e.g., precipitation), or put ellipses around points that belong to the same categorical variable (land-use category).
fort <- fortify(nmds) #Creates a data frame with 4 columns: score (sites, 70 rows; & Spp., 35 rows), label (site# or spp. code), and then the coordinates in NMDS ordination space. Gives lots of control
nrow(fort) # 1 row per site & 1 row per species
191 + 260

## Can visualize continuous predictor variables in ordination space in a number of ways 
# Base R
ordisurf(nmds ~ Site_covs2$Elev, plot = T, knots = 3) #Change plot to F to get ordisurf model info; see ?envifit and search ordisurf as well

## Ggplot
gg_ordisurf(nmds, env.var = Site_covs2$Elev, groups = Site_covs2$Habitat_cons, var.label = "Elevation") #+ labs(fill = "Shrubs") Help file says this should control legend but no luck

# Split sites & species dataframes
Sites_nmds <- fort %>% filter(score == "sites") %>% 
  rename(Id_muestreo = label) %>% 
  left_join(Site_covs2) %>% 
  tibble()
Spp_nmds <- fort %>% filter(score == "species") %>% 
  tibble()

# Sites colored by continuous variable (similar to Fig 2a and 2b in Karp et al. (2019))
ggplot() + geom_point(
  data = Sites_nmds, 
  aes(x = NMDS1, y = NMDS2, color = Elev, shape = as.factor(Habitat_cons)), 
  size = 4, alpha = .8
)

# 3-dimensional plot
plotly::plot_ly(Sites_nmds, 
                x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, 
                type = 'scatter3d', mode = 'markers',
                color = ~Habitat_cons)

# Ellipses + spp position plots ------------------------------------------
# Custom function leveraging gg_ordiplot to allow easy generation of centroids & ellipses for all axes combinations
ellipse_df <- function(choices_axes){
  gg_ordiplot(nmds, groups = Site_covs2$Habitat_cons, choices = choices_axes, plot = FALSE)
} 

# Create dataframes for all axes combinations
ellipses_axes12 <- ellipse_df(choices_axes = c(1,2))
ellipses_axes23 <- ellipse_df(choices_axes = c(2,3))
ellipses_axes13 <- ellipse_df(choices_axes = c(1,3))

# Custom function plot_nmds to plot ellipses, centroids, & species names
plot_nmds <- function(x, y, ellipse_axes_df){
  ggplot() +
    # species names
    ggrepel::geom_text_repel(
      data = Spp_nmds, aes(x = {{ x }}, y = {{  y }}, label = label)
      # Centroids
    ) + geom_point(data = ellipse_axes_df$df_mean.ord, 
                   aes(x = x, y = y, color = Group, size = 4)) + 
    # Ellipses
    geom_path(data = ellipse_axes_df$df_ellipse, 
              aes(x,y, color = Group), 
              linetype = "dashed", linewidth = 2, alpha = .8) + 
    theme(legend.position = "top") +
    #labs(title = "NMDS of 118 / 260 species") +
    guides(size = "none")
}

# Generate plots with 
p12 <- plot_nmds(x = NMDS1, y = NMDS2, ellipse_axes_df = ellipses_axes12)
p23 <- plot_nmds(x = NMDS2, y = NMDS3, ellipse_axes_df = ellipses_axes23)
p13 <- plot_nmds(x = NMDS1, y = NMDS3, ellipse_axes_df = ellipses_axes13)
nmds_plots_l <- list(axes_12 = p12, axes_23 = p23, axes_13 = p13)

# Save plots
if(FALSE){
  imap(nmds_plots_l, \(plot, name){
    print(plot)
    ggsave(paste0("NMDS_meta_", name, ".png"), 
           bg = "white", width = 15, height = 11)
  })
}
#ggarrange(p12, p23, p13, ncol = 1, common.legend = TRUE)

#gg_ordibubble(nmds, env.var = mite.env$SubsDens)

#Vegan::adonis calls a manova, essentially looking for differences in spp. (respone variable) by a categorical explanatory variable (e.g., shrubs: none, few, many).


# Statistical tests -------------------------------------------------------


#We conclude that bird communities differ significantly by the categorical var Habitat, but the model explains < 5% of the variation in community composition
aManova <- adonis2(Birds_mat ~ Habitat_cons, data = Site_covs2, permutations = 999, method = "bray")
aManova 

#Test whether mite communities vary along continuous gradients
envfit <- envfit(nmds ~ Elev + Avg_temp + Tot.prec, data = Site_covs2) #Both significant, and water content explains 70% of the variation in mite community
plot(envfit)
#We need to ensure that there is equal dispersion by shrub class to ensure that adonis results are statistically valid
dist.birds <- vegdist(Birds_mat, method = "bray") #bray curtis distances
anova(betadisper(dist.birds, Site_covs2$Habitat_cons)) #We conclude that variances are indeed equal