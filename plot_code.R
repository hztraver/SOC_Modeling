#### 01. K-MEANS CLUSTERS ####
## cluster transects based on landsat spectral indices and topographic variation

library(factoextra) ; library(FactoMineR) ; library(cluster); library(ggplot2); library(tidyr); library(dplyr) ; library(forcats); library(purrr) ; library(mgcv)

# Landsat and topographic variation aggregated within each transect
cluster.data = read.csv("Code & Data/01_transect_clusters.csv")

# scale variables for clustering
cvars = cluster.data[, 2:10]
cvars = scale(cvars)
row.names(cvars) = cluster.data$transect

# group transects into three clusters 
set.seed(123)
clusters = kmeans(cvars, centers = 3, nstart = 25)

# visualize clusters
fviz_cluster(clusters, data = cvars, ggtheme = theme_bw(), repel = T, 
             shape = 16, labelsize = 9, main = "", ellipse.alpha = 0.1, ellipse.type = 'confidence')+
                theme(axis.title = element_text(size = 14),
                   axis.text = element_text(size = 12),
                   legend.position = "bottom")

# get contributions of variables to principle components 
pca = PCA(cvars, graph=F)
contrib = data.frame(pca$var$contrib)
contrib$var = row.names(contrib)
contrib$order = c(8, 7, 9, 6, 3, 5, 1, 4, 1)

# plot contribution of variables along both dimensions
contrib %>% dplyr::select(Dim.1, Dim.2, var, order) %>% 
  pivot_longer(Dim.1:Dim.2, names_to = "Dimension", values_to = "Contribution") %>%
  ggplot(aes(x = reorder(var, order), y = Contribution, fill = Dimension)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  xlab("")+
  ylab("Contribtion %")+
  theme_bw()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.text = element_text(size = 14),
        legend.position = "bottom")

# print transect names and cluster types
factor(clusters$cluster) %>% forcats::fct_recode(., "Moist Tundra" = "1", "Mesic Woodland" = "2", "Submesic Woodland" = "3")

## boxplots of clustering variables extracted from sample plots
cluster.boxplots = read.csv("Code & Data/01b_transect_boxplots.csv")

#### Figure 3 ####
# Transform tree stems count (per 30m plot) to stems per ha
cluster.boxplots $`Stem Density` = (cluster.boxplots$tree_stem_count / 2827)*10000

long = cluster.boxplots %>% dplyr::select(c(cluster_names, EVI, TCW, TCG, DTM.range, Slope.mean, `Stem Density`)) %>% 
  pivot_longer(!cluster_names, names_to = "variable", values_to = "value")

# relevel factor to control order of plots
long$variable = factor(long$variable, levels = c('EVI', 'TCG', 'TCW', 'Slope.mean', 'DTM.range', 'Stem Density'))
long$variable = forcats::fct_recode(factor(long$variable), "Mean Slope" = "Slope.mean", "Elevation Range" = "DTM.range")
long$cluster_names = relevel(factor(long$cluster_names), ref = "Moist Tundra")

ggplot(aes(x = cluster_names, y = value, color = cluster_names), data = long) +
  facet_wrap(~variable, scales = "free")+
  scale_color_manual(values=c("#33608C","#F5B35E","#B81840"))+
  geom_boxplot()+
  xlab("Cluster")+
  ylab("")+
  theme_bw()+
  theme(aspect.ratio = 1, legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12))

#### 02. VEGETATION FUNCTIONAL GROUPS #### 

## percent cover of species within each vegetation plot (n=55)
veg_data = read.csv("Code & Data/02_veg_fg.csv")

# define functional types based on species
#### Figure 4 ####
functional.groups = data.table(transect = veg_data$Transect, plot = veg_data$plot,
           Shrubs = (veg_data$Betula + veg_data$Alder + veg_data$Salix + veg_data$Lab..Tea),
           Trees = veg_data$Spruce,
           Moss = veg_data$Moss,
           Lichen = veg_data$Lichen,
           Litter = (veg_data$Litter + veg_data$Dead.Wood),
           Graminoid = veg_data$Gramminoid,
           Herbs = (veg_data$Pedasites + veg_data$Rub..Cha + veg_data$Equisetum + veg_data$Pyrola + veg_data$Pedicularis + 
                      veg_data$Oxycocus),
           `Dwarf Shrubs` = (veg_data$Vacc..Vit + veg_data$Vacc..Uli + veg_data$Arctostaphylos +
                               veg_data$Dryas  + veg_data$Empetrum), 
           cluster_names = relevel(factor(veg_data$cluster), ref = "Moist Tundra"))
  
functional.groups %>% pivot_longer(Shrubs:`Dwarf Shrubs`) %>% 
  ggplot(aes(x = cluster_names, y = value, color = cluster_names)) +
  geom_boxplot() +
  scale_color_manual(values=c("#33608C","#F5B35E","#B81840"))+
  ylab("Percent Cover")+
  xlab("Cluster")+
  facet_wrap(~name, nrow = 2)+
  theme_bw()+
  theme(aspect.ratio = 1, legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12))
  
#### 03. SOC MODELLING ####

## Soil sample data joined to vegetation structure (3m plots) 
## and microtopography (pixel value at sample location) from the drone CHM and DTM
## tree stem density calculated from number of stems > 2.5m within 30m plot of each sample location

soc = read.csv("Code & Data/03_SOC.csv") %>% filter(mineral < 30) # soils with mixed mineral/organic not used for modeling
soc$cluster_names = relevel(factor(soc$cluster_names), ref = "Moist Tundra")
  
#### Figure 6 ####
# convert tree stem count to stem density per ha within plots 2827m2
ggplot(aes(x = (tree_stem_count/2827.433)*10000, y = SOC_kgCm2), data = soc) + 
  geom_smooth(colour = "#71A3CA", method = 'gam', formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3)) + # define gam smoothing parameters
  geom_point(alpha = 0.5, aes(size = org_soil_depth, color = density_kgm3))+
  scale_color_gradient2(low = "#33608C", mid = "#F5B35E", high = "#B81840", midpoint = 300) + 
  ylim(c(10,40))+
  scale_y_continuous(trans = "log10") +
  labs(y = bquote(paste("SOC (kgC/", m^2, ")")), x = 'Stems per hectare',
       size = "Organic Soil Depth (cm)", color = bquote(paste("Bulk Density (kg/", m^3, ")"))) +
  theme_bw()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

ggplot(aes(x = factor(stembins), y = SOC_kgCm2, color = cluster_names), data = soc) + 
  geom_boxplot()+
  scale_color_manual(values=c("#33608C","#F5B35E","#B81840"))+
  facet_wrap(~cluster_names, ncol = 1)+
  labs(y = bquote(paste("SOC (kgC/", m^2, ")")), x = 'Stems per hectare', color = "")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16))

#### linear models ####
## sample code for running stepwise variable selection for SOC storage within treeline types 
# cluster names
cluster.names = c("Moist Tundra", "Mesic Woodland", "Submesic Woodland")
# SOC storage, bulk density, carbon content (%), organic soil depth
variable.names = c("SOC_kgCm2", "density_kgm3", "C", "org_soil_depth") 

for (cluster in cluster.names) {

  for (variable in variable.names) {
    
    print("..........................")
    print(paste("Cluster: ", cluster))
    print(paste("Response Variable: ", variable))
    soc.subset = soc[soc$cluster_names == cluster,]
    
    # fit intercept only model & full model with all variables
    int_formula = as.formula(paste(variable, "~ 1"))
    all_formula = as.formula(paste(variable, "~ log(mean_height) + log(max_height) + sd_height + litter_layer + log(roughness) + slope + aspect + tpi"))
    
    int = lm(int_formula, data = soc.subset)
    all = lm(all_formula, data = soc.subset)
    # run stepwise selection 
    final = step(int, direction = 'both', scope = formula(all), trace = 0)
    
    # print p-values, coefficients and variable names
    pvals = summary(final)$coefficients[, "Pr(>|t|)"]
    estimate = summary(final)$coefficients[, "Estimate"]
    coef = names(pvals) 
    
    print("Stepwise Model Coefficients")
    print(estimate) 
    print("P-values")
    print(pvals)
  }
}

#### Figure 7 ####
## Soil carbon kgC/m2
soc %>% 
  ggplot(aes(x = mean_height, y = SOC_kgCm2, color = factor(cluster_names))) +
  facet_wrap(~cluster_names, scales = "fixed")+
  geom_point(alpha = 0.3) +
  scale_color_manual(values=c("#33608C","#F5B35E","#999999"))+
  geom_smooth(method = "lm", aes(group = cluster), linewidth = 1.2) +
  scale_x_continuous(trans = "log10")+
  labs(y = bquote(paste("SOC (kgC/", m^2, ")")), x = 'Mean Canopy Height')+
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(color = "cluster_names")

## Soil density kg/m3
soc %>% ggplot(aes(x = mean_height, y = density_kgm3, color = factor(cluster_names))) +
  facet_wrap(~cluster_names, scales = "fixed")+
  geom_point(alpha = 0.3) +
  scale_color_manual(values=c("#999999","#F5B35E","#999999"))+
  geom_smooth(method = "lm", aes(group = cluster), linewidth = 1.2) +
  labs(y = bquote(paste("Bulk Density (kg/", m^3, ")")), x = 'Mean Canopy Height')+
  scale_x_continuous(trans = "log10")+
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(color = "cluster_names")

## Carbon content (%)
soc %>% ggplot(aes(x = max_height, y = C, color = factor(cluster_names))) +
  facet_wrap(~cluster_names, scales = "fixed")+
  geom_point(alpha = 0.3) +
  scale_color_manual(values=c("#999999","#F5B35E","#B81840"))+
  geom_smooth(method = "lm", aes(group = cluster), linewidth = 1.2) +
  scale_x_continuous(trans = "log10")+
  xlab('Max Canopy Height')+
  ylab('C Content (%)')+
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(color = "cluster_names")

#### 04. PREDICT DRONE METRICS FROM LANDSAT ####

## Landsat pixels joined to drone metrics aggregated to 30m resolution 
##  mean ExG per pixel, maximum canopy height, tree stem count 
lsat = read.csv("Code & Data/04_landsat_predictors.csv")

#### Generalized Additive Models ####
landsat.predictors = c("EVI", "TCW", "TCG", "TCB")
drone.metrics = c("Excess Greenness", "Max Height", "Tree Stems")

for(drone.metric in drone.metrics) {
  
  print("..........................")
  print(paste("Drone Metric: ", drone.metric))
  subset = lsat[lsat$drone_metric == drone.metric,] %>% pivot_longer(., cols = c(EVI:TCG)) 
  
  for (landsat.predictor in landsat.predictors) {
    print(paste("Landsat Predictor: ", landsat.predictor))
    
    gam.subset = subset[subset$name == landsat.predictor,]
    gam_formula = as.formula("drone_value ~ s(value, bs = \"cs\", fx = TRUE, k = 3)")    
    
    # fit the gam and print r-squared and RMSE
    x = summary(mgcv::gam(formula = gam_formula, data = gam.subset))
    print(paste("r2 = ", x$r.sq))
    print(paste("RMSE = ", sqrt(x$scale)))
  }
}

#### Figure 8 ####
## Excess Greenness
lsat %>% filter(drone_metric == "Excess Greenness") %>% pivot_longer(., cols = c(EVI:TCG)) %>%
  ggplot(aes(x = value/100, y = drone_value, color = name), data = .) +
  geom_point(alpha = 0.1) + 
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3))+
  scale_color_manual(values=c("#999999", "#999999", "#999999","#D82632" ))+
  facet_wrap(~ name, scales = "free_x", ncol = 4)+
  ylab("Excess Greenness")+
  xlab("Index Value")+
  theme_bw()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        panel.spacing = unit(1.5, "lines"))

## Maximum Canopy Height
lsat %>% filter(drone_metric == "Max Height") %>% pivot_longer(., cols = c(EVI:TCG)) %>%
  ggplot(aes(x = value/100, y = drone_value, color = name), data = .) +
  geom_point(alpha = 0.1) + 
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3))+
  scale_color_manual(values=c("#999999", "#999999","#D82632", "#999999"))+
  facet_wrap(~ name, scales = "free_x", ncol = 4)+
  ylab("Maximum Canopy Height (m)")+
  xlab("Index Value")+
  theme_bw()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        panel.spacing = unit(1.5, "lines"))

## Tree Stems
lsat %>% filter(drone_metric == "Tree Stems") %>% pivot_longer(., cols = c(EVI:TCG)) %>%
  ggplot(aes(x = value/100, y = drone_value, color = name), data = .) +
  geom_point(alpha = 0.1) + 
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3))+
  scale_color_manual(values=c("#999999","#D82632", "#999999", "#999999"))+
  facet_wrap(~ name, scales = "free_x", ncol = 4)+
  scale_y_continuous(trans='log10')+
  ylab("Tree Stems")+
  xlab("Index Value")+
  theme_bw()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        panel.spacing = unit(1.5, "lines"))
