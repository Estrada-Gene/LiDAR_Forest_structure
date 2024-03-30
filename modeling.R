### LiDAR-derived Forest Structure Metrics as Predictors of Mammal Species Richness and Diversity

library(tidyverse)
library(vegan)
library(cowplot)

#mammal observation data - excludes all human observations
m <- read.csv(file = "data/pruned.csv", header = TRUE)
m <- m[,-1]
#ordering forest type
m$habitat <- factor(m$habitat, levels = c("Peat Swamp", "Freshwater Swamp", "Alluvial Bench", 
                                          "Lowland Sandstone", "Lowland Granite", "Upland Granite", "Montane"))

#mammal observation table
species_table <- as.data.frame(table(m$Species, m$habitat))
species_table <- pivot_wider(species_table, names_from = "Var2", values_from = "Freq")
colnames(species_table)[colnames(species_table) == 'Var1'] <- 'Species'
#write.csv(species_table, file = "species_observation_table_byFT.csv")

#LiDAR-derived metrics data
stand_mets <- read.csv(file = "data/select.stand.mets.csv", header = TRUE)
stand_mets <- stand_mets[,-1]
stand_mets$habitat <- factor(stand_mets$habitat, levels = c("Peat Swamp", "Freshwater Swamp", "Alluvial Bench", 
                                          "Lowland Sandstone", "Lowland Granite", "Upland Granite", "Montane"))

tree_mets <- read.csv(file = "data/all.tree.mets.csv", header = TRUE)
tree_mets <- tree_mets[,-1]
tree_mets$habitat <- factor(tree_mets$habitat, levels = c("Peat Swamp", "Freshwater Swamp", "Alluvial Bench", 
                                                            "Lowland Sandstone", "Lowland Granite", "Upland Granite", "Montane"))
#camera trap data
ct_elev <- read.csv(file = "data/cameradata_updatedZJ-ajm.csv", header = TRUE)
ct_elev <- ct_elev[!duplicated(ct_elev$locationID),]

#some observations noted 'Unid' - unidentified to species and/or genus level
table(m$Species)
table(m$Genus)
table(m$Family)


### SPECIES RICHNESS 

## Entire Study Area - excluding unid animals
#all camera trap locations: 55 species
m %>%
  filter(!(Species %in% c("Confused Squirrel or Treeshrew","Tragulus spp.","Tupai sp.","Unid bats",
                      "Unid civet","Unid rat","Unid squirrel"))) %>%
  summarise(SpeciesRichness = n_distinct(Species))

#only camera trap locations that have LiDAR scan data: 52 species
m %>%
  filter(locationID %in% stand_mets$locationID) %>%
  filter(!(Species %in% c("Confused Squirrel or Treeshrew","Tragulus spp.","Tupai sp.","Unid bats",
                          "Unid civet","Unid rat","Unid squirrel"))) %>%
  summarise(SpeciesRichness = n_distinct(Species))

## By Forest Type - excluding unid animals
#all ct locations
n_by_ft <- m %>%
  group_by(habitat) %>%
  filter(!(Species %in% c("Confused Squirrel or Treeshrew","Tragulus spp.","Tupai sp.","Unid bats",
                          "Unid civet","Unid rat","Unid squirrel"))) %>%
  summarise(n.all = n_distinct(Species)) %>%
  arrange(desc(habitat))

#only scanned locations
n_by_ft <- m %>%
  filter(locationID %in% stand_mets$locationID) %>%
  group_by(habitat) %>%
  filter(!(Species %in% c("Confused Squirrel or Treeshrew","Tragulus spp.","Tupai sp.","Unid bats",
                          "Unid civet","Unid rat","Unid squirrel"))) %>%
  summarise(n.scanned = n_distinct(Species)) %>%
  left_join(n_by_ft) %>%
  arrange(desc(habitat))

#the scanned sites consistently observe fewer species over the study period by 4 - 10 species
#more bias in the higher elevation FT's, especially the montane
n_by_ft$diff <- n_by_ft$n.all - n_by_ft$n.scanned
n_by_ft

## By Camera Trap Location - excluding unid animals
#all CT locations
n_by_ct <- m %>%
  group_by(locationID) %>%
  filter(!(Species %in% c("Confused Squirrel or Treeshrew","Tragulus spp.","Tupai sp.","Unid bats",
                          "Unid civet","Unid rat","Unid squirrel"))) %>%
  summarise(n.all = n_distinct(Species))

#compare 'all CT' and 'scanned CT' location richness distributions
#the 'scanned CT' locations seem biased towards higher species richness
hist(n_by_ct$n.all, breaks = 30)
hist(n_by_ct$n.all[n_by_ct$locationID %in% stand_mets$locationID], breaks = 30, xlim = c(0,30))

summary(n_by_ct$n.all)
summary(n_by_ct$n.all[n_by_ct$locationID %in% stand_mets$locationID])


### SPECIES DIVERSITY

## Shannon Diveristy
## by camera trap location
shannon_ct <- diversity(table(m$locationID, m$Species), index = "shannon")
shannon_ct <- rownames_to_column(as.data.frame(shannon_ct), var = "locationID")
shannon_ct$locationID <- as.integer(shannon_ct$locationID)

#compare 'all CT' and 'scanned CT' location Shannon Diversity value distributions
#the 'scanned CT' locations seem slightly biased towards higher diversity values
hist(shannon_ct$shannon_ct, breaks = 50, xlim = c(0,3))
hist(shannon_ct$shannon_ct[shannon_ct$locationID %in% stand_mets$locationID], breaks = 50, xlim = c(0,3))

summary(shannon_ct$shannon_ct)
summary(shannon_ct$shannon_ct[shannon_ct$locationID %in% stand_mets$locationID])

## by forest type
shannon_all <- diversity(table(m$habitat, m$Species), index = "shannon")
shannon_all <- rownames_to_column(as.data.frame(shannon_all), var = "habitat")

shannon_scan <- diversity(table(m$habitat[shannon_ct$locationID %in% stand_mets$locationID], 
                m$Species[shannon_ct$locationID %in% stand_mets$locationID]), index = "shannon")
shannon_scan <- rownames_to_column(as.data.frame(shannon_scan), var = "habitat")
shannon_ft <- left_join(shannon_all, shannon_scan)
rm(shannon_scan, shannon_all)

#again, using only the CT's that were scanned leads to most bias in the montane forest
shannon_ft$diff <- shannon_ft$shannon_all - shannon_ft$shannon_scan
shannon_ft


## Simpson Diversity
#by CT location
simp_ct <- diversity(table(m$locationID, m$Species), index = "simpson")
simp_ct <- rownames_to_column(as.data.frame(simp_ct), var = "locationID")
simp_ct$locationID <- as.integer(simp_ct$locationID)

#compare 'all CT' and 'scanned CT' location Simpson Diversity value distributions
#distributions look similar
hist(simp_ct$simp_ct, breaks = 50, xlim = c(0,1))
hist(simp_ct$simp_ct[simp_ct$locationID %in% stand_mets$locationID], breaks = 50, xlim = c(0,1))


#adding species richness and diversity data to stand metrics table
stand_mets <- left_join(stand_mets, n_by_ct)
stand_mets <- left_join(stand_mets, shannon_ct)
stand_mets <- left_join(stand_mets, simp_ct)

#calculating Pielou's diversity measure of species evenness
div_mets <- shannon_ct %>%
  left_join(., n_by_ct) %>%
  mutate(div_even = shannon_ct/log(n.all)) %>%
  left_join(simp_ct)
div_mets$div_even[div_mets$div_even == "NaN"] <- 0

#add CT survey effort data
ct <- read.csv(file = "data/ofp_deployments-2021-11-04.csv", header = TRUE)
ct_active_days <- ct %>%
  select(Deployment.Location.ID, Camera.Deployment.Begin.Date, Camera.Deployment.End.Date) %>%
  mutate(Camera.Deployment.Begin.Date = as.Date(Camera.Deployment.Begin.Date),
         Camera.Deployment.End.Date   = as.Date(Camera.Deployment.End.Date)) %>%
  mutate(act.per = as.numeric(difftime(Camera.Deployment.End.Date, 
                                       Camera.Deployment.Begin.Date, units = "days"))) %>%
  group_by(Deployment.Location.ID) %>%
  summarise(act.days = sum(act.per)) %>%
  rename(locationID = Deployment.Location.ID)

stand_mets <- left_join(stand_mets, ct_active_days)  
stand_mets <- left_join(stand_mets, div_mets[,c("locationID","div_even")])

#scale and center all numeric predictors 
stand_mets <- stand_mets %>%
  mutate(across(8:19, scale, .names = "{.col}.sc"))


## Richness and diversity visualization
ct <- ct %>%
  rename(locationID   = Deployment.Location.ID,
         habitat      = Treatment.Strata,
         location     = Location,
         on.off.trail = Feature.Type)

#add habitat data to diversity metrics table
div_mets <- left_join(div_mets, ct_elev[,c("locationID","habitat","altitude")])
div_mets$habitat <- factor(div_mets$habitat, levels = c("Peat Swamp","Freshwater Swamp","Alluvial Bench",
                                                        "Lowland Sandstone","Lowland Granite",
                                                        "Upland Granite","Montane"))
#copying for visualization use only
div_mets$dataset <- "all"
div_mets$dataset[div_mets$locationID %in% stand_mets$locationID] <- "scanned"
div_mets_viz <- div_mets
div_mets_viz_scan <- div_mets_viz[div_mets_viz$dataset == "scanned",]
div_mets_viz$dataset <- "all"
div_mets_viz <- rbind(div_mets_viz, div_mets_viz_scan)

#plotting species richness by CT location/forest type - all CT locations only
ggplot(div_mets, aes(x = habitat, y = n.all)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_flip() +
  theme_classic() +
  labs(title = "Species Richness by CT - all locations", x = "", y = "n species")

#plotting species richness by CT location/forest type - scanned sites only
ggplot(div_mets[div_mets$locationID %in% stand_mets$locationID,], aes(x = habitat, y = n.all)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_flip() +
  theme_classic() +
  labs(title = "Species Richness by CT - scanned sites", x = "", y = "n species")

#plotting both in same plot
p1 <- ggplot(div_mets_viz, aes(x = habitat, y = n.all, fill = dataset)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.5) +
  geom_jitter(position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("cornflowerblue","coral3")) +
  theme_classic() +
  labs(title = "Species Richness", x = "", y = "n species", fill = "") +
  coord_flip() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

#plotting Shannon diversity by CT location/forest type - all CT locations
shannon_ct %>%
  left_join(., ct[,c("locationID", "habitat")]) %>%
  mutate(habitat = factor(habitat, levels = c("Peat Swamp", "Freshwater Swamp", "Alluvial Bench", "Lowland Sandstone", 
                                              "Lowland Granite", "Upland Granite", "Montane"))) %>%
  ggplot(aes(x = habitat, y = shannon_ct)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_flip() +
  theme_classic() +
  ylim(0, 3) +
  labs(title = "Shannon Diversity by CT - all locations", x = "", y = "n species")

#plotting Shannon diversity by CT location/forest type - scanned sites only
ggplot(stand_mets, aes(x = habitat, y = shannon_ct)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_flip() +
  theme_classic() +
  ylim(0, 3) +
  labs(title = "Shannon Diversity by CT", x = "", y = "Shannon Diversity Index")

#plotting both in same plot
p2 <- ggplot(div_mets_viz, aes(x = habitat, y = shannon_ct, fill = dataset)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.5) +
  geom_jitter(position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("cornflowerblue","coral3")) +
  theme_classic() +
  labs(title = "Shannon Diversity", x = "", y = "Shannon Diversity Index", fill = "") +
  coord_flip() +
  theme(legend.position = "none", axis.text.y = element_blank(), 
        plot.title = element_text(hjust = 0.5))

#plotting Simpson Diversity by CT/forest type - all CT locations
simp_ct %>%
  left_join(., ct[,c("locationID", "habitat")]) %>%
  mutate(habitat = factor(habitat, levels = c("Peat Swamp", "Freshwater Swamp", "Alluvial Bench", "Lowland Sandstone", 
                                              "Lowland Granite", "Upland Granite", "Montane"))) %>%
  ggplot(aes(x = habitat, y = simp_ct)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_flip() +
  theme_classic() +
  labs(title = "Simpson Diversity by CT - all locations", x = "", y = "Simpson Diversity Index")

#plotting Simpson diversity by CT location/forest type - scanned sites only
ggplot(stand_mets, aes(x = habitat, y = simp_ct)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_flip() +
  theme_classic() +
  labs(title = "Simpson Diversity by CT - scanned locations", x = "", y = "Simpson Diversity Index")

#plotting both in same plot
ggplot(div_mets_viz, aes(x = habitat, y = simp_ct, fill = dataset)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.5) +
  geom_jitter(position = position_dodge(width = 0.75)) +
  theme_classic() +
  labs(title = "Simpson Diversity", x = "", y = "Simpson Diversity Index", fill = "") +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) +
  theme(legend.position = c(.875, .115))


#plotting species evenness by CT location/forest type - all CT locations
div_mets %>%
  left_join(., ct[,c("locationID", "habitat")]) %>%
  mutate(habitat = factor(habitat, levels = c("Peat Swamp", "Freshwater Swamp", "Alluvial Bench", "Lowland Sandstone", 
                                              "Lowland Granite", "Upland Granite", "Montane"))) %>%
  ggplot(aes(x = habitat, y = div_even)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_flip() +
  theme_classic() +
  labs(title = "Species Evenness by CT - all locations", x = "", y = "species evenness index")

#plotting both in same plot
p3 <- ggplot(div_mets_viz, aes(x = habitat, y = div_even, fill = dataset)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.5) +
  geom_jitter(position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("cornflowerblue","coral3")) +
  theme_classic() +
  labs(title = "Species Evenness", x = "", y = "Evenness Index", fill = "") +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) +
  theme(legend.position = "none", axis.text.y = element_blank(), 
        plot.title = element_text(hjust = 0.5))

plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(.4, .3, .3))


#species richness by elevation - 'all CT locations' and "scanned locations' in the same plot
p1 <- ggplot(div_mets, aes(x = altitude, y = n.all)) +
  geom_point(data = subset(div_mets, locationID %in% stand_mets$locationID), aes(color = "scanned")) +
  geom_point(data = subset(div_mets, !locationID %in% stand_mets$locationID), aes(color = "default")) +
  geom_smooth(data = subset(div_mets, locationID %in% stand_mets$locationID), aes(color = "scanned"), se = TRUE) +
  geom_smooth(data = div_mets, aes(color = "all"), se = TRUE) +  
  scale_color_manual(values = c("scanned" = "red", "all" = "blue", "default" = "black")) +  
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Species Richness", x = "elevation (m)", y = "n species")

#shannon diversity by elevation
p2 <- ggplot(div_mets, aes(x = altitude, y = shannon_ct)) +
  geom_point(data = subset(div_mets, locationID %in% stand_mets$locationID), aes(color = "scanned")) +
  geom_point(data = subset(div_mets, !locationID %in% stand_mets$locationID), aes(color = "default")) +
  geom_smooth(data = subset(div_mets, locationID %in% stand_mets$locationID), aes(color = "scanned"), se = TRUE) +
  geom_smooth(data = div_mets, aes(color = "all"), se = TRUE) +  
  scale_color_manual(values = c("scanned" = "red", "all" = "blue", "default" = "black")) +  
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Shannon Diversity", x = "elevation (m)", y = "Shannon Diversity Index")

#simpson diversity by elevation
ggplot(div_mets, aes(x = altitude, y = simp_ct)) +
  geom_point(data = subset(div_mets, locationID %in% stand_mets$locationID), aes(color = "scanned")) +
  geom_point(data = subset(div_mets, !locationID %in% stand_mets$locationID), aes(color = "default")) +
  geom_smooth(data = subset(div_mets, locationID %in% stand_mets$locationID), aes(color = "scanned"), se = TRUE) +
  geom_smooth(data = div_mets, aes(color = "all"), se = TRUE) +  
  scale_color_manual(values = c("scanned" = "red", "all" = "blue", "default" = "black")) +  
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Simpson Diversity", x = "elevation (m)", y = "Simpson Diversity Index")

#species evenness by elevation
p3 <- ggplot(div_mets, aes(x = altitude, y = div_even)) +
  geom_point(data = subset(div_mets, locationID %in% stand_mets$locationID), aes(color = "scanned")) +
  geom_point(data = subset(div_mets, !locationID %in% stand_mets$locationID), aes(color = "default")) +
  geom_smooth(data = subset(div_mets, locationID %in% stand_mets$locationID), aes(color = "scanned"), se = TRUE) +
  geom_smooth(data = div_mets, aes(color = "all"), se = TRUE) +  
  scale_color_manual(values = c("scanned" = "red", "all" = "blue", "default" = "black"),
                     name = NULL,
                     breaks = c("scanned", "all"),
                     labels = c("scanned locations", "all locations")) +  
  theme_classic() +
  theme(legend.position = c(.6, .2), plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Species Evenness", x = "elevation (m)", y = "Species Evenness Index")

plot_grid(p1,p2,p3, nrow = 1)


### MODELING
library(lme4)
library(glmmTMB)
library(coefplot)
library(coefplot2)
library(PerformanceAnalytics)
library(fitdistrplus)

#checking correlation between predictors first
chart.Correlation(stand_mets[,c("max.height.sc","sd.r.sc","CRR.rho.sc","perc.below.2m.sc","stand.dens.sc",
                                "basal.area.sc","stem.vol.sc","mean.tree.h.sc","mean.dbh.sc",
                                "n.all", "shannon_ct", "div_even")],
                  method = "pearson")

#checking distribution family for outcome vars
plotdist(stand_mets$n.all, histo = TRUE, demp = TRUE)
plotdist(stand_mets$shannon_ct, histo = TRUE, demp = TRUE)
plotdist(stand_mets$div_even, histo = TRUE, demp = TRUE)

descdist(stand_mets$n.all, boot = 100)
descdist(stand_mets$shannon_ct, boot = 100)
descdist(stand_mets$div_even, boot = 100)


## GLMM
#species richness as outcome variable - by CT location
richness_model <- glm(n.all ~ 
                 max.height.sc +
                 mean.tree.h.sc +
                 CRR.rho.sc +
                 sd.r.sc +
                 mean.dbh.sc +
                 basal.area.sc +
                 stand.dens.sc +
                 stem.vol.sc +
                 perc.below.2m.sc,
               family = Gamma(link = "log"), 
               offset = log(act.days),
               data = stand_mets)

summary(richness_model)
coefplot2(richness_model, top.axis = FALSE, main = "Species Richness",
          cex.pts = 1.5, lwd.1 = 4, lwd.2 = 2)

#species diversity as outcome variable - by CT location
shannon_model <- glm(shannon_ct ~ 
                         max.height.sc +
                         mean.tree.h.sc +
                         CRR.rho.sc +
                         sd.r.sc +
                         mean.dbh.sc +
                         basal.area.sc +
                         stand.dens.sc +
                         stem.vol.sc +
                         perc.below.2m.sc,
                        family = Gamma(link = "log"), 
                        offset = log(act.days),
                        data = stand_mets)

summary(shannon_model)
coefplot2(shannon_model, top.axis = FALSE, main = "Shannon Diversity",
          cex.pts = 1.5, lwd.1 = 4, lwd.2 = 2)

#species evenness as outcome variable - by CT location
even_model <- glm(div_even ~ 
                max.height.sc +
                mean.tree.h.sc +
                CRR.rho.sc +
                sd.r.sc +
                mean.dbh.sc +
                basal.area.sc +
                stand.dens.sc +
                stem.vol.sc +
                perc.below.2m.sc,
              family = Gamma(link = "log"), 
              offset = log(act.days),
              data = stand_mets)

summary(even_model)
coefplot2(even_model, top.axis = FALSE, main = "Species Evenness",
          cex.pts = 1.5, lwd.1 = 4, lwd.2 = 2)

#compare coefficients from all plots
multiplot(richness_model, shannon_model, even_model, intercept = FALSE) + 
  labs(title = "", x = "coefficient estimate", y = "") +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "horizontal")
