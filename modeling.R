### LiDAR-derived Forest Structure Metrics as Predictors of Mammal Species Richness and Diversity

library(tidyverse)
library(vegan)

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

#adding species richness and diversity data to stand metrics table
stand_mets <- left_join(stand_mets, n_by_ct)
stand_mets <- left_join(stand_mets, shannon_ct)

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

#scale and center all numeric predictors 
stand_mets <- stand_mets %>%
  mutate(across(8:19, scale, .names = "{.col}.sc"))


## Richness and diversity visualization
ct <- ct %>%
  rename(locationID   = Deployment.Location.ID,
         habitat      = Treatment.Strata,
         location     = Location,
         on.off.trail = Feature.Type)

#plotting species richness by CT location/forest type - all CT locations only
n_by_ct %>%
  left_join(., ct[,c("locationID", "habitat")]) %>%
  mutate(habitat = factor(habitat, levels = c("Peat Swamp", "Freshwater Swamp", "Alluvial Bench", "Lowland Sandstone", 
                                              "Lowland Granite", "Upland Granite", "Montane"))) %>%
  ggplot(aes(x = habitat, y = n.all)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_flip() +
  theme_classic() +
  labs(title = "Species Richness by CT - all locations", x = "", y = "n species")

#plotting species richness by CT location/forest type - scanned sites only
ggplot(stand_mets, aes(x = habitat, y = n.all)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_flip() +
  theme_classic() +
  labs(title = "Species Richness by CT - scanned sites", x = "", y = "n species")

#plotting species diversity by CT location/forest type - all CT locations
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
  labs(title = "Species Diversity by CT - all locations", x = "", y = "n species")

#plotting species diversity by CT location/forest type - scanned sites only
ggplot(stand_mets, aes(x = habitat, y = shannon_ct)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_flip() +
  theme_classic() +
  ylim(0, 3) +
  labs(title = "Species Diversity by CT", x = "", y = "Shannon Diversity Index")

#species richness by elevation
ct_elev <- read.csv(file = "data/cameradata_updatedZJ-ajm.csv", header = TRUE)
ct_elev <- ct_elev[!duplicated(ct_elev$locationID),]

n_by_ct %>%
  left_join(., ct_elev[,c("locationID", "altitude")]) %>%
  ggplot(aes(x = altitude, y = n.all)) +
  geom_point() +
  geom_smooth(se = TRUE) +
  theme_classic() +
  labs(title = "Species Richness by Elevation", x = "elevation (m)", y = "n species")

#species diversity by elevation
shannon_ct %>%
  left_join(., ct_elev[,c("locationID", "altitude")]) %>%
  ggplot(aes(x = altitude, y = shannon_ct)) +
  geom_point() +
  geom_smooth(se = TRUE) +
  theme_classic() +
  labs(title = "Species Diversity by Elevation", x = "elevation (m)", y = "Shannon Diversity Index")



### MODELING
library(lme4)
library(coefplot2)

## GLMM
#species richness as outcome variable - by CT location
richness_model <- glmer(n.all ~ 
                 max.height.sc +
                 sd.r.sc +
                 CRR.rho.sc +
                 perc.below.2m.sc +
                 basal.area.sc +
                 mean.tree.h.sc +
                 stem.vol.sc +
                 #habitat +
                 (1|habitat/locationID), 
               family = poisson(link = "log"), 
               offset = log(act.days),
               data = stand_mets)

summary(richness_model)
coefplot2(richness_model, top.axis = FALSE, main = "Species Richness",
          cex.pts = 1.5, lwd.1 = 4, lwd.2 = 2)

#species diversity as outcome variable - by CT location
shannon_model <- glmer(shannon_ct ~ 
                          max.height.sc +
                          sd.r.sc +
                          CRR.rho.sc +
                          perc.below.2m.sc +
                          basal.area.sc +
                          mean.tree.h.sc +
                          (1|habitat), 
                        family = gaussian(link = "identity"), 
                        offset = log(act.days),
                        data = stand_mets)

summary(shannon_model)
coefplot2(shannon_model, top.axis = FALSE, main = "Species Diversity (Shannon)",
          cex.pts = 1.5, lwd.1 = 4, lwd.2 = 2)

#maybe next I can try some GAM models and/or try to model individual species that are common?
