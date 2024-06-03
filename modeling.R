### LiDAR-derived Forest Structure Metrics as Predictors of Mammal Species Richness and Diversity

library(tidyverse)
library(vegan)
library(cowplot)

#mammal observation data - excludes all human observations
m <- read.csv(file = "data/pruned.csv", header = TRUE)
m <- m[,-1]

#adding updated CT location forest type and partition data
ct <- read.csv(file = "data/CTlocations_partitions_May2024.csv", header = TRUE)

#removing lat, long, and habitat designations and adding updated versions
m <- m[,-c(12:14)]
m <- left_join(m, ct[,c("locationID","latitude","longitude","habitat","partition")], by = "locationID")

#ordering forest type
m$habitat <- factor(m$habitat, levels = c("Peat Swamp", "Freshwater Swamp", "Alluvial Bench", 
                                          "Lowland Sandstone", "Lowland Granite", "Upland Granite", "Montane"),
                    ordered = TRUE)
#ordering partitions
m$partition <- factor(m$partition, levels = c("PS1", "FS1", "AB1", "AB2",
                                          "LS1", "LS2","LG1","LG2", "UG1","UG2", "MO1", "MO2"),
                    ordered = TRUE)

#mammal observation table
species_table <- as.data.frame(table(m$Species, m$habitat))
species_table <- pivot_wider(species_table, names_from = "Var2", values_from = "Freq")
colnames(species_table)[colnames(species_table) == 'Var1'] <- 'Species'
#write.csv(species_table, file = "species_observation_table_byFT.csv")

#LiDAR-derived metrics data
stand_mets <- read.csv(file = "data/select.stand.mets.csv", header = TRUE)
stand_mets$habitat <- factor(stand_mets$habitat, levels = c("Peat Swamp", "Freshwater Swamp", "Alluvial Bench", 
                                          "Lowland Sandstone", "Lowland Granite", "Upland Granite", "Montane"),
                             ordered = TRUE)
stand_mets$partition <- factor(stand_mets$partition, levels = c("PS1", "FS1", "AB1", "AB2",
                                              "LS1", "LS2","LG1","LG2", "UG1","UG2", "MO1", "MO2"),
                      ordered = TRUE)

#all tree metrics
tree_mets <- read.csv(file = "data/all.tree.mets.csv", header = TRUE)
tree_mets$habitat <- factor(tree_mets$habitat, levels = c("Peat Swamp", "Freshwater Swamp", "Alluvial Bench", 
                                                            "Lowland Sandstone", "Lowland Granite",
                                                          "Upland Granite", "Montane"),
                            ordered = TRUE)
tree_mets$partition <- factor(tree_mets$partition, levels = c("PS1", "FS1", "AB1", "AB2",
                                                                "LS1", "LS2","LG1","LG2", "UG1","UG2", "MO1", "MO2"),
                               ordered = TRUE)
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
  arrange(desc(habitat)) %>%
  print()

#only scanned locations
n_by_ft <- m %>%
  filter(locationID %in% stand_mets$locationID) %>%
  group_by(habitat) %>%
  filter(!(Species %in% c("Confused Squirrel or Treeshrew","Tragulus spp.","Tupai sp.","Unid bats",
                          "Unid civet","Unid rat","Unid squirrel"))) %>%
  summarise(n.scanned = n_distinct(Species)) %>%
  left_join(n_by_ft) %>%
  arrange(desc(habitat)) %>% 
  print()

#the scanned sites consistently observe fewer species over the study period by 4 - 10 species
#more bias in the higher elevation FT's, especially the montane
n_by_ft$diff <- n_by_ft$n.all - n_by_ft$n.scanned
n_by_ft

## By partition - excluding unid animals
#all ct locations
n_by_pt <- m %>%
  group_by(partition) %>%
  filter(!(Species %in% c("Confused Squirrel or Treeshrew","Tragulus spp.","Tupai sp.","Unid bats",
                          "Unid civet","Unid rat","Unid squirrel"))) %>%
  summarise(n.all = n_distinct(Species)) %>%
  arrange(desc(partition)) %>%
  print()

#only scanned locations
n_by_pt <- m %>%
  filter(locationID %in% stand_mets$locationID) %>%
  group_by(partition) %>%
  filter(!(Species %in% c("Confused Squirrel or Treeshrew","Tragulus spp.","Tupai sp.","Unid bats",
                          "Unid civet","Unid rat","Unid squirrel"))) %>%
  summarise(n.scanned = n_distinct(Species)) %>%
  left_join(n_by_pt) %>%
  arrange(desc(partition)) %>% 
  print()

#the scanned sites consistently observe fewer species but this is especially true in UG1 and LS2
n_by_pt$diff <- n_by_pt$n.all - n_by_pt$n.scanned
n_by_pt

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

## by partition
shannon_all <- diversity(table(m$partition, m$Species), index = "shannon")
shannon_all <- rownames_to_column(as.data.frame(shannon_all), var = "partition")

shannon_scan <- diversity(table(m$partition[shannon_ct$locationID %in% stand_mets$locationID], 
                                m$Species[shannon_ct$locationID %in% stand_mets$locationID]), index = "shannon")
shannon_scan <- rownames_to_column(as.data.frame(shannon_scan), var = "partition")
shannon_pt <- left_join(shannon_all, shannon_scan)
rm(shannon_scan, shannon_all)

#again, using only the CT's that were scanned leads to most bias in the montane forest
shannon_pt$diff <- shannon_pt$shannon_all - shannon_pt$shannon_scan
shannon_pt

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
ct2 <- read.csv(file = "data/ofp_deployments-2021-11-04.csv", header = TRUE)
ct_active_days <- ct2 %>%
  dplyr::select(Deployment.Location.ID, Camera.Deployment.Begin.Date, Camera.Deployment.End.Date) %>%
  mutate(Camera.Deployment.Begin.Date = as.Date(Camera.Deployment.Begin.Date),
         Camera.Deployment.End.Date   = as.Date(Camera.Deployment.End.Date)) %>%
  mutate(act.per = as.numeric(difftime(Camera.Deployment.End.Date, 
                                       Camera.Deployment.Begin.Date, units = "days"))) %>%
  group_by(Deployment.Location.ID) %>%
  summarise(act.days = sum(act.per)) %>%
  rename(locationID = Deployment.Location.ID)

stand_mets <- left_join(stand_mets, ct_active_days)  
stand_mets <- left_join(stand_mets, div_mets[,c("locationID","div_even")])
stand_mets$div_even <- as.numeric(stand_mets$div_even)


#scale and center all numeric predictors 
stand_mets <- stand_mets %>% mutate(across(c(5:16, 21:23, 25:27, 29), scale, .names = "{.col}.sc"))


## Richness and diversity visualization

#add habitat data to diversity metrics table
div_mets <- left_join(div_mets, ct[,c("locationID","habitat","partition","altitude")])
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

#plotting species richness by CT location/forest type - all locations
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
library(coefplot)
library(coefplot2)
library(PerformanceAnalytics)
library(fitdistrplus)
library(MuMIn)
library(sjPlot)

#checking correlation between predictors first
#basal area is highly correlated (>0.8) with stem vol and stand dens
chart.Correlation(stand_mets[,c("max.height.sc","sd.r.sc","CRR.rho.sc","perc.below.2m.sc","stand.dens.sc",
                                "basal.area.sc","stem.vol.sc","mean.tree.h.sc","mean.dbh.sc", "zentropy.sc",
                                "lad.max.sc","rumple.sc","vFRcanopy.sc","vzrumple.sc","ClosedGapSpace.sc")],
                  method = "spearman")

#checking distribution family for diversity outcome vars
plotdist(stand_mets$n.all, histo = TRUE, demp = TRUE)
plotdist(stand_mets$shannon_ct, histo = TRUE, demp = TRUE)
plotdist(stand_mets$div_even, histo = TRUE, demp = TRUE)

descdist(stand_mets$n.all, boot = 100)
descdist(stand_mets$shannon_ct, boot = 1000)
descdist(stand_mets$div_even, boot = 100)


## GLM
#species richness as outcome variable - by CT location
#fitting global model for use in dredge
rm1 <- glm(n.all ~  
             max.height.sc +
             mean.tree.h.sc +
             CRR.rho.sc +
             sd.r.sc +
             mean.dbh.sc +
             #basal.area.sc + #leaving out due to high corr - model doesn't converge
             stand.dens.sc +
             stem.vol.sc +
             pts.below.2m.sc +
             zentropy.sc +
             lad.max.sc +
             rumple.sc +
             vFRcanopy.sc +
             vzrumple.sc +
             ClosedGapSpace.sc,
           na.action = na.fail,
           family = Gamma(link = "log"), 
           offset = log(act.days),
           data = stand_mets)
#dredge
rich_models <- dredge(rm1)

#averaging models with dAIC<2 - this works but idk how to plot coefficients
#https://stackoverflow.com/questions/54962119/how-to-plot-from-mumin-model-avg-summary 
rich_models_avg <- model.avg(rich_models, subset = "delta" < 2)
summary(rich_models_avg)

#selecting top model by AIC
top_rich_model <- get.models(rich_models, subset = 1)

summary(top_rich_model[[1]])
coefplot2(top_rich_model[[1]], top.axis = FALSE, main = "Species Richness",
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
                         pts.below.2m.sc +
                         zentropy.sc +
                         lad.max.sc +
                         rumple.sc +
                         vFRcanopy.sc +
                         vzrumple.sc +
                         ClosedGapSpace.sc,
                        na.action = na.fail,
                        family = Gamma(link = "log"), 
                        offset = log(act.days),
                        data = stand_mets)

#dredge
shan_models <- dredge(shannon_model)

#averaging models with dAIC<2
shan_models_avg <- model.avg(shan_models, subset = "delta" < 2)
summary(shan_models_avg)

#selecting top model by AIC
top_shan_model <- get.models(shan_models, subset = 1)

summary(top_shan_model[[1]])
coefplot2(top_shan_model[[1]], top.axis = FALSE, main = "Shannon Diversity",
          cex.pts = 1.5, lwd.1 = 4, lwd.2 = 2)

#top 6 models show dAIC<2

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
                pts.below.2m.sc +
                zentropy.sc +
                lad.max.sc +
                rumple.sc +
                vFRcanopy.sc +
                vzrumple.sc +
                ClosedGapSpace.sc,
                na.action = na.fail,
              family = Gamma(link = "log"), 
              offset = log(act.days),
              data = stand_mets)

#dredge
even_models <- dredge(even_model)

#averaging models with dAIC<2
even_models_avg <- model.avg(even_models, subset = "delta" < 2)
summary(even_models_avg)

#selecting top model by AIC
top_even_model <- get.models(even_models, subset = 1)

summary(top_even_model[[1]])
coefplot2(top_even_model[[1]], top.axis = FALSE, main = "Species Evenness",
          cex.pts = 1.5, lwd.1 = 4, lwd.2 = 2)

#compare coefficients from all plots
multiplot(top_rich_model[[1]], top_shan_model[[1]], top_even_model[[1]], intercept = FALSE, 
          decreasing = TRUE, names = c("richness","diversity","evenness")) + 
  labs(title = "", x = "coefficient estimate", y = "") +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  scale_y_discrete(labels = rev(c("tree dens.","vert. complexity","leaf area dens.",
                            "rugosity","veg. vol.","max. h","gap vol.")))


#compare summary tables from all top models
tab_model(top_rich_model[[1]], top_shan_model[[1]], top_even_model[[1]], transform = NULL,
          #pred.labels = c("Intercept","Max. Tree h","SD tree h","Tree dens.","CRR"), #fix names
          dv.labels = c("Species Richness","Shannon Diversity","Species Evenness")
          #file = "table1.xls"                                                        #export as Excel file for manuscript submission
          )


############################################################################################################

## Census (line transect) data 

#importing census summary data - creating diversity metrics to compare to the CT observations

library(tidyverse)
library(vegan)
library(cowplot)

#reading in unique census sightings table
c <- read.csv(file = "data/all_uniquecensus sightings.csv", header = TRUE)
head(c)
#some animal ID's are in both upper and lower case. Also removing white spaces before and after characters
c$animal <- toupper(c$animal)
c$animal <- trimws(c$animal)

#3,976 unique CID's
table(c$cid)
length(table(c$cid))

#every id is unique to the row (22,404)
table(c$id)
length(table(c$id))

#reading in census segments table
ca <- read.csv(file = "data/census_analysis_all_segments.csv", header = TRUE)
head(ca)

#looking at 'meta_id' column - there are 4,081 unique ID's (is this the CID from the 'c' table?)
table(ca$meta_id)
length(table(ca$meta_id))

#every CID from the c table matches a value in the meta_id column
table(c$cid %in% ca$meta_id)

#how many census routes are there?
table(ca$census)
length(table(ca$census))  #34

#there are some segments that don't have a forest type or partition designation
table(ca$habitat)
table(ca$partishun)

#the routes w/o hab or part designations are 15, 16, and 17 A and B
ca[ca$habitat == "NULL",]
table(ca$census[ca$habitat == "NULL"])

#assigning ft and partition data to the observation table
p <- read.csv("data/cptrails_all-ajm.csv", header = TRUE)
#11 segments without habitat data - these are the segments labeled AB 17 - AB 999, for example
table(p$habitat)
table(p$partition)

#subsetting only segment, hab, and part columns
p <- p[,2:5]
#unique_values <- unique(as.data.frame(c(p$pointfrom, p$pointto)))

#creating list of trail names in main study area
trails <- separate(p, pointfrom, into = c("trail","num"), sep = "-", remove = FALSE) %>%
  select("trail") %>%
  unique() %>%
  pull(trail)

#adding two trails manually - not positive if all of AG and BR are freshwater swamp though
trails <- c(trails, "AG", "BR")
table(c$trail %in% trails)

#formatting the trail and trailmarker columns to match those in the partitions table
c$marker <- as.integer(c$trailmark)
c$marker <- formatC(c$marker, width = 3, format = "d", flag = "0")
c$pointfrom <-  paste(c$trail, c$marker, sep = "-")


#adding 1 to trail marker list
#splitting pointfrom column into trail and marker columns first
p <- p[,-2] %>%
  mutate(trail_name = substr(pointfrom, 1, 2),
         trail_number = as.numeric(substr(pointfrom, 4, 6)))

# Find the maximum trail number for each trail and create a new row with the incremented number
new_rows <- p %>%
  group_by(trail_name) %>%
  summarise(trail_number = max(trail_number) + 1,
            habitat = last(habitat),
            partition = last(partition)) %>%
  mutate(pointfrom = paste0(trail_name, "-", sprintf("%03d", trail_number)))

# Combine the new rows with the original data frame
p <- p %>%
  select(pointfrom, habitat, partition) %>%
  bind_rows(new_rows %>%
              select(pointfrom, habitat, partition)) %>%
  arrange(pointfrom)

rm(new_rows)

#joining and filtering for trails only in the study area
c <- left_join(c, p, by = "pointfrom")

#all of the non-matching trails are rangkong trails except AG and BR
#I might also think about imputing the forest type and partitions where the trail markers are NULL
table(c$trail[!c$trail %in% trails]) #some trails are listed as NULL

#filtering out rangkong observations
c <- c[c$trail %in% trails,]

#adding forest type and partition data for AG and BR
c$habitat[c$trail == "AG" | c$trail == "BR"] <- "FS"
c$partition[c$trail == "AG" | c$trail == "BR"] <- "FS.I"

#there are still 31 NA's - might need to look closer at these to fill them out
table(is.na(c$habitat))
table(c$trail[is.na(c$habitat)])
c[is.na(c$habitat),]

#adding species info from singkat file
singkat <- read.csv(file = "data/vert_singkat-2021-04-08.csv", header = TRUE)
c <- left_join(c, singkat[,c("New_Singkat","Latin_name","class")], join_by("animal" == "New_Singkat"))

#only including mammals - excluding humans and bats
#NULL or empty values in 'Latin_name': PR = Unid primate, TX = 'confused squirrel or treeshrew'
c <- c[which(c$class == "Mammal" & !c$Latin_name %in% c("Homo sapiens","","NULL","Unid bat")),]

#39 mammal species in the data set (maybe less, not sure if I should ignore "Tragulus spp")
length(unique(c$Latin_name))

#how many observations of each species?
as.data.frame(table(c$Latin_name)) %>%
  arrange(desc(Freq))


##########################################################################################################
#don't run this chunk again
#add arboreal/terrestrial data from PanTheria

#loading pantheria dataset (Jones et al. 2009)
#ref: https://esajournals.onlinelibrary.wiley.com/doi/10.1890/08-1494.1
library(traitdata)
data(pantheria)
terrest <- pantheria[c("Family", "Genus", "scientificNameStd","Terrestriality")]

#most species have same latin name between census and pantheria data, 12 don't
table(unique(c$Latin_name) %in% terrest$scientificNameStd)


#reading in existing trait data from Leech lab paper - species are from ct observations
traits <- read.csv(file = "data/mammal_list_ct_trait_202403.csv", header = TRUE)
traits$Terrestriality <- as.integer(traits$Terrestriality)
colnames(traits)[4] <- "scientificNameStd"

#10 species from the census data are not in the ct data
table(unique(c$Latin_name) %in% traits$Species)

new_sp <- as.data.frame(unique(c$Latin_name[!c$Latin_name %in% traits$Species]))
colnames(new_sp)[1] <- "Species"
new_sp$scientificNameStd <- new_sp$Species
new_sp$scientificNameStd[new_sp$scientificNameStd == "Cervus unicolor"] <- "Rusa unicolor"

new_sp <- left_join(new_sp, pantheria[,c("Family","Genus", "scientificNameStd", "ActivityCycle",
                                         "AdultBodyMass_g","Terrestriality", "TrophicLevel")])

traits <- bind_rows(traits, new_sp)

#still missing some terrestriality data for species, will have to fill it in on my own
traits[is.na(traits$Terrestriality) ,c("Species","Terrestriality")]

#write.csv(traits, file = "data/mammal_list_ct_trait_202405.csv")
############################################################################################################


#adding terrestriality data to observation table
traits <- read.csv(file = "data/mammal_list_ct_trait_202405.csv", header = TRUE)
c <- left_join(c, traits[,c("Species","Terrestriality")], join_by("Latin_name" == "Species"))

#no missing data here, and most observations (85%) are of arboreal mammals
table(is.na(c$Terrestriality))
table(c$Terrestriality)

#grouping observations by forest type and partition across study period
#can't do this by census tract since many cover multiple forest types, might not make sense for my analyses
obs_tab_ft <- c %>% 
  group_by(habitat, Latin_name) %>%
  summarise(n = sum(as.numeric(n_indiv))) %>%
  pivot_wider(names_from = Latin_name, values_from = n) %>%
  filter(habitat != "", !is.na(habitat))  #removing the rows associated with no hab designation 

obs_tab_ft$habitat <- factor(obs_tab_ft$habitat, 
                           levels = c("PS", "FS","AB","LS","LG","UG", "MO"),
                             ordered = TRUE)

obs_tab_pt <- c %>% 
  group_by(partition, Latin_name) %>%
  summarise(n = sum(as.numeric(n_indiv))) %>%
  pivot_wider(names_from = 'Latin_name', values_from = 'n') %>%
  filter(partition != "", !is.na(partition))

obs_tab_pt$partition <- factor(obs_tab_pt$partition, 
                           levels = c("PS.I","FS.I","AB.I","AB.II","LS.I","LS.II",
                                      "LG.I","LG.II","UG.I","UG.II","MO.I","MO.II","MO.III"),
                           ordered = TRUE)

#########################################################################################################
### species richness

#by entire study area - 39 unique mammal species observed
length(obs_tab_ft[,-1])

#by forest type
#adding survey effort - STILL NEED TO FIX HOW 15A & B SURVEY EFFORT IS FACTORED IN HERE
obs_tab_ft <- ca %>%
  group_by(habitat) %>%
  summarize(effort = sum(actual_effort)) %>%
  left_join(obs_tab_ft, ., by = "habitat")

div_tab_ft <- obs_tab_ft %>%
  rowwise() %>%
  mutate(richness = sum(!is.na(c_across(-1)))) %>%
  select(habitat = 1, richness) %>%
  left_join(., obs_tab_ft[,c("habitat","effort")]) %>%
  mutate(rich_std = richness/effort) %>%
  arrange(desc(habitat))


#by partition
#adding survey effort
obs_tab_pt <- ca %>%
  group_by(partishun) %>%
  summarize(effort = sum(actual_effort)) %>%
  left_join(obs_tab_pt, ., join_by("partition" == "partishun"))

div_tab_pt <- obs_tab_pt %>%
  rowwise() %>%
  mutate(richness = sum(!is.na(c_across(-1)))) %>%
  select(partition = 1, richness) %>%
  left_join(., obs_tab_pt[,c("partition","effort")]) %>%
  mutate(rich_std = richness/effort) %>%
  arrange(desc(partition))

div_tab_pt_terr <- obs_tab_pt %>%
  select(partition, which(colnames(obs_tab_pt) %in% traits$Species[traits$Terrestriality == 1])) %>%
  rowwise() %>%
  mutate(richness = sum(!is.na(c_across(-1)))) %>%
  select(partition = 1, richness) %>%
  left_join(., obs_tab_pt[,c("partition","effort")]) %>%
  mutate(rich_std = richness/effort) %>%
  arrange(desc(partition))

div_tab_pt_arb <- obs_tab_pt %>%
  select(partition, which(colnames(obs_tab_pt) %in% traits$Species[traits$Terrestriality == 2])) %>%
  rowwise() %>%
  mutate(richness = sum(!is.na(c_across(-1)))) %>%
  select(partition = 1, richness) %>%
  left_join(., obs_tab_pt[,c("partition","effort")]) %>%
  mutate(rich_std = richness/effort) %>%
  arrange(desc(partition))

### shannon diversity

#by entire study area - not controlling for survey effort
diversity(colSums(obs_tab_ft[,setdiff(names(obs_tab_ft), c("habitat", "effort"))], na.rm = TRUE), 
          index = "shannon")

#by forest type
obs_tab_ft_std <- sweep(obs_tab_ft[,-1], 1, obs_tab_ft$effort, FUN = "/")
rownames(obs_tab_ft_std) <- obs_tab_ft$habitat
obs_tab_ft_std[is.na(obs_tab_ft_std)] <- 0

shan_by_ft <- as.data.frame(diversity(obs_tab_ft_std[,-40], index = "shannon"))
colnames(shan_by_ft)[1] <- 'shannon'
shan_by_ft <- rownames_to_column(shan_by_ft, var = "habitat")

div_tab_ft <- left_join(div_tab_ft, shan_by_ft, by = "habitat")
rm(shan_by_ft)

#by partition
obs_tab_pt_std <- sweep(obs_tab_pt[,-1], 1, obs_tab_pt$effort, FUN = "/")
rownames(obs_tab_pt_std) <- obs_tab_pt$partition
obs_tab_pt_std[is.na(obs_tab_pt_std)] <- 0

shan_by_pt <- as.data.frame(diversity(obs_tab_pt_std, index = "shannon"))
colnames(shan_by_pt)[colnames(shan_by_pt) == 'diversity(obs_tab_pt_std, index = "shannon")'] <- 'shannon'
shan_by_pt <- rownames_to_column(shan_by_pt, var = "partition")

div_tab_pt <- left_join(div_tab_pt, shan_by_pt, by = "partition")
rm(shan_by_pt)

#by partition - terrestrial species only
shan_by_pt <- as.data.frame(diversity(obs_tab_pt_std[which(colnames(obs_tab_pt) %in% traits$Species[traits$Terrestriality == 1])], index = "shannon"))
colnames(shan_by_pt)[1] <- 'shannon'
shan_by_pt <- rownames_to_column(shan_by_pt, var = "partition")

div_tab_pt_terr <- left_join(div_tab_pt_terr, shan_by_pt, by = "partition")
rm(shan_by_pt)

#by partition - arboreal species only
shan_by_pt <- as.data.frame(diversity(obs_tab_pt_std[which(colnames(obs_tab_pt) %in% traits$Species[traits$Terrestriality == 2])], index = "shannon"))
colnames(shan_by_pt)[1] <- 'shannon'
shan_by_pt <- rownames_to_column(shan_by_pt, var = "partition")

div_tab_pt_arb <- left_join(div_tab_pt_arb, shan_by_pt, by = "partition")
rm(shan_by_pt)


### Pielou's diversity measure of species evenness

#by forest type
div_tab_ft$evenness <- div_tab_ft$shannon / log(div_tab_ft$richness)

#by partition - all 
div_tab_pt$evenness <- div_tab_pt$shannon / log(div_tab_pt$richness)

#by partition - terrestrial
div_tab_pt_terr$evenness <- div_tab_pt_terr$shannon / log(div_tab_pt_terr$richness)

#by partition - arboreal
div_tab_pt_arb$evenness <- div_tab_pt_arb$shannon / log(div_tab_pt_arb$richness)

#ordering forest types and partitions
div_tab_ft$habitat <- factor(div_tab_ft$habitat, 
                             levels = c("PS", "FS","AB","LS","LG","UG", "MO"),
                             ordered = TRUE)

div_tab_pt$partition <- factor(div_tab_pt$partition, 
                               levels = c("PS.I","FS.I","AB.I","AB.II","LS.I","LS.II",
                                          "LG.I","LG.II","UG.I","UG.II","MO.I","MO.II","MO.III"),
                               ordered = TRUE)

div_tab_pt_terr$partition <- factor(div_tab_pt_terr$partition, 
                               levels = c("PS.I","FS.I","AB.I","AB.II","LS.I","LS.II",
                                          "LG.I","LG.II","UG.I","UG.II","MO.I","MO.II","MO.III"),
                               ordered = TRUE)

div_tab_pt_arb$partition <- factor(div_tab_pt_arb$partition, 
                                    levels = c("PS.I","FS.I","AB.I","AB.II","LS.I","LS.II",
                                               "LG.I","LG.II","UG.I","UG.II","MO.I","MO.II","MO.III"),
                                    ordered = TRUE)

### Diversity metrics visualization 

#species richness by ft
#standardized
ggplot(div_tab_ft, aes(habitat, rich_std)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Mammal Richness by Forest Type", subtitle = "all census data", 
       x = "", y = "n species / survey effort")

#raw counts
ggplot(div_tab_ft, aes(habitat, richness)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Mammal Richness by Forest Type - Raw Counts", subtitle = "all census data", 
       x = "", y = "n species / survey effort")

#species richness by pt
ggplot(div_tab_pt, aes(partition, rich_std)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Mammal Richness by Partition", subtitle = "all census data", 
       x = "", y = "n species / survey effort")

#species richness by pt - terrestrial only
ggplot(div_tab_pt_terr, aes(partition, rich_std)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Terrestrial Mammal Richness by Partition", 
       x = "", y = "n species / survey effort")

#species richness by pt - arboreal only
ggplot(div_tab_pt_arb, aes(partition, rich_std)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Arboreal Mammal Richness by Partition", 
       x = "", y = "n species / survey effort")


#species diversity by ft
ggplot(div_tab_ft, aes(habitat, shannon)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Mammal Diversity by Forest Type", subtitle = "all census data", 
       x = "", y = "Shannon diversity index / survey effort")

#species diversity by pt
ggplot(div_tab_pt, aes(partition, shannon)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Mammal Diversity by Partition", subtitle = "all census data", 
       x = "", y = "Shannon diversity index / survey effort")

#species diversity by pt - terrestrial
ggplot(div_tab_pt_terr, aes(partition, shannon)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Terrestrial Mammal Diversity by Partition", 
       x = "", y = "Shannon diversity index / survey effort")

#species diversity by pt - arboreal - THIS DOESN'T LOOK RIGHT
ggplot(div_tab_pt_arb, aes(partition, shannon)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Arboreal Mammal Diversity by Partition", 
       x = "", y = "Shannon diversity index / survey effort")

########## LEFT OFF HERE

#species evenness by ft
ggplot(div_tab_ft, aes(habitat, evenness)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Mammal Evenness by Forest Type", subtitle = "all census data", 
       x = "", y = "Pielou index / survey effort")

#species evenness by pt
ggplot(div_tab_pt, aes(partition, evenness)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Mammal Evenness by Partition", subtitle = "all census data", 
       x = "", y = "Pielou index / survey effort")

