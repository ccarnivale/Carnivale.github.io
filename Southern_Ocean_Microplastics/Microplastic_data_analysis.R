#This is the script to analyze the Plastics Polymers that have been enumerated
#Throughout the water column

#First env setup for data cleaning and analysis

# ****NOTE 1 ABS was found in PS15 airblank*****
# NEEDS TO BE REMOVED

library(ggplot2)
library(magrittr)
library(tidyverse)
library(ggeffects)
library(gridExtra)
library(grid)
library(vegan)
library(ggpubr)
library(lme4)
library(ggrepel)
library(reshape)
#Data Prepping--------
setwd("/Users/christophercarnivale/Desktop/Dissertation_data/Raman Spectral Data copy/")

set.seed(444)

plastic_df <- read.csv("Env_Data_CC_wPlastic_newCF_blankCor.csv", header = TRUE)



plastic_df_calc <- mutate(plastic_df,
                          ABS_L = ABS * Correction_factor/Volume_filt,
                          PET_L = PET.1 * Correction_factor/Volume_filt,
                          PE_L = PE * Correction_factor/Volume_filt,
                          PP_L = PP * Correction_factor/Volume_filt,
                          PA_L = PA * Correction_factor/Volume_filt,
                          total_L = Total_plastic * Correction_factor/Volume_filt)


plastic_df_calc$depth <- factor(plastic_df_calc$depth,
                                levels =c ("Surface","14m","15m","23m","40m","70m","110m","120m","150m","205m","250m","Deep"))


plastic_df_calc$Size_fraction <- factor(plastic_df_calc$Size_fraction,
                                levels =c ("<20",">20",">100"))

plastic_df_calc$station <- factor(plastic_df_calc$station,
                                        levels =c ("R","PS15"))

ggplot(plastic_df_calc)+
  geom_col(aes(x = station, y = total_L, fill = station))


ggplot(plastic_df_calc)+
  geom_col(aes(x = depth_m, y = total_L, fill = station))

plastic_df_calc_grp <- group_by(plastic_df_calc, station, depth)

plastic_df_calc_grp_sum <- summarise(plastic_df_calc_grp,
                                     ABS_sum = sum(ABS_L),
                                     PET_sum = sum(PET_L),
                                     PE_sum = sum(PE_L),
                                     PP_sum = sum(PP_L),
                                     PA_sum = sum(PA_L),
                                     total_sum = sum(total_L))

#I need to split this for ease of calculation
plastic_df_calc_StnR <- dplyr::filter(plastic_df_calc, station == "R")

#plastic_df_calc_StnR$depth <- factor(plastic_df_calc_StnR$depth,
#                                     levels = c("Surface","14m","23m","70m","110m","Deep"))

plastic_df_calc_StnPS15 <- dplyr::filter(plastic_df_calc, station == "PS15")

#plastic_df_calc_StnPS15$depth <- factor(plastic_df_calc_StnPS15$depth,
#                                     levels =c ("Surface","15m","40m","120m","150m","205m","250m","Deep"))

#Initial EDA ----------
#Set theme for publication and cleaner presentations
theme_set(theme_bw())
theme_update(element_line(linewidth = 0), panel.grid.major = element_blank())


ggplot(plastic_df_calc_StnR)+
  geom_col(aes(x = depth, y = total_L, fill = station))


ggplot(plastic_df_calc_StnPS15)+
  geom_col(aes(x = depth, y = total_L, fill = station))

ggplot(plastic_df_calc)+
  geom_col(aes(x = depth, y = total_L, fill = station), position = position_dodge())

ggplot(plastic_df_calc)+
  geom_col(aes(x = depth, y = PE_L, fill = station), position = position_dodge())

ggplot(plastic_df_calc)+
  geom_col(aes(x = depth, y = ABS_L, fill = station), position = position_dodge())

ggplot(plastic_df_calc)+
  geom_col(aes(x = depth, y = PET_L, fill = station), position = position_dodge())

ggplot(plastic_df_calc)+
  geom_col(aes(x = depth, y = PP_L, fill = station), position = position_dodge())

ggplot(plastic_df_calc)+
  geom_col(aes(x = depth, y = PA_L, fill = station), position = position_dodge())

ggplot(plastic_df_calc)+
  geom_col(aes(x = Size_fraction, y = total_L, fill = station), position = position_dodge())

ggplot(plastic_df_calc)+
  geom_col(aes(x = Size_fraction, y = PP_L, fill = station), position = position_dodge())

ggplot(dplyr::filter(plastic_df_calc, Size_fraction == "<20"))+
  geom_point(aes(x = CStarTr0, y = total_L, color = station))

ggplot(dplyr::filter(plastic_df_calc, Size_fraction == "<20"))+
  geom_point(aes(x = flECO.AFL, y = total_L, color = station))

ggplot(dplyr::filter(plastic_df_calc, Size_fraction == "<20"))+
  geom_point(aes(x = Density, y = total_L, color = station))

ggplot(dplyr::filter(plastic_df_calc, Size_fraction == "<20"))+
  geom_point(aes(x = sal00, y = total_L, color = station))

ggplot(dplyr::filter(plastic_df_calc, Size_fraction == "<20"))+
  geom_point(aes(x = depth_m, y = total_L, color = station))

ggplot(plastic_df_calc)+
  geom_point(aes(x = depth_m, y = total_L, color = station))

ggplot(plastic_df_calc)+
  geom_point(aes(x = flECO.AFL, y = total_L, color = station))

ggplot(plastic_df_calc)+
  geom_point(aes(x =  par, y = total_L, color = station))

ggplot(plastic_df_calc)+
  geom_col(aes(x = station, y = PE_L, fill = Size_fraction))

ggplot(plastic_df_calc)+
  geom_col(aes(x = station, y = ABS_L, fill = Size_fraction))

ggplot(plastic_df_calc)+
  geom_col(aes(x = station, y = PET_L, fill = Size_fraction))

ggplot(plastic_df_calc)+
  geom_col(aes(x = station, y = PP_L, fill = Size_fraction))

ggplot(plastic_df_calc)+
  geom_col(aes(x = station, y = PA_L, fill = Size_fraction))

ggplot(plastic_df_calc)+
  geom_col(aes(x = Size_fraction, y = PE_L))

ggplot(plastic_df_calc)+
  geom_col(aes(x = Size_fraction, y = ABS_L))

ggplot(plastic_df_calc)+
  geom_col(aes(x = Size_fraction, y = PET_L))

ggplot(plastic_df_calc)+
  geom_col(aes(x = Size_fraction, y = PP_L))

ggplot(plastic_df_calc)+
  geom_col(aes(x = Size_fraction, y = PA_L))

#Some initial notes 
#Total plastics
#There is are a lot of 0s in this dataset, especially in the year 2022.
#There are very few estimated plastics estimated by either of the larger size fractions
#There is an order of magnitude difference as we move from largest to smallest size
#Almost all of the estimated plastics come in the smallest size fraction.
#Stn_R has a nearly triple of the amount of plastic estimated...BUT this is not
# a proper water column integration.
#There is some signal (I haven't tested this yet) with a variety of different
# env vars such as:
# salinity, density, pressure, depth, size fraction, %transmission, flourescence

#Polymer specific notes
# All identified polymers occurs at both sampling periods.
# There are major differences between size fractions for each polymer and 
# some but less important between stations

#ABS only occurs at the smallest size fraction. This is the only size specific 
# polymer in the dataset.

# By concentration, ABS is the second largest in amount of microplastics

# PE is the largest by total concentration and the majority of them is in the 
# smallest size fraction.

# PP was the smallest of the 5 plastic polymers identified.
# almost entirely absent from station R but quantified to be a few hundred per L

# PET is the only polymer to show any real signal at the largest size fraction
# again largest amount in the smallest size fraction.

# PA doesn't show anything different than previously mentioned

# I need to do composition plots with depth bar/frequency plots
#First I need to change from wide format to long format to properly address this

plastic_df_calc_forMELT <- plastic_df_calc[,c(1:3,32,41:45)]

plastic_df_calc_melt <- melt(plastic_df_calc_forMELT, 
                             id.vars = c("station","depth", "depth_m","Size_fraction"),
                             measure.vars = c("ABS_L","PET_L","PE_L","PP_L","PA_L"))

ggplot(plastic_df_calc_melt)+
  geom_col(aes(x = Size_fraction, y = value, fill = variable))

ggplot(plastic_df_calc_melt)+
  geom_line(aes(x = depth_m, y = value, color = variable, linetype = station))+
  coord_flip()+
  scale_y_continuous(position = "right")+
  scale_x_reverse()

ggplot(plastic_df_calc_melt)+
  geom_line(aes(x = depth_m, y = value, color = variable))+
  coord_flip()+
  scale_y_continuous(position = "right")+
  scale_x_reverse()+facet_wrap(station~Size_fraction)

#This approach works for lineplot but not for a stacked barplot...I need to fix
# the data structure to be able to visualize the data

plastic_df_calc_forMELT_stacked <- plastic_df_calc[,c(1:3,32,41:46)] %>% 
                              mutate(ABS_prop = ABS_L/total_L,
                                     PET_prop = PET_L/total_L,
                                     PE_prop = PE_L/total_L,
                                     PP_prop = PP_L/total_L,
                                     PA_prop = PA_L/total_L) %>% dplyr::filter(total_L != 0)

plastic_df_calc_melt_stacked <- melt(plastic_df_calc_forMELT_stacked, 
                             id.vars = c("station","depth", "depth_m","Size_fraction"),
                             measure.vars = c("ABS_prop","PET_prop","PE_prop","PP_prop","PA_prop"))

ggplot(plastic_df_calc_melt_stacked)+
  geom_col(aes(x = Size_fraction, y = value, fill = variable))

ggplot(plastic_df_calc_melt_stacked)+
  geom_col(aes(x = depth, y = value, fill = variable))+facet_wrap(station~Size_fraction, nrow = 1)

ggplot(plastic_df_calc_melt_stacked)+
  geom_line(aes(x = depth_m, y = value, color = variable))+
  coord_flip()+
  scale_y_continuous(position = "right")+
  scale_x_reverse()+facet_grid(cols = vars(Size_fraction))

ggplot(plastic_df_calc_melt_stacked)+
  geom_line(aes(x = depth_m, y = value, color = variable))+
  coord_flip()+
  scale_y_continuous(position = "right", label = c("0","0.25", "0.5", "0.75","1"))+
  scale_x_reverse()+facet_grid(cols = vars(station))

ggplot(plastic_df_calc_melt_stacked)+
  geom_line(aes(x = depth_m, y = value, color = variable))+
  coord_flip()+
  scale_y_continuous(position = "right", label = c("0","0.25", "0.5", "0.75","1"))+
  scale_x_reverse()+facet_wrap(station~Size_fraction)

#Some notes from additional visualizations:
# I don't like the line plots until I separate it into separate station or size.
#It's too busy and jumpy to make any obvious conclusions.

#When line plots are subdivided by station:


 #Making of the NMDS --------
#Subset and transform the data into relative "abundances" for each of the plastic
# "species"

#There are several site at which there are no plastics found at all. With these
# I need to remove the sites from all analysis

#NMDS_df <- plastic_df_calc[, 41:46]
#
#NMDS_df_transformed <- decostand(NMDS_df, "hellinger")
#
#NMDS_df_transformed_distance <- vegdist(NMDS_df, "bray")
#
#all_NMDS <- metaMDS(NMDS_df_transformed_distance)
#Previous chunk did not work

#Transforming data set to filter out sites with total_L = 0

plastic_df_calc_no0 <- dplyr::filter(plastic_df_calc, total_L != 0)

group_by(plastic_df_calc_no0, Size_fraction) %>% summarise(counts = n())

NMDS_df <- plastic_df_calc_no0[, 41:45]

env_df <- plastic_df_calc_no0[, 1:32]

NMDS_df_transformed <- decostand(NMDS_df, "hellinger")

NMDS_df_transformed_distance <- vegdist(NMDS_df, "bray")

all_NMDS <- metaMDS(NMDS_df_transformed)

all_NMDS_untran <- metaMDS(NMDS_df)

all_NMDS_fit <- envfit(all_NMDS ~ ., env_df, perm = 999)

all_NMDS_fit_untran <- envfit(all_NMDS_untran ~ ., env_df, perm = 999)

#Extracting scores for ggploting
all_NMDS_scores <- as.data.frame(scores(all_NMDS, display = "sites"))
all_NMDS_species <- as.data.frame(scores(all_NMDS, display = "species"))
all_NMDS_scores$depth <- plastic_df_calc_no0$depth
all_NMDS_scores$depth_m <- plastic_df_calc_no0$depth_m
all_NMDS_scores$station <- plastic_df_calc_no0$station
all_NMDS_scores$Size_fraction <- plastic_df_calc_no0$Size_fraction
#scores of vectors from the envfit used to then plot on top of original PCAs they were fit to
all_NMDS_fit_scores_vectors <- as.data.frame(scores(all_NMDS_fit, "vectors"))

all_NMDS_fit_scores_factors <- as.data.frame(scores(all_NMDS_fit, "factors"))

ggplot(data = all_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = depth))

ggplot(data = all_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = station))

ggplot(data = all_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = depth_m, shape = Size_fraction))

ggplot(data = all_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(size = Size_fraction, shape = station, color = depth_m))+
  #geom_segment(data = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),], aes(x = 0, y = 0, xend = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),1], yend = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),2]))+
  #geom_text_repel(data = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),], aes(x = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),1], y = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),2]-0.04), fontface = "bold", label = rownames(all_NMDS_fit_scores_vectors[c(3,5,19,20,23),]))+
  geom_segment(data = all_NMDS_fit_scores_factors[c(12:14),], aes(x = 0, y = 0, xend = all_NMDS_fit_scores_factors[c(12:14),1], yend = all_NMDS_fit_scores_factors[c(12:14),2]))+
  geom_text_repel(data = all_NMDS_fit_scores_factors[c(12:14),], aes(x = all_NMDS_fit_scores_factors[c(12:14),1], y = all_NMDS_fit_scores_factors[c(12:14),2]-0.04), fontface = "bold", label = rownames(all_NMDS_fit_scores_factors[c(12:14),]))+
  geom_segment(data = all_NMDS_species, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2/2), color = "red")+
  geom_text_repel(data = all_NMDS_species, aes(x = NMDS1, y = NMDS2/2), label = rownames(all_NMDS_species))#+
#geom_text(size = 3,label = rownames(all_NMDS_scores))

#I need to repeat with untransformed data to get and display the species scores
#Extracting scores for ggploting
all_NMDS_scores_untran <- as.data.frame(scores(all_NMDS_untran, display = "sites"))
all_NMDS_species_untran <- as.data.frame(scores(all_NMDS_untran, display = "species"))
all_NMDS_scores_untran$depth <- plastic_df_calc_no0$depth
all_NMDS_scores_untran$depth_m <- plastic_df_calc_no0$depth_m
all_NMDS_scores_untran$station <- plastic_df_calc_no0$station
all_NMDS_scores_untran$Size_fraction <- plastic_df_calc_no0$Size_fraction

all_NMDS_fit_scores_vectors_untran <- as.data.frame(scores(all_NMDS_fit_untran, "vectors"))

all_NMDS_fit_scores_factors_untran <- as.data.frame(scores(all_NMDS_fit_untran, "factors"))

ggplot(data = all_NMDS_scores_untran, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(size = station, shape = Size_fraction, color = depth_m))+
  geom_segment(data = all_NMDS_fit_scores_vectors[c(17,19,20),], aes(x = 0, y = 0, xend = all_NMDS_fit_scores_vectors[c(17,19,20),1], yend = all_NMDS_fit_scores_vectors[c(17,19,20),2]))+
  geom_text_repel(data = all_NMDS_fit_scores_vectors[c(17,19,20),], aes(x = all_NMDS_fit_scores_vectors[c(17,19,20),1], y = all_NMDS_fit_scores_vectors[c(17,19,20),2]-0.04), fontface = "bold", label = rownames(all_NMDS_fit_scores_vectors[c(17,19,20),]))+
  geom_segment(data = all_NMDS_species_untran, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "red")+
  geom_text_repel(data = all_NMDS_species_untran, aes(x = NMDS1, y = NMDS2), label = rownames(all_NMDS_species_untran))#+
  #geom_text(size = 3,label = rownames(all_NMDS_scores_untran))

#Some initial notes from the NMDS analysis
#There are very few significant environmental factors even relatively significant.
#This may be due to a lot of autocorrelation within the dataset and may not be true
#If I redo the analysis with just the <20 size fraction

#There seems to be a difference of polymers estimated between and within the size
# fractions and size fraction is the most obvious difference within the dataset

# NMDS w/0s **USELESS data** -----------
#Transforming data set to filter out sites with total_L = 0

#plastic_df_calc_no20 <- dplyr::filter(plastic_df_calc, Size_fraction != "<20")
#
#no20_NMDS_df <- plastic_df_calc_no20[, 41:45]
#
#no20_env_df <- plastic_df_calc_no20[, 1:32]
#
#NMDS_df_transformed <- decostand(no20_NMDS_df, "hellinger")
#
#NMDS_df_transformed_distance <- vegdist(no20_NMDS_df, "bray")
#
#no20_NMDS <- metaMDS(NMDS_df_transformed)
#
#no20_NMDS_untran <- metaMDS(no20_NMDS_df, distance = "bray" )
#
#no20_NMDS_fit <- envfit(all_NMDS ~ ., no20_env_df, perm = 999)
#
#no20_NMDS_fit_untran <- envfit(no20_NMDS_untran ~ ., no20_env_df, perm = 999)
#
##Extracting scores for ggploting
#all_NMDS_scores <- as.data.frame(scores(all_NMDS, display = "sites"))
#all_NMDS_species <- as.data.frame(scores(all_NMDS, display = "species"))
#all_NMDS_scores$depth <- plastic_df_calc_no0$depth
#all_NMDS_scores$depth_m <- plastic_df_calc_no0$depth_m
#all_NMDS_scores$station <- plastic_df_calc_no0$station
#all_NMDS_scores$Size_fraction <- plastic_df_calc_no0$Size_fraction
##scores of vectors from the envfit used to then plot on top of original PCAs they were fit to
#all_NMDS_fit_scores_vectors <- as.data.frame(scores(all_NMDS_fit, "vectors"))
#
#all_NMDS_fit_scores_factors <- as.data.frame(scores(all_NMDS_fit, "factors"))
#
#ggplot(data = all_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
#  geom_point(aes(color = depth))
#
#ggplot(data = all_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
#  geom_point(aes(color = station))
#
#ggplot(data = all_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
#  geom_point(aes(color = depth_m, shape = Size_fraction))
#
#ggplot(data = all_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
#  geom_point(aes(color = Size_fraction, shape = station))+
#  geom_segment(data = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),], aes(x = 0, y = 0, xend = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),1], yend = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),2]))+
#  geom_text_repel(data = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),], aes(x = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),1], y = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),2]-0.04), fontface = "bold", label = rownames(all_NMDS_fit_scores_vectors[c(3,5,19,20,23),]))+
#  geom_segment(data = all_NMDS_species, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2/2), color = "red")+
#  geom_text_repel(data = all_NMDS_species, aes(x = NMDS1, y = NMDS2/2), label = rownames(all_NMDS_species))#+
##geom_text(size = 3,label = rownames(all_NMDS_scores))
#
##I need to repeat with untransformed data to get and display the species scores
##Extracting scores for ggploting
#no20_NMDS_scores_untran <- as.data.frame(scores(no20_NMDS_untran, display = "sites"))
#no20_NMDS_species_untran <- as.data.frame(scores(no20_NMDS_untran, display = "species"))
#no20_NMDS_scores_untran$depth <- plastic_df_calc_no20$depth
#no20_NMDS_scores_untran$depth_m <- plastic_df_calc_no20$depth_m
#no20_NMDS_scores_untran$station <- plastic_df_calc_no20$station
#no20_NMDS_scores_untran$Size_fraction <- plastic_df_calc_no20$Size_fraction
#
#no20_NMDS_fit_scores_vectors_untran <- as.data.frame(scores(no20_NMDS_fit_untran, "vectors"))
#
#no20_NMDS_fit_scores_factors_untran <- as.data.frame(scores(no20_NMDS_fit_untran, "factors"))
#
#ggplot(data = no20_NMDS_scores_untran, aes(x = NMDS1, y = NMDS2))+
#  geom_point(aes(color = depth_m, shape = station))+
#  #geom_segment(data = no20_NMDS_fit_scores_vectors_untran[c(17,19,20),], aes(x = 0, y = 0, xend = no20_NMDS_fit_scores_vectors_untran[c(17,19,20),1], yend = all_NMDS_fit_scores_vectors[c(17,19,20),2]))+
#  #geom_text_repel(data = no20_NMDS_fit_scores_vectors_untran[c(17,19,20),], aes(x = no20_NMDS_fit_scores_vectors_untran[c(17,19,20),1], y = all_NMDS_fit_scores_vectors[c(17,19,20),2]-0.04), fontface = "bold", label = rownames(all_NMDS_fit_scores_vectors[c(17,19,20),]))+
#  geom_segment(data = no20_NMDS_species_untran, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "red")+
#  geom_text_repel(data = no20_NMDS_species_untran, aes(x = NMDS1, y = NMDS2), label = rownames(no20_NMDS_species_untran))
#  #geom_text(size = 3,label = rownames(no20_NMDS_scores_untran))
#
# G 20 only NMDS -----------

plastic_df_calc_G20 <- dplyr::filter(plastic_df_calc, Size_fraction == ">20", total_L != 0)

G20_NMDS_df <- plastic_df_calc_G20[, 41:45]

G20_env_df <- plastic_df_calc_G20[, 1:32]

G20_NMDS_df_transformed <- decostand(G20_NMDS_df, "hellinger")

G20_NMDS_df_transformed_distance <- vegdist(G20_NMDS_df, "bray")

G20_NMDS <- metaMDS(G20_NMDS_df_transformed)

G20_NMDS_untran <- metaMDS(G20_NMDS_df, distance = "bray" )

G20_NMDS_fit <- envfit(G20_NMDS ~ ., G20_env_df, perm = 999)

G20_NMDS_fit_untran <- envfit(G20_NMDS_untran ~ ., G20_env_df, perm = 999)

#Extracting scores for ggploting
G20_NMDS_scores <- as.data.frame(scores(G20_NMDS, display = "sites"))
G20_NMDS_species <- as.data.frame(scores(G20_NMDS, display = "species"))
G20_NMDS_scores$depth <- plastic_df_calc_G20$depth
G20_NMDS_scores$depth_m <- plastic_df_calc_G20$depth_m
G20_NMDS_scores$station <- plastic_df_calc_G20$station
G20_NMDS_scores$Size_fraction <- plastic_df_calc_G20$Size_fraction
#scores of vectors from the envfit used to then plot on top of original PCAs they were fit to
G20_NMDS_fit_scores_vectors <- as.data.frame(scores(G20_NMDS_fit, "vectors"))

G20_NMDS_fit_scores_factors <- as.data.frame(scores(G20_NMDS_fit, "factors"))

ggplot(data = G20_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = depth))

ggplot(data = G20_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = station))

ggplot(data = G20_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = depth_m, shape = Size_fraction))

ggplot(data = G20_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = depth_m, shape = station), position = position_jitter(0.1))+
  geom_segment(data = G20_NMDS_fit_scores_vectors[c(9,10,13),], aes(x = 0, y = 0, xend = G20_NMDS_fit_scores_vectors[c(9,10,13),1], yend = G20_NMDS_fit_scores_vectors[c(9,10,13),2]))+
  geom_text_repel(data = G20_NMDS_fit_scores_vectors[c(9,10,13),], aes(x = G20_NMDS_fit_scores_vectors[c(9,10,13),1], y = G20_NMDS_fit_scores_vectors[c(9,10,13),2]-0.04), fontface = "bold", label = rownames(all_NMDS_fit_scores_vectors[c(9,10,13),]))+
  geom_segment(data = G20_NMDS_species, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2/2), color = "red")+
  geom_text_repel(data = G20_NMDS_species, aes(x = NMDS1, y = NMDS2/2), label = rownames(G20_NMDS_species))#+
#geom_text(size = 3,label = rownames(all_NMDS_scores))

#I need to repeat with untransformed data to get and display the species scores
#Extracting scores for ggploting
G20_NMDS_scores_untran <- as.data.frame(scores(G20_NMDS_untran, display = "sites"))
G20_NMDS_species_untran <- as.data.frame(scores(G20_NMDS_untran, display = "species"))
G20_NMDS_scores_untran$depth <- plastic_df_calc_G20$depth
G20_NMDS_scores_untran$depth_m <- plastic_df_calc_G20$depth_m
G20_NMDS_scores_untran$station <- plastic_df_calc_G20$station
G20_NMDS_scores_untran$Size_fraction <- plastic_df_calc_G20$Size_fraction

G20_NMDS_fit_scores_vectors_untran <- as.data.frame(scores(G20_NMDS_fit_untran, "vectors"))

G20_NMDS_fit_scores_factors_untran <- as.data.frame(scores(G20_NMDS_fit_untran, "factors"))

ggplot(data = G20_NMDS_scores_untran, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = depth_m, shape = station), position = position_jitter(0.1))+
  #geom_segment(data = no20_NMDS_fit_scores_vectors_untran[c(17,19,20),], aes(x = 0, y = 0, xend = no20_NMDS_fit_scores_vectors_untran[c(17,19,20),1], yend = all_NMDS_fit_scores_vectors[c(17,19,20),2]))+
  #geom_text_repel(data = no20_NMDS_fit_scores_vectors_untran[c(17,19,20),], aes(x = no20_NMDS_fit_scores_vectors_untran[c(17,19,20),1], y = all_NMDS_fit_scores_vectors[c(17,19,20),2]-0.04), fontface = "bold", label = rownames(all_NMDS_fit_scores_vectors[c(17,19,20),]))+
  geom_segment(data = G20_NMDS_species_untran, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "red")+
  geom_text_repel(data = G20_NMDS_species_untran, aes(x = NMDS1, y = NMDS2), label = rownames(G20_NMDS_species_untran))
#geom_text(size = 3,label = rownames(no20_NMDS_scores_untran))

# <20 µm NMDS --------
plastic_df_calc_less20 <- dplyr::filter(plastic_df_calc, Size_fraction == "<20")



NMDS_df_less20 <- plastic_df_calc_less20[, 41:45]

env_df_less20 <- plastic_df_calc_less20[, 1:32]

NMDS_df_transformed_less20 <- decostand(NMDS_df_less20, "hellinger")

#NMDS_df_transformed_distance <- vegdist(NMDS_df, "bray")

all_NMDS_less20 <- metaMDS(NMDS_df_transformed_less20)

all_NMDS_untran_less20 <- metaMDS(NMDS_df)

all_NMDS_fit_less20 <- envfit(all_NMDS_less20 ~ ., env_df_less20, perm = 999)

all_NMDS_fit_untran_less20 <- envfit(all_NMDS_less20 ~ ., env_df_less20, perm = 999)

#Extracting scores for ggploting
all_NMDS_scores_less20 <- as.data.frame(scores(all_NMDS_less20, display = "sites"))
all_NMDS_species_less20 <- as.data.frame(scores(all_NMDS_less20, display = "species"))
all_NMDS_scores_less20$depth <- plastic_df_calc_less20$depth
all_NMDS_scores_less20$depth_m <- plastic_df_calc_less20$depth_m
all_NMDS_scores_less20$station <- plastic_df_calc_less20$station
all_NMDS_scores_less20$Size_fraction <- plastic_df_calc_less20$Size_fraction
#scores of vectors from the envfit used to then plot on top of original PCAs they were fit to
all_NMDS_fit_scores_vectors_less20 <- as.data.frame(scores(all_NMDS_fit_less20, "vectors"))

all_NMDS_fit_scores_factors_less20 <- as.data.frame(scores(all_NMDS_fit_less20, "factors"))

ggplot(data = all_NMDS_scores_less20, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = depth), position = position_jitter(.1))

ggplot(data = all_NMDS_scores_less20, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = station), position = position_jitter(.1))

ggplot(data = all_NMDS_scores_less20, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = depth, shape = station), position = position_jitter(0.1))+
  geom_segment(data = all_NMDS_species_less20, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2/2), color = "red")+
  geom_text_repel(data = all_NMDS_species_less20, aes(x = NMDS1, y = NMDS2/2), label = rownames(all_NMDS_species_less20))+
  geom_text_repel(size = 3,label = rownames(all_NMDS_scores_less20))


# Hard to really say much with the extremely low stress but we see a slight
# separation of depth, Break between 15 m and 23 m sample so "deeper" and "Shallower"
# Could say 3 breaks, "Surface", "Mid", "Deep"
# nothing of significance due to the low amount of sample size in this size fraction
# 

# There is little distincts separation 


#Correlation Matrix ---------
#Need to create a correlation matrix AND a p-value matrix for said correlations
#Old code repeated this step but I since this is a mirrored matrix 1 set of names 
# is needed
#I need to trim the df to remove the non-numeric/non-informative columns
plastic_df_calc_forMatrix <- plastic_df_calc[, -c(1,2,26,27,32:40)] %>% 
                                                  mutate(ABS_prop = ABS_L/total_L,
                                                          PET_prop = PET_L/total_L,
                                                          PE_prop = PE_L/total_L,
                                                          PP_prop = PP_L/total_L,
                                                          PA_prop = PA_L/total_L,
                                                         Light_L = PE_L + PP_L,
                                                         Heavy_L = total_L-Light_L,
                                                         Light_prop = PE_prop + PP_prop,
                                                         Heavy_prop = 1-Light_L)

cor_names <- colnames(plastic_df_calc_forMatrix)

cor_matrix <- matrix(nrow = length(colnames(plastic_df_calc_forMatrix)),
                     ncol = length(colnames(plastic_df_calc_forMatrix)), dimnames = list(cor_names, cor_names))
cor_p_matrix <- matrix(nrow = length(colnames(plastic_df_calc_forMatrix)),
                     ncol = length(colnames(plastic_df_calc_forMatrix)), dimnames = list(cor_names, cor_names))

#str(ant_community_final_calc_grp_remStations[,-(1:3)])
temp_i_vector <- vector()
temp_j_vector <- vector()
for(i in colnames(cor_matrix)){
  for(j in rownames(cor_matrix)){
    temp_i_vector <- as.numeric(unlist(plastic_df_calc_forMatrix[,i]))
    temp_j_vector <- as.numeric(unlist(plastic_df_calc_forMatrix[,j]))
    temp_cor <- cor.test(temp_i_vector,temp_j_vector)
    cor_matrix[j,i] <- temp_cor[["estimate"]]
    cor_p_matrix[j,i] <- temp_cor[["p.value"]]
  }
}

cor_matrix_melt <- melt(cor_matrix)
cor_p_matrix_melt <- melt(cor_p_matrix)

ggplot(cor_matrix_melt, aes(X1,X2))+
  geom_tile(aes(fill = value))+scale_fill_gradient2()

cor_matrix_melt$pvalue <- cor_p_matrix_melt$value

#This is a plot that contains all of the variables
# I would need to split this dataset into env only and response only for dissertation
ggplot(cor_matrix_melt, aes(X1,X2))+
  geom_tile(aes(fill = value, color = pvalue<0.05), lwd = 0.5, height = 0.75, width = 0.95)+scale_color_manual(name ='pvalue<0.05', values = setNames(c('black','White'), c(T,F)))+scale_fill_gradient2()+
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 8))+
  geom_text(label = round(cor_matrix_melt$value, 3), size = 1.5)+
  scale_x_discrete(guide = guide_axis(angle = 90))

# Initial notes from this plots....
# Very little correlation with plastic concentrations and the polymers themselves

# Only 2 plastic polymers shared a significant correlation with % transmission
# PET and ABS
# They both also shared correlation with PAR suggesting a connection between PAR 
# and % trans

# ABS is also associated with salinity and oxygen levels

# One mechanistic possibility is the smallest fraction gets caught in the 
# microbial biomass, which prevents the sinking

# Otherwise there is no signiicant signal with any other environmental vars

# Info on the Proportional data:
# Mostly the same for the raw concetration data. There are some differences:

# For ABS, signal for %transmission is lost however conductivity is now significant

# For PA, density, depth, and TimeS (correlated with depth)

# For PET, TimeJ, ZML, ZCM, Ze are all correlated and signal for differences for
# Entire station

# PE and PP both show no significant signal with env vars
# Even though PE is the largest estimated plastic polymer
# Maybe because it is the most common is doesn't show any significance
# Would this still be true if I subset by size fraction?

plastic_df_calc_forMELT_hell <- plastic_df_calc

plastic_df_calc_forMELT_hell[, 41:45] <- decostand(plastic_df_calc[, 41:45], "hellinger")

plastic_df_calc_forMatrix_hell <- plastic_df_calc_forMELT_hell[, -c(1,2,26,27,32:40)]%>% 
  mutate(ABS_prop = ABS_L/total_L,
         PET_prop = PET_L/total_L,
         PE_prop = PE_L/total_L,
         PP_prop = PP_L/total_L,
         PA_prop = PA_L/total_L,
         Light_L = PE_L + PP_L,
         Heavy_L = total_L-Light_L,
         Light_prop = PE_prop + PP_prop,
         Heavy_prop = 1-Light_L)

cor_names <- colnames(plastic_df_calc_forMatrix_hell)

cor_matrix_hell <- matrix(nrow = length(colnames(plastic_df_calc_forMatrix_hell)),
                     ncol = length(colnames(plastic_df_calc_forMatrix_hell)), dimnames = list(cor_names, cor_names))
cor_p_matrix_hell <- matrix(nrow = length(colnames(plastic_df_calc_forMatrix_hell)),
                       ncol = length(colnames(plastic_df_calc_forMatrix_hell)), dimnames = list(cor_names, cor_names))

#str(ant_community_final_calc_grp_remStations[,-(1:3)])
temp_i_vector <- vector()
temp_j_vector <- vector()
for(i in colnames(cor_matrix_hell)){
  for(j in rownames(cor_matrix_hell)){
    temp_i_vector <- as.numeric(unlist(plastic_df_calc_forMatrix_hell[,i]))
    temp_j_vector <- as.numeric(unlist(plastic_df_calc_forMatrix_hell[,j]))
    temp_cor <- cor.test(temp_i_vector,temp_j_vector)
    cor_matrix_hell[j,i] <- temp_cor[["estimate"]]
    cor_p_matrix_hell[j,i] <- temp_cor[["p.value"]]
  }
}

cor_matrix_melt_hell <- melt(cor_matrix_hell)
cor_p_matrix_melt_hell <- melt(cor_p_matrix_hell)

ggplot(cor_matrix_melt_hell, aes(X1,X2))+
  geom_tile(aes(fill = value))+scale_fill_gradient2()

cor_matrix_melt_hell$pvalue <- cor_p_matrix_melt_hell$value

ggplot(cor_matrix_melt_hell, aes(X1,X2))+
  geom_tile(aes(fill = value, color = pvalue<0.05), lwd = 0.5, height = 0.75, width = 0.95)+scale_color_manual(name ='pvalue<0.05', values = setNames(c('black','White'), c(T,F)))+scale_fill_gradient2()+
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 8))+
  geom_text(label = round(cor_matrix_melt_hell$value, 3), size = 1.5)+
  scale_x_discrete(guide = guide_axis(angle = 45))

# Need to repeat this step for both the <20 and >20 size fractions
# Filter from previous chunk: plastic_df_calc_forMatrix_hell
plastic_df_calc_forMELT_hell_L20 <- filter(plastic_df_calc_forMELT_hell, Size_fraction == "<20")

plastic_df_calc_forMatrix_hell_L20 <- plastic_df_calc_forMELT_hell_L20[, -c(1,2,26,27,32:40)]%>% 
  mutate(ABS_prop = ABS_L/total_L,
         PET_prop = PET_L/total_L,
         PE_prop = PE_L/total_L,
         PP_prop = PP_L/total_L,
         PA_prop = PA_L/total_L,
         Light_L = PE_L + PP_L,
         Heavy_L = total_L-Light_L,
         Light_prop = PE_prop + PP_prop,
         Heavy_prop = 1-Light_L)

cor_names <- colnames(plastic_df_calc_forMatrix_hell_L20)

cor_matrix_hell_L20 <- matrix(nrow = length(colnames(plastic_df_calc_forMatrix_hell_L20)),
                          ncol = length(colnames(plastic_df_calc_forMatrix_hell_L20)), dimnames = list(cor_names, cor_names))
cor_p_matrix_hell_L20 <- matrix(nrow = length(colnames(plastic_df_calc_forMatrix_hell_L20)),
                            ncol = length(colnames(plastic_df_calc_forMatrix_hell_L20)), dimnames = list(cor_names, cor_names))

#str(ant_community_final_calc_grp_remStations[,-(1:3)])
temp_i_vector <- vector()
temp_j_vector <- vector()
for(i in colnames(cor_matrix_hell_L20)){
  for(j in rownames(cor_matrix_hell_L20)){
    temp_i_vector <- as.numeric(unlist(plastic_df_calc_forMatrix_hell_L20[,i]))
    temp_j_vector <- as.numeric(unlist(plastic_df_calc_forMatrix_hell_L20[,j]))
    temp_cor <- cor.test(temp_i_vector,temp_j_vector)
    cor_matrix_hell_L20[j,i] <- temp_cor[["estimate"]]
    cor_p_matrix_hell_L20[j,i] <- temp_cor[["p.value"]]
  }
}

cor_matrix_melt_hell_L20 <- melt(cor_matrix_hell_L20)
cor_p_matrix_melt_hell_L20 <- melt(cor_p_matrix_hell_L20)

ggplot(cor_matrix_melt_hell_L20, aes(X1,X2))+
  geom_tile(aes(fill = value))+scale_fill_gradient2()

cor_matrix_melt_hell_L20$pvalue <- cor_p_matrix_melt_hell_L20$value

ggplot(cor_matrix_melt_hell_L20, aes(X1,X2))+
  geom_tile(aes(fill = value, color = pvalue<0.05), lwd = 0.5, height = 0.75, width = 0.95)+scale_color_manual(name ='pvalue<0.05', values = setNames(c('black','White'), c(T,F)))+scale_fill_gradient2()+
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 8))+
  geom_text(label = round(cor_matrix_melt_hell_L20$value, 3), size = 1.5)+
  scale_x_discrete(guide = guide_axis(angle = 45))

# Again for G20
# Filter from previous chunk: plastic_df_calc_forMatrix_hell
plastic_df_calc_forMELT_hell_G20 <- filter(plastic_df_calc_forMELT_hell, Size_fraction == ">20")

plastic_df_calc_forMatrix_hell_G20 <- plastic_df_calc_forMELT_hell_G20[, -c(1,2,26,27,32:40)]%>% 
  mutate(ABS_prop = ABS_L/total_L,
         PET_prop = PET_L/total_L,
         PE_prop = PE_L/total_L,
         PP_prop = PP_L/total_L,
         PA_prop = PA_L/total_L,
         Light_L = PE_L + PP_L,
         Heavy_L = total_L-Light_L,
         Light_prop = PE_prop + PP_prop,
         Heavy_prop = 1-Light_L)

cor_names <- colnames(plastic_df_calc_forMatrix_hell_G20)

cor_matrix_hell_G20 <- matrix(nrow = length(colnames(plastic_df_calc_forMatrix_hell_G20)),
                              ncol = length(colnames(plastic_df_calc_forMatrix_hell_G20)), dimnames = list(cor_names, cor_names))
cor_p_matrix_hell_G20 <- matrix(nrow = length(colnames(plastic_df_calc_forMatrix_hell_G20)),
                                ncol = length(colnames(plastic_df_calc_forMatrix_hell_G20)), dimnames = list(cor_names, cor_names))

#str(ant_community_final_calc_grp_remStations[,-(1:3)])
temp_i_vector <- vector()
temp_j_vector <- vector()
for(i in colnames(cor_matrix_hell_G20)){
  for(j in rownames(cor_matrix_hell_G20)){
    temp_i_vector <- as.numeric(unlist(plastic_df_calc_forMatrix_hell_G20[,i]))
    temp_j_vector <- as.numeric(unlist(plastic_df_calc_forMatrix_hell_G20[,j]))
    temp_cor <- cor.test(temp_i_vector,temp_j_vector)
    cor_matrix_hell_G20[j,i] <- temp_cor[["estimate"]]
    cor_p_matrix_hell_G20[j,i] <- temp_cor[["p.value"]]
  }
}

cor_matrix_melt_hell_G20 <- melt(cor_matrix_hell_G20)
cor_p_matrix_melt_hell_G20 <- melt(cor_p_matrix_hell_G20)

ggplot(cor_matrix_melt_hell_G20, aes(X1,X2))+
  geom_tile(aes(fill = value))+scale_fill_gradient2()

cor_matrix_melt_hell_G20$pvalue <- cor_p_matrix_melt_hell_G20$value

ggplot(cor_matrix_melt_hell_G20, aes(X1,X2))+
  geom_tile(aes(fill = value, color = pvalue<0.05), lwd = 0.5, height = 0.75, width = 0.95)+scale_color_manual(name ='pvalue<0.05', values = setNames(c('black','White'), c(T,F)))+scale_fill_gradient2()+
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 8))+
  geom_text(label = round(cor_matrix_melt_hell_G20$value, 3), size = 1.5)+
  scale_x_discrete(guide = guide_axis(angle = 45))

# Again for G100
# Filter from previous chunk: plastic_df_calc_forMatrix_hell
plastic_df_calc_forMELT_hell_G100 <- filter(plastic_df_calc_forMELT_hell, Size_fraction == ">100")

plastic_df_calc_forMatrix_hell_G100 <- plastic_df_calc_forMELT_hell_G100[, -c(1,2,26,27,32:40)]%>% 
  mutate(ABS_prop = ABS_L/total_L,
         PET_prop = PET_L/total_L,
         PE_prop = PE_L/total_L,
         PP_prop = PP_L/total_L,
         PA_prop = PA_L/total_L,
         Light_L = PE_L + PP_L,
         Heavy_L = total_L-Light_L,
         Light_prop = PE_prop + PP_prop,
         Heavy_prop = 1-Light_L)

cor_names <- colnames(plastic_df_calc_forMatrix_hell_G100)

cor_matrix_hell_G100 <- matrix(nrow = length(colnames(plastic_df_calc_forMatrix_hell_G100)),
                              ncol = length(colnames(plastic_df_calc_forMatrix_hell_G100)), dimnames = list(cor_names, cor_names))
cor_p_matrix_hell_G100 <- matrix(nrow = length(colnames(plastic_df_calc_forMatrix_hell_G100)),
                                ncol = length(colnames(plastic_df_calc_forMatrix_hell_G100)), dimnames = list(cor_names, cor_names))

#str(ant_community_final_calc_grp_remStations[,-(1:3)])
temp_i_vector <- vector()
temp_j_vector <- vector()
for(i in colnames(cor_matrix_hell_G100)){
  for(j in rownames(cor_matrix_hell_G100)){
    temp_i_vector <- as.numeric(unlist(plastic_df_calc_forMatrix_hell_G100[,i]))
    temp_j_vector <- as.numeric(unlist(plastic_df_calc_forMatrix_hell_G100[,j]))
    temp_cor <- cor.test(temp_i_vector,temp_j_vector)
    cor_matrix_hell_G100[j,i] <- temp_cor[["estimate"]]
    cor_p_matrix_hell_G100[j,i] <- temp_cor[["p.value"]]
  }
}

cor_matrix_melt_hell_G100 <- melt(cor_matrix_hell_G100)
cor_p_matrix_melt_hell_G100 <- melt(cor_p_matrix_hell_G100)

ggplot(cor_matrix_melt_hell_G100, aes(X1,X2))+
  geom_tile(aes(fill = value))+scale_fill_gradient2()

cor_matrix_melt_hell_G100$pvalue <- cor_p_matrix_melt_hell_G100$value

ggplot(cor_matrix_melt_hell_G100, aes(X1,X2))+
  geom_tile(aes(fill = value, color = pvalue<0.05), lwd = 0.5, height = 0.75, width = 0.95)+scale_color_manual(name ='pvalue<0.05', values = setNames(c('black','White'), c(T,F)))+scale_fill_gradient2()+
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 8))+
  geom_text(label = round(cor_matrix_melt_hell_G20$value, 3), size = 1.5)+
  scale_x_discrete(guide = guide_axis(angle = 45))
# I will likely do a spearmen's rank for categorical vars to look for size 
# implications

# Now I need to make some summary pivot tables for proper visualization
# of specific sub categories with the data in long format

plastic_df_calc_forMELT_stacked_0 <- plastic_df_calc[,c(1:3,32,41:46)] %>% 
  mutate(ABS_prop = ABS_L/total_L,
         PET_prop = PET_L/total_L,
         PE_prop = PE_L/total_L,
         PP_prop = PP_L/total_L,
         PA_prop = PA_L/total_L)

plastic_df_calc_forMELT_all_0 <- plastic_df_calc %>% 
  mutate(ABS_prop = ABS_L/total_L,
         PET_prop = PET_L/total_L,
         PE_prop = PE_L/total_L,
         PP_prop = PP_L/total_L,
         PA_prop = PA_L/total_L)

plastic_df_calc_melt_all_0 <- melt(plastic_df_calc_forMELT_all_0, 
                                     id.vars = colnames(plastic_df_calc[-(35:length(plastic_df_calc))]),
                                     measure.vars = c("ABS_L","PET_L","PE_L","PP_L","PA_L"))

plastic_df_calc_melt_all_0_prop <- melt(plastic_df_calc_forMELT_all_0, 
                                   id.vars = colnames(plastic_df_calc[-(35:length(plastic_df_calc))]),
                                   measure.vars = c("ABS_prop","PET_prop","PE_prop","PP_prop","PA_prop"))

plastic_df_calc_melt_all_0_prop$value[which(plastic_df_calc_melt_all_0_prop$value == "NaN")] <- 0

plastic_df_calc_melt_all_0_prop_test <- melt(plastic_df_calc_forMELT_all_0, 
                                        id.vars = c("station","depth", "depth_m","Size_fraction", "total_L"),
                                        measure.vars = c("ABS_prop","PET_prop","PE_prop","PP_prop","PA_prop"))

plastic_df_calc_melt_all_0_prop$end_use <- with(plastic_df_calc_melt_all_0_prop, 
                                                ifelse(variable == "ABS_prop", "Auto/construction",
                                                       ifelse(variable == "PET_prop", "Packaging",
                                                              ifelse(variable == "PA_prop", "Textile/fishing","Generalist"))))

plastic_df_calc_melt_all_0$end_use <- with(plastic_df_calc_melt_all_0, 
                                                ifelse(variable == "ABS_L", "Auto/construction",
                                                       ifelse(variable == "PET_L", "Packaging",
                                                              ifelse(variable == "PA_L", "Textile/fishing","Generalist"))))

# I need to redo the line plots with all of the 0s 

ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == "<20"))+
  geom_col(aes(x = depth, y = value, fill = end_use))+
  facet_wrap(station~Size_fraction, nrow = 1)+
  theme(axis.text.x = element_text(angle = 305))
  

ggplot(plastic_df_calc_melt_all_0_prop)+
  geom_col(aes(x = depth, y = value, fill = end_use))+facet_wrap(station~Size_fraction, nrow = 2)

ggplot(plastic_df_calc_melt_all_0)+
  geom_line(aes(x = depth_m, y = value, color = variable))+
  coord_flip()+
  scale_y_continuous(position = "right")+
  scale_x_reverse()+facet_grid(cols = vars(Size_fraction))

ggplot(plastic_df_calc_melt_all_0_prop)+
  geom_line(aes(x = depth_m, y = value, color = variable))+
  coord_flip()+
  scale_y_continuous(position = "right", label = c("0","0.25", "0.5", "0.75","1"))+
  scale_x_reverse()+facet_grid(cols = vars(station))

ggplot(plastic_df_calc_melt_all_0_prop)+
  geom_line(aes(x = depth_m, y = value, color = variable), linewidth = 1)+
  coord_flip()+
  scale_y_continuous(position = "right", label = c("0",".25", ".5", ".75","1"))+
  scale_x_reverse()+facet_wrap(station~Size_fraction)

ggplot(plastic_df_calc_melt_all_0)+
  geom_line(aes(x = depth_m, y = value, color = variable))+
  coord_flip()+
  scale_y_continuous(position = "right")+
  scale_x_reverse()+facet_wrap(station~Size_fraction)

# PCA of env vars
env_df_all <- plastic_df_calc[, 1:32]
env_df_forPCA <- env_df_all[,3:31]
env_df_forPCA <- env_df_forPCA[,colnames(env_df_forPCA) !="flag"]
env_df_forPCA <- env_df_forPCA[,c(1,3,4,7,9,10,12,13,17,18,21,23)]

env_df_PCA <- prcomp(env_df_forPCA, center = TRUE, scale = TRUE)

env_df_PCA_scores <- as.data.frame(scores(env_df_PCA, display = "sites"))
env_df_PCA_species <- as.data.frame(scores(env_df_PCA, display = "species"))*10
env_df_PCA_scores$depth <- plastic_df_calc$depth
env_df_PCA_scores$depth_m <- plastic_df_calc$depth_m
env_df_PCA_scores$station <- plastic_df_calc$station
env_df_PCA_scores$Size_fraction <- plastic_df_calc$Size_fraction


ggplot(data = env_df_PCA_scores, aes(x = PC1, y = PC2))+
  geom_point(aes(color = depth_m, shape = station), size = 3)+
  geom_segment(data = env_df_PCA_species, aes(x = 0, y = 0, xend = PC1, yend = PC2), color = "red")+
  geom_text_repel(data = env_df_PCA_species, aes(x = PC1, y = PC2), label = rownames(env_df_PCA_species))

#Results Figures needed
# Figure 1 -----------
# Total abundance plots for the size fractions across depth for each station
# I need to split more
figure_1A <- ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == "<20", station == "R"))+
  geom_col(aes(x = depth, y = value, fill = variable))+facet_wrap(~station)+coord_flip()+
  scale_x_discrete(limits = rev)

figure_1B <- ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == "<20", station == "PS15"))+
  geom_col(aes(x = depth, y = value, fill = variable), position = "stack")+facet_wrap(~station)+facet_wrap(~station)+coord_flip()+
  scale_x_discrete(limits = rev)

figure_1C <- ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == ">20", station == "R"))+
  geom_col(aes(x = depth, y = value, fill = variable))+facet_wrap(~station)+coord_flip()+
  scale_x_discrete(limits = rev)

figure_1D <- ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == ">20", station == "PS15"))+
  geom_col(aes(x = depth, y = value, fill = variable))+facet_wrap(~station)+coord_flip()+
  scale_x_discrete(limits = rev)


figure_1E <- ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == ">100" & station == "R"))+
  geom_col(aes(x = depth, y = value, fill = variable))+facet_wrap(~station)+coord_flip()+
  scale_x_discrete(limits = rev)

#ggarrange(figure_1A,figure_1B,figure_1C, ncol = 1, common.legend = TRUE, legend = "bottom")

ggsave("figure_1A.pdf",figure_1A)
ggsave("figure_1B.pdf",figure_1B)
ggsave("figure_1C.pdf",figure_1C)
ggsave("figure_1D.pdf",figure_1D)
ggsave("figure_1E.pdf",figure_1E)

# Figure 2 -------------
# Each polymer abundance with depth and size fraction
ggplot(plastic_df_calc_melt_all_0_prop)+
  geom_line(aes(x = depth_m, y = value, color = variable), linewidth = 1)+
  coord_flip()+
  scale_y_continuous(position = "right", label = c("0",".25", ".5", ".75","1"))+
  scale_x_reverse()+facet_wrap(station~Size_fraction)

figure_2 <- ggplot(plastic_df_calc_melt_all_0_prop)+
  geom_col(aes(x = depth, y = value, fill = variable))+
  coord_flip()+
  facet_wrap(station~Size_fraction, nrow = 2)+
  theme(axis.text.x = element_text(angle = 305))+
  scale_x_discrete(limits = rev)

ggsave("figure_2.pdf",figure_2)

figure_2A <- ggplot(filter(plastic_df_calc_melt_all_0_prop, Size_fraction == "<20", station == "R"))+
  geom_col(aes(x = depth, y = value, fill = variable))+facet_wrap(~station)+coord_flip()+
  scale_x_discrete(limits = rev)

figure_2B <- ggplot(filter(plastic_df_calc_melt_all_0_prop, Size_fraction == "<20", station == "PS15"))+
  geom_col(aes(x = depth, y = value, fill = variable), position = "stack")+facet_wrap(~station)+facet_wrap(~station)+coord_flip()+
  scale_x_discrete(limits = rev)

figure_2C <- ggplot(filter(plastic_df_calc_melt_all_0_prop, Size_fraction == ">20", station == "R"))+
  geom_col(aes(x = depth, y = value, fill = variable))+facet_wrap(~station)+coord_flip()+
  scale_x_discrete(limits = rev)

figure_2D <- ggplot(filter(plastic_df_calc_melt_all_0_prop, Size_fraction == ">20", station == "PS15"))+
  geom_col(aes(x = depth, y = value, fill = variable))+facet_wrap(~station)+coord_flip()+
  scale_x_discrete(limits = rev)


figure_2E <- ggplot(filter(plastic_df_calc_melt_all_0_prop, Size_fraction == ">100" & station == "R"))+
  geom_col(aes(x = depth, y = value, fill = variable))+facet_wrap(~station)+coord_flip()+
  scale_x_discrete(limits = rev)

#ggarrange(figure_1A,figure_1B,figure_1C, ncol = 1, common.legend = TRUE, legend = "bottom")

ggsave("figure_2A.pdf",figure_2A)
ggsave("figure_2B.pdf",figure_2B)
ggsave("figure_2C.pdf",figure_2C)
ggsave("figure_2D.pdf",figure_2D)
ggsave("figure_2E.pdf",figure_2E)

# Figure 3 ------------
# NMDS of the community with A (All size fractions) *DONE*
figure_3A <- ggplot(data = all_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(size = depth_m, shape = station, color = Size_fraction))+
  #geom_segment(data = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),], aes(x = 0, y = 0, xend = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),1], yend = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),2]))+
  #geom_text_repel(data = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),], aes(x = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),1], y = all_NMDS_fit_scores_vectors[c(3,5,19,20,23),2]-0.04), fontface = "bold", label = rownames(all_NMDS_fit_scores_vectors[c(3,5,19,20,23),]))+
  geom_segment(data = all_NMDS_fit_scores_factors[c(12:14),], aes(x = 0, y = 0, xend = all_NMDS_fit_scores_factors[c(12:14),1], yend = all_NMDS_fit_scores_factors[c(12:14),2]))+
  geom_text_repel(data = all_NMDS_fit_scores_factors[c(12:14),], aes(x = all_NMDS_fit_scores_factors[c(12:14),1], y = all_NMDS_fit_scores_factors[c(12:14),2]-0.04), fontface = "bold", label = rownames(all_NMDS_fit_scores_factors[c(12:14),]))+
  geom_segment(data = all_NMDS_species, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2/2), color = "red")+
  geom_text_repel(data = all_NMDS_species, aes(x = NMDS1, y = NMDS2/2), label = rownames(all_NMDS_species))
# B (G 20 size fraction) *DONE*;
figure_3B <- ggplot(data = G20_NMDS_scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = depth_m, shape = station), position = position_jitter(0.1))+
  geom_segment(data = G20_NMDS_fit_scores_vectors[c(9,10,13),], aes(x = 0, y = 0, xend = G20_NMDS_fit_scores_vectors[c(9,10,13),1], yend = G20_NMDS_fit_scores_vectors[c(9,10,13),2]))+
  geom_text_repel(data = G20_NMDS_fit_scores_vectors[c(9,10,13),], aes(x = G20_NMDS_fit_scores_vectors[c(9,10,13),1], y = G20_NMDS_fit_scores_vectors[c(9,10,13),2]-0.04), fontface = "bold", label = rownames(all_NMDS_fit_scores_vectors[c(9,10,13),]))+
  geom_segment(data = G20_NMDS_species, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2/2), color = "red")+
  geom_text_repel(data = G20_NMDS_species, aes(x = NMDS1, y = NMDS2/2), label = rownames(G20_NMDS_species))#+
# C) NMDS of the smallest size fraction *DONE* 
figure_3C <- ggplot(data = all_NMDS_scores_less20, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = depth_m, shape = station), position = position_jitter(0.1))+
  geom_segment(data = all_NMDS_species_less20, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2/2), color = "red")+
  geom_text_repel(data = all_NMDS_species_less20, aes(x = NMDS1, y = NMDS2/2), label = rownames(all_NMDS_species_less20))+
  geom_text_repel(size = 3,label = rownames(all_NMDS_scores_less20))

ggsave("figure_3A.pdf", figure_3A)
ggsave("figure_3B.pdf", figure_3B)
ggsave("figure_3C.pdf", figure_3C)
# No significant impacts of env vars --- power is too low to detect impacts would need more samples

# G100 is not usable as there are too few samples
# I will use the hellinger transformed data as this reduces the influence of zero inflated data

# No PERMANOVA will be done as there isn't enough variation explained via these nMDS tests

# Figure 4 --------
# Abundance totals by industry
figure_4A <- ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == "<20"))+
  geom_col(aes(x = depth, y = value, fill = end_use))+
  facet_wrap(station~Size_fraction, nrow = 1)+
  theme(axis.text.x = element_text(angle = 305))+
  scale_x_discrete(limits = rev)+coord_flip()+ylab(expression("Particles L"^-1))
figure_4B <- ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == ">20"))+
  geom_col(aes(x = depth, y = value, fill = end_use))+
  facet_wrap(station~Size_fraction, nrow = 1)+
  theme(axis.text.x = element_text(angle = 305))+
  scale_x_discrete(limits = rev)+coord_flip()+ylab(expression("Particles L"^-1))
figure_4C <- ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == ">100" & station == "R"))+
  geom_col(aes(x = depth, y = value, fill = end_use))+
  facet_wrap(station~Size_fraction, nrow = 1)+
  theme(axis.text.x = element_text(angle = 305))+
  scale_x_discrete(limits = rev)+coord_flip()+ylab(expression("Particles L"^-1))

figure_4A <- ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == "<20"))+
  geom_col(aes(x = depth, y = value, fill = Den_class_new))+
  facet_wrap(station~Size_fraction, nrow = 1)+
  theme(axis.text.x = element_text(angle = 305))+
  scale_x_discrete(limits = rev)+coord_flip()+ylab(expression("Particles L"^-1))

figure_4A <- ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == "<20"))+
  geom_col(aes(x = depth, y = value, fill = variable))+
  facet_wrap(station~Size_fraction, nrow = 1)+
  theme(axis.text.x = element_text(angle = 305))+
  scale_x_discrete(limits = rev)+coord_flip()+ylab(expression("Particles L"^-1))

ggsave("figure_4A.pdf", figure_4A)
ggsave("figure_4B.pdf", figure_4B)
ggsave("figure_4C.pdf", figure_4C)

# Figure 5 **OPTIONAL** could move to supplemental --------
# Abundance proportions by industry   ***NEED TO FIX ***

#This doesn't work as nice
#ggplot(plastic_df_calc_melt_all_0_prop)+
#  geom_line(aes(x = depth_m, y = value, color = end_use), linewidth = 1)+
#  coord_flip()+
#  scale_y_continuous(position = "right", label = c("0",".25", ".5", ".75","1"))+
#  scale_x_reverse()+facet_wrap(station~Size_fraction)

figure_5 <- ggplot(plastic_df_calc_melt_all_0_prop)+
  geom_col(aes(x = depth, y = value, fill = end_use))+
  coord_flip()+
  scale_y_continuous(position = "right", label = c("0",".25", ".5", ".75","1"))+
  scale_x_discrete(limits = rev)+facet_wrap(station~Size_fraction)

ggsave("figure_5.pdf", figure_5)

# Total estimated plastics heavily dependent upon size
# Order of magnitude differences between each size fraction - strongest determinant 
# Composition is similar between the two smaller size fraction but significantly different 
# at the largest size fraction
# When we delve deeper into particular size fractions
#   % trans and par are signifant vars for predicting compositional differences between sites

# There are significant differences with year and season.
# Size fraction is the largest determinant
# There are certain env parameters associated with specific polymers


# Summary tables for plastics (A) and env data and end_use (B) categories

# Supplemental 
# Figure S1 - env data PCA --------
figure_S1 <- ggplot(data = env_df_PCA_scores, aes(x = PC1, y = PC2))+
  geom_point(aes(color = depth_m, shape = station), size = 3)+
  geom_segment(data = env_df_PCA_species, aes(x = 0, y = 0, xend = PC1, yend = PC2), color = "red", arrow = arrow(length = unit(0.1,"inches")))+
  geom_text_repel(data = env_df_PCA_species, aes(x = PC1, y = PC2), label = rownames(env_df_PCA_species))+
  xlab("PC1 (48%)")+ylab("PC2 (33%)")+scale_color_gradient(low = "#56B1F7", high = "#132B43")

ggsave("figure_S1.pdf", figure_S1)

# Figure S2 Correlation matrix ------------
# I need a cor matrix for each of my NMDS plots
# I need 2 more correlation matrices < 20 and > 20 **DONE**
# Trim the correlation matrix to be more readable and informed
# x axis to be resonses, y axis to be explanatory vars

#vector of reponse var names
response_name <- c("ABS_L","ABS_prop","PET_L","PET_prop","PE_L","PE_prop",
                   "PP_L","PP_prop","PA_L","PA_prop", "Heavy_L","Heavy_prop","Light_L","Light_prop", "total_L")
#vector of explanatory var names
explanatory_name <- c(colnames(env_df_forPCA),"Ze","ZCM","timeJ",
                      "ZML.Dens..0.023Kg.m3.","UML..based.on.buyancy.frequency.N.")

explanator_name_final <- c("depth_m", "CStarTr0", "flECO.AFL","Density","ZML.Dens..0.023Kg.m3.")

cor_matrix_melt_allfilt <- filter(cor_matrix_melt, X1 %in% response_name, X2 %in% explanator_name_final)
cor_matrix_melt_allfilt$X1 <- factor(cor_matrix_melt_allfilt$X1, levels = response_name)

figure_S2A <- ggplot(cor_matrix_melt_allfilt, aes(X1,X2))+
  geom_tile(aes(fill = value, color = pvalue<0.05), lwd = 0.5, height = 0.75, width = 0.95)+scale_color_manual(name ='pvalue<0.05', values = setNames(c('black','White'), c(T,F)))+scale_fill_gradient2()+
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 6))+
  geom_text(label = round(cor_matrix_melt_allfilt$value, 3), size = 1.5)+
  scale_y_discrete(guide = guide_axis(angle = 45))+theme(axis.title.x = element_blank(),
                                                         axis.title.y = element_blank())+
  coord_flip()
# G20 cor matrix
cor_matrix_melt_G20filt <- filter(cor_matrix_melt_hell_G20, X1 %in% response_name, X2 %in% explanator_name_final)
cor_matrix_melt_G20filt$X1 <- factor(cor_matrix_melt_G20filt$X1, levels = response_name)

figure_S2B <- ggplot(cor_matrix_melt_G20filt, aes(X1,X2))+
  geom_tile(aes(fill = value, color = pvalue<0.05), lwd = 0.5, height = 0.75, width = 0.95)+scale_color_manual(name ='pvalue<0.05', values = setNames(c('black','White'), c(T,F)))+scale_fill_gradient2()+
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 6))+
  geom_text(label = round(cor_matrix_melt_G20filt$value, 3), size = 1.5)+
  scale_y_discrete(guide = guide_axis(angle = 45))+theme(axis.title.x = element_blank(),
                                                         axis.title.y = element_blank())+
  coord_flip()
# L20 cor matrix
cor_matrix_melt_L20filt <- filter(cor_matrix_melt_hell_L20, X1 %in% response_name, X2 %in% explanator_name_final)

cor_matrix_melt_L20filt$X1 <- factor(cor_matrix_melt_L20filt$X1, levels = response_name)

figure_S2C <- ggplot(cor_matrix_melt_L20filt, aes(X1,X2))+
  geom_tile(aes(fill = value, color = pvalue<0.05), lwd = 0.5, height = 0.75, width = 0.95)+scale_color_manual(name ='pvalue<0.05', values = setNames(c('black','White'), c(T,F)))+scale_fill_gradient2()+
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 6))+
  geom_text(label = round(cor_matrix_melt_L20filt$value, 3), size = 1.5)+
  scale_y_discrete(guide = guide_axis(angle = 45))+theme(axis.title.x = element_blank(),
                                                         axis.title.y = element_blank())+
  coord_flip()

# G100 cor matrix
#Not all polymers were ID'd for this size class
response_name_G100 <- c("PET_L","PET_prop","PE_L","PE_prop","PA_L","PA_prop", "Heavy_L","Heavy_prop","Light_L","Light_prop", "total_L")
explanatory_name_G100 <- c(colnames(env_df_forPCA))
cor_matrix_melt_G100filt <- filter(cor_matrix_melt_hell_G100, X1 %in% response_name_G100, X2 %in% explanator_name_final)

cor_matrix_melt_G100filt$X1 <- factor(cor_matrix_melt_G100filt$X1, levels = response_name_G100)

figure_S2D <- ggplot(cor_matrix_melt_G100filt, aes(X1,X2))+
  geom_tile(aes(fill = value, color = pvalue<0.05), lwd = 0.5, height = 0.75, width = 0.95)+scale_color_manual(name ='pvalue<0.05', values = setNames(c('black','White'), c(T,F)))+scale_fill_gradient2()+
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 6))+
  geom_text(label = round(cor_matrix_melt_G100filt$value, 3), size = 1.5)+
  scale_y_discrete(guide = guide_axis(angle = 45))+theme(axis.title.x = element_blank(),
                                                         axis.title.y = element_blank())+
  coord_flip()
#**NOTE** There are greys in the proportions because there is only 1 station
#This means there is no seasonal information to glean from this analysis


ggsave("figure_S2A.pdf", figure_S2A)
ggsave("figure_S2B.pdf", figure_S2B)
ggsave("figure_S2C.pdf", figure_S2C)
ggsave("figure_S2D.pdf", figure_S2D)

#all
figure_S2A
# G20
figure_S2B
# L20
figure_S2C
# G 100
figure_S2D

# Figure S3 Abundance by polymer type

figure_S3A <- ggplot(plastic_df_calc_melt_all_0)+
  geom_col(aes(x = variable, y = value, fill = station), position = "stack")+
  xlab("Polymer Type")+ylab(expression("Particles L"^-1))

figure_S3B <- ggplot(plastic_df_calc_melt_all_0)+
  geom_col(aes(x = end_use, y = value, fill = station), position = "stack")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  xlab("Plastic End Use")+ylab(expression("Particles L"^-1))

ggsave("figure_S3A.pdf", figure_S3A)
ggsave("figure_S3B.pdf", figure_S3B)

# What I need to change for paper:

# add in density classification scheme
# merge plots of 4 and 5 to 1 and 2
# remove nMDS and add in linear models and scatter plots of most important 
# relationships 
# Trim correlation matrix for just the vars needed

# New figures configuration -----------
# Figure 1
plastic_df_calc_melt_all_0_prop$Den_class <- with(plastic_df_calc_melt_all_0_prop, 
                                                ifelse(variable == "ABS_prop" | variable == "PET_prop" | variable == "PA_prop", "Denser","Lighter"))

plastic_df_calc_melt_all_0$Den_class <- with(plastic_df_calc_melt_all_0, 
                                             ifelse(variable == "ABS_L" | variable == "PET_L" | variable == "PA_L", "Denser","Lighter"))


ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == "<20"))+
  geom_col(aes(x = depth, y = value, fill = Den_class))+
  facet_wrap(station~Size_fraction, nrow = 1)+
  theme(axis.text.x = element_text(angle = 305))+
  scale_x_discrete(limits = rev)+coord_flip()+ylab(expression("Particles L"^-1))

ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == ">20"))+
  geom_col(aes(x = depth, y = value, fill = Den_class))+
  facet_wrap(station~Size_fraction, nrow = 1)+
  theme(axis.text.x = element_text(angle = 305))+
  scale_x_discrete(limits = rev)+coord_flip()+ylab(expression("Particles L"^-1))

 ggplot(filter(plastic_df_calc_melt_all_0, Size_fraction == ">100" & station == "R"))+
  geom_col(aes(x = depth, y = value, fill = Den_class))+
  facet_wrap(station~Size_fraction, nrow = 1)+
  theme(axis.text.x = element_text(angle = 305))+
  scale_x_discrete(limits = rev)+coord_flip()+ylab(expression("Particles L"^-1))
 
 
 ggplot(plastic_df_calc_melt_all_0_prop)+
   geom_col(aes(x = depth, y = value, fill = Den_class))+
   coord_flip()+
   scale_y_continuous(position = "right", label = c("0",".25", ".5", ".75","1"))+
   scale_x_discrete(limits = rev)+facet_wrap(station~Size_fraction)
 
ggplot(filter(plastic_df_calc_melt_all_0_prop, station == "R"))+
   geom_col(aes(x = depth, y = value, fill = Den_class))+
   coord_flip()+
   scale_y_continuous(position = "right", label = c("0",".25", ".5", ".75","1"))+
   scale_x_discrete(limits = rev)+facet_wrap(station~Size_fraction)

ggplot(filter(plastic_df_calc_melt_all_0_prop, station == "PS15", Size_fraction != ">100"))+
  geom_col(aes(x = depth, y = value, fill = Den_class))+
  coord_flip()+
  scale_y_continuous(position = "right", label = c("0",".25", ".5", ".75","1"))+
  scale_x_discrete(limits = rev)+facet_wrap(station~Size_fraction)

# Distinct seasonal/year differences
# spring saw much more in the total plastic abundance
# Size fraction saw significant differences
# No strong env signal
# Ice coverage not added but would be a great addition

ggplot(filter(plastic_df_calc_melt_all_0, station == "R" & Size_fraction == "<20"))+
  geom_point(aes(x = Density, y = Den_class, color = value))

# Size matters - important for influence of "community composition"
#   - the < 100 µm and < 20 µm communities are similar and suggest the larger feeds the smaller
#   - the largest size fraction is significantly different suggesting a different source/origin
# Generalists are the most abundant polymer type
#   - the types of plastic are the most abundant with the most generalized use plastics
#   - However, in the largest size fraction - packaging is the most common end use
#       - Which is also true for the density class
# Seasonality matters but what drives them is unclear
#   - Potential drivers being density and % transmission
#   - Could the main difference be currents? mixing? Ice formation? more research is needed
#       - several explanations remain but this is speculative...some lit review needed here
#       - maybe larger sizes fall out while smaller can remain suspended...reynolds # Q


# Notes from JD to fix
# Merge Figures from end use and total plastic production
# merge plots of 4 and 5 to 1 and 2

# is there autocorrelation between end use and density? 
# *NEED to tease that out before settling on this use*
# All of the light density particles are generalists, while the heavy are spred throughout
# the other send_use classes
# **DONE**

# I need to calculate detection limits for each size fraction
#       **DONE**

# Need a statistical test to compare the 0s in the dataset
#       **DONE**

# Maybe I need to rethink how I classified the density classes
# I only used the values and compared to general saltwater values but I didn't use the
# the density at that specific depth (which I have) to create the class.
# Directly calculating at the specific depth is more accurate and would be more denfendable
#       **DONE**

# Trim vars in the environmental PCA
#       **DONE**

# Trim correlation matrix for just the vars needed
#       **DONE**

# remove nMDS and add in linear models with scatter plots of most important 
# relationships 

# current 2 main issues with the study have to do with the ID process
# 1) the threshold I used is lower than todays standard
# 2) I did not acid digest my samples prior to ID

# I'm okay with explaining the lower threshold as I'm using an older technique
# I need to confer with Bob and JD about possible issues but Blanks should remove any 
# contaminations AND there are organic molecules in the opensource DB I used.

# Fixes from JD's notes --------------
# Dynamic calculation of density class

#first need to add density of each polymer type and then perform a simple logical
# P densities are
# ABS - 1.05 g/cm^3
# PE - 0.94 g/cm^3
# PP - 0.9 g/cm^3
# PET - 1.35 g/cm^3
# PA - 1.14 g/cm^3

plastic_df_calc_melt_all_0_prop$poly_density <- with(plastic_df_calc_melt_all_0_prop, 
                                                  ifelse(variable == "ABS_prop", 1.05,
                                                         ifelse(variable == "PET_prop", 1.35,
                                                                ifelse(variable == "PA_prop", 1.14,
                                                                       ifelse(variable == "PE_prop", 0.94, 0.9)))))

plastic_df_calc_melt_all_0$poly_density <- with(plastic_df_calc_melt_all_0, 
                                                ifelse(variable == "ABS_L", 1.05,
                                                       ifelse(variable == "PET_L", 1.35,
                                                              ifelse(variable == "PA_L", 1.14,
                                                                     ifelse(variable == "PE_L", 0.94, 0.9)))))


plastic_df_calc_melt_all_0_prop$Den_class_new <- with(plastic_df_calc_melt_all_0_prop, 
                                                     ifelse(Density < 1000*poly_density, "heavy", "light"))

plastic_df_calc_melt_all_0$Den_class_new <- with(plastic_df_calc_melt_all_0, 
                                                      ifelse(Density < 1000*poly_density, "heavy", "light"))

# Chi-squared comparison for the number of 0s between sites
# important for microplastic total number discussion
# Can handle differing number of samples within sites
# expected lower

# Num of 0s for each group
# Stn R
# 0s: 1            non-0s: 14
# Stn PS15
# 0s: 13            non-0s: 7

contig_tb <- matrix(c(1, 14, 13, 7), nrow = 2, byrow = TRUE)
rownames(contig_tb) <- c("Stn R", "Stn PS15")
colnames(contig_tb) <-  c("Zeros", "Non_Zero")

chisq.test(contig_tb)

fisher.test(contig_tb)

# Repeat for size fraction
contig_tb_Size <- matrix(c(0, 7, 5, 9, 9, 5), nrow = 3, byrow = TRUE)
rownames(contig_tb_Size) <- c("<20", ">20",">100")
colnames(contig_tb_Size) <-  c("Zeros", "Non_Zero")

chisq.test(contig_tb_Size)

fisher.test(contig_tb_Size)


# In both tests, we reject the null hypothesis and there is a significant difference
# in the ratio of 0s in the dataset between the 2 sites.

# Detection limit calculation for each size fraction
# how: 1/ fraction of filter * vol filtered
# <20 µm - 1/ 0.00178 filter * 4L = 140.25 particles per L
# >20 µm - 1/ 0.0222 filter * 20L = 2.25 particles per L
# >100 µm - 1/0.5 filter*20L = 0.1 particles per L

# Checking Autocorrelation between end-use and density
ggplot()+
  geom_point(data = filter(plastic_df_calc_melt_all_0_prop, Size_fraction != ">100"), aes(x = Density, y = value, color = Den_class_new))

ggplot()+
  geom_point(data = filter(plastic_df_calc_melt_all_0_prop, Size_fraction == ">100"), aes(x = Density, y = value, color = Den_class_new))

ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0_prop, aes(x = Density, y = value, color = Den_class_new))

ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0_prop, aes(x = depth_m, y = value, color = Den_class_new))

# Using the most important vars for linear regression

Figure_3new <- ggplot()+
  geom_boxplot(data = filter(plastic_df_calc_melt_all_0, value != 0), aes(x = Size_fraction, y = log(value+1), fill = station))+
  ylab("Microplastic Abundance (log)")+xlab("Size Fraction")+scale_x_discrete(labels = c("<20", "20 - 100", ">100"))

ggsave("Figure_3new.pdf", Figure_3new)

#Not useful
#ggplot()+
#  geom_boxplot(data = filter(plastic_df_calc_melt_all_0, value != 0), aes(x = variable, y = log(value+1), fill = station))+
#  ylab("Microplastic Abundance (log)")+xlab("Size Fraction")+facet_grid(cols = vars(Size_fraction))+theme(axis.text.x = element_text(angle = 45))

#Doesn't work as intended...scale is the same when using a wrap
#ggplot()+
#  geom_boxplot(data = filter(plastic_df_calc_melt_all_0, value != 0), aes(x = variable, y = log(value+1)))+
#  ylab("Microplastic Abundance (log)")+xlab("Size Fraction")+facet_wrap(~Size_fraction, ncol = 1)

Figure_4new <- ggplot()+
  geom_boxplot(data = filter(plastic_df_calc_melt_all_0_prop, value != 0), aes(x = variable, y = value, fill = station))+
  ylab("Microplastic Proportional Composition")+xlab("Size Fraction")+facet_grid(cols = vars(Size_fraction))+theme(axis.text.x = element_text(angle = 45))

ggsave("Figure_4new.pdf", Figure_4new)

ggplot()+
  geom_boxplot(data = filter(plastic_df_calc_melt_all_0_prop, value != 0), aes(x = Size_fraction, y = log(value+1), fill = variable))+
  ylab("Microplastic Abundance (log)")+xlab("Size Fraction")

ggplot()+
  geom_boxplot(data = filter(plastic_df_calc_melt_all_0_prop, value != 0 & Size_fraction == ">100"), aes(x = variable, y = log(value+1)))+
  ylab("Microplastic Abundance (log)")+xlab("Size Fraction")

ggplot()+
  geom_boxplot(data = filter(plastic_df_calc_melt_all_0_prop, value != 0 & Size_fraction == ">20"), aes(x = variable, y = log(value+1)))+
  ylab("Microplastic Abundance (log)")+xlab("Size Fraction")

ggplot()+
  geom_boxplot(data = filter(plastic_df_calc_melt_all_0_prop, value != 0 & Size_fraction == "<20"), aes(x = variable, y = log(value+1)))+
  ylab("Microplastic Abundance (log)")+xlab("Size Fraction")

ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0_prop, aes(x = Size_fraction, y = value))

summary(lm(log(value+1) ~ Size_fraction, plastic_df_calc_melt_all_0))
summary(aov(log(value+1) ~ Size_fraction, plastic_df_calc_melt_all_0))
summary(aov(log(value+1) ~ depth_m, plastic_df_calc_melt_all_0))

summary(lm(log(value+1) ~ depth_m, filter(plastic_df_calc_melt_all_0, Size_fraction == ">100")))

ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0, aes(x = end_use, y = log(value+1)))

ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0, aes(x = Den_class_new, y = log(value+1), color = variable))

ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0, aes(x = CStarTr0, y = log(value+1)))

ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0, aes(x = depth_m, y = log(value+1), color = Size_fraction))

ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0, aes(x = depth_m, y = log(value+1), color = variable))+facet_wrap(~Size_fraction)

Figure_S3new <- ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0, aes(x = depth_m, y = log(value+1), color = station))+
  xlab("Depth (m)")+ylab("Microplastic Abundance (log)")+
  facet_grid(rows = vars(variable), cols = vars(Size_fraction))+coord_flip()+scale_x_reverse()

ggsave("Figure_S3new.pdf", Figure_S3new)

ggplot()+
  geom_point(data = filter(plastic_df_calc_melt_all_0, Size_fraction != "<20"), aes(x = depth_m, y = log(value+1)))+facet_grid(rows = vars(variable), cols = vars(Size_fraction))

Figure_S4Anew <- ggplot()+
  geom_point(data = filter(plastic_df_calc_melt_all_0, value != 0), aes(x = depth_m, y = log(value+1)))+facet_grid(cols = vars(Size_fraction, station))

ggsave("Figure_S4Anew.pdf", Figure_S4Anew)

Figure_S4Bnew <- ggplot()+
  geom_point(data = filter(plastic_df_calc_melt_all_0, value != 0), aes(x = depth_m, y = log(value+1)))

ggsave("Figure_S4Bnew.pdf", Figure_S4Bnew)

ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0, aes(x = depth_m, y = log(value+1), color = variable, shape = Size_fraction))

ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0_prop, aes(x = depth_m, y = value, color = variable, shape = Size_fraction))

ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0_prop, aes(x = Size_fraction, y = value, color = variable))

# Summary stats by abundance - not needed as they can't be compared ----------

plastic_df_calc_melt_all_0_grp <- group_by(plastic_df_calc_melt_all_0, station, depth, Size_fraction, variable)

#plastic_df_calc_melt_all_0_grp <- group_by(plastic_df_calc_melt_all_0, depth, Size_fraction, variable)

plastic_df_calc_melt_all_0_grp_sum_poly <- summarise(plastic_df_calc_melt_all_0_grp, sum1 = sum(value), avg = mean(value),avg_log = mean(log(value+1)),
                                                     min = min(value), max = max(value), depth_m = mean(depth_m))

plastic_df_calc_melt_all_0_grp_sum_SizeFraction <- summarise(plastic_df_calc_melt_all_0_grp_sum_poly, sum1 = sum(sum1), avg = mean(avg), avg_log = mean(avg_log),
                                                             min = min(min), max = max(max), depth_m = mean(depth_m))

plastic_df_calc_melt_all_0_grp_sum_Depth <- summarise(plastic_df_calc_melt_all_0_grp_sum_SizeFraction, sum1 = sum(sum1),avg = mean(avg),avg_log = mean(avg_log),
                                                      min = min(min), max = max(max))

# I need to recalculate these numbers on a per sample basis without attributing to a polymer
# Need to Sum the values for each sample first and then perform the summary statistics
#This first is an avg across the polymer estimate types per station and depth
plastic_df_calc_melt_all_0_grp_SF_bystnDepth <- group_by(plastic_df_calc_melt_all_0, Size_fraction,station, depth, variable)

plastic_df_calc_melt_all_0_grp_sum_SF_bystnDepth <- summarise(plastic_df_calc_melt_all_0_grp_SF_bystnDepth, sum = sum(value))

plastic_df_calc_melt_all_0_grp_sum_SF_bystnDepth1 <- summarise(plastic_df_calc_melt_all_0_grp_sum_SF_bystnDepth, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                               min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_SF_bystn <- group_by(plastic_df_calc_melt_all_0, Size_fraction,station, depth)

plastic_df_calc_melt_all_0_grp_sum_SF_bystn <- summarise(plastic_df_calc_melt_all_0_grp_SF_bystn, sum = sum(value))

plastic_df_calc_melt_all_0_grp_sum_SF_bystn1 <- summarise(plastic_df_calc_melt_all_0_grp_sum_SF_bystn, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                          min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, station, variable)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_grp_sum_total_SF <- summarise(plastic_df_calc_melt_all_0_grp_total_SF, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_SF1 <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_SF, avg = sum(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                          min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_sum_total_SF2 <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_SF1, avg1 = mean(avg), sd = sum(sd), avg_log = mean(avg_log),
                                                          min = min(avg), max = max(avg))

plastic_df_calc_melt_all_0_grp_total_Stn <- group_by(plastic_df_calc_melt_all_0,station, Size_fraction,variable)

plastic_df_calc_melt_all_0_grp_sum_total_Stn <- summarise(plastic_df_calc_melt_all_0_grp_total_Stn, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_Stn1 <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_Stn, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                           min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_sum_total_Stn2 <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_Stn1, avg1 = sum(avg), sd = sum(sd), avg_log = mean(avg_log),
                                                           min = min(avg), max = max(avg))
#Repeat for each polymer
# ABS
plastic_df_calc_melt_all_0_grp_SF_bystn_ABS <- group_by(filter(plastic_df_calc_melt_all_0, variable == "ABS_L"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_grp_sum_SF_bystn_ABS <- summarise(plastic_df_calc_melt_all_0_grp_SF_bystn_ABS, sum = sum(value))

plastic_df_calc_melt_all_0_grp_sum_SF_bystn1_ABS <- summarise(plastic_df_calc_melt_all_0_grp_sum_SF_bystn_ABS, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                              min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_SF_ABS <- group_by(filter(plastic_df_calc_melt_all_0, variable == "ABS_L"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_grp_sum_total_SF_ABS <- summarise(plastic_df_calc_melt_all_0_grp_total_SF_ABS, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_SF1_ABS <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_SF_ABS, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                              min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_Stn_ABS <- group_by(filter(plastic_df_calc_melt_all_0, variable == "ABS_L"),station, depth)

plastic_df_calc_melt_all_0_grp_sum_total_Stn_ABS <- summarise(plastic_df_calc_melt_all_0_grp_total_Stn_ABS, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_Stn1_ABS <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_Stn_ABS, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                               min = min(sum), max = max(sum))

# PET
plastic_df_calc_melt_all_0_grp_SF_bystn_PET <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PET_L"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_grp_sum_SF_bystn_PET <- summarise(plastic_df_calc_melt_all_0_grp_SF_bystn_PET, sum = sum(value))

plastic_df_calc_melt_all_0_grp_sum_SF_bystn1_PET <- summarise(plastic_df_calc_melt_all_0_grp_sum_SF_bystn_PET, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                              min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_SF_PET <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PET_L"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_grp_sum_total_SF_PET <- summarise(plastic_df_calc_melt_all_0_grp_total_SF_PET, sum = sum(value))

plastic_df_calc_melt_all_0_grp_sum_total_SF1_PET <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_SF_PET, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                              min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_Stn_PET <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PET_L"),station, depth)

plastic_df_calc_melt_all_0_grp_sum_total_Stn_PET <- summarise(plastic_df_calc_melt_all_0_grp_total_Stn_PET, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_Stn1_PET <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_Stn_PET, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                               min = min(sum), max = max(sum))
# PA
plastic_df_calc_melt_all_0_grp_SF_bystn_PA <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PA_L"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_grp_sum_SF_bystn_PA <- summarise(plastic_df_calc_melt_all_0_grp_SF_bystn_PA, sum = sum(value))

plastic_df_calc_melt_all_0_grp_sum_SF_bystn1_PA <- summarise(plastic_df_calc_melt_all_0_grp_sum_SF_bystn_PA, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                             min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_SF_PA <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PA_L"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_grp_sum_total_SF_PA <- summarise(plastic_df_calc_melt_all_0_grp_total_SF_PA, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_SF1_PA <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_SF_PA, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                             min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_Stn_PA <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PA_L"),station, depth)

plastic_df_calc_melt_all_0_grp_sum_total_Stn_PA <- summarise(plastic_df_calc_melt_all_0_grp_total_Stn_PA, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_Stn1_PA <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_Stn_PA, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                              min = min(sum), max = max(sum))
# PE
plastic_df_calc_melt_all_0_grp_SF_bystn_PE <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PE_L"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_grp_sum_SF_bystn_PE <- summarise(plastic_df_calc_melt_all_0_grp_SF_bystn_PE, sum = sum(value))

plastic_df_calc_melt_all_0_grp_sum_SF_bystn1_PE <- summarise(plastic_df_calc_melt_all_0_grp_sum_SF_bystn_PE, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                             min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_SF_PE <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PE_L"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_grp_sum_total_SF_PE <- summarise(plastic_df_calc_melt_all_0_grp_total_SF_PE, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_SF1_PE <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_SF_PE, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                             min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_Stn_PE <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PE_L"),station, depth)

plastic_df_calc_melt_all_0_grp_sum_total_Stn_PE <- summarise(plastic_df_calc_melt_all_0_grp_total_Stn_PE, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_Stn1_PE <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_Stn_PE, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                              min = min(sum), max = max(sum))

# PP
plastic_df_calc_melt_all_0_grp_SF_bystn_PP <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PP_L"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_grp_sum_SF_bystn_PP <- summarise(plastic_df_calc_melt_all_0_grp_SF_bystn_PP, sum = sum(value))

plastic_df_calc_melt_all_0_grp_sum_SF_bystn1_PP <- summarise(plastic_df_calc_melt_all_0_grp_sum_SF_bystn_PP, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                             min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_SF_PP <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PP_L"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_grp_sum_total_SF_PP <- summarise(plastic_df_calc_melt_all_0_grp_total_SF_PP, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_SF1_PP <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_SF_PP, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                             min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_Stn_PP <- group_by(filter(plastic_df_calc_melt_all_0, variable == "PP_L"),station, depth)

plastic_df_calc_melt_all_0_grp_sum_total_Stn_PP <- summarise(plastic_df_calc_melt_all_0_grp_total_Stn_PP, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_Stn1_PP <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_Stn_PP, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                              min = min(sum), max = max(sum))

# Density class Heavy
plastic_df_calc_melt_all_0_grp_SF_bystn_H <- group_by(filter(plastic_df_calc_melt_all_0, Den_class_new == "heavy"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_grp_sum_SF_bystn_H <- summarise(plastic_df_calc_melt_all_0_grp_SF_bystn_H, sum = sum(value))

plastic_df_calc_melt_all_0_grp_sum_SF_bystn1_H <- summarise(plastic_df_calc_melt_all_0_grp_sum_SF_bystn_H, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                            min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_SF_H <- group_by(filter(plastic_df_calc_melt_all_0, Den_class_new == "heavy"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_grp_sum_total_SF_H <- summarise(plastic_df_calc_melt_all_0_grp_total_SF_H, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_SF1_H <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_SF_H, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                            min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_Stn_H <- group_by(filter(plastic_df_calc_melt_all_0, Den_class_new == "heavy"),station, depth)

plastic_df_calc_melt_all_0_grp_sum_total_Stn_H <- summarise(plastic_df_calc_melt_all_0_grp_total_Stn_H, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_Stn1_H <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_Stn_H, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                             min = min(sum), max = max(sum))

# Density class light
plastic_df_calc_melt_all_0_grp_SF_bystn_L <- group_by(filter(plastic_df_calc_melt_all_0, Den_class_new == "light"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_grp_sum_SF_bystn_L <- summarise(plastic_df_calc_melt_all_0_grp_SF_bystn_L, sum = sum(value))

plastic_df_calc_melt_all_0_grp_sum_SF_bystn1_L <- summarise(plastic_df_calc_melt_all_0_grp_sum_SF_bystn_L, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                            min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_SF_L <- group_by(filter(plastic_df_calc_melt_all_0, Den_class_new == "light"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_grp_sum_total_SF_L <- summarise(plastic_df_calc_melt_all_0_grp_total_SF_L, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_SF1_L <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_SF_L, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                            min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_grp_total_Stn_L <- group_by(filter(plastic_df_calc_melt_all_0, Den_class_new == "light"),station, depth)

plastic_df_calc_melt_all_0_grp_sum_total_Stn_L <- summarise(plastic_df_calc_melt_all_0_grp_total_Stn_L, sum = mean(value))

plastic_df_calc_melt_all_0_grp_sum_total_Stn1_L <- summarise(plastic_df_calc_melt_all_0_grp_sum_total_Stn_L, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                             min = min(sum), max = max(sum))

# Summary stats By proportion - not needed as they can't be compared -------
plastic_df_calc_melt_all_0_prop_grp <- group_by(plastic_df_calc_melt_all_0_prop, station, depth, Size_fraction, variable)

plastic_df_calc_melt_all_0_prop_grp <- group_by(plastic_df_calc_melt_all_0_prop, depth, Size_fraction, variable)

plastic_df_calc_melt_all_0_prop_grp_sum_poly <- summarise(plastic_df_calc_melt_all_0_prop_grp, avg = mean(value),
                                                          min = min(value), max = max(value))

plastic_df_calc_melt_all_0_prop_grp_sum_SizeFraction <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_poly, avg = mean(avg),
                                                                  min = min(min), max = max(max))

plastic_df_calc_melt_all_0_prop_grp_sum_Depth <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_SizeFraction, avg = mean(avg),
                                                           min = min(min), max = max(max))

# I need to recalculate these numbers on a per sample basis without attributing to a polymer
# Need to Sum the values for each sample first and then perform the summary statistics
#This first is an avg across the polymer estimate types per station and depth
#For proportion usng this to analyze is incorrect and unneccessary and incorrect approach
plastic_df_calc_melt_all_0_prop_grp_SF_bystnDepth <- group_by(plastic_df_calc_melt_all_0_prop, Size_fraction,station, depth, variable)

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystnDepth <- summarise(plastic_df_calc_melt_all_0_prop_grp_SF_bystnDepth, sum = sum(value))

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystnDepth1 <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystnDepth, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                    min = min(sum), max = max(sum))

#plastic_df_calc_melt_all_0_prop_grp_SF_bystn <- group_by(plastic_df_calc_melt_all_0_prop, Size_fraction,station, depth)

#plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn <- summarise(plastic_df_calc_melt_all_0_prop_grp_SF_bystn, sum = sum(value))

#plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn1 <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
#                                                          min = min(sum), max = max(sum))

#plastic_df_calc_melt_all_0_prop_grp_total_SF <- group_by(plastic_df_calc_melt_all_0_prop, Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
#plastic_df_calc_melt_all_0_prop_grp_sum_total_SF <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_SF, sum = sum(value))

#plastic_df_calc_melt_all_0_prop_grp_sum_total_SF1 <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_SF, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
#                                                          min = min(sum), max = max(sum))

#plastic_df_calc_melt_all_0_prop_grp_total_Stn <- group_by(plastic_df_calc_melt_all_0_prop,station, Size_fraction)

#plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_Stn, sum = sum(value))

#plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn1 <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
#                                                           min = min(sum), max = max(sum))
#Repeat for each polymer
# ABS
plastic_df_calc_melt_all_0_prop_grp_SF_bystn_ABS <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "ABS_prop"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_ABS <- summarise(plastic_df_calc_melt_all_0_prop_grp_SF_bystn_ABS, sum = sum(value))

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn1_ABS <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_ABS, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                   min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_SF_ABS <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "ABS_prop"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_ABS <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_SF_ABS, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_SF1_ABS <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_ABS, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                   min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_Stn_ABS <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "ABS_prop"),station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_ABS <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_Stn_ABS, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn1_ABS <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_ABS, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                    min = min(sum), max = max(sum))

# PET
plastic_df_calc_melt_all_0_prop_grp_SF_bystn_PET <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PET_prop"),station, Size_fraction, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_PET <- summarise(plastic_df_calc_melt_all_0_prop_grp_SF_bystn_PET, sum = sum(value))

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn1_PET <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_PET, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                   min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_SF_PET <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PET_prop"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_PET <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_SF_PET, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_SF1_PET <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_PET, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                   min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_Stn_PET <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PET_prop"),station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_PET <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_Stn_PET, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn1_PET <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_PET, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                    min = min(sum), max = max(sum))
# PA
plastic_df_calc_melt_all_0_prop_grp_SF_bystn_PA <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PA_prop"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_PA <- summarise(plastic_df_calc_melt_all_0_prop_grp_SF_bystn_PA, sum = sum(value))

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn1_PA <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_PA, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                  min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_SF_PA <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PA_prop"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_PA <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_SF_PA, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_SF1_PA <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_PA, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                  min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_Stn_PA <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PA_prop"),station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_PA <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_Stn_PA, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn1_PA <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_PA, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                   min = min(sum), max = max(sum))
# PE
plastic_df_calc_melt_all_0_prop_grp_SF_bystn_PE <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PE_prop"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_PE <- summarise(plastic_df_calc_melt_all_0_prop_grp_SF_bystn_PE, sum = sum(value))

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn1_PE <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_PE, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                  min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_SF_PE <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PE_prop"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_PE <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_SF_PE, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_SF1_PE <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_PE, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                  min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_Stn_PE <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PE_prop"),station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_PE <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_Stn_PE, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn1_PE <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_PE, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                   min = min(sum), max = max(sum))

# PP
plastic_df_calc_melt_all_0_prop_grp_SF_bystn_PP <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PP_prop"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_PP <- summarise(plastic_df_calc_melt_all_0_prop_grp_SF_bystn_PP, sum = sum(value))

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn1_PP <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_PP, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                  min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_SF_PP <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PP_prop"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_PP <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_SF_PP, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_SF1_PP <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_PP, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                  min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_Stn_PP <- group_by(filter(plastic_df_calc_melt_all_0_prop, variable == "PP_prop"),station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_PP <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_Stn_PP, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn1_PP <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_PP, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                   min = min(sum), max = max(sum))

# Density Class Heavy
plastic_df_calc_melt_all_0_prop_grp_SF_bystn_H <- group_by(filter(plastic_df_calc_melt_all_0_prop, Den_class_new == "heavy"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_H <- summarise(plastic_df_calc_melt_all_0_prop_grp_SF_bystn_H, sum = sum(value))

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn1_H <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_H, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                 min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_SF_H <- group_by(filter(plastic_df_calc_melt_all_0_prop, Den_class_new == "heavy"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_H <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_SF_H, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_SF1_H <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_H, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                 min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_Stn_H <- group_by(filter(plastic_df_calc_melt_all_0_prop,  Den_class_new == "heavy"),station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_H <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_Stn_H, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn1_H <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_H, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                  min = min(sum), max = max(sum))

# Density Class Light
plastic_df_calc_melt_all_0_prop_grp_SF_bystn_L <- group_by(filter(plastic_df_calc_melt_all_0_prop,  Den_class_new == "light"), Size_fraction,station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_L <- summarise(plastic_df_calc_melt_all_0_prop_grp_SF_bystn_L, sum = sum(value))

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn1_L <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_L, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                 min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_SF_L <- group_by(filter(plastic_df_calc_melt_all_0_prop,  Den_class_new == "light"), Size_fraction, depth)
# Only needed for calculation for the largest size fraction 
#plastic_df_calc_melt_all_0_grp_total_SF <- group_by(plastic_df_calc_melt_all_0, Size_fraction, depth)
plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_L <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_SF_L, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_SF1_L <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_SF_L, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                 min = min(sum), max = max(sum))

plastic_df_calc_melt_all_0_prop_grp_total_Stn_L <- group_by(filter(plastic_df_calc_melt_all_0_prop,  Den_class_new == "light"),station, depth)

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_L <- summarise(plastic_df_calc_melt_all_0_prop_grp_total_Stn_L, sum = mean(value))

plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn1_L <- summarise(plastic_df_calc_melt_all_0_prop_grp_sum_total_Stn_L, avg = mean(sum), sd = sd(sum), avg_log = mean(log(sum+1)),
                                                                  min = min(sum), max = max(sum))


plastic_df_calc_melt_all_0_prop_grp_SF_bystn_both <- group_by(plastic_df_calc_melt_all_0_prop, Size_fraction,station, depth, Den_class_new)

plastic_df_calc_melt_all_0_prop_grp_sum_SF_bystn_both <- summarise(plastic_df_calc_melt_all_0_prop_grp_SF_bystn_both, sum = sum(value))

plastic_df_calc_melt_all_0_grp_SF_bystn_both <- group_by(plastic_df_calc_melt_all_0, Size_fraction,station, depth, Den_class_new)

plastic_df_calc_melt_all_0_grp_sum_SF_bystn_both <- summarise(plastic_df_calc_melt_all_0_grp_SF_bystn_both, sum = sum(value))



write.table(as.data.frame(plastic_df_calc_melt_all_0_prop_grp_total_Stn_L), "plastic_df_calc_melt_all_0_prop_grp_total_Stn_L.csv", row.names = TRUE, sep = ",")
write.table(as.data.frame(plastic_df_calc_melt_all_0_prop), "plastic_df_calc_melt_all_0_prop.csv", row.names = TRUE, sep = ",")
# Summary tables for environmental vars

env_df_all_grp <- group_by(env_df_all,station)

env_df_all_grp_sum_avg <- summarise_if(env_df_all_grp, is.numeric, mean, na.rm = TRUE)

env_df_all_grp_sum_sd <- summarise_if(env_df_all_grp, is.numeric, sd, na.rm = TRUE)

write.table(as.data.frame(t(env_df_all_grp_sum_avg)), "env_df_all_grp_sum_avg.csv", row.names = TRUE, sep = ",")
write.table(as.data.frame(t(env_df_all_grp_sum_sd)), "env_df_all_grp_sum_sd.csv", row.names = TRUE, sep = ",")

# Figures post meeting with JD and BOB---------
# Making figure 1  - Merging the information from figure 1 with figure 2
# Make a point plot over top of the figure 2. Figure 2 remains mostly unchanged
# Make a continuous point plot and overlay the bar plot

Figure_S4Anew <- ggplot()+
  geom_point(data = filter(plastic_df_calc_melt_all_0, value != 0), aes(x = depth_m, y = log(value+1)))+facet_wrap(Size_fraction ~station)

#This is the summary dataset to use for Figure 1 with the proper information
Figure_1Anew <- ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0_grp_sum_SizeFraction, aes(x = depth_m, y = log(sum1)))+facet_wrap(station ~ Size_fraction)+
  coord_flip()+scale_x_reverse()

Figure_1Anew <- ggplot()+
  geom_point(data = plastic_df_calc_melt_all_0_grp_sum_SizeFraction, aes(x = depth_m, y = log(sum1+1)))+facet_grid(station ~ Size_fraction)+
  coord_flip()+scale_x_reverse()+
  annotate("text", x = -40, y = 0.25, label = "0", size= 3.5)+
  annotate("text", x = -40, y = 4, label = "0.5", size= 3.5)+
  annotate("text", x = -40, y = 8, label = "1", size = 3.5)

Figure_1Anew

ggsave("Figure_1Anew.pdf", Figure_1Anew)

# Making additional Figure 2 figures
# Reading in the parsed CTD data for the profiles
# Use depDM as depth vs Density calculations
stnR_ctd <- read.csv("/Users/christophercarnivale/Desktop/Dissertation_data/NBP_1910_R_ctd_data_edited.csv", header = TRUE)

stnPS15_ctd <- read.csv("/Users/christophercarnivale/Desktop/Dissertation_data/NBP_2205_PS15_ctd_data_edited.csv", header = TRUE)

stnPS15_ctd_NewD <- read.csv("/Users/christophercarnivale/Desktop/Dissertation_data/Raman Spectral Data copy/Density_PB_CC_edit_PS15.csv", header = TRUE)

stnR_ctd_filt <- stnR_ctd[,c("depSM","Density","t090C","CStarTr0")]
stnPS15_ctd_filt <- stnPS15_ctd[,c("depSM","density11","t090C","CStarTr0")]

stnR_ctd_filt$station <- "R"
stnPS15_ctd_filt$station <- "PS15"
stnPS15_ctd_NewD$station <- "PS15"
colnames(stnPS15_ctd_filt) <- c("depSM","Density","t090C","CStarTr0","station")
colnames(stnPS15_ctd_NewD) <- c("depSM","Density","station")

stnPS15_ctd_filt_new_D <- stnPS15_ctd_filt[-c(1:9),]

stnPS15_ctd_filt_new_D$Density <- stnPS15_ctd_NewD$Density

CTD_both <- rbind(stnR_ctd_filt,stnPS15_ctd_filt)


CTD_both_newD <- rbind(stnR_ctd_filt,stnPS15_ctd_filt_new_D)

Figure_2A_newAgain <- ggplot()+
  geom_line(data=CTD_both_newD, aes(x = depSM, y = Density, color = station), linewidth = 2)+
  coord_flip()+scale_x_reverse()+scale_y_continuous(position = "right")+xlab("Depth (m)")+
  ylab(expression(paste("Density (", rho, ")")))+labs(color = "Season")+
  scale_color_manual(values = c("#F8766D","#00BFC4"),labels = c("Early Winter","Late Spring"))+
  theme(legend.position = "none")+
  geom_vline(data = plastic_df_calc_melt_all_0, aes(xintercept = depth_m, color = station), linetype = "dashed")

ggsave("Figure_2A_newAgain.pdf", Figure_2A_newAgain)

Figure_2B_newAgain <- ggplot()+
  geom_line(data=CTD_both, aes(x = depSM, y = t090C, color = station), linewidth = 2)+
  coord_flip()+scale_x_reverse()+scale_y_continuous(position = "right")+xlab("Depth (m)")+
  ylab("Temperature (Celsius)")+labs(color = "Season")+
  scale_color_manual(values = c("#F8766D","#00BFC4"),labels = c("Early Winter","Late Spring"))+
  theme(legend.position = c(0.8,0.2))+
  geom_vline(data = plastic_df_calc_melt_all_0, aes(xintercept = depth_m, color = station), linetype = "dashed")

Figure_2C_newAgain <- ggplot()+
  geom_line(data=CTD_both, aes(x = depSM, y = CStarTr0, color = station), linewidth = 2)+
  coord_flip()+scale_x_reverse()+scale_y_continuous(position = "right")+xlab("Depth (m)")+
  ylab("Beam Transmission (%)")+labs(color = "Season")+
  scale_color_manual(values = c("#F8766D","#00BFC4"),labels = c("Early Winter","Late Spring"))+
  geom_vline(data = plastic_df_calc_melt_all_0, aes(xintercept = depth_m, color = station), linetype = "dashed")+
  theme(legend.position = c(0.2,0.2))

ggsave("Figure_2B_newAgain.pdf", Figure_2B_newAgain)
ggsave("Figure_2C_newAgain.pdf", Figure_2C_newAgain)

# Final Table Calculation: Range of % particles ID'd as Plastic-------
plastic_df_percentID <- read.csv("Env_Data_CC_wPlastic_newCF_blankCor_particlesampled.csv", header = TRUE)

plastic_df_percentID$depth <- factor(plastic_df_percentID$depth,
                                levels =c ("Surface","14m","15m","23m","40m","70m","110m","120m","150m","205m","250m","Deep"))


plastic_df_percentID$Size_fraction <- factor(plastic_df_percentID$Size_fraction,
                                        levels =c ("<20",">20",">100"))

plastic_df_percentID$station <- factor(plastic_df_percentID$station,
                                  levels =c ("R","PS15"))

plastic_df_percentID <- mutate(plastic_df_percentID,
                               percent_particle_ID = Total_plastic/Plastic.Sampled*100)

# Calculating all of the Mins and Maxs
# >100
min(filter(plastic_df_percentID, station == "R", Size_fraction == ">100")$percent_particle_ID)

max(filter(plastic_df_percentID, station == "R", Size_fraction == ">100")$percent_particle_ID)

min(filter(plastic_df_percentID, station == "PS15", Size_fraction == ">100")$percent_particle_ID)

max(filter(plastic_df_percentID, station == "PS15", Size_fraction == ">100")$percent_particle_ID)

# >20
min(filter(plastic_df_percentID, station == "R", Size_fraction == ">20")$percent_particle_ID)

max(filter(plastic_df_percentID, station == "R", Size_fraction == ">20")$percent_particle_ID)

min(filter(plastic_df_percentID, station == "PS15", Size_fraction == ">20")$percent_particle_ID)

max(filter(plastic_df_percentID, station == "PS15", Size_fraction == ">20")$percent_particle_ID)

#<20
min(filter(plastic_df_percentID, station == "R", Size_fraction == "<20")$percent_particle_ID)

max(filter(plastic_df_percentID, station == "R", Size_fraction == "<20")$percent_particle_ID)

min(filter(plastic_df_percentID, station == "PS15", Size_fraction == "<20")$percent_particle_ID)

max(filter(plastic_df_percentID, station == "PS15", Size_fraction == "<20")$percent_particle_ID)


