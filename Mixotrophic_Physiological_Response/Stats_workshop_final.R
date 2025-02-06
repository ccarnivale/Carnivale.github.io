# Final assignment: Mixotrophic physiological respose to microplastics
# Chris Carnivale
# Performing mixed models on feeding and primary production responses
# Data: mixo_mp_pp.csv (PP and growth data), Mixo_MP_feeding_csv.csv (bacterivory data),
#       mixo_mp_adundance (growth curve data)

#QUESTION: What are the physiological responses of mixotrophic organisms to 
#increasing amounts of microplastic pollution?

#H1:  Mixotrophic photosynthetic efficiency will be reduced and reduce growth.
#H2:  Mixotrophic organisms will exhibit false satiation or malnutrition causing 
# increased energy expenditure for bacterivory and reduce growth.

#Loading Packages----------
#require(tidyverse) not needed for this course
require(seacarb)
require(modelr)
require(phytotools)
require(magrittr)
require(ggpubr)
require(lme4)
require(lmerTest)
require(plotrix)
require(mgcv)
require(tidymv)
require(ggplot2)
require(tidyverse)

setwd("/Users/christophercarnivale/Desktop/Dissertation_data/Microplastic_physiological_response/")
#Reading in datasets (3) ---------
#PP
mixo_phys_pp <- read.csv("mixo_mp_pp.csv")
#Abundance
mixo_mp_abund <- read.csv("mixo_mp_abundance.csv")
#Bacterivory
mixo_mp_feed <- read.csv("Mixo_MP_feeding_csv.csv")

#Primary Production Tidying for amixo_phys_pp_baseR_calc--------
#PP
str(mixo_phys_pp)

summary(mixo_phys_pp)

mixo_phys_pp$Time_days <- factor(mixo_phys_pp$Time_days)

mixo_phys_pp$Rep <- factor(mixo_phys_pp$Rep)

mixo_phys_pp$Treatment_plastic <- factor(mixo_phys_pp$Treatment_plastic,
                                         levels = c("control", "10","1000", "100000"))

#Creating calculation columns needed from raw data w/o mutate
mixo_phys_pp_baseR_calc <- mixo_phys_pp
mixo_phys_pp_baseR_calc$pp_DPM <- mixo_phys_pp_baseR_calc$pp_CPM/mixo_phys_pp_baseR_calc$C_14_efficiency
mixo_phys_pp_baseR_calc$bkgr_DPM <- mixo_phys_pp_baseR_calc$background_CPM/mixo_phys_pp_baseR_calc$C_14_efficiency/0.1
mixo_phys_pp_baseR_calc$ugram_Cfixed_per_L_hr <- mixo_phys_pp_baseR_calc$pp_DPM*mixo_phys_pp_baseR_calc$DIC.umol.kg./mixo_phys_pp_baseR_calc$bkgr_DPM*12*1.029/2
mixo_phys_pp_baseR_calc$cell_per_L <- mixo_phys_pp_baseR_calc$Cell_counts/1.35*mixo_phys_pp_baseR_calc$Zeiss_CF*1000*1.08
mixo_phys_pp_baseR_calc$ugram_Cfixed_perCELL_perHR <- mixo_phys_pp_baseR_calc$ugram_Cfixed_per_L_hr/mixo_phys_pp_baseR_calc$cell_per_L
mixo_phys_pp_baseR_calc$pg_perCELL_perHR <- mixo_phys_pp_baseR_calc$ugram_Cfixed_perCELL_perHR*1000000

#0.1 for bckgr_CPM because 0.1 ml was processed
#12 is the molar mass of carbon o convert from umol to ugrams
#1.029 is to convert from kg to L of Seawater
#2 is to get an hourly uptake since it was exposed for 2 hours

#Can't use mutate()
#mixo_phy_pp_calc <- mixo_phys_pp %>% 
#  mutate(pp_DPM = pp_CPM/C_14_efficiency, 
#          bkgr_DPM = background_CPM/C_14_efficiency/0.1, 
#          ugram_Cfixed_per_L_hr = pp_DPM*DIC.umol.kg./bkgr_DPM*12*1.029/2)

#Abundance Tidying ---------
# Its now part of the PP dataset BUT not all 7 days
#Abundance
str(mixo_mp_abund)

mixo_mp_abund <- as.data.frame(mixo_mp_abund)

mixo_mp_abund$Rep <- factor(mixo_phys_pp$Rep)

mixo_mp_abund$Treatment <- factor(mixo_mp_abund$Treatment,
                                   levels = c("Control", "10","1000", "100000"))

mixo_mp_abund$cell_per_L <- mixo_mp_abund$Cell_counts/1.35*mixo_mp_abund$Zeiss_CF*1000*1.02

#Can't use mutate()
#mixo_mp_gc <- mixo_mp_abund %>% 
# mutate(cell_per_L = Cell_counts/1.35*Zeiss_CF*1000)

#mixo_mp_gc_sum <- mixo_mp_gc %>% 
#  group_by(Time_days, Treatment) %>% 
#  summarise(avg_cell_per_L = mean(cell_per_L))

#Bacterivory Tidying ----------
str(mixo_mp_feed)
mixo_mp_feed$treatment_rep <- factor(mixo_mp_feed$treatment_rep)

mixo_mp_feed$Time_days <- factor(mixo_mp_feed$Time_days)

mixo_mp_feed$Treatment_plastic <- factor(mixo_mp_feed$Treatment_plastic,
                                  levels = c("control", "10","1000", "100000"))

mixo_mp_feed$bac_conc_perML <- mixo_mp_feed$bac_avg/mixo_mp_feed$bac_micro_filtered*mixo_mp_feed$dil_factor*mixo_mp_feed$cor_factor_grid
mixo_mp_feed$micro_conc_perML <- mixo_mp_feed$micro_avg/mixo_mp_feed$bac_micro_filtered*mixo_mp_feed$dil_factor*mixo_mp_feed$cor_factor_grid
mixo_mp_feed$BM_ratio <- mixo_mp_feed$bac_conc_perML/mixo_mp_feed$micro_conc_perML
mixo_mp_feed$BSP_ingested <- mixo_mp_feed$conc_bkgrd_per*mixo_mp_feed$BM_ratio + mixo_mp_feed$conc_bkgrd_per
mixo_mp_feed$ing_percell_perhour <- mixo_mp_feed$BSP_ingested*2 #30 mins feeding time

mixo_mp_feed$cell_L <- mixo_mp_feed$cell_avg*10000*1000*mixo_mp_feed$dil_factor

#mixo_mp_feed_rates <- mixo_mp_feed %>% 
#  mutate(bac_conc_perML = bac_avg/bac_micro_filtered*dil_factor*cor_factor_grid,
#         micro_conc_perML = micro_avg/bac_micro_filtered*dil_factor*cor_factor_grid,
#         BM_ratio = bac_conc_perML/micro_conc_perML,
#         BSP_ingested = (conc_bkgrd_per*BM_ratio) + conc_bkgrd_per,
#         ing_percell_perhour = BSP_ingested*2)


#Data exploration: Growth Results ----------
#Growth curves
plot(cell_per_L ~ Time_days, data = mixo_mp_abund, type = "p", 
     col = mixo_mp_abund$Treatment, xlab = "Time (Days)", ylab = "Cells/mL",
     main = "Ochromonas abundance over 
    one week exposure to Microplastics")



#plot(mixo_mp_abund$Time_days, mixo_mp_abund$cell_per_L, type = "p",col = mixo_mp_abund$Treatment)

#lines(mixo_mp_abund$Time_days,mixo_mp_abund$cell_per_L, col = mixo_mp_abund$Treatment)
#This lines() doesn't work
#Need to make a line plot of the average values

#NOW this is a convoluted fix...BUT NOW I will convert Time_days into a factor
# in the dataset so that I can loop and calculate the means I need for this plot.

mixo_mp_abund$Time_days <- as.factor(mixo_mp_abund$Time_days)
mixo_mp_abund_avg_colnames <- c("Time_days","Treatment", "cell_per_L")
mixo_mp_gcPlot_means_fixed <- aggregate(mixo_mp_abund$cell_per_L, list(mixo_mp_abund$Time_days,mixo_mp_abund$Treatment), mean)
colnames(mixo_mp_gcPlot_means_fixed) <- mixo_mp_abund_avg_colnames
mixo_mp_gcPlot_means_fixed$Time_days <- as.numeric(as.character(mixo_mp_gcPlot_means_fixed$Time_days))

#Now subset for each line plot
mixo_mp_gcPlot_means_fixed_CONTROL <- subset(mixo_mp_gcPlot_means_fixed, Treatment == "Control")
mixo_mp_gcPlot_means_fixed_10 <- subset(mixo_mp_gcPlot_means_fixed, Treatment == "10")
mixo_mp_gcPlot_means_fixed_1000 <- subset(mixo_mp_gcPlot_means_fixed, Treatment == "1000")
mixo_mp_gcPlot_means_fixed_100000 <- subset(mixo_mp_gcPlot_means_fixed, Treatment == "100000")
#Checking default palette to match my new lines with the points plot
palette()
plot.new()

lines(mixo_mp_gcPlot_means_fixed_CONTROL$Time_days, mixo_mp_gcPlot_means_fixed_CONTROL$cell_per_L)
lines(mixo_mp_gcPlot_means_fixed_10$Time_days, mixo_mp_gcPlot_means_fixed_10$cell_per_L, col = "#DF536B")
lines(mixo_mp_gcPlot_means_fixed_1000$Time_days, mixo_mp_gcPlot_means_fixed_1000$cell_per_L, col = "#61D04F")
lines(mixo_mp_gcPlot_means_fixed_100000$Time_days, mixo_mp_gcPlot_means_fixed_100000$cell_per_L, col = "#2297E6")
legend("topleft", levels(mixo_mp_abund$Treatment), fill = c("black","#DF536B","#61D04F","#2297E6"))


#Growth Rate and doubling time Calculation and tabulation ------
#From abundance plots:  Days 2-4 is where largest growth occurs for all treatments
#mixo_mp_gc <- subset(mixo_mp_abund, Time_days == "2"| Time_days == "4")
mixo_mp_gcd2 <- subset(mixo_mp_abund, Time_days == "2")
mixo_mp_gcd4 <- subset(mixo_mp_abund, Time_days == "4")

mixo_mp_gc_2_4 <- cbind(mixo_mp_gcd2$cell_per_L,mixo_mp_gcd4$cell_per_L)
colnames(mixo_mp_gc_2_4) <- c("Day_2_Cells/L", "Day_4_Cells/L")
mixo_mp_gc_2_4 <- as.data.frame(mixo_mp_gc_2_4)
mixo_mp_gc_2_4$Treatment <- mixo_mp_gcd2$Treatment
mixo_mp_gc_2_4$Growth_Rate_perday <-log(mixo_mp_gc_2_4$`Day_4_Cells/L`/mixo_mp_gc_2_4$`Day_2_Cells/L`)/2
mixo_mp_gc_2_4$Growth_Rate_perday[16] <- 0
mixo_mp_gc_2_4$Doubling_time_perday <- log(2)/log(1+mixo_mp_gc_2_4$Growth_Rate_perday)
mixo_mp_gc_2_4$Doubling_time_perday[16] <- 0


Growth_Rate <- tapply(mixo_mp_gc_2_4$Growth_Rate_perday, mixo_mp_gc_2_4$Treatment, mean)

Doubling_time_adj <- log(2)/log(1+Growth_Rate)

Doubling_time_adj_se <- tapply(mixo_mp_gc_2_4$Doubling_time_perday, mixo_mp_gc_2_4$Treatment, function(x){sd(x)/sqrt(4)})

Doubling_Time <- tapply(mixo_mp_gc_2_4$Doubling_time_perday, mixo_mp_gc_2_4$Treatment, mean)

gc_sd <- tapply(mixo_mp_gc_2_4$Growth_Rate_perday, mixo_mp_gc_2_4$Treatment, sd)

gc_se <- tapply(mixo_mp_gc_2_4$Growth_Rate_perday, mixo_mp_gc_2_4$Treatment, function(x){sd(x)/sqrt(4)})

Plastic_factor_lvl_vector <- c("Control", "10", "1000","100000")

mixo_mp_gc_table <- t(rbind(Plastic_factor_lvl_vector, round(Growth_Rate,2),round(Doubling_time_adj,2)))
colnames(mixo_mp_gc_table) <- c("Treatment","Growth Rate", "Doubling Time")

#Plot table using base R graphics
plot(0,type='n',axes=FALSE,ann=FALSE)

addtable2plot(0.5,-0.20, mixo_mp_gc_table, bg = "grey",bty = "o", hlines = T, 
              vlines = T, cex = 1)

# Now I have negative values due to 1 Rep...Rep D of 1e5 Treatment.
#NEED to consider this later on in the project.
#NEGATIVE values indicate population decline and thus follow the trend of which
#I was expecting.
#Data exploration: Primary Production ----------
boxplot(pg_perCELL_perHR ~ Treatment_plastic*Time_days, mixo_phys_pp_baseR_calc,
        col = c("black","#DF536B","#61D04F","#2297E6"), xaxt = "n",
        xlab = "Time (Days)", ylab = "pg/(cell*hr)",
        main = "Ochromonas per cell Carbon fixation rates
        over one week")
axis(side = 1, at = c(2.5,6.5,10.5,14.5), labels = c("Day 0", "Day 2", "Day 4","Day 7"))
legend("topright", levels(mixo_phys_pp_baseR_calc$Treatment), fill = c("black","#DF536B","#61D04F","#2297E6"))


#Data exploration: Bacterivory ---------
boxplot(ing_percell_perhour ~ Treatment_plastic*Time_days, mixo_mp_feed,
        col = c("black","#DF536B","#61D04F","#2297E6"), xaxt = "n",
        xlab = "Time (Days)", ylab = "BSP ingestion/(cell*hr)",
        main = "Ochromonas bacterial sized ingestion (BSP)
        rates over one week")
axis(side = 1, at = c(2.5,6.5,10.5), labels = c("Day 0", "Day 3","Day 7"))
legend("topright", levels(mixo_phys_pp_baseR_calc$Treatment), fill = c("black","#DF536B","#61D04F","#2297E6"), cex = 0.7)

#Reasoning behind exploratory Graphs/Tables -------

#For abundance, I used a points/lines graph to explore the feeding over the 
#course of the week. In the analysis for this physiological response I'm not 
#as interested in the distribution of the abundances but finding the days
#at which growth looks to be at its maximum (max slope) for all treatments. 
#In this experiment, max growth is between 2-4 days.From there I calculate the 
# specific growth rates and create a table which is closest to being a reportable
#table I am capable in making in base r code (so for this project a results 
#visualization.

#For Primary Production and bacterivory, the graphs both shared a similar purporse.
#In these visualizations, I wanted to see what the relationship of each response
# to increasing microplastic contamination is to get an idea of the appropriate
#model to use to describe my results. In both, the control had a linear decreasing
#while all microplastic samples had a nonlinear response. Thus, changing my original
#plan os using a GLMM to a GAMM (which I'm hoping I get close enough to right on this,
#since we didn't get to this during the semester).










#Abundance Analysis -------
#I know these will be be pretty similar trends but I'm going to create 2 models:
# 1 for growth rate and 1 for doubling time

GR_aov <- aov(Growth_Rate_perday ~ Treatment, mixo_mp_gc_2_4)
summary(GR_aov)
model.tables(GR_aov)
summary(lm(GR_aov))
DT_aov <- aov(Doubling_time_perday ~ Treatment, mixo_mp_gc_2_4)
summary(DT_aov)
model.tables(DT_aov)
summary(lm(DT_aov))

#Assumptions Testing
#Assumptions of ANOVA analysis:
#     Homogeneity of variance
#     Normally distributed
#     Independence
#     Fixed X
plot(lm(GR_aov),1)
plot(lm(GR_aov),2)
plot(lm(GR_aov),3)
plot(lm(GR_aov),4)

#From these diagnostic plots, I would argue that my data do not violate the 
#normality assumption as shown by the QQ plot. From the Residuals vs. Fitted plot
# and the Scale-location it seems as though the assumption of homogeneity of
#variance may be violated. #There does seem to potentially be some points of
# high leverage but none that are influential to the model. Assumptions
#independence and Fixed X are not violated, as well.

#Primary Production Analysis --------
#Need to make time a numeric again
mixo_phys_pp_baseR_calc$Time_days <- as.numeric(as.character(mixo_phys_pp_baseR_calc$Time_days))
#base model with no grouping
mixo_mp_pp_gam <- gam(pg_perCELL_perHR ~ s(Time_days, k = 3), data = mixo_phys_pp_baseR_calc)
summary(mixo_mp_pp_gam)

mixo_mp_pp_gam_byPlastic <- gam(pg_perCELL_perHR ~ s(Time_days, by = Treatment_plastic, k = 2), data = mixo_phys_pp_baseR_calc)
summary(mixo_mp_pp_gam_byPlastic)

#Here, by setting k = 2, I am manually setting the maximum number of knots to be 2.
#Since the intercept take up 1 the model becomes an intercept model and thus the significant result
# for the control smoother.

#Now lets add additional turns
mixo_mp_pp_gam_byPlastic_k3 <- gam(pg_perCELL_perHR ~ s(Time_days, by = Treatment_plastic, k = 3), data = mixo_phys_pp_baseR_calc)
summary(mixo_mp_pp_gam_byPlastic_k3)

#Here, by setting k = 3, I am manually setting the maximum number of knots to be 2.
#Since the intercept take up 1 the model becomes a linear model and thus the significant result
# for the control smoother.

mixo_mp_pp_gam_byPlastic_k4 <- gam(pg_perCELL_perHR ~ s(Time_days, by = Treatment_plastic, k = 4), data = mixo_phys_pp_baseR_calc)
summary(mixo_mp_pp_gam_byPlastic_k4)

#Now by allowing for a more polynomial like structure the smoother is now significant
#For multiple plastic treatments. This is also, the maximum value of knots allowed
#in this model as I would then the term would have fewer covariate combinations
#than the maximum degrees of freedom.

#Assumption Testing
plot(mixo_phys_pp_baseR_calc$Time_days, residuals(mixo_mp_pp_gam_byPlastic_k4),
     xlab = "Time (Days)", ylab = "PP Residuals", main = "GAM Model residuals vs Time")

#gam.check is better
gam.check(mixo_mp_pp_gam_byPlastic_k4)

#plot.gam(mixo_mp_pp_gam_byPlastic_k4, residuals = TRUE)

#Assumptions of GAM analysis:
#     Homogeneity of variance
#     Normally distributed
#     Independence
#     Fixed X

#Just by looking at the residuals plot I can see that the homogeneity of variance
#could be violated. It seems as though there is a trend that I am not capturing in 
#this model. Now I have been trying but I have not been able to get the GAMM working
#so I am running with this for now.
#As part of gam.check(), I can also confirm that the K value of 4 used in this model
#is appropriate. Too high or too low of a k value can cause overfitting.


#Bacterivory Analysis --------
#I will be using the same model for the feeding results but with adjusted k values
# I need to convert into numeric for every model estimation
mixo_mp_feed$Time_days <- as.numeric(as.character(mixo_mp_feed$Time_days))
mixo_mp_feed_gam <- gam(ing_percell_perhour ~ s(Time_days, by = Treatment_plastic, k = 2), data = mixo_mp_feed)
#summary(mixo_mp_feed_gam)

mixo_mp_feed_gam3 <- gam(ing_percell_perhour ~ s(Time_days, by = Treatment_plastic, k = 3), data = mixo_mp_feed)
#summary(mixo_mp_feed_gam3)

#Assumption Testing
gam.check(mixo_mp_feed_gam3)

#Assumptions of GAM analysis:
#     Homogeneity of variance
#     Normally distributed
#     Independence
#     Fixed X

#The bacterivory dataset is very different from the PP dataset. It adheres to the
#Normality and homogeneity of variance assumptions. And the k value selected
# is confirmed to be appropriate for the model.

#Making of the GAM no longer works with this dataset for some reason I can't figure out.
# There is no model that it can fit and converge or there isn't enough data points for
# the model to estimate from

#Results Write Up ------------
#Abundance Assessment
#In both specific growth rate and doubling time, there were not statistically
#significant differences between group means as determined by the one-way ANOVA,
#(F(3,12) = 0.541, p = 0.663, and F(3,12) = 0.919, p = 0.461, respectively).
#Therefor no post hoc tests were supplied.

#NOT part of the write up but I want to get my thoughs down for later...
#It seems as though there is a linear trend BUT the variance is too high to lend
#to significant results, as suggested by the broken assumption of homogeneity of variance.
#Such is the life of doing science with live cultures.

#Primary Production
#Carbon fixation rates were significantly affected by contamination of 
#microplastics at levels of 10 per ml and 100000 per ml(p = 0.01, f = 4.3; p <0.04,
#f = 3.1). The response curve of the microplastic treatments exhibited a polynomial?
# shape, where as the control treatment exhibited a negative linear response (p < 0.001, F = 13.5).

#Bacterivory
#Bacterial sized particle ingestion rates were significantly affected by all
#microplastic exposure levels (p <0.001, F = 37.3; p < 0.001, F = 26.7; p = 0.006, F =  6.5).
#Similar to Primary Production, the microplastic contamination exhibited a 
#statistically significant non-response curve shape (p < 0.001) but only at the 
#higher levels of contamination (1000 per ml, 100000 per ml). The lowest treatment
#level saw a negative linear decrease, whereas the control responded in a 
#positive linear fashion.








#Vizualization of Results ----------
#Plot table using base R graphics
#plot(NULL)

#addtable2plot(0,0.1, mixo_mp_gc_table, bg = "grey",bty = "o", hlines = T, 
#              vlines = T, cex = 1)
#plot(predict(mixo_mp_pp_gam_byPlastic_k4), type = "l")

#ewdata_pp <- mixo_phys_pp_baseR_calc[,c("Treatment_plastic", "Time_days")]
#mixo_phys_pp_predict <- predict(mixo_mp_pp_gam_byPlastic_k4,newdata = newdata_pp)
#I can't get predict to work properly for me to be in the same axis as the original
#data.

#Primary Production
#plot(mixo_mp_pp_gam_byPlastic_k4)

pp_gam_predict <- predict_gam(mixo_mp_pp_gam_byPlastic_k4)

plot(pp_gam_predict)
#Bacterivory
#plot(mixo_mp_feed_gam3)

feed_game_predict <- predict_gam(mixo_mp_feed_gam3)

print(feed_game_predict$Time_days == 3)

#Interpreting Results -------------
#While nothing turned out as I expected, there are still significant results to report.
#However, microplastics seem to have no significant impact on the specific growth rate
#and doubling time, even though max population was much lower for the microplastic
#treatments. Thus I have to reject my hypothesis. But primary production and
#bacterivory tell a different story. Both physiological responses saw signficant, yet
#nonlinear responses to microplastic treatment over the course of a week. In both
# parts of the experiment, the control had a negative linear response, whereas
#microplastics treatments expressed a sinusoidal(polynomial?) like response.
#These difference may come down to the suppression of the max population reached
#in the plastic treatments allowing for less nutrient limitation and increased prey
#abundance over the week relative to the control. However, without futher investigation
#it would be impossible to do more than speculation. Irregardless of shape, in primary
#production increased relative to control after 7 days and bacterivory decreased
#relative to control, confirming my second hypothesis.



#Redoing Code to use in GGPLOT ---------

#GGPLOT is much aestheticically pleasing and more pub ready AND easier to use...

#Setting a new theme for all plots that fit better with publications
theme_set(theme_bw())
theme_update(element_line(linewidth = 0), panel.grid.major = element_blank())

ggplot(mixo_phys_pp_baseR_calc)+
    geom_boxplot(aes(x = as.factor(Time_days), y = pg_perCELL_perHR, fill = Treatment_plastic))+
    geom_line(data = pp_gam_predict, aes(x = Time_days, y = fit, color = Treatment_plastic))
#This almost works

#rework boxplot to have on continuous scale
ggplot(mixo_phys_pp_baseR_calc)+
    geom_boxplot(aes(x = Time_days, y = pg_perCELL_perHR, group = interaction(Time_days,Treatment_plastic), fill = Treatment_plastic), alpha = 0.8)+
    geom_line(data = pp_gam_predict, aes(x = Time_days, y = fit, color = Treatment_plastic))+
    geom_ribbon(data = pp_gam_predict, aes(x = Time_days,ymin = fit-se.fit, ymax = fit+se.fit, fill = Treatment_plastic), alpha = 0.2)

#want to get them both on the same scale - will try with geom points

ggplot(mixo_phys_pp_baseR_calc)+
    geom_point(aes(x = Time_days, y = pg_perCELL_perHR, color = Treatment_plastic))+
    geom_line(data = pp_gam_predict, aes(x = Time_days, y = fit, color = Treatment_plastic))+
    geom_ribbon(data = pp_gam_predict, aes(x = Time_days,ymin = fit-se.fit, ymax = fit+se.fit, fill = Treatment_plastic), alpha = 0.2)

#THIS WORKS GREAT

#Now to repeat this with bacterivory

ggplot(mixo_mp_feed)+
    geom_point(aes(x = Time_days, y = ing_percell_perhour, color = Treatment_plastic))+
    geom_line(data = feed_game_predict, aes(x = Time_days, y = fit, color = Treatment_plastic))+
    geom_ribbon(data = feed_game_predict, aes(x = Time_days,ymin = fit-se.fit, ymax = fit+se.fit, fill = Treatment_plastic), alpha = 0.2)

#Looking at bacteria concentrations for all conditions over time

ggplot(mixo_mp_feed, aes(x = Time_days, y = bac_conc_perML))+
    geom_point(aes(color = Treatment_plastic, size = 1.2))+
    theme(axis.text = element_text(size = 12),
          legend.text = element_text(size = 12))
#I didn't count bacteria for each Rep on for Rep A in each condition.
#May need to look at this.

#I want to keep all of the data and see what is looks like from the perspective
#of change from control because we had a drastic drop in population from day 5
# to day 6. Thus, any data after that is not reliable. HOWEVER, all treatments
# responded the same way. Thus I'm doing this to retain reliance of all the data.

#This creates a vector of the same length as the dataframe with group means for 
# each observation's corresponding group

mixo_mp_feed_grpmean <- ave(mixo_mp_feed$ing_percell_perhour, mixo_mp_feed$Treatment_plastic, mixo_mp_feed$Time_days, FUN = mean)

mixo_mp_feed$ing_percell_perhour_scaledbyCntlmean <- mixo_mp_feed$ing_percell_perhour - mixo_mp_feed_grpmean


ggplot(mixo_mp_feed)+
    geom_point(aes(x = Time_days, y = ing_percell_perhour_scaledbyCntlmean, color = Treatment_plastic))

ggplot(mixo_mp_feed)+
    geom_boxplot(aes(x = as.factor(Time_days), y = ing_percell_perhour_scaledbyCntlmean, fill = Treatment_plastic))

#Did NOT work...I've included all group means so I scaled all treatments to 0
# I need to only apply the control group means for each time point and apply
# to all data within that time point.

#I'll use aggregate to calculate each group mean and then create a vector of 
#only the control group means.

mixo_mp_group_means_df <- aggregate(mixo_mp_feed, by = list(mixo_mp_feed$Treatment_plastic, mixo_mp_feed$Time_days), FUN = "mean")

mixo_mp_feed_ctlgrpmean <- c(rep(mixo_mp_group_means_df$ing_percell_perhour[1], 16), rep(mixo_mp_group_means_df$ing_percell_perhour[5], 16),
      rep(mixo_mp_group_means_df$ing_percell_perhour[9], 16))

mixo_mp_feed$ing_percell_perhour_scaledbyCntlmean <- mixo_mp_feed$ing_percell_perhour - mixo_mp_feed_ctlgrpmean

ggplot(mixo_mp_feed)+
    geom_point(aes(x = Time_days, y = ing_percell_perhour_scaledbyCntlmean, color = Treatment_plastic))

ggplot(mixo_mp_feed)+
    geom_boxplot(aes(x = as.factor(Time_days), y = ing_percell_perhour_scaledbyCntlmean, fill = Treatment_plastic))

#This gives me the results that I need. Now to repeat for abundance and Primary Production Data


mixo_phys_pp_baseR_calc_grpmeansdf <- aggregate(mixo_phys_pp_baseR_calc, by = list(mixo_phys_pp_baseR_calc$Treatment_plastic, mixo_phys_pp_baseR_calc$Time_days), FUN = "mean")

mixo_phys_pp_baseR_calc_ctlgrpmeans <- c(rep(mixo_phys_pp_baseR_calc_grpmeansdf$pg_perCELL_perHR[1], 16), rep(mixo_phys_pp_baseR_calc_grpmeansdf$pg_perCELL_perHR[5], 16),
                                         rep(mixo_phys_pp_baseR_calc_grpmeansdf$pg_perCELL_perHR[9], 16),rep(mixo_phys_pp_baseR_calc_grpmeansdf$pg_perCELL_perHR[13], 16))

mixo_phys_pp_baseR_calc$pg_perCELL_perHR_scaledbyCntlmean <- mixo_phys_pp_baseR_calc$pg_perCELL_perHR - mixo_phys_pp_baseR_calc_ctlgrpmeans

ggplot(mixo_phys_pp_baseR_calc)+
    geom_boxplot(aes(x = as.factor(Time_days), y = pg_perCELL_perHR_scaledbyCntlmean, fill = Treatment_plastic))

#These are the new data for both PP and Feeding to use all 7 days

#Revisualize abundance data

#NOW I need to check the PP model after removing the outliers in the control in days 0 and 2.

mixo_mp_abund_grpmeans <- aggregate(mixo_mp_abund, by = list(mixo_mp_abund$Treatment, mixo_mp_abund$Time_days), FUN = "mean")

mixo_mp_abund_grpmeans_tidyverse <- group_by(mixo_mp_abund, Time_days, Treatment) %>% summarise(cell_per_L_means = mean(cell_per_L))

mixo_mp_abund_cntlgrpmeans <- c(rep(mixo_mp_abund_grpmeans$cell_per_L[1], 16),
                                rep(mixo_mp_abund_grpmeans$cell_per_L[5], 16),
                                rep(mixo_mp_abund_grpmeans$cell_per_L[9], 16),
                                rep(mixo_mp_abund_grpmeans$cell_per_L[13], 16),
                                rep(mixo_mp_abund_grpmeans$cell_per_L[17], 16),
                                rep(mixo_mp_abund_grpmeans$cell_per_L[21], 16),
                                rep(mixo_mp_abund_grpmeans$cell_per_L[25], 16),
                                rep(mixo_mp_abund_grpmeans$cell_per_L[29], 16))

mixo_mp_abund$cell_per_L_scaledbtCntlmean <- mixo_mp_abund$cell_per_L - mixo_mp_abund_cntlgrpmeans

ggplot(mixo_mp_abund)+
    geom_boxplot(aes(x = as.factor(Time_days), y = cell_per_L_scaledbtCntlmean, fill = Treatment))

#Need to check the overall PP and bacterivory for the culture by biovolume of the treatments
#as in check the PP of the entire culture and not on a per cell basis

mixo_mp_abund_feeding_days <- filter(mixo_mp_abund, Time_days == "0"| Time_days == "3"| Time_days =="7")

mixo_mp_abund_pp_days <- filter(mixo_mp_abund, Time_days == "0"| Time_days == "2" | Time_days == "4" | Time_days =="7")

mixo_mp_abund_feeding_days$overall_feeding <- mixo_mp_abund_feeding_days$cell_per_L * mixo_mp_feed$ing_percell_perhour

mixo_mp_abund_pp_days$overall_pp <- mixo_mp_abund_pp_days$cell_per_L * mixo_phys_pp_baseR_calc$pg_perCELL_perHR

ggplot(mixo_mp_abund_feeding_days, aes(x = Time_days, y = overall_feeding))+
    geom_boxplot(aes(fill = Treatment))

ggplot(mixo_mp_abund_pp_days, aes(x = Time_days, y = overall_pp))+
    geom_boxplot(aes(fill = Treatment))

#Final Considerations for project --------------

# need to consider which are the outliers in which measured response...
# then decide whether to remove them or not. ***Probably won't***
# Re-do growth curve calculations using a different time frame. 
    # Iterative derivative to find max growth for each treatment
# need to consider how to calculate clearance rates for feeding and compare from usphere technique
# look at the bacterial concentrations of the treatments and how they differ (may need to count more)
        # limitations to bacterial counting due to stain quick burnoff
# need to consider how the group means are calculated and if a different grp mean approch is better (next code chunk)
#DONE
# need to check if plastics had an impact on PP estimates from scint counter
# should make a chart for current research known for different plankton trophic groups
# LOOK through notes to check for any notions of cellular size differences.
    # could go back and test for biomass differences by using average cellular size
# Light and temperature regimes are OKAY no need to worry about those for now.


#Checking a second way of calculating the group means for a more "true mean" approch
#This involves adding the totals of each PP and estimated grazing. THEN divide
#by the total population and see if there's a difference in estimates and how that 
#impacts the model  
#  *** this may NOT be possible with feeding with how it's currently calculated ***

#Testing alt group mean calc -------------

# extract only needed columns
mixo_mp_pp_alt <- mixo_phys_pp_baseR_calc[,c(1:3,12,13)]
#calculate total sums of both PP and abundance
mixo_mp_pp_alt_grp <- aggregate(cbind(ugram_Cfixed_per_L_hr, cell_per_L) ~ Treatment_plastic*Time_days, data = mixo_mp_pp_alt, FUN = "sum")
#add in new per cell calculations as done before
mixo_mp_pp_alt_grp$ugram_Cfixed_perCELL_perHR <- mixo_mp_pp_alt_grp$ugram_Cfixed_per_L_hr/mixo_mp_pp_alt_grp$cell_per_L
mixo_mp_pp_alt_grp$pg_perCELL_perHR <- mixo_mp_pp_alt_grp$ugram_Cfixed_perCELL_perHR*1000000

cbind(mixo_phys_pp_baseR_calc_grpmeansdf$pg_perCELL_perHR,mixo_mp_pp_alt_grp$pg_perCELL_perHR)

#NOW this does change some things but not much for most treatments but does for some controls

#lets see how different the PP curves look with the new scaling

mixo_mp_pp_alt_grpmeans <- c(rep(mixo_mp_pp_alt_grp$pg_perCELL_perHR[1], 16), rep(mixo_mp_pp_alt_grp$pg_perCELL_perHR[5], 16),
                             rep(mixo_mp_pp_alt_grp$pg_perCELL_perHR[9], 16),rep(mixo_mp_pp_alt_grp$pg_perCELL_perHR[13], 16))

mixo_phys_pp_baseR_calc$pg_perCELL_perHR_scaledbyCntlmean <- mixo_phys_pp_baseR_calc$pg_perCELL_perHR - mixo_mp_pp_alt_grpmeans

ggplot(mixo_phys_pp_baseR_calc)+
    geom_boxplot(aes(x = as.factor(Time_days), y = pg_perCELL_perHR_scaledbyCntlmean, fill = Treatment_plastic))

#Not too much visual difference except at time 0. The plotted control means seem to be closer to 0.

# Next I will test this with abundance counts

mixo_mp_abund_alt <- mixo_mp_abund[,c(1:3, 6)]

mixo_mp_abund_alt_grpmeans <- aggregate(cell_per_L~Treatment*Time_days, mixo_mp_abund_alt, FUN = "sum")

#Checking to see if there is a difference "4" is the mean since there are 4 reps and a sum was calculated
cbind(mixo_mp_abund_grpmeans$cell_per_L, mixo_mp_abund_alt_grpmeans$cell_per_L/4)

# NO difference thus is stops here

#Lastly, check to see if there is a way to change it for feeding...

# Theoretically, the calculation is dependent only on the ratio counted of uspheres counted
# therefore, only the rate is averaged and not normalize by the raw population values as for PP.
#UNLESS I take a total microsphere to bacteria ratio and 

#Outlier checking ----------------
#Primary Production & abundance, respectively
ggplot(mixo_phys_pp_baseR_calc)+
    geom_boxplot(aes(x = as.factor(Time_days), y = pg_perCELL_perHR, fill = Treatment_plastic))

ggplot(filter(mixo_phys_pp_baseR_calc, Treatment_plastic == "control"))+
    geom_point(aes(x = Time_days, y = pg_perCELL_perHR, col = Rep))+
    geom_text(aes(x = Time_days, y = pg_perCELL_perHR,label = Rep), nudge_x = 0.25, size = 10)

ggplot(filter(mixo_mp_abund, Treatment == "Control"))+
    geom_point(aes(x = Time_days, y = cell_per_L, col = Rep))+
    geom_text(aes(x = Time_days, y = cell_per_L,label = Rep), nudge_x = 0.25, size = 10)

ggplot(filter(mixo_phys_pp_baseR_calc, Treatment_plastic == "10"))+
    geom_point(aes(x = Time_days, y = pg_perCELL_perHR, col = Treatment_plastic))+
    geom_text(aes(x = Time_days, y = pg_perCELL_perHR,label = Rep, check_overlap = TRUE))

ggplot(filter(mixo_mp_abund, Treatment == "10"))+
    geom_point(aes(x = Time_days, y = cell_per_L, col = Treatment))+
    geom_text(aes(x = Time_days, y = cell_per_L,label = Rep))

ggplot(filter(mixo_phys_pp_baseR_calc, Treatment_plastic == "1000"))+
    geom_point(aes(x = Time_days, y = pg_perCELL_perHR, col = Treatment_plastic))+
    geom_text(aes(x = Time_days, y = pg_perCELL_perHR,label = Rep, check_overlap = TRUE))

ggplot(filter(mixo_mp_abund, Treatment == "1000"))+
    geom_point(aes(x = Time_days, y = cell_per_L, col = Treatment))+
    geom_text(aes(x = Time_days, y = cell_per_L,label = Rep))

ggplot(filter(mixo_phys_pp_baseR_calc, Treatment_plastic == "100000"))+
    geom_point(aes(x = Time_days, y = pg_perCELL_perHR, col = Treatment_plastic))+
    geom_text(aes(x = Time_days, y = pg_perCELL_perHR,label = Rep, check_overlap = TRUE))

ggplot(filter(mixo_mp_abund, Treatment == "100000"))+
    geom_point(aes(x = Time_days, y = cell_per_L, col = Treatment))+
    geom_text(aes(x = Time_days, y = cell_per_L,label = Rep))

#Data visualization without Rep B of control -----------
#PP
#Helps remove a lot of variation in the control, slightly changes interpretation
raw_PP_directExp <- ggplot(filter(mixo_phys_pp_baseR_calc, !(Rep == "B" & Treatment_plastic == "control")))+
    geom_boxplot(aes(x = as.factor(Time_days), y = pg_perCELL_perHR, fill = Treatment_plastic))+
  xlab("Time (days)")+ylab(expression("Primary Production (pg C"~"cell"^-1~"hr"^-1~")"))+labs(fill = "Treatment")
#Bacterivory
#No change to variation or interpretation
raw_ing_directExp <- ggplot(filter(mixo_mp_feed, !(treatment_rep == "B" & Treatment_plastic == "control")))+
    geom_boxplot(aes(x = as.factor(Time_days), y = ing_percell_perhour, fill = Treatment_plastic))+
  xlab("Time (days)")+ylab(expression("Bacterial Sized Particle Ing (BSP"~"cell"^-1~"hr"^-1~")"))+labs(fill = "Treatment")
#Growth Curve
#Still consistently higher but not as signifigantly. The peak abundance is signifigant though
ggplot(filter(mixo_mp_abund, !(Rep == "B" & Treatment == "Control")))+
    geom_boxplot(aes(x = as.factor(Time_days), y = cell_per_L, fill = Treatment))

#Need to remove control Rep B from all datasets before redo of the mean scaling
mixo_mp_abund_noB <- filter(mixo_mp_abund, !(Rep == "B" & Treatment == "Control"))

mixo_phys_pp_baseR_calc_noB <- filter(mixo_phys_pp_baseR_calc, !(Rep == "B" & Treatment_plastic == "control"))

mixo_mp_feed_noB <- filter(mixo_mp_feed, !(treatment_rep == "B" & Treatment_plastic == "control"))

# ggsaving the raw data plots
ggsave("raw_ing_directExp.pdf",raw_ing_directExp)
ggsave("raw_PP_directExp.pdf",raw_PP_directExp)
#Now to do the mean scaling

#PP scaling

mixo_phys_pp_baseR_calc_grpmeansdf_noB <-  aggregate(pg_perCELL_perHR ~ Treatment_plastic*Time_days, data = mixo_phys_pp_baseR_calc_noB, FUN = "mean")

mixo_phys_pp_baseR_calc_ctlgrpmeans_noB <- c(rep(mixo_phys_pp_baseR_calc_grpmeansdf_noB$pg_perCELL_perHR[1], 15), rep(mixo_phys_pp_baseR_calc_grpmeansdf_noB$pg_perCELL_perHR[5], 15),
                             rep(mixo_phys_pp_baseR_calc_grpmeansdf_noB$pg_perCELL_perHR[9], 15),rep(mixo_phys_pp_baseR_calc_grpmeansdf_noB$pg_perCELL_perHR[13], 15))

mixo_phys_pp_baseR_calc_noB$pg_perCELL_perHR_scaledbyCntlmean <- mixo_phys_pp_baseR_calc_noB$pg_perCELL_perHR - mixo_phys_pp_baseR_calc_ctlgrpmeans_noB

primary_prod_scaled_direct <- ggplot(mixo_phys_pp_baseR_calc_noB)+
    geom_boxplot(aes(x = as.factor(Time_days), y = pg_perCELL_perHR_scaledbyCntlmean, fill = Treatment_plastic))+
  xlab("Time (days)")+ylab(expression("Primary Production (pg C"~"cell"^-1~"hr"^-1~")"))+labs(fill = "Treatment")

ggsave("primary_prod_scaled_direct.pdf", primary_prod_scaled_direct)

#Abundance scaling .... this is where I expect to see the most difference

mixo_mp_abund_grpmeans_noB <- aggregate(cell_per_L ~ Treatment*Time_days, data = mixo_mp_abund_noB, FUN = "mean")

mixo_mp_abund_cntlgrpmeans_noB <- c(rep(mixo_mp_abund_grpmeans_noB$cell_per_L[1], 15),
                                rep(mixo_mp_abund_grpmeans_noB$cell_per_L[5], 15),
                                rep(mixo_mp_abund_grpmeans_noB$cell_per_L[9], 15),
                                rep(mixo_mp_abund_grpmeans_noB$cell_per_L[13], 15),
                                rep(mixo_mp_abund_grpmeans_noB$cell_per_L[17], 15),
                                rep(mixo_mp_abund_grpmeans_noB$cell_per_L[21], 15),
                                rep(mixo_mp_abund_grpmeans_noB$cell_per_L[25], 15),
                                rep(mixo_mp_abund_grpmeans_noB$cell_per_L[29], 15))

mixo_mp_abund_noB$cell_per_L_scaledbtCntlmean <- mixo_mp_abund_noB$cell_per_L - mixo_mp_abund_cntlgrpmeans_noB

abund_scaled_plot_direct <- ggplot(mixo_mp_abund_noB)+
    geom_boxplot(aes(x = as.factor(Time_days), y = cell_per_L_scaledbtCntlmean, fill = Treatment))+
  xlab("Time (days)")+ ylab(expression("Abundance"~"(cell L"^-1~")"))+labs(fill = "Treatment")


ggsave("abund_scaled_plot_direct.pdf", abund_scaled_plot_direct)
#Bacterivory scaling

mixo_mp_group_means_df_noB <- aggregate(ing_percell_perhour ~ Treatment_plastic*Time_days, data = mixo_mp_feed_noB, FUN = "mean")

mixo_mp_feed_ctlgrpmean_noB <- c(rep(mixo_mp_group_means_df_noB$ing_percell_perhour[1], 15), rep(mixo_mp_group_means_df_noB$ing_percell_perhour[5], 15),
                             rep(mixo_mp_group_means_df_noB$ing_percell_perhour[9], 15))

mixo_mp_feed_noB$ing_percell_perhour_scaledbyCntlmean <- mixo_mp_feed_noB$ing_percell_perhour - mixo_mp_feed_ctlgrpmean_noB

ggplot(mixo_mp_feed_noB)+
    geom_boxplot(aes(x = as.factor(Time_days), y = ing_percell_perhour_scaledbyCntlmean, fill = Treatment_plastic))

#Overall PP and Bacterivory without B

#Need to check the overall PP and bacterivory for the culture by biovolume of the treatments
#as in check the PP of the entire culture and not on a per cell basis

mixo_mp_abund_feeding_days_noB <- filter(mixo_mp_abund_noB, Time_days == "0"| Time_days == "3"| Time_days =="7")

mixo_mp_abund_pp_days_noB <- filter(mixo_mp_abund_noB, Time_days == "0"| Time_days == "2" | Time_days == "4" | Time_days =="7")

mixo_mp_abund_feeding_days_noB$overall_feeding <- mixo_mp_abund_feeding_days_noB$cell_per_L * mixo_mp_feed_noB$ing_percell_perhour

mixo_mp_abund_pp_days_noB$overall_pp <- mixo_mp_abund_pp_days_noB$cell_per_L * mixo_phys_pp_baseR_calc_noB$pg_perCELL_perHR

ggplot(mixo_mp_abund_feeding_days_noB, aes(x = Time_days, y = overall_feeding))+
    geom_boxplot(aes(fill = Treatment))

ggplot(mixo_mp_abund_pp_days_noB, aes(x = Time_days, y = overall_pp))+
    geom_boxplot(aes(fill = Treatment))




# General Results notes and some discussion points ---------------
#Abundance trends
#   Control began to separate and was consistently higher after days 3-4
#   only 1 day to be considered signifigantly higher and that's at day 5 - peak abundace
#   ALL microplastic treatments were reduced in reference to the control
#   and in descending order as we would expect
#this suggests an impact on cell health and a decrease in growth rate or an increase in mortality rate
#No way to really identify which is happening or if both are happening. In the end,
# I believe the growth rate estimated parameters will not be signifigantly different from one another.
#NO overt claims can be made here.

#PP results
#Over week trends for per cell
#Control was stable until day 4...begins a reduction in per cell PP
#All microplastic treatments had a nonlinear relationship over the course of the week
# Day 0: starts about the same or a little lower
# Day 2: microplastic treatments are higher than control
# Day 4: microplastic treatment then drops lower
# Day 7: microplastic treatment much higher
# overall weekly trends different shaped response between control and plastic treatments
# Last Day see biggest difference between treatments and control.

#Overall culture PP
# ***Have not scaled by control means yet...BUT***
# General trend Control Peaks at day 2-4 dips at start and end of week. (potential carrying capacity limit)
# Plastics peak at day 2 and all other days seem a little lower
# Day 4 is the only signifigant difference from control and was lower
# most days the plastic was not significantly different from contol

# Bacterivory Results
#Per cell Trends
# Control has net positive over the course of the week
#10^1 has a negative trend over the weekly period while the others exhibit a more
#non-linear relationship.
#low then a peak at day 3 and a dip at day 7. Signoficantly different from control on day 7

#Culture Trends
# Control sees increase across the whole week
# all 3 microplastic treatments see a non-linear trend similar to the per cell
# Signifigantly lower on day 7

#Removing Rep B helped resolve most of the variance from the dataset and fixed some 
#of the potential interpretability issues.

#Overall thoughts for discussion --------
#There's seems to be an issue during a 24 hour period where all treatments and control
# dropped. This makes it harder to see a true impact, however, since all of them
# reacted in the same way by scaling the impact to control we can get a comparison
# for the entire week. May need to mention this as a caveat.

#In regards to abundance, there is not much significance except at the day of the peak abundance, Day 5
#Treatments are lower suggesting some suppression of peak growth on the microplastic treatment.

#What could be causing the difference at the end of the experiment. Well, on a per cell basis,
# while control saw a linear negative correlation over the course of the week, due to a boost in abundance
# larger than a boost in overall PP, the microplastic treatments response differed over time. It seems
# over a 1 week period there is a difference adjustment to microplastic stress. It seems that initially
# there may be a an adjustment for both PP and bacterivory but over the course of 7 days the effects from
# microplastics start to settle in. As seen if we compare total carbon intake that the control has the most
# carbon incorporated which is indicative of a higher C biomass out of all the cultures.

#To do list -------------
# 1) I need to double check the doubling time calculation is correct
# 2) Scale the PP and feeding to the control group means as it has been applied to abundance
# 3) make figures for analysis and export in proper theme formating theme_bw() as pdfs
# 4) Perform statistical analysis GAMs or repeated measures regression
# 5) make model summary tables and report everything in a results section

# Rerun Analysis ------------

RR_df <- read.csv("Mixo_MP_feeding_Rerun.csv")

RR_df_calc <- mutate(RR_df,
                     pp_DPM = pp_CPM/C_14_efficiency,
                     bg_DPM = bg_CPM/C_14_efficiency/0.1,
                     dark_pp_DPM = dark_CPM/C_14_efficiency,
                     dark_bg_DPM = dark_bg_CPM/C_14_efficiency/0.1, 
                     ugram_Cfixed_per_L_hr = pp_DPM*DIC..umol.kg./bg_DPM*12*1.029/2,
                     dark_ugram_Cfixed_per_L_hr = dark_pp_DPM*DIC..umol.kg./dark_bg_DPM*12*1.029/2,
                     cell_per_mL = cell_avg/1.35*cor_factor_grid,
                     cell_per_L = cell_avg/1.35*cor_factor_grid*1000,
                     ugram_Cfixed_per_L_hr_cell = ugram_Cfixed_per_L_hr/(cell_per_L/2),
                     pg_perCELL_perHR = ugram_Cfixed_per_L_hr_cell*1000000,
                     bac_conc_perML = bac_avg/bac_micro_filtered*dil_factor*cor_factor_grid,
                     micro_conc_perML = micro_avg/bac_micro_filtered*dil_factor*cor_factor_grid,
                     BM_ratio = bac_conc_perML/micro_conc_perML,
                     BSP_ingested = conc_bkgrd_per*BM_ratio+conc_bkgrd_per,
                     ing_percell_perhour = BSP_ingested*2,
                     overall_feeding = ing_percell_perhour*cell_per_L,
                     overall_pp = pg_perCELL_perHR*cell_per_L,
                     ugram_Cfeed_percell_perhour_L = ing_percell_perhour*(2.2*10^-14)*1000000,
                     pg_Cfeed_percell_perhour_L = ugram_Cfeed_percell_perhour_L*1000000,
                     ochromonas_biomass_pg_carbon_percell = (10^-1.694)*(268.0826^1.218),
                     #These numbers are for radii
                     #I need 4um radius for 8 um protist
                     # 4 um = 4/3*pi*4^3 = 268.0826
                     feed_carbon_percell = pg_Cfeed_percell_perhour_L/ochromonas_biomass_pg_carbon_percell,
                     PP_carbon_percell = pg_perCELL_perHR/ochromonas_biomass_pg_carbon_percell,
                     total_carbon_incorporated_pg_cell = feed_carbon_percell+PP_carbon_percell, 
                     ochromonas_biomass_total_pg = ochromonas_biomass_pg_carbon_percell*cell_per_L
                     )


RR_df_CF_test <- RR_df
RR_df_CF_test$cor_factor_grid <- 10000

mixo_phys_pp_baseR_calc$pp_DPM <- mixo_phys_pp_baseR_calc$pp_CPM/mixo_phys_pp_baseR_calc$C_14_efficiency
mixo_phys_pp_baseR_calc$bkgr_DPM <- mixo_phys_pp_baseR_calc$background_CPM/mixo_phys_pp_baseR_calc$C_14_efficiency/0.1
mixo_phys_pp_baseR_calc$ugram_Cfixed_per_L_hr <- mixo_phys_pp_baseR_calc$pp_DPM*mixo_phys_pp_baseR_calc$DIC.umol.kg./mixo_phys_pp_baseR_calc$bkgr_DPM*12*1.029/2
mixo_phys_pp_baseR_calc$cell_per_L <- mixo_phys_pp_baseR_calc$Cell_counts/1.35*mixo_phys_pp_baseR_calc$Zeiss_CF*1000*1.08
mixo_phys_pp_baseR_calc$ugram_Cfixed_perCELL_perHR <- mixo_phys_pp_baseR_calc$ugram_Cfixed_per_L_hr/mixo_phys_pp_baseR_calc$cell_per_L
mixo_phys_pp_baseR_calc$pg_perCELL_perHR <- mixo_phys_pp_baseR_calc$ugram_Cfixed_perCELL_perHR*1000000



RR_df_CF_test_calc <- mutate(RR_df_CF_test,
                     pp_DPM = pp_CPM/C_14_efficiency,
                     bg_DPM = bg_CPM/C_14_efficiency/0.1,
                     dark_pp_DPM = dark_CPM/C_14_efficiency,
                     dark_bg_DPM = dark_bg_CPM/C_14_efficiency/0.1, 
                     ugram_Cfixed_per_L_hr = pp_DPM*DIC..umol.kg./bg_DPM*12*1.029/2,
                     dark_ugram_Cfixed_per_L_hr = dark_pp_DPM*DIC..umol.kg./dark_bg_DPM*12*1.029/2,
                     cell_per_mL = cell_avg/1*cor_factor_grid,
                     cell_per_L = cell_avg/1*cor_factor_grid*1000,
                     ugram_Cfixed_per_L_hr_cell = ugram_Cfixed_per_L_hr/cell_per_L,
                     pg_perCELL_perHR = ugram_Cfixed_per_L_hr_cell*1000000,
                     bac_conc_perML = bac_avg/bac_micro_filtered*dil_factor*cor_factor_grid,
                     micro_conc_perML = micro_avg/bac_micro_filtered*dil_factor*cor_factor_grid,
                     BM_ratio = bac_conc_perML/micro_conc_perML,
                     BSP_ingested = conc_bkgrd_per*BM_ratio+conc_bkgrd_per,
                     ing_percell_perhour = BSP_ingested*2,
                     overall_feeding = ing_percell_perhour*cell_per_L,
                     overall_pp = pg_perCELL_perHR*cell_per_L,
                     ugram_Cfeed_percell_perhour = ing_percell_perhour*(2.2*10^-14)*1000000,
                     pg_Cfeed_percell_perhour = ugram_Cfeed_percell_perhour*1000000,
                     ochromonas_biomass_pg_carbon_percell = (10^-1.694)*(904.7787^1.218),
                     #These numbers are for radii
                     #I need 4um radius for 8 um protist
                     # 4 um = 4/3*pi*4^3 = 268.0826
                     #for 5 um = 523.5988
                     #for 6 um = 904.7787 PP fails to account for enough carbon here
                     feed_carbon_percell = pg_Cfeed_percell_perhour/ochromonas_biomass_pg_carbon_percell*100,
                     PP_carbon_percell = pg_perCELL_perHR/ochromonas_biomass_pg_carbon_percell*100,
                     total_carbon_incorporated_pg_cell = feed_carbon_percell+PP_carbon_percell, 
                     ochromonas_biomass_total_pg = ochromonas_biomass_pg_carbon_percell*cell_per_L
)

mixo_phys_pp_baseR_calc$pp_DPM/mixo_phys_pp_baseR_calc$bkgr_DPM
RR_df_CF_test_calc$pp_DPM/RR_df_CF_test_calc$bg_DPM

mixo_phys_pp_baseR_calc$pp_DPM
RR_df_CF_test_calc$pp_DPM
mixo_phys_pp_baseR_calc$bkgr_DPM
RR_df_CF_test_calc$bg_DPM

RR_df_CF_test_calc$pg_Cfeed_percell_perhour

#I need to fix the categorical vars to be factors


RR_df_calc$Time_days <- factor(RR_df_calc$Time_days)

RR_df_calc$treatment_rep <- factor(RR_df_calc$treatment_rep)

RR_df_calc$Treatment_plastic <- factor(RR_df_calc$Treatment_plastic,
                                         levels = c("control", "100000"))

RR_df_calc$Dialysis_or_Not <- factor(RR_df_calc$Dialysis_or_Not,
                                     levels = c("No", "Dialysis Bag"))

RR_df_calc <- drop_na(RR_df_calc)

# Same for test dataset

RR_df_CF_test_calc$Time_days <- factor(RR_df_CF_test_calc$Time_days)

RR_df_CF_test_calc$treatment_rep <- factor(RR_df_CF_test_calc$treatment_rep)

RR_df_CF_test_calc$Treatment_plastic <- factor(RR_df_CF_test_calc$Treatment_plastic,
                                       levels = c("control","100000"))

RR_df_CF_test_calc$Dialysis_or_Not <- factor(RR_df_CF_test_calc$Dialysis_or_Not,
                                     levels = c("No", "Dialysis Bag"))

RR_df_CF_test_calc <- drop_na(RR_df_CF_test_calc)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = Time_days, y = cell_per_L, fill = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = Time_days, y = cell_per_mL/2, fill = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = Time_days, y = pg_perCELL_perHR, fill = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = Time_days, y = ing_percell_perhour, fill = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_point(aes(x = Time_days, y = bac_conc_perML, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)


# I need to compare the numbers between the two runs.
# I will first look at the raw values and compare them, but I think it would be best
# to compare on a per cell basis

#Time was converted by into numerics for model estimation and will be converted in line
# here for boxplots as I do not want to ruin it for point and lines plots to follow

ggplot(mixo_mp_feed)+
  geom_boxplot(aes(x = as.factor(Time_days), y = ing_percell_perhour, fill = Treatment_plastic), position = position_dodge())+
  geom_boxplot(data = RR_df_calc, aes(x = Time_days, y = ing_percell_perhour, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)

ggplot(mixo_mp_feed)+
  geom_boxplot(aes(x = as.factor(Time_days), y = bac_conc_perML, fill = Treatment_plastic), position = position_dodge())+
  geom_boxplot(data = RR_df_calc, aes(x = Time_days, y = bac_conc_perML, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)

ggplot(mixo_mp_abund)+
  geom_boxplot(aes(x = Time_days, y = cell_per_L, fill = Treatment), position = position_dodge())+
  geom_boxplot(data = mixo_mp_feed, aes(x = as.factor(Time_days), y = cell_L, color = Treatment_plastic), position = position_dodge())

ggplot(mixo_mp_abund)+
  geom_boxplot(aes(x = Time_days, y = cell_per_L, fill = Treatment), position = position_dodge())+
  geom_boxplot(data = RR_df_calc, aes(x = Time_days, y = cell_per_L/2, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)

ggplot(mixo_mp_feed)+
  geom_boxplot(aes(x = as.factor(Time_days), y = cell_L, fill = Treatment_plastic), position = position_dodge())+
  geom_boxplot(data = filter(RR_df_calc, treatment_rep != "A"), aes(x = Time_days, y = cell_per_L, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)

ggplot(mixo_mp_feed)+
  geom_boxplot(aes(x = as.factor(Time_days), y = cell_L, fill = Treatment_plastic), position = position_dodge())+
  geom_boxplot(data = RR_df_calc, aes(x = Time_days, y = cell_per_L, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)

ggplot(mixo_phys_pp_baseR_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = ugram_Cfixed_per_L_hr, fill = Treatment_plastic), position = position_dodge())+
  geom_boxplot(data = RR_df_calc, aes(x = Time_days, y = ugram_Cfixed_per_L_hr, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)

ggplot(mixo_phys_pp_baseR_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = pg_perCELL_perHR, fill = Treatment_plastic), position = position_dodge())+
  geom_boxplot(data = RR_df_calc, aes(x = Time_days, y = pg_perCELL_perHR*2, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)


ggplot(mixo_mp_abund)+
  geom_point(aes(x = Time_days, y = cell_per_L, color = Treatment))

# I need to make clearer figures for BOB to see but I may have already made them earlier
# in the script I will search for them first and then recalculate 

# Filtering only the last 2 days for direct comparison between benchtop and incubator.
mixo_phys_pp_baseR_calc_forCOMP <- filter(mixo_phys_pp_baseR_calc, Time_days == 0 | Time_days == 7)
mixo_mp_feed_forCOMP <- filter(mixo_mp_feed, Time_days == 0 | Time_days == 7)
mixo_mp_abund_forCOMP <- filter(mixo_mp_abund, Time_days == 0 | Time_days == 7)

RR_PP_comp <- ggplot(mixo_phys_pp_baseR_calc_forCOMP)+
  geom_boxplot(aes(x = as.factor(Time_days), y = pg_perCELL_perHR, fill = Treatment_plastic), position = position_dodge())+
  geom_boxplot(data = RR_df_calc, aes(x = Time_days, y = pg_perCELL_perHR*2, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)+scale_color_manual("Incubator", values = scales::hue_pal(h = c(50,250))(2))+
  scale_fill_manual("Benchtop", values = scales::hue_pal()(4))

RR_feed_comp <-ggplot(mixo_mp_feed_forCOMP)+
  geom_boxplot(aes(x = as.factor(Time_days), y = ing_percell_perhour, fill = Treatment_plastic), position = position_dodge())+
  geom_boxplot(data = RR_df_calc, aes(x = Time_days, y = ing_percell_perhour, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)+scale_color_manual("Incubator", values = scales::hue_pal(h = c(50,250))(2))+
  scale_fill_manual("Benchtop", values = scales::hue_pal()(4))

RR_bac_comp <-ggplot(mixo_mp_feed_forCOMP)+
  geom_boxplot(aes(x = as.factor(Time_days), y = bac_conc_perML, fill = Treatment_plastic), position = position_dodge())+
  geom_boxplot(data = RR_df_calc, aes(x = Time_days, y = bac_conc_perML, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)+scale_color_manual("Incubator", values = scales::hue_pal(h = c(50,250))(2))+
  scale_fill_manual("Benchtop", values = scales::hue_pal()(4))

RR_abund_Olympus_comp <-ggplot(mixo_mp_feed_forCOMP)+
  geom_boxplot(aes(x = as.factor(Time_days), y = cell_L, fill = Treatment_plastic), position = position_dodge())+
  geom_boxplot(data = RR_df_calc, aes(x = Time_days, y = cell_per_L, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)+scale_color_manual("Incubator", values = scales::hue_pal(h = c(50,250))(2))+
  scale_fill_manual("Benchtop", values = scales::hue_pal()(4))

RR_abund_comp <-ggplot(mixo_mp_abund_forCOMP)+
  geom_boxplot(aes(x = as.factor(Time_days), y = cell_per_L, fill = Treatment), position = position_dodge())+
  geom_boxplot(data = RR_df_calc, aes(x = Time_days, y = cell_per_L, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)+scale_color_manual("Incubator", values = scales::hue_pal(h = c(50,250))(2))+
  scale_fill_manual("Benchtop", values = scales::hue_pal()(4))

ggsave("RR_PP_comp.pdf", RR_PP_comp)

ggsave("RR_feed_comp.pdf", RR_feed_comp)

ggsave("RR_bac_comp.pdf", RR_bac_comp)

ggsave("RR_abund_Olympus_comp.pdf", RR_abund_Olympus_comp)

ggsave("RR_abund_comp.pdf", RR_abund_comp)

# Raw value plots for indirect experiment
primary_production_raw_indirect <- ggplot(RR_df_CF_test_calc)+
  geom_boxplot(aes(x = Time_days, y = pg_perCELL_perHR, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)+
  xlab("Time (days)")+ ylab(expression("Primary Production"~"(pg C cell"^-1~")"))+labs(fill = "Treatment")

bacterivory_raw_indirect <- ggplot(RR_df_CF_test_calc)+
  geom_boxplot(aes(x = Time_days, y = pg_Cfeed_percell_perhour, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)+
  xlab("Time (days)")+ ylab(expression("Bacterivory"~"(pg C cell"^-1~")"))+labs(fill = "Treatment")

ggsave("primary_production_raw_indirect.pdf", primary_production_raw_indirect)
ggsave("bacterivory_raw_indirect.pdf", bacterivory_raw_indirect)
# Taking the mean controlled data from the earlier graphs

RR_df_calc_grp <- group_by(RR_df_CF_test_calc, Treatment_plastic, Dialysis_or_Not, Time_days)

RR_df_calc_grpsum <- summarise(RR_df_calc_grp,
                               ing_percell_perhour_MEAN = mean(ing_percell_perhour),
                               ing_percell_perhour_SE = sd(ing_percell_perhour)/sqrt(4),
                               pg_Cfeed_percell_perhour_MEAN = mean(pg_Cfeed_percell_perhour),
                               pg_Cfeed_percell_perhour_SE = sd(pg_Cfeed_percell_perhour)/sqrt(4),
                               bac_conc_perML_MEAN = mean(bac_conc_perML),
                               bac_conc_perML_SE = sd(bac_conc_perML)/sqrt(4),
                               pg_perCELL_perHR_MEAN = mean(pg_perCELL_perHR),
                               pg_perCELL_perHR_SE = sd(pg_perCELL_perHR)/sqrt(4),
                               cell_per_L_MEAN = mean(cell_per_L),
                               cell_per_L_SE = sd(cell_per_L)/sqrt(4),
                               ugram_Cfixed_per_L_hr_MEAN = mean(ugram_Cfixed_per_L_hr),
                               ugram_Cfixed_per_L_hr_SE = sd(ugram_Cfixed_per_L_hr)/sqrt(4)
                               )

ggplot(RR_df_calc_grpsum)+
  geom_point(aes(x = Time_days, y = ing_percell_perhour_MEAN, color = Treatment_plastic, shape = Dialysis_or_Not), size = 2.5)+
  geom_errorbar(aes(x = Time_days,
                    ymin = ing_percell_perhour_MEAN - ing_percell_perhour_SE,
                    ymax = ing_percell_perhour_MEAN + ing_percell_perhour_SE, color = Treatment_plastic))

ggplot(RR_df_calc_grpsum)+
  geom_point(aes(x = Time_days, y = bac_conc_perML_MEAN, color = Treatment_plastic, shape = Dialysis_or_Not), size = 2.5)+
  geom_errorbar(aes(x = Time_days,
                    ymin = bac_conc_perML_MEAN - bac_conc_perML_SE,
                    ymax = bac_conc_perML_MEAN + bac_conc_perML_SE, color = Treatment_plastic))

ggplot(RR_df_calc_grpsum)+
  geom_point(aes(x = Time_days, y = pg_perCELL_perHR_MEAN, color = Treatment_plastic, shape = Dialysis_or_Not), size = 2.5)+
  geom_errorbar(aes(x = Time_days,
                    ymin = pg_perCELL_perHR_MEAN - pg_perCELL_perHR_SE,
                    ymax = pg_perCELL_perHR_MEAN + pg_perCELL_perHR_SE, color = Treatment_plastic))

ggplot(RR_df_calc_grpsum)+
  geom_point(aes(x = Time_days, y = ugram_Cfixed_per_L_hr_MEAN, color = Treatment_plastic, shape = Dialysis_or_Not), size = 2.5)+
  geom_errorbar(aes(x = Time_days,
                    ymin = ugram_Cfixed_per_L_hr_MEAN - ugram_Cfixed_per_L_hr_SE,
                    ymax = ugram_Cfixed_per_L_hr_MEAN + ugram_Cfixed_per_L_hr_SE, color = Treatment_plastic))

ggplot(RR_df_calc_grpsum)+
  geom_point(aes(x = Time_days, y = cell_per_L_MEAN, color = Treatment_plastic, shape = Dialysis_or_Not), size = 2.5)+
  geom_errorbar(aes(x = Time_days,
                    ymin = cell_per_L_MEAN - cell_per_L_SE,
                    ymax = cell_per_L_MEAN + cell_per_L_SE, color = Treatment_plastic))


RR_means_for_scaling_ing <- c(rep(RR_df_calc_grpsum[[1, 'ing_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[3, 'ing_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[1, 'ing_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[3, 'ing_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[2, 'ing_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[4, 'ing_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[2, 'ing_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[4, 'ing_percell_perhour_MEAN']], 4))

RR_means_for_scaling_ing_pg <- c(rep(RR_df_calc_grpsum[[1, 'pg_Cfeed_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[3, 'pg_Cfeed_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[1, 'pg_Cfeed_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[3, 'pg_Cfeed_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[2, 'pg_Cfeed_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[4, 'pg_Cfeed_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[2, 'pg_Cfeed_percell_perhour_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[4, 'pg_Cfeed_percell_perhour_MEAN']], 4))

RR_means_for_scaling_PP <- c(rep(RR_df_calc_grpsum[[1, 'pg_perCELL_perHR_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[3, 'pg_perCELL_perHR_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[1, 'pg_perCELL_perHR_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[3, 'pg_perCELL_perHR_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[2, 'pg_perCELL_perHR_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[4, 'pg_perCELL_perHR_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[2, 'pg_perCELL_perHR_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[4, 'pg_perCELL_perHR_MEAN']], 4))

RR_means_for_scaling_abund <- c(rep(RR_df_calc_grpsum[[1, 'cell_per_L_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[3, 'cell_per_L_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[1, 'cell_per_L_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[3, 'cell_per_L_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[2, 'cell_per_L_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[4, 'cell_per_L_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[2, 'cell_per_L_MEAN']], 4),
                              rep(RR_df_calc_grpsum[[4, 'cell_per_L_MEAN']], 4))


ggplot(RR_df_calc)+
  geom_point(aes(x = Time_days, y = ing_percell_perhour - RR_means_for_scaling_ing, color = Treatment_plastic, shape = Dialysis_or_Not), size = 2.5)

ggplot(RR_df_calc)+
  geom_point(aes(x = Time_days, y = pg_perCELL_perHR - RR_means_for_scaling_PP, color = Treatment_plastic, shape = Dialysis_or_Not), size = 2.5)

ggplot(RR_df_calc)+
  geom_point(aes(x = Time_days, y = cell_per_L - RR_means_for_scaling_abund, color = Treatment_plastic, shape = Dialysis_or_Not), size = 2.5)

RR_abund_scaled <- ggplot(RR_df_CF_test_calc)+
  geom_boxplot(aes(x = Time_days, y = cell_per_L - RR_means_for_scaling_abund, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)+
  xlab("Time (days)")+ ylab(expression("Abundance"~"(cell L"^-1~")"))+labs(fill = "Treatment")

RR_PP_scaled <- ggplot(RR_df_calc)+
  geom_boxplot(aes(x = Time_days, y = pg_perCELL_perHR - RR_means_for_scaling_PP, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)+
  xlab("Time (days)")+ ylab(expression("Primary Pruduction"~"(pg Ccell L"^-1~")"))+labs(fill = "Treatment")

RR_feed_scaled <- ggplot(RR_df_CF_test_calc)+
  geom_boxplot(aes(x = Time_days, y = pg_Cfeed_percell_perhour - RR_means_for_scaling_ing_pg, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)+
  xlab("Time (days)")+ ylab(expression("Bacterivory"~"(pg Ccell L"^-1~")"))+labs(fill = "Treatment")

#Need it calculated for the means and sd
RR_df_CF_test_calc$pg_perCELL_perHR_scaled <- RR_df_CF_test_calc$pg_perCELL_perHR - RR_means_for_scaling_PP
RR_df_CF_test_calc$pg_Cfeed_percell_perhour_scaled <- RR_df_CF_test_calc$pg_Cfeed_percell_perhour- RR_means_for_scaling_ing_pg
RR_df_CF_test_calc$cell_per_L_scaled <- RR_df_CF_test_calc$cell_per_L - RR_means_for_scaling_abund

RR_df_calc_grp_scaled <- group_by(RR_df_CF_test_calc, Treatment_plastic, Dialysis_or_Not, Time_days)

RR_df_calc_grpsum_scaled <- summarise(RR_df_calc_grp_scaled,
                               pg_Cfeed_percell_perhour_scaled_MEAN = mean(pg_Cfeed_percell_perhour_scaled),
                               pg_Cfeed_percell_perhour_scaled_SE = sd(pg_Cfeed_percell_perhour_scaled)/sqrt(4),
                               pg_perCELL_perHR_scaled_MEAN = mean(pg_perCELL_perHR_scaled),
                               pg_perCELL_perHR_scaled_SE = sd(pg_perCELL_perHR_scaled)/sqrt(4),
                               cell_per_L_scaled_MEAN = mean(cell_per_L_scaled),
                               cell_per_L_scaled_SE = sd(cell_per_L_scaled)/sqrt(4),
)

ggsave("RR_abund_scaled.pdf", RR_abund_scaled)

ggsave("RR_PP_scaled.pdf", RR_PP_scaled)

ggsave("RR_feed_scaled.pdf", RR_feed_scaled)



# Overall population
ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = overall_feeding, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = overall_pp, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = cell_per_L, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = ugram_Cfeed_percell_perhour_L, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = ugram_Cfixed_per_L_hr, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = ugram_Cfixed_per_L_hr_cell, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

#Carbon Budget Analysis ------------
#Biomass Validation test will be used using the estimated biovolumes of C content per cell from Menden-Deuer and Lessard 2000
#Equation for C based on Volume for each group
# log pgC per cell = log a + b * Vol (um^3)
#OR y = a vol^b --- log y = b log vol + log a, with y = pgC cell', b = slope, and a = y-intercept of the regression equation.
# Using chrysophytes (*Ochromonas sp.*) um^3 ~ 9-15 um in radius
#I will compare using both the less than 9 and greater than 9 considering much of the biomass was smaller and around 10-12 um in radius.
# Chrysophytes: log a = -1.694 and b = 1.218

# Chrysophytes: a = 10^-1.694 and b = 1.218
#final equation is 10^-1.694 * Vol^1.218
#For volumes of each I could create a table with each to create a dataframe with avg radius and estimates of volume and pg C per cell for each one and then calculate based on each...but I will just start with the first one.
#12 um = 7234.56
#8 um = 2143.573

#Done in the code chunk above and I really need to double check this calculation


ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = feed_carbon_percell, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = PP_carbon_percell, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = ugram_Cfixed_per_L_hr, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(mixo_phys_pp_baseR_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = ugram_Cfixed_per_L_hr, fill = Treatment_plastic))

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = cell_per_L, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = total_carbon_incorporated_pg_cell, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = pg_Cfeed_percell_perhour_L/ochromonas_biomass_pg_carbon_percell, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = (pg_perCELL_perHR), fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(mixo_phys_pp_baseR_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = (pg_perCELL_perHR), fill = Treatment_plastic))

mixo_mp_feed$ing_percell_perhour <- mixo_mp_feed$BSP_ingested*2


ggplot(RR_df_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = ochromonas_biomass_pg_carbon_percell, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)
# Problems I am currently facing....
# 1) there isn't enough carbon being consumed to sustain sizes larger than ~8 um
# 2) the estimated feeding is well below those estiamted by the caron group using HKBs
# 3) only 10 - 30% ~ of carbon is from PP and the "rest" should be from feeding but this
#     isn't the case as there isn't enough feeding from the estimated carbon


ggplot(RR_df_CF_test_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = cell_per_L, fill = Treatment_plastic))+
  geom_boxplot(data = mixo_mp_abund,aes(x = as.factor(Time_days), y = cell_per_L,fill = Treatment))
  #facet_grid(~Dialysis_or_Not)

ggplot(mixo_mp_feed_forCOMP)+
  geom_boxplot(aes(x = as.factor(Time_days), y = cell_L, fill = Treatment_plastic), position = position_dodge())+
  geom_boxplot(data = RR_df_CF_test_calc, aes(x = as.factor(Time_days), y = cell_per_L, color = Treatment_plastic), position = position_dodge())+
  scale_color_manual("Incubator", values = scales::hue_pal(h = c(50,250))(5))+
  geom_boxplot(data = mixo_mp_abund_forCOMP, aes(x = as.factor(Time_days), y = cell_per_L, color = Treatment), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)+
  scale_fill_manual("Benchtop", values = scales::hue_pal()(4))

ggplot(mixo_mp_abund_forCOMP)+
  geom_boxplot(aes(x = as.factor(Time_days), y = cell_per_L, fill = Treatment), position = position_dodge())+
  geom_boxplot(data = RR_df_CF_test_calc, aes(x = as.factor(Time_days), y = cell_per_L, color = Treatment_plastic), position = position_dodge())+
  facet_grid(~Dialysis_or_Not)+scale_color_manual("Incubator", values = scales::hue_pal(h = c(50,250))(2))+
  scale_fill_manual("Benchtop", values = scales::hue_pal()(4))


ggplot(RR_df_CF_test_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = pg_Cfeed_percell_perhour/2, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_CF_test_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = pg_perCELL_perHR, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_CF_test_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = feed_carbon_percell/2, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggplot(RR_df_CF_test_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = PP_carbon_percell, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

test_plot <- ggplot(RR_df_CF_test_calc)+
  geom_boxplot(aes(x = as.factor(Time_days), y = PP_carbon_percell, fill = Treatment_plastic))+
  facet_grid(~Dialysis_or_Not)

ggsave("test_plot.pdf",test_plot, path = "/Users/christophercarnivale/Desktop/Dissertation_data/Microplastic_physiological_response")

# Now going to calculate the amount of carbon needed for doubling time
# I need the amount of carbon per cell...I don't have the measured average BUT I can calculate
# this on a range of sizes to see where the size range for growth is limited.
# Doubling time is in days so I need to convverst all of the carbon to day and compare
# For this I need the prior dataset and to accumulate all of the data into figures and
# a story. One from the initial experiment and one from the new experiment

#The doubling time was calculated over a preiod of 2-4 days so rates at day 3
# will be used to estimate carbon needed.

# Feeding df:     mixo_mp_feed
# PP df:          mixo_phys_pp_baseR_calc
# PP GAM inter:   pp_gam_predict


# calculating carbon needed for a single cell to double
carbon_needed_fordoubling_perCELL <- RR_df_CF_test_calc$ochromonas_biomass_pg_carbon_percell[1] * as.numeric(mixo_mp_gc_table[1,3])

carbon_needed_fordoubling_perCELL <- RR_df_CF_test_calc$ochromonas_biomass_pg_carbon_percell[1] * as.numeric(mixo_mp_gc_table[4,3])

# pulling the interpolation data for PP rates
PP_D3_rate_control_pg_percell <- dplyr::filter(pp_gam_predict, Treatment_plastic == "control" & Time_days == 3)[3]

PP_D3_rate_100000_pg_percell <- dplyr::filter(pp_gam_predict, Treatment_plastic == "100000" & Time_days == 3)[3]

PP_D3_rate_1000_pg_percell <- dplyr::filter(pp_gam_predict, Treatment_plastic == "1000" & Time_days == 3)[3]

PP_D3_rate_10_pg_percell <- dplyr::filter(pp_gam_predict, Treatment_plastic == "10" & Time_days == 3)[3]

# converting rates from first exp into Carbon and subsetting as with PP
mixo_mp_feed$ing_percell_perhour_pg <- mixo_mp_feed$ing_percell_perhour*(2.2*10^-13)*1000000000*1000

# I need to calculate the means for control and 10000

# then apply that to the amount of carbon needed for 

mixo_mp_feed$ing_percell_perhour_pg*1.988*24
PP_D3_rate_control_pg_percell*1.988*24

mixo_mp_feed$ing_percell_perhour_pg*2.728*24
PP_D3_rate_100000_pg_percell*2.728*24

mixo_mp_gc_table_df <- as.data.frame(mixo_mp_gc_table)
colnames(mixo_mp_gc_table_df) <- c("Treatment", "Growth_Rate","Doubling_Time") 
mixo_mp_gc_table_df <- mutate(mixo_mp_gc_table_df, 
                              Carbon_to_Double = round(RR_df_CF_test_calc$ochromonas_biomass_pg_carbon_percell[1]*2, digits = 2),
                              )
colnames(mixo_mp_gc_table_df) <- c("Treatment", "Growth Rate (µ)","Doubling Time (Days)", "Carbon to Double") 
ggtexttable(mixo_mp_gc_table_df, rows = NULL, theme = ttheme("light",base_size = 8))

#I want to create the table with all of the data for the Carbon Budget analysis
# I need to summarise and average all of the PP and feeding data

# Carbon Needed to Double is 161.5
RR_df_CF_test_calc$ochromonas_biomass_pg_carbon_percell[1]*2

#Feeding data first
mixo_mp_feed_grp <- group_by(mixo_mp_feed, Treatment_plastic, Time_days)

mixo_mp_feed_grp_summary <- summarise(mixo_mp_feed_grp, 
                                      ing_mean = mean(ing_percell_perhour_pg), 
                                      ing_sd = sd(ing_percell_perhour_pg))

day_3_feeding_rate_percellHr <- mixo_mp_feed_grp_summary[mixo_mp_feed_grp_summary[['Time_days']]==3,]

day_7_feeding_rate_percellHr <- mixo_mp_feed_grp_summary[mixo_mp_feed_grp_summary[['Time_days']]==7,]
# Not needed as the data is interpolated using the GAMs
#mixo_phys_pp_baseR_calc_grp <- group_by(mixo_phys_pp_baseR_calc, Treatment_plastic, Time_days)

#mixo_phys_pp_baseR_calc_grp_sum <- summarise(mixo_phys_pp_baseR_calc_grp,
#                                             PP_mean = mean(ing_percell_perhour_pg), 
#                                             PP_sd = sd(ing_percell_perhour_pg))

# Primary Production data
#Already pulled from interpolation in previous section
day_3_PP_rate_percellHr <- data.frame(Treatment = c("control", "10","1000","100000"),Time_days = 3)
day_3_PP_tibble <- rbind(PP_D3_rate_control_pg_percell,PP_D3_rate_10_pg_percell, PP_D3_rate_1000_pg_percell, PP_D3_rate_100000_pg_percell)
day_3_PP_rate_percellHr$PP_perhrpercell <- round(day_3_PP_tibble$fit,digits = 2)
day_7_PP_rate_percellHr <- mixo_phys_pp_baseR_calc_noB_grp_summary[mixo_phys_pp_baseR_calc_noB_grp_summary[['Time_days']]==7,]
# Add all of the info into the table
mixo_mp_gc_table_df$Ingested_Carbon <- round(day_7_feeding_rate_percellHr$ing_mean*24*.5,digits = 2)
mixo_mp_gc_table_df$PP_Carbon <- round(day_7_PP_rate_percellHr$pp_mean*24, digits = 2)
mixo_mp_gc_table_df <- mutate(mixo_mp_gc_table_df,
                              total_C_uptake = Ingested_Carbon+PP_Carbon)
colnames(mixo_mp_gc_table_df) <- c("Treatment", "Growth Rate (µ)","Doubling Time (days)", "Carbon to Double",
                                   "Ingested Carbon","PP Carbon","Total Carbon Uptake") 

CB_name_row <- c("Treatment", "Growth Rate (µ)","Doubling Time (days)", "Carbon to Double",
                             "Ingested Carbon","PP Carbon","Total Carbon Uptake")

mixo_mp_gc_table_df <- rbind(CB_name_row,mixo_mp_gc_table_df)

rownames(mixo_mp_gc_table_df) <- NULL
colnames(mixo_mp_gc_table_df) <- NULL
mixo_mp_gc_table_df_noCdoub <- mixo_mp_gc_table_df[,-c(3,4,7)]
Carbon_budget_table <- ggtexttable(mixo_mp_gc_table_df_noCdoub, rows = NULL, theme = ttheme("blank")) %>% 
  tab_add_footnote(text = "*Acquired carbon reported as daily pg of Carbon per cell; Ingested Carbon received a 50% penalty due to respiration",
                   size = 6, vjust = -.5, hjust = 0.95)%>% 
  tab_add_hline(at.row = 1:2, row.side = c("top"), linewidth = 3) %>% 
  tab_add_hline(at.row = 5, row.side = "bottom", linewidth = 2) %>% table_cell_font(row = 1, column = 1:4,face = "bold", size = 10)

ggsave("Carbon_budget_table.pdf",Carbon_budget_table, width = 8, height = 4.5)

# Summary statistics and model testiong for differences between means -------
# Repeated Measures ANOVA --------
summary(aov(cell_per_L ~ Treatment, data = filter(mixo_mp_abund_noB, Time_days == 7)))

summary(aov(log(cell_per_L) ~ Treatment + Error(Rep), data = mixo_mp_abund))
#Days 5 and 6 are the only days of significant difference between abundances
# Repeated measures ANOVA indicated significant reductions by plastic treatment

# Growth rate Analysis -------
#ANOVA
#Need to fix df for ANOVA
mixo_mp_gc_table_df <- mixo_mp_gc_table_df[-1,]
colnames(mixo_mp_gc_table_df) <- c("Treatment", "Growth_Rate","Doubling_Time", "Carbon_to_Double",
                                   "Ingested_Carbon","PP_Carbon","Total_Carbon_Uptake") 
summary(aov(Doubling_Time ~ Treatment, data = mixo_mp_gc_table_df))

gc_se

#Bacterivory tests -------
summary(aov(ing_percell_perhour_pg ~ Treatment_plastic, data = filter(mixo_mp_feed, Time_days == 0)))

summary(aov(log(ing_percell_perhour_pg) ~ Treatment_plastic + Error(treatment_rep), data = mixo_mp_feed))

summary(mixo_mp_feed_gam3)

anova.gam(mixo_mp_feed_gam3)

#Primary Production tests ---------
summary(aov(pg_perCELL_perHR_scaledbyCntlmean ~ Treatment_plastic, data = filter(mixo_phys_pp_baseR_calc_noB, Time_days == 4)))
summary(aov(pg_perCELL_perHR ~ Treatment_plastic, data = filter(mixo_phys_pp_baseR_calc_noB, Time_days == 7)))

summary(aov(pg_perCELL_perHR ~ Treatment_plastic + Error(Rep), data = mixo_phys_pp_baseR_calc_noB))

summary(mixo_mp_pp_gam_byPlastic_k4)

anova.gam(mixo_mp_pp_gam_byPlastic_k4)

# Abundance Summary statistics -------
mixo_mp_abund_noB_grp <- group_by(mixo_mp_abund_noB, Treatment, Time_days)

mixo_mp_abund_noB_grp_summary <- summarise(mixo_mp_abund_noB_grp, 
                                      abund_mean = mean(cell_per_L), 
                                      abund_sd = sd(cell_per_L))

raw_abund_plot <- ggplot(mixo_mp_abund_noB_grp_summary)+
  geom_point(data = mixo_mp_abund_noB, aes(as.numeric(Time_days), y = cell_per_L, color = Treatment), size = 0.5)+
  geom_line(aes(x = as.numeric(Time_days), y = abund_mean, color = Treatment))+
  geom_errorbar(aes(x = as.numeric(Time_days),ymin = abund_mean-abund_sd, ymax = abund_mean+abund_sd, color = Treatment), linewidth = 0.2)+
  xlab("Time (days)")+ ylab(expression("Abundance"~"(cell L"^-1~")"))

ggsave("raw_abund_plot.pdf", raw_abund_plot)
# Bacterivory Summary statistics --------
mixo_mp_feed_noB$ing_percell_perhour_pg <- mixo_mp_feed_noB$ing_percell_perhour*(2.2*10^-13)*1000000000*1000
mixo_mp_feed_noB_grp <- group_by(mixo_mp_feed_noB, Treatment_plastic, Time_days)

mixo_mp_feed_noB_grp_summary <- summarise(mixo_mp_feed_noB_grp, 
                                      ing_mean = mean(ing_percell_perhour_pg), 
                                      ing_sd = sd(ing_percell_perhour_pg))

# Primary Production Summary statistics ----------
mixo_phys_pp_baseR_calc_noB_grp <- group_by(mixo_phys_pp_baseR_calc_noB, Treatment_plastic, Time_days)

mixo_phys_pp_baseR_calc_noB_grp_summary <- summarise(mixo_phys_pp_baseR_calc_noB_grp, 
                                      pp_mean = mean(pg_perCELL_perHR), 
                                      pp_sd = sd(pg_perCELL_perHR))

# Repeat but now for the scaled from mean values ---------
# Abundance Summary statistics -------

mixo_mp_abund_noB_grp <- group_by(mixo_mp_abund_noB, Treatment, Time_days)

mixo_mp_abund_noB_grp_summary_scale <- summarise(mixo_mp_abund_noB_grp, 
                                           abund_mean = mean(cell_per_L_scaledbtCntlmean), 
                                           abund_sd = sd(cell_per_L_scaledbtCntlmean))

# Bacterivory Summary statistics --------
mixo_mp_feed_ctlgrpmean_noB_pg <- c(rep(mixo_mp_feed_noB_grp_summary$ing_mean[1], 15), rep(mixo_mp_feed_noB_grp_summary$ing_mean[2], 15),
                                    rep(mixo_mp_feed_noB_grp_summary$ing_mean[3], 15))

mixo_mp_feed_noB$ing_percell_perhour_pg_scaled <- mixo_mp_feed_noB$ing_percell_perhour_pg-mixo_mp_feed_ctlgrpmean_noB_pg

ing_pg_scaled_plot <- ggplot(mixo_mp_feed_noB)+
  geom_boxplot(aes(x = as.factor(Time_days), y = ing_percell_perhour_pg_scaled, fill = Treatment_plastic))+
  xlab("Time (days)")+ ylab(expression("Bacterivory (pg C"~"cell"^-1~"hr"^-1~")"))+labs(fill = "Treatment")


raw_ing_pg_plot <- ggplot(mixo_mp_feed_noB)+
  geom_boxplot(aes(x = as.factor(Time_days), y = ing_percell_perhour_pg, fill = Treatment_plastic))+
  xlab("Time (days)")+ylab(expression("Bacterivory (pg C"~"cell"^-1~"hr"^-1~")"))+labs(fill = "Treatment")

ggsave("ing_pg_scaled_plot.pdf",ing_pg_scaled_plot)
ggsave("raw_ing_pg_plot.pdf",raw_ing_pg_plot)

mixo_mp_feed_noB_grp <- group_by(mixo_mp_feed_noB, Treatment_plastic, Time_days)

mixo_mp_feed_noB_grp_summary_scaled <- summarise(mixo_mp_feed_noB_grp, 
                                          ing_mean = mean(ing_percell_perhour_pg_scaled), 
                                          ing_sd = sd(ing_percell_perhour_pg_scaled))

# Primary Production Summary statistics ---------
mixo_phys_pp_baseR_calc_noB_grp <- group_by(mixo_phys_pp_baseR_calc_noB, Treatment_plastic, Time_days)

mixo_phys_pp_baseR_calc_noB_grp_summary_scaled <- summarise(mixo_phys_pp_baseR_calc_noB_grp, 
                                                     pp_mean = mean(pg_perCELL_perHR_scaledbyCntlmean), 
                                                     pp_sd = sd(pg_perCELL_perHR_scaledbyCntlmean))

# Extra coding ideas if asked for -----------
# Repeat analysis for second experiment to make sure it's consistent....IF NEEDED
RR_df_CF_test_calc_CB <- mutate(RR_df_CF_test_calc, 
                                Ingested_Carbon = pg_Cfeed_percell_perhour*24*0.5,
                                PP_Carbon = pg_perCELL_perHR*24)
RR_df_CF_test_calc_CB_start <- RR_df_CF_test_calc_CB[RR_df_CF_test_calc_CB$Time_days == 0,]
RR_df_CF_test_calc_CB_end <- RR_df_CF_test_calc_CB[RR_df_CF_test_calc_CB$Time_days == 7,]

RR_df_CF_test_calc_CB_end$Growth_Rate_perday <-log(RR_df_CF_test_calc_CB_end$`cell_per_L`/RR_df_CF_test_calc_CB_start$`cell_per_L`)/7
RR_df_CF_test_calc_CB_end$Doubling_time_perday <- log(2)/log(1+RR_df_CF_test_calc_CB_end$Growth_Rate_perday)

RR_df_CF_test_calc_CB_end$Ingested_Carbon <- RR_df_CF_test_calc_CB_end$Ingested_Carbon*RR_df_CF_test_calc_CB_end$Doubling_time_perday
RR_df_CF_test_calc_CB_end$PP_Carbon <- RR_df_CF_test_calc_CB_end$PP_Carbon*RR_df_CF_test_calc_CB_end$Doubling_time_perday


#RR_df_CF_test_calc_CB_start_sum <- group_by(RR_df_CF_test_calc_CB_start[,c(2:5,48:51)], Treatment_plastic, Time_days, Dialysis_or_Not) %>% 
#  summarise(Ingested_Carbon = mean(Ingested_Carbon),
#            PP_Carbon = mean(PP_Carbon))

RR_df_CF_test_calc_CB_end_sum <- group_by(RR_df_CF_test_calc_CB_end[,c(2:5, 48:51)], Treatment_plastic, Time_days, Dialysis_or_Not) %>% 
  summarise(Ingested_Carbon = mean(Ingested_Carbon),
            PP_Carbon = mean(PP_Carbon),
            Doubling_Time = mean(Doubling_time_perday),
            Growth_Rate = mean(Growth_Rate_perday)) %>% 
  mutate(Total_C_Intake = Ingested_Carbon+PP_Carbon)

RR_df_CF_test_calc_CB_end_sum_order <- RR_df_CF_test_calc_CB_end_sum[,c('Treatment_plastic','Dialysis_or_Not',
                                 'Growth_Rate','Doubling_Time','Ingested_Carbon','PP_Carbon','Total_C_Intake')]

Carbon_demand_RR_table <- ggtexttable(RR_df_CF_test_calc_CB_end_sum_order, rows = NULL, theme = ttheme("blank",base_size = 9)) %>% 
  tab_add_footnote(text = "*All measures are report as pg of Carbon per cell; Ingested Carbon is reported as a 50% penalty due to respiration",
                   size = 8, vjust = -.5, hjust = 0.945)%>% 
  tab_add_hline(at.row = 1:2, row.side = c("top"), linewidth = 3) %>% 
  tab_add_hline(at.row = 5, row.side = "bottom", linewidth = 2) 
#%>% table_cell_font(row = 1, column = 5,face = "bold", size = 10)

ggsave("Carbon_demand_RR_table.pdf",Carbon_demand_RR_table, width = 8, height = 4.5)


#mixo_mp_gc_2_4$Growth_Rate_perday <-log(RR_df_CF_test_calc_CB_end$`Day_4_Cells/L`/RR_df_CF_test_calc_CB_end$`Day_2_Cells/L`)/2
#mixo_mp_gc_2_4$Growth_Rate_perday[16] <- 0
#mixo_mp_gc_2_4$Doubling_time_perday <- log(2)/log(1+mixo_mp_gc_2_4$Growth_Rate_perday)
#mixo_mp_gc_2_4$Doubling_time_perday[16] <- 0


#Growth_Rate <- tapply(mixo_mp_gc_2_4$Growth_Rate_perday, mixo_mp_gc_2_4$Treatment, mean)

#Doubling_time_adj <- log(2)/log(1+Growth_Rate)

# Converting all plots to point plots -------------
# Using the dataframes from the summary statistics would be usable for the point plots
#mixo_phys_pp_baseR_calc
#mixo_mp_feed
#mixo_mp_abund

# Direct experiment---------
#These plots are good and fixed
Mphys_Figure_1A <- ggplot(filter(mixo_mp_abund_noB_grp_summary_scale, Treatment != "Control", Time_days == "0" | Time_days == "7"))+
  geom_col(aes(x = Time_days, y = abund_mean, fill = Treatment), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = abund_mean - abund_sd, ymax = abund_mean + abund_sd, color = Treatment), position = position_dodge())+
  scale_color_manual(values = c("black","black","black"))+scale_fill_manual(values = c("#7CAE00","#00BFC4","#C77CFF"))+
  xlab("Time (days)")

Mphys_Figure_1B <-ggplot(filter(mixo_phys_pp_baseR_calc_noB_grp_summary_scaled, Treatment_plastic != "control", Time_days == "0" | Time_days == "7"))+
  geom_col(aes(x = Time_days, y = pp_mean, fill = Treatment_plastic), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = pp_mean - pp_sd, ymax = pp_mean + pp_sd, color = Treatment_plastic), position = position_dodge())+
  scale_color_manual(values = c("black","black","black"))+scale_fill_manual(values = c("#7CAE00","#00BFC4","#C77CFF"))+
  xlab("Time (days)")

Mphys_Figure_1C <-ggplot(filter(mixo_mp_feed_noB_grp_summary_scaled, Treatment_plastic != "control", Time_days == "0" | Time_days == "7"))+
  geom_col(aes(x = Time_days, y = ing_mean, fill = Treatment_plastic), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = ing_mean - ing_sd, 
                    ymax = ing_mean + ing_sd, color = Treatment_plastic), position = position_dodge())+
  scale_color_manual(values = c("black","black","black"))+scale_fill_manual(values = c("#7CAE00","#00BFC4","#C77CFF"))+
  xlab("Time (days)")


# Repeating but with the raw data for supplementary figures
mixo_mp_abund_noB_grp_summary_se_log <- summarise(mixo_mp_abund_noB_grp, 
                                           abund_mean = mean(log(cell_per_L)), 
                                           abund_se = sd(log(cell_per_L))/sqrt(3))

mixo_mp_abund_noB_grp_summary_se <- summarise(mixo_mp_abund_noB_grp, 
                                                  abund_mean = mean(cell_per_L), 
                                                  abund_se = sd(cell_per_L)/sqrt(3))

mixo_mp_feed_noB_grp_summary_se <- summarise(mixo_mp_feed_noB_grp, 
                                          ing_mean = mean(ing_percell_perhour_pg), 
                                          ing_se = sd(ing_percell_perhour_pg)/sqrt(3))

mixo_phys_pp_baseR_calc_noB_grp_summary_se <- summarise(mixo_phys_pp_baseR_calc_noB_grp, 
                                                     pp_mean = mean(pg_perCELL_perHR), 
                                                     pp_se = sd(pg_perCELL_perHR)/sqrt(3))

Mphys_Figure_2B <-ggplot(mixo_phys_pp_baseR_calc_noB_grp_summary_se)+
  geom_col(aes(x = Time_days, y = pp_mean, fill = Treatment_plastic), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = pp_mean - pp_se, ymax = pp_mean + pp_se, color = Treatment_plastic), position = position_dodge())+
  scale_color_manual(values = c("black","black","black","black"))+
  xlab("Time (days)")

Mphys_Figure_2A <-ggplot(mixo_mp_abund_noB_grp_summary_se)+
  geom_col(aes(x = Time_days, y = abund_mean, fill = Treatment), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = abund_mean - abund_se, ymax = abund_mean + abund_se, color = Treatment), position = position_dodge())+
  scale_color_manual(values = c("black","black","black","black"))+xlab("Time (days)")

Mphys_Figure_2C <-ggplot(mixo_mp_feed_noB_grp_summary_se)+
  geom_col(aes(x = Time_days, y = ing_mean, fill = Treatment_plastic), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = ing_mean - ing_se, 
                    ymax = ing_mean + ing_se, color = Treatment_plastic), position = position_dodge())+
  scale_color_manual(values = c("black","black","black","black"))+
  xlab("Time (days)")

# Indirect experiment----------
# Data for bar plots
# RR_df_calc_grp_scaled <- group_by(RR_df_CF_test_calc, Treatment_plastic, Dialysis_or_Not, Time_days)

RR_df_calc_grpsum_final <- summarise(RR_df_calc_grp_scaled,
                                      ing_scaled_mean = mean(pg_Cfeed_percell_perhour_scaled),
                                      ing_scaled_se = sd(pg_Cfeed_percell_perhour_scaled)/sqrt(4),
                                      pp_scaled_mean = mean(pg_perCELL_perHR_scaled),
                                      pp_scaled_se = sd(pg_perCELL_perHR_scaled)/sqrt(4),
                                      abund_scaled_mean = mean(cell_per_L_scaled),
                                      abund_scaled_se = sd(cell_per_L_scaled)/sqrt(4),
                                      ing_mean = mean(pg_Cfeed_percell_perhour),
                                      ing_se = sd(pg_Cfeed_percell_perhour/sqrt(4)),
                                      pp_mean = mean(pg_perCELL_perHR),
                                      pp_se = sd(pg_perCELL_perHR)/sqrt(4),
                                      abund_mean = mean(cell_per_L),
                                      abund_se = sd(cell_per_L)/sqrt(4),
)

# Effect size
Mphys_Figure_S1B <-ggplot(filter(RR_df_calc_grpsum_final, Treatment_plastic != "control"))+
  geom_col(aes(x = Time_days, y = pp_scaled_mean, fill = Treatment_plastic), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = pp_scaled_mean - pp_scaled_se, ymax = pp_scaled_mean + pp_scaled_se, color = Treatment_plastic), position = position_dodge())+
  xlab("Time (days)")+scale_color_manual(values = "black")+scale_fill_manual(values = "#C77CFF")+facet_wrap(~Dialysis_or_Not)

Mphys_Figure_S1A <-ggplot(filter(RR_df_calc_grpsum_final, Treatment_plastic != "control"))+
  geom_col(aes(x = Time_days, y = abund_scaled_mean, fill = Treatment_plastic), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = abund_scaled_mean - abund_scaled_se, ymax = abund_scaled_mean + abund_scaled_se, color = Treatment_plastic), position = position_dodge())+
  scale_color_manual(values = "black")+scale_fill_manual(values = "#C77CFF")+
  xlab("Time (days)")+facet_wrap(~Dialysis_or_Not)

Mphys_Figure_S1C <-ggplot(filter(RR_df_calc_grpsum_final, Treatment_plastic != "control"))+
  geom_col(aes(x = Time_days, y = ing_scaled_mean, fill = Treatment_plastic), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = ing_scaled_mean - ing_scaled_se, 
                    ymax = ing_scaled_mean + ing_scaled_se, color = Treatment_plastic), position = position_dodge())+
  scale_color_manual(values = "black")+scale_fill_manual(values = "#C77CFF")+
  xlab("Time (days)")+facet_wrap(~Dialysis_or_Not)

# Raw data
Mphys_Figure_S2B <-ggplot(RR_df_calc_grpsum_final)+
  geom_col(aes(x = Time_days, y = pp_mean, fill = Treatment_plastic), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = pp_mean - pp_se, ymax = pp_mean + pp_se, color = Treatment_plastic), position = position_dodge())+
  scale_color_manual(values = c("black","black","black","black"))+
  xlab("Time (days)")+facet_wrap(~Dialysis_or_Not)+scale_fill_manual(values = c("#F8766D","#C77CFF"))

Mphys_Figure_S2A <-ggplot(RR_df_calc_grpsum_final)+
  geom_col(aes(x = Time_days, y = abund_mean, fill = Treatment_plastic), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = abund_mean - abund_se, ymax = abund_mean + abund_se, color = Treatment_plastic), position = position_dodge())+
  scale_color_manual(values = c("black","black","black","black"))+xlab("Time (days)")+
  facet_wrap(~Dialysis_or_Not)+scale_fill_manual(values = c("#F8766D","#C77CFF"))

Mphys_Figure_S2C <-ggplot(RR_df_calc_grpsum_final)+
  geom_col(aes(x = Time_days, y = ing_mean, fill = Treatment_plastic), position = position_dodge())+
  geom_errorbar(aes(x = Time_days, ymin = ing_mean - ing_se, 
                    ymax = ing_mean + ing_se, color = Treatment_plastic), position = position_dodge())+
  scale_color_manual(values = c("black","black","black","black"))+
  xlab("Time (days)")+facet_wrap(~Dialysis_or_Not)+scale_fill_manual(values = c("#F8766D","#C77CFF"))

# Saving all of the new figures
ggsave("Mphys_Figure_1A.pdf",Mphys_Figure_1A)
ggsave("Mphys_Figure_1B.pdf",Mphys_Figure_1B)
ggsave("Mphys_Figure_1C.pdf",Mphys_Figure_1C)
ggsave("Mphys_Figure_2A.pdf",Mphys_Figure_2A)
ggsave("Mphys_Figure_2B.pdf",Mphys_Figure_2B)
ggsave("Mphys_Figure_2C.pdf",Mphys_Figure_2C)
ggsave("Mphys_Figure_S1A.pdf",Mphys_Figure_S1A)
ggsave("Mphys_Figure_S1B.pdf",Mphys_Figure_S1B)
ggsave("Mphys_Figure_S1C.pdf",Mphys_Figure_S1C)
ggsave("Mphys_Figure_S2A.pdf",Mphys_Figure_S2A)
ggsave("Mphys_Figure_S2B.pdf",Mphys_Figure_S2B)
ggsave("Mphys_Figure_S2C.pdf",Mphys_Figure_S2C)

# I need to redo the Carbon demand analysis again with the indirect experiment
# I'm going to need this to confirm the findings from the redo experiment.

ggsave("Mphys_Figure_1A_D7.pdf",Mphys_Figure_1A)
ggsave("Mphys_Figure_1B_D7.pdf",Mphys_Figure_1B)
ggsave("Mphys_Figure_1C_D7.pdf",Mphys_Figure_1C)

