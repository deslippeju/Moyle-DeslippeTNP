### R Script for Moyle and Deslippe 2023: Invasion Alters Plant and Mycorrhizal Communities in an Alpine Tussock Grassland######
##This script includes analyses for C. rubra morphometric data; root colonisation; soil nutrients; plant community composition
# Victoria University of Wellington 
# 2023

#Libraries
library(tidyverse)
library(readxl)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(broom)
library(ggfortify)
library(markdown)
library(lme4)
library(lmerTest)
library(car)

#########C. rubra morphometrics ###########

### Tussock Diameter #####

library(readxl)
Tussock_Data <- read_excel("TNP_Manuscript_Data.xlsx", 
                           sheet = "Tussock Info")
View(Tussock_Data)

library(grid)
# How to alter the x- and y-axis labels
ggplot(Tussock_Data,aes(x=Heather_Density,y=Tussock_Diameter))+
  geom_point(size=2)+
  geom_smooth(method = 'lm')+
  xlab(~italic(C.vulgaris)*" density (%)")+
  ylab(~italic(C.rubra)* " diameter (cm)")+
  theme_bw()+
  theme(axis.text.x= element_text(colour="black", size=14)) +
  theme(axis.text.y= element_text(colour="black", size=14)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15))

# ADD this code on to make ticks inwards
theme(
  axis.ticks.length.x = unit(-.25, "cm"),
  axis.ticks.length.y = unit(-.25, "cm"),
  axis.text.x = element_text(margin = margin(t = .3, unit = "cm")),
  axis.text.y = element_text(margin = margin(r = .3, unit = "cm"))
)

# Save and export the graph 
ggsave("T1(Diameter).png", width = 10, height = 8, dpi = 1000)

#Simple Linear Regression
model_tussock<-lm(Tussock_Diameter ~ Heather_Density,
                  data=Tussock_Data)

anova(model_tussock)
summary(model_tussock)

# Checking assumptions
plot(model_tussock, 1)

model.diag.metrics <- augment(model_tussock)
head(model.diag.metrics)

ggplot(model.diag.metrics, aes(Heather_Density, Tussock_Diameter)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = Heather_Density, yend = .fitted), color = "red", size = 0.3)

par(mfrow = c(2, 2))
plot(model_tussock)


#lmer to account for study design

TDLmer1 <- lmer(Tussock_Diameter ~ Heather_Density + (1|Geog_Area), data=Tussock_Data) 

anova(TDLmer1)
summary(TDLmer1)

### Tussock Cover #####

library(grid)
# How to alter the x- and y-axis labels
ggplot(Tussock_Data,aes(x=Heather_Density,y=Tussock_Cover))+
  geom_point(size=2)+
  geom_smooth(method = 'lm')+
  xlab(~italic(C.vulgaris)*" density (%)")+
  ylab(~italic(C.rubra)* " cover")+
  theme_bw()+
  theme(axis.text.x= element_text(colour="black", size=14)) +
  theme(axis.text.y= element_text(colour="black", size=14)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15))


# Save and export the graph 
ggsave("Tussock dimater.png", width = 10, height = 8, dpi = 1000)

#Simple Linear Regression
mod_tus_cov<-lm(Tussock_Cover ~ Heather_Density,
                  data=Tussock_Data)

anova(mod_tus_cov)
summary(mod_tus_cov)

# Checking assumptions
plot(mod_tus_cov, 1)

model.diag.metrics <- augment(mod_tus_cov)
head(model.diag.metrics)

ggplot(model.diag.metrics, aes(Heather_Density, Tussock_Cover)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = Heather_Density, yend = .fitted), color = "red", size = 0.3)

par(mfrow = c(2, 2))
plot(mod_tus_cov)


#lmer to account for study design

TCLmer1 <- lmer(Tussock_Cover ~ Heather_Density + (1|Geog_Area), data=Tussock_Data) 

anova(TCLmer1)
summary(TCLmer1)


### Tussock Height #####

library(grid)
# How to alter the x- and y-axis labels
ggplot(Tussock_Data,aes(x=Heather_Density,y=Tussock_Height))+
  geom_point(size=2)+
  geom_smooth(method = 'lm')+
  xlab(~italic(C.vulgaris)*" density (%)")+
  ylab(~italic(C.rubra)* " cover")+
  theme_bw()+
  theme(axis.text.x= element_text(colour="black", size=14)) +
  theme(axis.text.y= element_text(colour="black", size=14)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15))


# Save and export the graph 
ggsave("Tussock height.png", width = 10, height = 8, dpi = 1000)

#Simple Linear Regression
mod_height<-lm(Tussock_Height ~ Heather_Density,
                data=Tussock_Data)

anova(mod_height)
summary(mod_height)

# Checking assumptions
plot(mod_height, 1)

model.diag.metrics <- augment(mod_height)
head(model.diag.metrics)

ggplot(model.diag.metrics, aes(Heather_Density, Tussock_Height)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = Heather_Density, yend = .fitted), color = "red", size = 0.3)

par(mfrow = c(2, 2))
plot(mod_height)


#lmer to account for study design

THLmer1 <- lmer(Tussock_Height ~ Heather_Density + (1|Geog_Area), data=Tussock_Data) 

anova(THLmer1)
summary(THLmer1)

######### C. rubra Root colonisation data###########

#Set working directory
setwd("/Users/deslipju/Library/CloudStorage/OneDrive-VictoriaUniversityofWellington-STAFF/Manuscripts/Darby_heather invasion/DATA/Rcode+data")

#Import Data Sheet
library(readxl)

Root_Colonisation_2 <- read_excel("TNP_Manuscript_Data.xlsx", 
                                  sheet = "Root_Colonisation_2")

View(Root_Colonisation_2)

Root_Colonisation_2$Region <- factor(Root_Colonisation_2$Region)
Root_Colonisation_2$Site <- factor(Root_Colonisation_2$Site)
Root_Colonisation_2$Area <- factor(Root_Colonisation_2$Area)
Root_Colonisation_2$CV_Density_Cat <- factor(Root_Colonisation_2$CV_Density_Cat)


##### Plotting Root colonisation ~ heather density####

ggplot(Root_Colonisation_2,aes(x=Heather_Density,y=Fungal_Colonisation))+
  geom_point(size=2)+
  geom_smooth(method = 'lm')+
  xlab(~italic(C.vulgaris)*" density (%)")+
  ylab("Fungal colonisation (%)")+
  theme_bw()+
  theme(axis.text.x= element_text(colour="black", size=14)) +
  theme(axis.text.y= element_text(colour="black", size=14)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15)) 

# Save and export the graph 
ggsave("T1(RootColonisation(.png", width = 10, height = 8, dpi = 1000)

####Simple Linear Regression 
m1<-lm(Fungal_Colonisation ~ Heather_Density,
                 data=Root_Colonisation_2)

# anova table 
anova(m1)
summary(m1)

# Checking assumptions
plot(m1, 1)

model.diag.metrics <- augment(m1)
head(model.diag.metrics)

ggplot(model.diag.metrics, aes(Heather_Density, Fungal_Colonisation)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = Heather_Density, yend = .fitted), color = "red", size = 0.3)

par(mfrow = c(2, 2))
plot(m1)

##log transforming data##
m2<-lm(log(Fungal_Colonisation) ~ Heather_Density,
       data=Root_Colonisation_2)

# anova table 
anova(m2)
summary(m2)

# Checking assumptions
plot(m2, 1)

model.diag.metrics <- augment(m2)
head(model.diag.metrics)

ggplot(model.diag.metrics, aes(Heather_Density, log(Fungal_Colonisation))) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = Heather_Density, yend = .fitted), color = "red", size = 0.3)

par(mfrow = c(2, 2))
plot(m2)

####lmer models accounting for root colonisation with nutrients as co-variates and random effects structure for study design ####


####Deciding the best random effects structure

Lmer1 <- lmer(log(Fungal_Colonisation) ~ Heather_Density+ (1|Region), data=Root_Colonisation_2)
Lmer2 <- lmer(log(Fungal_Colonisation) ~ Heather_Density+ (1|Area), data=Root_Colonisation_2)


#Compare models
AIC(Lmer1,Lmer2)

summary(Lmer2)
anova(Lmer2)

###Considering nutrient as possible co-variates

Lmer6 <- lmer(log(Fungal_Colonisation) ~ Heather_Density+ Nitrogen + Olsen_P+ Carbon+ CN+ (1|Area), data=Root_Colonisation_2)
Lmer7 <- lmer(log(Fungal_Colonisation) ~ Heather_Density+ Nitrogen + Olsen_P+ Carbon+  (1|Area), data=Root_Colonisation_2)
Lmer8 <- lmer(log(Fungal_Colonisation) ~ Heather_Density+  Olsen_P+ CN  + (1|Area), data=Root_Colonisation_2)

AIC(Lmer6,Lmer7,Lmer8)

#best model
step_result <- step(Lmer7)
step_result
final_model7 <- get_model(step_result)
#this final_model7 <- lmer(log(Fungal_Colonisation) ~ Heather_Density + (1 | Area)) 

summary(final_model7)
anova(final_model7)


########Soil Nutrients #########

Soil_Nut <- read_excel("TNP_Manuscript_Data.xlsx", 
                       sheet = "Soil Nutrients")
View(Soil_Nut)

#Olsen P

nut_Lmer1 <- lmer(Olsen_P ~ Avg_heather_density+ (1|Area), data=Soil_Nut) 
nut_Lmer2<- lmer(Olsen_P ~ Avg_heather_density+ (1|Region), data=Soil_Nut)

#Compare random effects structures
AIC(nut_Lmer1,nut_Lmer2)

#test-results best model
summary(nut_Lmer2)
anova(nut_Lmer2)

#Anova-like table for random-effect term using likelihood ratio tests
ranova(nut_Lmer2)

#Carbon

nut_Lmer3 <- lmer(Carbon ~ Avg_heather_density+ (1|Area), data=Soil_Nut) 
nut_Lmer4<- lmer(Carbon ~ Avg_heather_density+ (1|Region), data=Soil_Nut)

#Singular

#test-results best model
summary(nut_Lmer3)
anova(nut_Lmer3)

#Nitrogen

nut_Lmer5 <- lmer(Nitrogen ~ Avg_heather_density+ (1|Area), data=Soil_Nut) 
nut_Lmer6<- lmer(Nitrogen ~ Avg_heather_density+ (1|Region), data=Soil_Nut)

#Compare random effects structures
AIC(nut_Lmer5,nut_Lmer6)

#test-results best model
summary(nut_Lmer5)
anova(nut_Lmer5)

#Anova-like table for random-effect term using likelihood ratio tests
ranova(nut_Lmer5)

#C:N

nut_Lmer7 <- lmer(CN ~ Avg_heather_density+ (1|Area), data=Soil_Nut) 
nut_Lmer8<- lmer(CN ~ Avg_heather_density+ (1|Region), data=Soil_Nut)

#Singular

#test-results 
summary(nut_Lmer7)
anova(nut_Lmer7)


### heather density had no effect on any soil nutrient, therefore...
###Simple lm for study design

#Olsen_P

nut_Lm1 <- lm(Olsen_P ~ Area, data=Soil_Nut) 
nut_Lm2<- lm(Olsen_P ~ Region, data=Soil_Nut)

#Compare random effects structures
AIC(nut_Lm1,nut_Lm2)

#test-results best model
summary(nut_Lm2)
anova(nut_Lm2)

#Carbon

nut_Lm3 <- lm(Carbon ~ Area, data=Soil_Nut) 
nut_Lm4<- lm(Carbon ~ Region, data=Soil_Nut)

#Compare random effects structures
AIC(nut_Lm3,nut_Lm4)

#test-results best model
summary(nut_Lm4)
anova(nut_Lm4)

#Nitrogen

nut_Lm5 <- lm(Nitrogen ~ Area, data=Soil_Nut) 
nut_Lm6<- lm(Nitrogen ~ Region, data=Soil_Nut)

#Compare random effects structures
AIC(nut_Lm5,nut_Lm6)

#test-results best model
summary(nut_Lm6)
anova(nut_Lm6)

#CN

nut_Lm7 <- lm(CN ~ Area, data=Soil_Nut) 
nut_Lm8<- lm(CN ~ Region, data=Soil_Nut)

#Compare random effects structures
AIC(nut_Lm7,nut_Lm8)

#test-results best model
summary(nut_Lm8)
anova(nut_Lm8)


##########Plant community ##########
# Code for Figure 1 B: Percent cover of subordinate plant community by heather density category (high, medium and low), by plant structural class
# excludes the study species C. vulgaris and C. rubra) 

library(readxl)
structuralClass <- read_excel("TNP_Manuscript_Data.xlsx", 
                              sheet = "Plantcommunity(NOcalvul_chirub)")
View(structuralClass)

# Stacked bargraph 
ggplot(data=structuralClass,aes(x = Heather_Category, y=Plant_Species_cover/25, fill=Structural_class)) +
  geom_bar(stat = "identity") +
  labs(x=~italic(C.vulgaris)*" density (%)",
       y="Plant species cover (%)",
       fill="Structural Class") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.2, 0.8)) +
  theme(axis.text.x= element_text(colour="black", size=14)) +
  theme(axis.text.y= element_text(colour="black", size=14)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15)) +
  scale_x_discrete( limits = c("High", "Medium", "Low"))


