### R Script for Moyle and Deslippe 2023: Invasion Alters Plant and Mycorrhizal Communities in an Alpine Tussock Grassland######
# Victoria University of Wellington 
# 2023

######### TNP PLFA & NLFA DATA ############

###Set working directory ###
setwd("/Users/deslipju/Library/CloudStorage/OneDrive-VictoriaUniversityofWellington-STAFF/Manuscripts/Darby_heather invasion/DATA/Rcode+data")

#### Load libraries ####
library(tidyverse)
library(readxl)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(broom)
library(ggfortify)
library(markdown)

# Import the PLFA data from GCMS output
GCMS_PLFA <- read_excel("TNP_Manuscript_Data.xlsx", 
                        sheet = "GCMS PLFA")
View(GCMS_PLFA)

# Import the FAME biomarker excel sheet 

library(readr)
FAME_biomarker <- read_csv("FAME_RFF_biomarker.csv")
View(FAME_biomarker)

# Merge the FAME sheet with the TNP PLFA information sheet

PLFA_FAME<-merge(GCMS_PLFA, FAME_biomarker, by="Name")

# Import the sheet with the Information about each of the TNP sites, e.g heather density

library(readxl)
plfameta <- read_excel("TNP_Manuscript_Data.xlsx", 
                       sheet = "PLFA_META")
View(plfameta)

# Merge the TNP site information with both the PLFA data and the FAME data
TNP_Data<-merge(PLFA_FAME, plfameta, by.x="Sample",by.y =  "Sample_No")

##### FAME correction ######
# Use Fame Sheet to correct and gain the concentration of each of the lipids
#Calculate correction factor for each sample
correction <- TNP_Data %>%
  select(c(Sample,Area, Name))%>%
  filter(Name=="19:0")%>%
  mutate(IS_area=1*Area) %>%  
  select(c(Sample,IS_area))

head(correction)  

#Calculate the concentration of each FAME
TNP_Concentrations<-  left_join(TNP_Data, correction, by ="Sample")%>%
  mutate(conc_temp=Area*(0.266667/IS_area)*(1/RFF)) %>%  ###the conc is now in e-03 mol/l    
  mutate(conc_FA= 75*conc_temp/Subsample_weight)  # The unit is now in nmol  g-1 DW soil

# Exporting this sheet
library(xlsx)
write.xlsx(TNP_Concentrations,"TNP_Concentrations.xlsx")

###### Creating proportions#########
#Creating a table with the different biomarker groups
func_groups <- 
  TNP_Concentrations %>% 
  group_by(Sample, Biom2Soil) %>%
  summarise(sum_groups=sum(conc_FA)) %>% 
  spread(key=Biom2Soil, value = sum_groups) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate (sum_all2 = bacteria + Actinomycetales + gram.neg + Fungi + AMF) %>%   ## Gram.pos is not included here as it was not found in samples
  mutate(Actin_per2 = 100*Actinomycetales/sum_all2, 
         gram.neg_per2 = 100*gram.neg/sum_all2, 
         Fungi_per2 = 100*Fungi/sum_all2,
         bacteria_per2 = 100*bacteria/sum_all2,
         AMF_per2 = 100*AMF/sum_all2,
         fb=Fungi/bacteria+gram.neg) %>%
  add_column(method="PLFA2")

#Here are the final results 
TNP_Proportions<- func_groups %>% select(Sample, bacteria_per2, Actin_per2, gram.neg_per2, Fungi_per2, AMF_per2, sum_all2) 

# Exporting this sheet
library(xlsx)
write.csv(TNP_Proportions,"TNP_Proportions.csv")


##### Graph Total Fungi (%) by heather density ########

Fungi_graph<-merge(TNP_Proportions, plfameta, by.x="Sample",by.y =  "Sample_No")

Fungiplot <- ggplot(Fungi_graph,aes(x=Heather_Density,
                                    y=Fungi_per2))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab(~italic(C.vulgaris)*" density (%)")+
  ylab(expression("Proportional abundance of PLFA fungi (%)"))+
  theme_bw()+
  theme(text = element_text(size = 15))

FUNGI <- Fungiplot + theme_bw() + 
  theme(legend.position = c(0.2, 0.9)) +
  theme(axis.text.x= element_text(colour="black", size=14)) +
  theme(axis.text.y= element_text(colour="black", size=14)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15))

FUNGI

##### LM Proportional abundance of Fungi  ~ heather density ########

# Simple Linear model
model<-lm(Fungi_per2 ~ Heather_Density,
          data=Fungi_graph)

# anova table 
anova(model)
summary(model)

# Checking assumptions
plot(model, 1)

model.diag.metrics <- augment(model)
head(model.diag.metrics)

ggplot(model.diag.metrics, aes(Heather_Density, Fungi_per2)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = Heather_Density, yend = .fitted), color = "red", size = 0.3)

par(mfrow = c(2, 2))
plot(model)


######## Log transforming Fungi data ######


First.model <- lm(Fungi_per2 ~ Heather_Density,
                  data=Fungi_graph)

Transformed.model <- lm(log(Fungi_per2) ~ Heather_Density, data=Fungi_graph)

# anova table 
anova(Transformed.model)
summary(Transformed.model)

# Checking assumptions
plot(Transformed.model, 3)

par(mfrow = c(2, 2))
plot(Transformed.model)

####################fungi - lmer to account for study design#######
library(lme4)
library(lmerTest)

lmer1 <-lmer(log(Fungi_per2) ~ Heather_Density + (1|Geog_Area), data=Fungi_graph)
lmer2 <-lmer(log(Fungi_per2) ~ Heather_Density + (1|Region), data=Fungi_graph)

AIC (lmer1,lmer2)
summary(lmer1)
anova(lmer1)

#### NLFA DATA ###########

#### Create new sheets with the new units, this is NLFA TNP DATA 

library(readxl)
GCMS_NLFA <- read_excel("TNP_Manuscript_Data.xlsx", 
                        sheet = "GCMS NLFA")
View(GCMS_NLFA)

# Import the FAME biomarker sheet 

library(readr)
FAME_biomarker <- read_csv("FAME_RFF_biomarker.csv")
View(FAME_biomarker)

# Merge the FAME sheet with the TNP PLFA information sheet

TNP_NLFA <-merge(GCMS_NLFA, FAME_biomarker, by="Name")

# Import the sheet with the Information about each of the TNP sites, e.g heather density

library(readxl)
NLmeta <- read_excel("TNP_Manuscript_Data.xlsx", 
                     sheet = "NLFA_META")
View(NLmeta)


# Merge the TNP site information with both the PLFA data and the FAME data
TNP_NLFA_data <-merge(TNP_NLFA, NLmeta, by.x="Sample",by.y =  "Sample_No")

# Use Fame Sheet to correct and gain the concentration of each of the lipids

#Calculate correction factor for each sample
correction <- TNP_NLFA_data %>%
  select(c(Sample,Area, Name))%>%
  filter(Name=="19:0")%>%
  mutate(IS_area=1*Area) %>%  

head(correction)  

#Calculate the concentration of each FAME
TNP_NLFA<-  left_join(TNP_NLFA_data, correction, by ="Sample")%>%
  mutate(conc_temp=Area*(0.266667/IS_area)*(1/RFF)) %>%  ###the conc is now in e-03 mol/l    
  mutate(conc_FA= 75*conc_temp/Subsample_weight)   # The unit is now in nmol  g-1 DW soil

# Exporting this sheet
library(xlsx)
write.xlsx(TNP_NLFA,"TNP_NLFA.xlsx")

############ NLFA AMF graph ######################

# Importing the TNP_NLFA_SF sheet, with only the AMF lipid in it (16:1ω5).

library(readxl)
NLFA_AMF <- read_excel("TNP_Manuscript_Data.xlsx", 
                       sheet = "NLFA AMF")
View(NLFA_AMF)


AMFGraph <-ggplot(NLFA_AMF,aes(x=Heather_Density,
                               y=conc_FA))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab(~italic(C.vulgaris)*" density (%)")+
  ylab(expression("Total NLFA AMF (nmol g"^-1* "DW soil)"))+
  theme_bw()+
  theme(text = element_text(size = 15))

# Getting rid of the gridlines on the scatter plot 
AMF <- AMFGraph + theme_bw() + 
  theme(legend.position = c(0.2, 0.9)) +
  theme(axis.text.x= element_text(colour="black", size=14)) +
  theme(axis.text.y= element_text(colour="black", size=14)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15))

AMF


# Simple Linear Regression
model_AMF<-lm(conc_FA ~ Heather_Density,
              data=NLFA_AMF)

#anaylsis of varience 
summary(model_AMF)
anova(model_AMF)

##Lmer 
model_AMFlmer<-lmer(conc_FA ~ Heather_Density + (1|Area),
                  data=NLFA_AMF)

#anaylsis of varience 
summary(model_AMFlmer)
anova(model_AMFlmer)

# Checking assumptions
plot(model_AMF, 1)

model.diag.metrics <- augment(model_AMF)
head(model.diag.metrics)

ggplot(model.diag.metrics, aes(Heather_Density, conc_FA)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = Heather_Density, yend = .fitted), color = "red", size = 0.3)

par(mfrow = c(2, 2))
plot(model_AMF)

### NLFA Binominal graph & stats #####

library(readxl)
Bionominal <- read_excel("TNP_Manuscript_Data.xlsx", 
                         sheet = "Binominal")
View(Bionominal)

set.seed(42)

AMF_Binominal <- ggplot(Bionominal, aes(x=Heather_Density, y=AMF_Presence)) + geom_point()+
  stat_smooth(method = "glm", method.args=list(family="binomial"), se=TRUE)+
  theme_classic() + 
  labs(x="Heather Density (%)",
       y="Soil AMF presence/absence")

AMF_Binominal
ggsave("AMF_Binominal.png", width = 10, height = 8, dpi = 1000)

#The stats 
#Regression

library(glm.predict)
library(fmsb)

glmtest <- glm(AMF_Presence ~ Heather_Density, data=Bionominal, family = "binomial")

summary(glmtest)
NagelkerkeR2(glmtest)



library(sjPlot)
LM.1 <- lm(AMF_Presence ~ Heather_Density, data = Bionominal)

plot_model(LM.1, type = "eff", axis.lim = c(0, 100),
           themes=c("Heather_Density")) + theme_sjplot2()

plot_model(LM.1, type = "diag")

Bionominal$Predicted.Value <-predict(LM.1)
summary(Bionominal$Predicted.Value)

LR.1 <- glm(AMF_Presence ~ Heather_Density, data = Bionominal, family = binomial(link = "logit"))
plot_model(LR.1, type = "eff", 
           terms=c("Heather_Density"), axis.title = "Predicted probability of soil AMF") + theme_sjplot2() +
  theme(text = element_text(size = 15))

summary(LR.1)
anova(LR.1, test = "Chi")

ggsave("BinominalGraph.png", width = 10, height = 8, dpi = 1000)
##### Present and absence 

head(Glm_YesNoSheet)

library(xlsx)
write.xlsx(Glm_YesNoSheet,"Glm_YesNoSheet_PV.xlsx")



##### Microbial Biomass across heather densities ######

Sum_graph<-merge(TNP_Proportions, plfameta, by.x="Sample",by.y =  "Sample_No")

# Exporting this sheet
library(xlsx)
write.xlsx(Sum_graph,"Sum_graph.xlsx")

library(readxl)
Sum_graph <- read_excel("Sum_graph.xlsx")
View(Sum_graph)

ggplot(Sum_graph,aes(x=Heather_Density,
                     y=sum_all2))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab("C.vulgaris Density (%)")+
  ylab(expression("Microbial Biomass (nmol g"^-1* "DW soil)"))+
  theme_bw()+
  theme(text = element_text(size = 15))

# Microbial biomass - Simple Linear Regression 
model<-lm(sum_all2 ~ Heather_Density,
          data=Sum_graph)

# anova table 
anova(model)
summary(model)

####lmer model####
library(lmerTest)
library(lme4)
library(MuMIn)

Lmer3 <- lmer(sum_all2 ~ Heather_Density + (1|Geog_Area), data=Sum_graph) 
Lmer4 <- lmer(sum_all2 ~ Heather_Density + (1|Region), data=Sum_graph) 

AIC(Lmer3,Lmer4)
anova(Lmer4)
summary(Lmer4)
ranova(Lmer4)


ggsave("MicrobialBiomass_vs_heather.png", width = 10, height = 8, dpi = 1000)

# Checking assumptions
plot(model, 1)

model.diag.metrics <- augment(model)
head(model.diag.metrics)

ggplot(model.diag.metrics, aes(Heather_Density, sum_all2)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = Heather_Density, yend = .fitted), color = "red", size = 0.3)

par(mfrow = c(2, 2))
plot(model)

######## Making the multiple scatter plot - Fig 3 #######
# Includes all the FAME data, Fungi, NLFA AMF and microbial biomass ##
# log transforming the data so that they all fit on the same y-axis 

# Import a sheet with all the data we need 

library(readxl)
MultipleRegression <- read_excel("TNP_Manuscript_Data.xlsx", 
                                 sheet = "MultipleRegression")
View(MultipleRegression)

# Exporting this sheet
library(xlsx)
write.xlsx(MultipleRegression,"MultipleRegression.xlsx")

# Change point shapes by the levels of FAME
ggplot(MultipleRegression, aes(x=Heather_Density, y=conc_FA, shape=FAME)) +
  geom_point()
# Change point shapes and colors
ggplot(MultipleRegression, aes(x=Heather_Density, y=conc_FA, shape=FAME, color=FAME)) +
  geom_point()

# Add regression lines
ggplot(MultipleRegression, aes(x=Heather_Density, y=conc_FA, color=FAME, shape=FAME)) +
  geom_point() + 
  geom_smooth(method=lm)

#THIS IS THE CODE to set line types manually
myplot <- ggplot(MultipleRegression, aes(x = Heather_Density, y = conc_FA, color = FAME, shape = FAME)) +
  geom_point() +
  geom_smooth(method = "lm", aes(fill = FAME, linetype = FAME), se = TRUE) +
  xlab(expression(paste(italic("C. vulgaris"), " density (%)")))+
  ylab(expression("FAME (nmol g"^-1* "DW soil)")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_linetype_manual(values = c("AMF" = "dashed", "Microbial Biomass" = "dashed", "Fungi" = "solid")) +
  scale_y_log10() +
  theme(axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.2, 0.9),  # Adjust position to float within the axes bounds
        legend.margin = margin(0.1, 0.1, 0.1, 0.1),  # Set margin to zero
        legend.direction = "vertical",
        legend.box = "horizontal")
# Print the plot
print(myplot)

ggsave("LOGscaleGraphlines.png", width = 10, height = 8, dpi = 1000)

