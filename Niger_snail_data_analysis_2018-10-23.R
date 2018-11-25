#####################
### 2018-08-03 ##REBOOT ANALYSIS
##cleaned site coll and survey data
##remerge param and survey data
####################

library(tidyverse)
library(reshape2)
library(lubridate)
library(ggfortify)
library(glmmTMB)
library(mgcv)
library(boot)
library(pscl)
library(aod)

setwd("C:\\Users\\murr\\Documents\\R_FILES-2015")
##merge water chem and survey data
new1 <- read.csv("Niger_survey_clean_2018-10-13.csv")
new2 <- read.csv("Niger_param_clean_2018-10-13.csv")
new3 <- merge(new1, new2, by = "comp_key", all = T)
write.csv(new3, "survey_para_merge_ALL_2018-10-13.csv")

new1 <- read.csv("survey_para_merge_ALL_2018-10-13.csv")
new2 <- read.csv("FTA_all_2018-11-24.csv")
new3 <- merge(new1, new2, by = "comp_key", all = T)
write.csv(new3, "survey_FTA_merge.csv")

## Cleanup- removed 7 rows from 2011- pre july- para data only
## filter min and max- odd site records for min, max- at least 2 of other filter =y
### HERE removed lines with para data only- 973 rows. should keep- site visit no snails..
##redo w select
##survey <- read.csv("survey_para_merge_2018-10-13.csv")

##METADATA
newer <- read.csv("survey_coll_headings_wide.csv")
names(newer)
newerer <-
  newer %>%
  gather(rec_surv:notes)
write.csv(newerer, "survey_long.csv")
############## 

######## SUMMARY DATA
#species, month, year, site type, village
## do long and wide version- spread

new <- survey%>%
  group_by(month_code) %>%
  summarise(
    BT_tot = sum(BT_tot, na.rm = T),
    BT_pos = sum(BT_pos, na.rm = T),
    B.for_tot = sum(BF_tot, na.rm = T),
    B.for_pos = sum(BF_pos, na.rm = T),
    BG_tot = sum(BG_tot, na.rm = T),
    BS_tot = sum(BS_tot, na.rm = T),
    BP_tot = sum(BP_tot, na.rm = T),
    BP_pos = sum(BP_pos, na.rm = T),
    L_tot = sum(L_tot, na.rm = T))
#   mean_BT = mean(Bulinus_tot),
#   sd_BT = sd(Bulinus_tot)
write.csv(new, "counts_by_month_code.csv") 

#### do summary- locality- count month code
## error- can't use for factors
new <- newer%>%
  group_by(locality) %>%
  count(newer$site_type_clean, na.rm = T)
    
##PREVALENCE
survey3 <-
  survey %>%
  group_by(month, year) %>%
  summarise(
    BT_tots = sum(BT_tot, na.rm = T),
    BT_pos = sum(BT_pos, na.rm = T))
write.csv(survey3, "test_group_survey_month_year.csv")

##### COUNTS ## no site visits, count by month year etc
surveyor <-
  survey %>%
  count(month_code)
write.csv(surveyor, "month_code.csv")
surveyor <- read.csv("month_code.csv")
pie(surveyor$n, labels=surveyor$month_code)

range(survey$Bulinus_tot, na.rm = T)

#### CLEAN UP THIS BIT PLUS ADD IN SURVEY.NOZERO PLOTS
### EXPLORE DISTRIBUTIONS of continuous variables code fr R for data science Hadley Wickham
###coord_cartesian, xlim, ylim= optimise display

ggplot(data = survey)+ geom_histogram(mapping = aes(x = Bulinus_tot), binwidth = 0.6) +
  coord_cartesian((ylim = c(0, 50)))
ggplot(data = survey)+ geom_histogram(mapping = aes(x = BP_tot), binwidth = 2)
ggplot(data = survey)+ geom_histogram(mapping = aes(x = PPM), binwidth = 5)
ggplot(data = survey)+ geom_histogram(mapping = aes(x = Conductivite), binwidth = 5)
ggplot(data = survey)+ geom_histogram(mapping = aes(x = Temp_Air), binwidth = 0.1)
ggplot(data = survey, mapping = aes(x = Bulinus_tot, colour = site_type_eng)) +
  geom_freqpoly(binwidth = 100)

##pp 99 PLOTTING COUNTS
ggplot(survey) + geom_count(mapping = aes(x = site_code, y = site_type)) +coord_flip()
survey %>%
  count(site_code, site_type_eng) %>%
  ggplot(mapping = aes(x = site_code, y = site_type_eng)) +
  geom_tile(mapping = aes(fill = n))

#######
## GENERAL PLOTTING
## trends per species by year, month, by site type
## B forkalii remove outlier (value at 800, Lymnaea and Biomphalaria removed zero values)
#ggplot(data = survey, mapping = aes(x = passage.code, y = Bulinus_tot)) + geom_boxplot() #year
ggplot(data = survey, mapping = aes(x = month_code, y = Bulinus_tot, 
                                    fill = month_code)) + geom_boxplot() + coord_flip()
ggplot(data = survey, mapping = aes(x = site_code, y = Bulinus_tot,  
                                    fill = site_code)) + geom_boxplot()
ggplot(data = survey, mapping = aes(x = year.code, y = Bulinus_tot)) + geom_boxplot() + coord_flip()
ggplot(data = survey, mapping = aes(x = site_type_eng, y = L_tot, 
                                    fill = site_type_eng)) + geom_boxplot() + coord_flip()

###PARAM DATA ONLY
#CLEANING

para <- read.csv("surv.para_merge_PARA_version_2018-08-13.csv")
para <-
  para1 %>%
  filter(filter_year != "y") ## remove pre july 2011 records

range(para$pH, na.rm = T)
new <- 
  para %>%
  filter(pH < 4 | pH > 14) %>%
  arrange(pH)
new$pH ##2.88 3.00 3.00 3.80 3.97

##plots by month - do by passage, village?
plot(para$Temp_Eau, para$Temp_Air)
plot(para$month, para$Temp_Eau) 
plot(para$month, para$Temp_Air)## still m shape but less clear
plot(para$month, para$water_depth)
plot(para$month, para$water_level) ##regig coding for water level
plot(para$month, para$pH)## fluctuate
plot(para$month, para$Conductivite) ## low in april- or v narrow range, high july
plot(para$month, para$PPM) ## again low in april- less data points fr this month

para2 <-
  para %>%
  filter(month ==  8)%>%
  count(month, year) %>%
  write.csv(surveyor, "survey_month_check.csv")
range(survey$   , na.rm = T)

#####
##### ADD IN WMO DATA
##### Lubridate for fixing dates

data <- read.csv("WMO_clean1_2018-10-26.csv")
data2 <- 
  data %>%
  mutate(date_clean2 = dmy(conc))#%>%
# select(-date_orig) ## if wanted to remove org date column
write.csv(data2, "WMO_clean2_2018-10-26.csv")
new1 <- read.csv("WMO_clean2_2018-10-26.csv")
new2 <- read.csv("survey_para_merge_ALL_2018-10-13.csv")
new3 <- merge(new1, new2, by = "date_val", all.y = T)
write.csv(new3, "wmo_merge.csv")

## what site visits inaccessible- and when
##"flooded"
new2 <- 
  newer %>%
  filter(surv_note_cleaned == "flooded")
new.sum <-
  new2 %>%
  count(site_type_clean)
write.csv(new.sum, "site_count_inacc.csv")

new2 <-
  newer %>%
  filter(BT_tot >= 1)
################
## WMO data- do average precipiation

wmo <- read.csv("wmo_2018-10-30.csv")
wmo.new <-
  newer%>%
  group_by(site_type_code) %>%
  summarise(
    air_temp = mean(Temp_Air, na.rm = T),
    water_temp = mean(Temp_Eau, na.rm = T),
    mean_prec = mean(wmo_prec,  na.rm = T),
#   mean_wmo.min_temp = mean(wmo_min_temp,  na.rm = T), #join by date so not valid
#   mean_wmo.max_temp = mean(wmo_max_temp,  na.rm = T), #sum by month fine as
    mean_pH = mean(pH, na.rm = T),
    mean_PPM = mean(PPM, na.rm = T),
    mean_cond = mean(Conductivite, na.rm = T))
#mean_wat.lev = mean(water_level, na.rm = T)) ## words
range(newer$water_depth, na.rm = T)
range(newer$water_speed_ms,  na.rm = T)
range(newer$water_depth,  na.rm = T)

write.csv(wmo.new, "wmo2.csv")


#######
#######
#### STATS

newer <- read.csv("survey_para_merge_ALL_2018-10-13.csv")

## response variable= Bulinus pres- 0/1
## Bulinus pos 0/1
#Error in eval(expr, envir, enclos) : y values must be 0 <= y <= 1
#changed to Bulinus presence
#BT_pos 
#bf_pos
#Bulinus_pres 
#BT_pres 
#bf_pres

## predictor variables 
#Temp_Air
#Temp_Eau
#water_depth
#WMO_prec
#WMO_max_temp # WMO
#pH
#PPM
#Conductivite
#Stagnante
#water_speed_ms
#snails.present##other snail presense- interaction term?
#month - month as a variable
#w_level_ed
#access
#season_WMO
#seaon
#site_type_code ## merge of para and survey- if na in para add fr para
#canal
#site_code
#site_CAT

## mutate- average niger data temp by date- wmo vals will be same 
## for all sites for given date
## WMO max min av prec- diff results- some pos some neg
# see what site types recorded as dry or flooded

mylogit <- glm(Bul_pres ~ Temp_Air + Temp_Eau + wmo_max_temp + wmo_prec + PPM + Cond,
               data = newerer, family = "binomial")
mylogit <- glm(bf_pos ~ site_code, data = newer, family = "binomial")
mylogit <- glm(Bul_pres ~ season, data = newer, family = "binomial")

summary(mylogit)

## DO FOR POS

confint(mylogit)

wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 1:19)
## WALD AOD LIBRARY, TERMS = terms in model e.g predictor ~ term value 1 + term 2 etc

## odds ratios- only
exp(coef(mylogit))
## odds ratios plus CI
exp(cbind(OR = coef(mylogit), confint(mylogit)))

## do for BT, BF, BP nad LN

### ANOVA

mylogit <- glm(Bul_pres ~ seas_wmo, data = newer, family = "binomial")
autoplot(mylogit)
anova(mylogit)
##FOR TUKEY
mylogit.aov <- aov(mylogit) ##repackage model for tukey test
tukey.out <- TukeyHSD(mylogit.aov)
tukey.out

#OTHER ANOV ## do site type code

anov <- aov(Bul_pres~seas_wmo, data = newer) #runs the ANOVA test.
ls(anov) #lists the items stored by the test.
summary(anov) #give the basic ANOVA output.

### Kruskall Wallis test

#zero inflated negative binomial GLMM - Anouk paper for Bulinus pos count

glm <- glm(Bulinus_pos~seas_wmo, data = newer) 
glm_nb <- glm.nb(Bulinus_pos~seas_wmo, data = newer) 
gam <- gam(Bulinus_pos~seas_wmo, data = newer) 
glmmTMB <- glmmTMB(Bulinus_pos~seas_wmo, data = newer) 
zeroinfla <- zeroinfl(Bulinus_pos~seas_wmo, data = newer) 
zeroid <- update(zeroinfla, . ~ 1)
pchisq(2 * (logLik(zeroinfla) - logLik(zeroid)), df = 3, lower.tail=FALSE)

glm <- glm(Bulinus_pos~Bulinus_tot, data = newer) 

newer <- read.csv("survey_para_merge_ALL_2018-10-13.csv")
mylogit <- glm(Bul_pres ~ Temp_Air + wmo_prec + Cond, data = newerer, family = "binomial")

newerer <-
  newer%>%
  filter(filter.max == "y")  

