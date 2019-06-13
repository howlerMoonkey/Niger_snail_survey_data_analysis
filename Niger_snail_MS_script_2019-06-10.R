### Prepare the environment

# Remove all items in environment, for a clean start
rm(list = ls(all=TRUE))  

# Load the necessary packages
library(tidyverse)
library(lubridate)
library(glmmTMB)
library(ggthemes)
library(DHARMa) # for model diagnostics
library(emmeans)  # for post hoc test
library(DataExplorer)   # for data exploration
library(gridExtra)
library(lme4)
library(car)
library(reshape2)
library(ggfortify)
library(boot)
library(aod)
library(mgcv)
library(pscl)

# Change font size and ggthemes globally
theme_set(ggthemes::theme_few(base_size = 8))

# Change globally how numbers are displayed
options("scipen"=100, "digits"=4)

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  

# Import some helpful functions
source("HighstatLibV10.R")
setwd("..")

################################################################################################################################
################################################################################################################################
##Data cleaning and exploratory data analysis- 
##commented out as not part final analysis for manuscript
##merge water chem and survey data

#new1 <- read.csv("Niger_survey_clean_2018-10-13.csv")
#new2 <- read.csv("Niger_param_clean_2018-10-13.csv")
#new3 <- merge(new1, new2, by = "comp_key", all = T)
#write.csv(new3, "survey_para_merge_ALL_2018-10-13.csv")
#new1 <- read.csv("survey_para_merge_ALL_2018-10-13.csv")
#new2 <- read.csv("FTA_all_2018-11-24.csv")
#new3 <- merge(new1, new2, by = "comp_key", all = T)
#write.csv(new3, "survey_FTA_merge.csv")

##METADATA
#newer <- read.csv("survey_coll_headings_wide.csv")
#names(newer)
#newerer <-
#  newer %>%
#  gather(names, values, rec_surv:notes)
#write.csv(newerer, "survey_long.csv")

### EXPLORE DISTRIBUTIONS of continuous variables
###coord_cartesian, xlim, ylim= optimise display
#ggplot(data = survey)+ geom_histogram(mapping = aes(x = Bulinus_tot), binwidth = 0.6) +
#  coord_cartesian((ylim = c(0, 50)))
#ggplot(survey) + geom_count(mapping = aes(x = site_code, y = site_type)) +coord_flip()
#survey %>%
#  count(site_code, site_type_eng) %>%
#  ggplot(mapping = aes(x = site_code, y = site_type_eng)) +
#  geom_tile(mapping = aes(fill = n))

#######
## GENERAL PLOTTING
## trends per species by year, month, by site type
## B forkalii remove outlier (value at 800, Lymnaea and Biomphalaria removed zero values)
#ggplot(data = survey, mapping = aes(x = month_code, y = Bulinus_tot, 
#                                    fill = month_code)) + geom_boxplot() + coord_flip()
#ggplot(data = survey, mapping = aes(x = site_code, y = Bulinus_tot,  
#                                    fill = site_code)) + geom_boxplot()
#ggplot(data = survey, mapping = aes(x = year.code, y = Bulinus_tot)) + geom_boxplot() + coord_flip()
#ggplot(data = survey, mapping = aes(x = site_type, y = L_tot, 
##                                    fill = site_type)) + geom_boxplot() + coord_flip()

###PARAM DATA ONLY
#CLEANING
#para <- read.csv("surv.para_merge_PARA_version_2018-08-13.csv")
#para <-
#  para1 %>%
#  filter(filter_year != "y") ## remove pre july 2011 records
#
#range(para$pH, na.rm = T)
#new <- 
#  para %>%
#  filter(pH < 4 | pH > 14) %>%
#  arrange(pH)
#new$pH ##2.88 3.00 3.00 3.80 3.97

#plot(para$month, para$Temp_Eau) 
#plot(para$month, para$Temp_Air)
#plot(para$month, para$water_depth)
#plot(para$month, para$water_level)
#plot(para$month, para$pH
#plot(para$month, para$Conductivite) 
#plot(para$month, para$PPM)

##### ADD IN WMO DATA
#data <- read.csv("WMO_clean1_2018-10-26.csv")
#data2 <- 
#  data %>%
#  mutate(date_clean2 = dmy(conc))#%>%
# select(-date_orig) ## if wanted to remove org date column
#write.csv(data2, "WMO_clean2_2018-10-26.csv")
#new1 <- read.csv("WMO_clean2_2018-10-26.csv")
#new2 <- read.csv("survey_para_merge_ALL_2018-10-13.csv")
#new3 <- merge(new1, new2, by = "date_val", all.y = T)
#write.csv(new3, "wmo_merge.csv")
################################################################################################################################
################################################################################################################################

##FINAL DATASET FOR ANALYSIS

#' ### Import and prepare the data set
fulldf <- read.csv("Niger_snail_survey_clean_2019_v3.csv") %>% 
  dplyr::select(filter.min, coll_date, month, locality, site_irn, BP_tot, BP_pos_tot, BF_tot, BF_pos_tot, 
                BT_tot, BT_pos_tot, BT_prev, BF_prev, Bulinus_tot, Bulinus_pos_tot,
                bp_pres, bt_pres, bf_pres, visit_no, site_type, L_tot, 
                Temp_Air,Temp_Water, water_speed_ms, water_depth, pH, Cond, PPM, Latitude.v, Longitude.v, Stagnante,
                wmo_min_temp, wmo_max_temp, wmo_av_temp, wmo_prec, seas_wmo, duration, Heure,
                water_level.v, water_level3) %>% 
  # Remove NA values in predictor variables
  dplyr::filter(Bulinus_tot <= 1000,
                !is.na(site_irn),
                site_type != "stream") %>% 
  mutate(coll_date = lubridate::dmy(coll_date), tz = "Africa/Niger",
         duration = as.numeric(as.character(duration)),
         month = as.factor(month),
         year = lubridate::year(coll_date),
         site_irn = as.factor(site_irn), 
         visit_no = as.factor(visit_no),
         bp_pres = as.factor(as.character(bp_pres)),
         bf_pres = as.factor(as.character(bf_pres)),
         Heure = as.factor(as.character(Heure))) %>% 
  # Re-arrange the columns
  dplyr::select(locality, site_irn, visit_no, Bulinus_tot, Bulinus_pos_tot, coll_date, everything()) %>% 
  arrange(locality, site_irn) %>% 
  # This site has no measurements in the wet season, and a lot in the dry season. Later in the models, this leads to 
  # very high variance, when we're looking at interactions containing the season variable.
  # For now, remove this site from the analysis.
  filter(filter.min != "y") %>% 
  # Rescale variables with large values, change water depth unit to meters
  mutate(Cond = scale(Cond, center = 0),
         PPM = scale(PPM, center = 0),
         wmo_prec = scale(wmo_prec, center = 0),
         Temp_Water = scale(Temp_Water, center = 0),
         water_speed_ms = scale(water_speed_ms, center = 0),
         water_depth = water_depth/100,
         water_depth = scale(water_depth, center = 0),
         pH = scale(pH, center = 0)) %>% 
  group_by(site_irn) %>% 
  mutate(Latitude = mean(Latitude.v, na.rm = TRUE),
         Longitude = mean(Longitude.v, na.rm = TRUE))
fulldf$Heure[fulldf$Heure == ""] <- NA

#' Overview over missing values
plot_missing(fulldf) + theme_few()
# For 33.5% of the data, information on the duration of sampling is missing. How is sampling duration distributed?
# To decide how to deal with the missing values for duration, we need to know if the values are missing at complete random 
# (in that case, removing them would not be an issue)

# Is the duration more likely to be missing, given certain other values?
fulldf$dur_missing <- ifelse(is.na(fulldf$duration), 1, 0)
miss_model <- glm(dur_missing ~ Temp_Air + BP_tot + BT_tot + BF_tot +
                    site_type, 
                  family = binomial,
                  data = fulldf)
summary(miss_model) # Not really any strong effects

################################################################################################################################
# Remove NA values in duration, as well as predictors
subdf <- fulldf %>% 
  filter(!is.na(duration))
################################################################################################################################

#' Overview over the data structure
plot_str(subdf) + theme_few()
introduce(subdf)

#' Overview over predictor variables
plot_histogram(subdf)  # Outliers in Cond, PPM, water_depth, wmo_prec, water_speed. These outliers will drive results. Look at them critically

#' Look closer at the outliers

# 1: Cond
plot(subdf$Cond, subdf$Bulinus_tot)

# Where do these high values occur?
subdf$site_irn[subdf$Cond > 6]   
View(subdf[subdf$site_irn == 382877,] %>% dplyr::select(Cond, site_irn, water_depth, pH, PPM, everything()))
View(subdf[subdf$site_irn == 487402,] %>% dplyr::select(Cond, site_irn, water_depth, pH, PPM, everything()))

# High Cond seem to coincide with high PPM, so it's not very likely that these are measurement outliers

# 2: water_speed_ms
plot(subdf$water_speed_ms, subdf$Bulinus_tot)  # One strong outlier

# Where do these high values occur?
subdf$site_irn[subdf$water_speed_ms > 3]   
View(subdf[subdf$site_irn == 382867,] %>% dplyr::select(water_speed_ms, Cond, site_irn, water_depth, pH, PPM, everything()))

# Is it generally an outlier, looking at all site types?
plot(subdf$site_type, subdf$water_speed_ms)   # General outlier

# The measurement of 5 is possible, but a general outlier. It was the single measurement higher than 2, throughout
# the whole study. Even though it is most likely a real measurement, the complete lack of other measurements
# in that range means that it is hard to say that this is an actual representation of counts at high water speeds,
# yet the measurement looks like it's a highly influential outliers, due to its position and high leverage.
# It might be best to exclude this observation from the data.

# 2: wmo_prec
plot(subdf$wmo_prec, subdf$Bulinus_tot)  # Many observations at most extreme precipitation

# Where do these high values occur?
View(subdf[subdf$wmo_prec > 60,] %>% dplyr::select(wmo_prec, seas_wmo, site_irn, water_depth, pH, PPM, everything()))  

# All outliers here come from the same, very rainy day at Lata Kabia (2014-08-02)
ggplot(subdf[subdf$locality == "Lata Kabia",]) +
  geom_point(aes(x = as.factor(month), y = wmo_prec))

#' Make a data frame with outliers removed
subdf_out <- subdf %>% 
  dplyr::filter(water_speed_ms < 2,
                Cond < 8)

#' Are there any variance inflation factors (multicollinearity)? Check using a function from Zuur et al. 2010

pairs(subdf[,c("Temp_Air", "Temp_Water", "water_speed_ms", "water_depth", "pH", "Cond", "PPM",
               "wmo_av_temp", "Bulinus_tot", "wmo_prec")],
      lower.panel = panel.cor)

corvif(data.table::as.data.table(subdf)[, c("Temp_Air", "Temp_Water", "water_speed_ms", "water_depth", "pH", "Cond", "PPM",
                                            "wmo_av_temp", "wmo_prec"), with=FALSE])

# Cond and PPM have high GVIF values (10). For values of higher than 4, only one of the two variables should be used in models, 
# to avoid multicollinearity. 


#' Look at general correlation matrix
plot_correlation(na.omit(subdfc), maxcat = 5L)   


# Study design:
# I have count data of snails per date, counted over many dates at sites, nested in localities.
# So, in each locality the snail counts come from several different sites, repeatedly sampled on different dates.

# Goal:
# Test if snail counts differ between localities, and test influence of environmental factors (e.g. water pH)

# Things to account for:
# A: Sites in localities might show variation in intercepts due to higher initial snail abundance
# B: Sampling duration differed (5-33 minutes), which will most likely influence counts

# How to account for it:
# A: Include site as a random intercept, to account for variation in counts between the sites
# B: Include sampling duration as an offset, to account for differences in sampling effort

#########################################################################################
############################      Bulinus total : To investigate water attribute effects specifically
#########################################################################################

# There are generally many zeroes. This doesn't mean that the data is zero-inflated, but it might be worth checking for zero-inflation anyways

#' ### Make a GLMM

#' Make a maximum model. glmmTMB is a new package by Ben Bolker, that fits models faster, and allows to include arguments
# to account for zero-inflation, if needed.
# Include the sampling duration as an offset. It needs to be specified as log(), since we are using a family distribution
# with log link (nbinom2)

# Decide for a family. It's count count data, looking very overdispersed, so negative binomial is most likely appropriate
Bulinus_poiss <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn) + (1|coll_date) +
                           Temp_Water + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                           locality + site_type + Bulinus_pos_tot + 
                           site_type*Temp_Water + site_type*pH + site_type*Cond + site_type*wmo_prec + 
                           offset(log(duration)),
                         data=subdf,
                         family=poisson)

Bulinus_m <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn) + (1|coll_date) +
                       Temp_Water + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                       locality + site_type + Bulinus_pos_tot + 
                       site_type*Temp_Water + site_type*pH + site_type*Cond + site_type*wmo_prec + 
                       offset(log(duration)),
                     data=subdf,
                     family=nbinom1)

Bulinus_m1 <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn) + (1|coll_date) +
                        Temp_Water + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                        locality + site_type + Bulinus_pos_tot + 
                        site_type*Temp_Water + site_type*pH + site_type*Cond + site_type*wmo_prec + 
                        offset(log(duration)),
                      data=subdf,
                      family=nbinom2)

#' Model diagnostics
# Visually check the model fit using the DHARMa package, a package for model diagnostics
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#formal-goodness-of-fit-tests-on-the-scaled-residuals
sim_residualsPoiss <- DHARMa::simulateResiduals(Bulinus_poiss, 1000)  # Ignore warnings
# Plot the residuals to visually test for over-/underdispersion
plot(sim_residualsPoiss)  # significant deviation

sim_residuals1 <- DHARMa::simulateResiduals(Bulinus_m, 1000)  # Ignore warnings
# Plot the residuals to visually test for over-/underdispersion
plot(sim_residuals1)  # Good fit

sim_residuals2 <- DHARMa::simulateResiduals(Bulinus_m1, 1000) 
plot(sim_residuals2)  # Virtually the same. Difference between nbinom1 and nbinom2 seems neglible

testZeroInflation(sim_residuals2)  # There is no evidence for zero-inflation

# Model results
summary(Bulinus_m1)  # Next to no variance is coming from locality, so technically we could exclude this from the model.
# Aesthetically, we can keep it in, to highlight the nested design of the study
Anova.glmmTMB(Bulinus_m1)

######################################
Bulinus_m1 <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn) + (1|coll_date) +
                        Temp_Water + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                        site_type + Bulinus_pos_tot + month + locality + 
                        site_type*Temp_Water + site_type*pH + site_type*Cond + site_type*wmo_prec + 
                        offset(log(duration)),
                      data=subdf,
                      family=nbinom2)

#########################################################################################
########################      Bulinus truncatus total      ##############################  
#########################################################################################
BT_m1 <- glmmTMB(BT_tot ~ (1|locality/site_irn) + (1|coll_date) + wmo_prec +
                   site_type + month + BT_pos_tot + locality +
                   offset(log(duration)), 
                 data=subdf[subdf$locality != "Gantchi Bassarou",],
                 data=subdf,
                 family=nbinom2)


sim_res_BT_m1 <- DHARMa::simulateResiduals(test, 1000)
plot(sim_res_BT_m1)  # Not overdispersed.
testZeroInflation(sim_res_BT_m1)  # No zero-inflation. 

# Model results
summary(BT_m1)  
Anova.glmmTMB(BT_m1)


#########################################################################################
########################      Bulinus forskalii total      ##############################
#########################################################################################
BF_m1 <- glmmTMB(BF_tot ~ (1|locality/site_irn) + (1|coll_date) + wmo_prec +
                   site_type + month + BF_pos_tot + locality +
                   offset(log(duration)),
                 data=subdf,
                 family=nbinom2)

sim_res_BF_m1 <- DHARMa::simulateResiduals(BF_m1, 1000)
plot(sim_res_BF_m1)  # Not overdispersed.
testZeroInflation(sim_res_BF_m1)  # No zero-inflation. 

# Model results
summary(BF_m1)
Anova.glmmTMB(BF_m1)


#########################################################################################
##############################      Radix total      ###################################
#########################################################################################
subdf %>% group_by(site_type) %>% summarize(n = sum(L_tot))

L_m1 <- glmmTMB(L_tot ~ (1|locality/site_irn) + (1|coll_date) + wmo_prec +
                  site_type + month +  locality +
                  offset(log(duration)),
                data= subdf[subdf$locality != "Koutoukale Zeno" & subdf$locality != "Say" &
                              subdf$locality != "Seberi" & subdf$site_type != "spillway", ],
                family=nbinom2)

sim_res_L_m1 <- DHARMa::simulateResiduals(L_m1, 1000)
plot(sim_res_L_m1)  # Not overdispersed.
testZeroInflation(sim_res_L_m1)  # No zero-inflation. 

# Model results
summary(L_m1)
Anova.glmmTMB(L_m1)


#########################################################################################
##############################      B truncatus positive total      #####################
#########################################################################################
btpos_m1 <- glmmTMB(BT_pos_tot ~ (1|locality/site_irn) + (1|coll_date) + 
                      site_type + locality + month +
                      offset(log(duration)),
                    data= fulldf[fulldf$locality != "Gantchi Bassarou" 
                                 & fulldf$locality != "Tagabati" 
                                 & fulldf$locality != "Tiaguirire" 
                                 & fulldf$locality != "Kohan Garantche" 
                                 & fulldf$locality != "Yoreize Koira",],
                    family=nbinom2)

sim_res_btpos_m1 <- DHARMa::simulateResiduals(btpos_m1, 1000)
plot(sim_res_btpos_m1) 
testZeroInflation(sim_res_btpos_m1)  

# Model results
summary(btpos_m1)
Anova.glmmTMB(btpos_m1)


#########################################################################################
##############################      B. forskalii positive total      ####################
#########################################################################################
## NOT CONVERGING WITH ALL 3 TERMS- RUN SEPARATELY

BFpos_m1 <- glmmTMB(BF_pos_tot ~ (1|locality/site_irn) + (1|coll_date) + 
                      site_type + 
                      offset(log(duration)),
                    data= all, 
                    family=nbinom2)

BFpos_m1 <- glmmTMB(BF_pos_tot ~ (1|locality/site_irn) + (1|coll_date) + 
                      site_type + locality + month +
                      offset(log(duration)),
                    data= fulldf[fulldf$locality != "Tagabati" 
                                 & fulldf$locality != "Tiaguirire" 
                                 & fulldf$locality != "Kohan Garantche" 
                                 & fulldf$locality != "Karma" 
                                 & fulldf$locality != "Bangou Koirey" 
                                 & fulldf$locality != "Diambala" 
                                 & fulldf$locality != "Doguel Kaina" 
                                 & fulldf$locality != "Koutoukale Zeno" 
                                 & fulldf$locality != "Say" 
                                 & fulldf$locality != "Yoreize Koira",],
                    family=nbinom2)


sim_res_L_m1 <- DHARMa::simulateResiduals(BFpos_m1, 1000)
plot(sim_res_L_m1)  # Not overdispersed.
testZeroInflation(sim_res_L_m1)  # No zero-inflation. 

# Model results
summary(BFpos_m1)
Anova.glmmTMB(BFpos_m1)

#########################################################################################
##############################      Biomphalaria positive total      ####################
#########################################################################################

BP_pos_m1 <- glmmTMB(BP_pos_tot ~ (1|locality/site_irn) + (1|coll_date) + 
                       locality +
                       offset(log(duration)),
                     data=subset(fulldf, locality %in% c("Namari Goungou", "Diambala")),
                     family=nbinom2)

sim_res_L_m1 <- DHARMa::simulateResiduals(BP_pos_m1, 1000)
plot(sim_res_L_m1)  # Not overdispersed.
testZeroInflation(sim_res_L_m1)  # No zero-inflation. 

# Model results
summary(BP_pos_m1)
Anova.glmmTMB(BP_pos_m1)



#########################################################################################
###############################      Prevalence      ####################################
#########################################################################################
df_monthly <- subdf %>% group_by(month) %>% 
  mutate(av_prec = mean(wmo_prec)) %>% ungroup() %>% 
  group_by(month, year, locality, site_type, seas_wmo) %>% 
  summarize(BT_pos_tot = as.numeric(sum(BT_pos_tot)),
            BT_tot = as.numeric(sum(BT_tot)),
            BP_pos_tot = as.numeric(sum(BP_pos_tot)),
            BP_tot = as.numeric(sum(BP_tot)),
            BF_pos_tot = as.numeric(sum(BF_pos_tot)),
            BF_tot = as.numeric(sum(BF_tot)),
            av_prec = mean(av_prec))

# test <- df_monthly %>% 
#       gather("species", "value", -c(month, locality, site_type, seas_wmo, av_prec)) %>% 
#       separate(col = "species", into = c("species", "var"), sep = "_") %>% 
#       spread(key = var, value = value)

BT_prev_m1 <- glm(BT_pos_tot/BT_tot ~ site_type + locality + month,
                  weights = BT_tot,
                  data= df_monthly[df_monthly$locality != "Gantchi Bassarou" 
                                   & df_monthly$locality != "Kohan Garantche" 
                                   & df_monthly$locality != "Tagabati" 
                                   & df_monthly$locality != "Tiaguirire" 
                                   & df_monthly$locality != "Yoreize Koira",],
                  family= binomial)


sim_residuals_BT_prev <- DHARMa::simulateResiduals(BT_prev_m1, 1000)  
plot(sim_residuals_BT_prev) 
DHARMa::testDispersion(sim_residuals_BT_prev)
car::Anova(BT_prev_m1)
summary(BT_prev_m1)
DHARMa::testZeroInflation(sim_residuals_BT_prev)

BF_prev_m1 <- glm(BF_pos_tot/BF_tot ~ site_type + locality + month,
                  weights = BF_tot,
                  data= df_monthly[df_monthly$locality != "Tagabati" 
                                   & df_monthly$locality != "Tiaguirire" 
                                   & df_monthly$locality != "Kohan Garantche" 
                                   & df_monthly$locality != "Karma" 
                                   & df_monthly$locality != "Bangou Koirey" 
                                   & df_monthly$locality != "Diambala" 
                                   & df_monthly$locality != "Doguel Kaina" 
                                   & df_monthly$locality != "Koutoukale Zeno" 
                                   & df_monthly$locality != "Say" 
                                   & df_monthly$locality != "Yoreize Koira",],
                  family= binomial)
car::Anova(BF_prev_m1)
summary(BF_prev_m1)


BP_prev_m1 <- glm(BP_pos_tot/BP_tot ~ site_type + locality + month,
                  weights = BP_tot,
                  data=subset(df_monthly, locality %in% c("Namari Goungou", "Diambala")),
                  family= binomial)

sim_residuals_BP_prev <- DHARMa::simulateResiduals(BP_prev_m1, 1000)  
plot(sim_residuals_BP_prev) 
DHARMa::testDispersion(sim_residuals_BP_prev)

car::Anova(BP_prev_m1)
summary(BP_prev_m1)


#########################################################################################
############################      Plot of emmeans      ##################################
#########################################################################################

#' FOR SUBMISSION: site type abundance figure
#' 
site_type_BF  <- data.frame(emmeans(BF_m1,  ~ site_type, type = "response")) %>% mutate(species = "Bulinus forskalii ")
site_type_BT  <- data.frame(emmeans(BT_m1,  ~ site_type, type = "response")) %>% mutate(species = "Bulinus truncatus")

pairs(emmeans(Bulinus_m1,  ~ site_type), type = "response")

p1 <- ggplot(site_type_BF) +
  geom_bar(aes(x = site_type, y = rate), col = "black", fill = "lightgrey",
           position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = site_type, ymin = ifelse((rate-SE) < 0, 0, (rate-SE)), ymax = rate+SE), 
                position = position_dodge(), width = .2) +
  ggtitle("Bulinus forskalii") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank()) +
  ylim(0,13)

p2 <- ggplot(site_type_BT) +
  geom_bar(aes(x = site_type, y = rate), col = "black", fill = "lightgrey",
           position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = site_type, ymin = ifelse((rate-SE) < 0, 0, (rate-SE)), ymax = rate+SE), 
                position = position_dodge(), width = .2) +
  ggtitle("Bulinus truncatus") +
  ylim(0,13) +
  ylab("Mean abundance") +
  theme(axis.title.x=element_blank())

#jpeg(file="Figures/count_site_type.jpg", width = 170, height = 85, units = "mm", res = 600)
grid.arrange(p2, p1, ncol = 2)
#dev.off()

#' FOR SUBMISSION: Locality abundance figure 
#' 
locality_Bulinus  <- data.frame(emmeans(Bulinus_m1,  ~ locality, type = "response")) %>% mutate(species = "Bulinus")
locality_BF  <- data.frame(emmeans(BF_m1,  ~ locality, type = "response")) %>% mutate(species = "Bulinus forskalii ")
locality_BT  <- data.frame(emmeans(BT_m1,  ~ locality, type = "response")) %>% mutate(species = "Bulinus truncatus")
locality_Rad  <- data.frame(emmeans(L_m1,  ~ locality, type = "response")) %>% mutate(species = "Radix natalensis")
temp_Rad <- data.frame(c("Koutoukale Zeno", "Lata Kabia", "Seberi"), c(0,0,0), c(0,0,0), c(0,0,0), c(0,0,0), c(0,0,0), c(0,0,0))
names(temp_Rad) <- names(locality_Rad); locality_Rad <- rbind(temp_Rad, locality_Rad)

temp_BT <- data.frame("Gantchi Bassarou", 0, 0, 0, 0, 0, 0)
names(temp_BT) <- names(locality_BT); locality_BT <- rbind(temp_BT, locality_BT)

locality_Bulinus <- locality_Bulinus %>% arrange(as.character(locality))
locality_BF <- locality_BF %>% arrange(as.character(locality))
locality_BT <- locality_BT %>% arrange(as.character(locality))
locality_Rad <- locality_Rad %>% arrange(as.character(locality))

locality_Bulinus$locality <- as.factor(as.character(locality_Bulinus$locality))
locality_BF$locality <- as.factor(as.character(locality_BF$locality))
locality_BT$locality <- as.factor(as.character(locality_BT$locality))
locality_Rad$locality <- as.factor(as.character(locality_Rad$locality))

p1 <- ggplot(locality_Bulinus) +
  geom_bar(aes(x = locality, y = rate), col = "black", fill = "lightgrey",
           position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = locality, ymin = ifelse((rate-SE) < 0, 0, (rate-SE)), ymax = rate+SE), 
                position = position_dodge(), width = .2) +
  ggtitle("Bulinus") +
  ylab("Mean abundance") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())
p2 <- ggplot(locality_BF) +
  geom_bar(aes(x = locality, y = rate), col = "black", fill = "lightgrey",
           position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = locality, ymin = ifelse((rate-SE) < 0, 0, (rate-SE)), ymax = rate+SE), 
                position = position_dodge(), width = .2) +
  ggtitle("Bulinus forskalii") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) 
p3 <- ggplot(locality_BT) +
  geom_bar(aes(x = locality, y = rate), col = "black", fill = "lightgrey",
           position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = locality, ymin = ifelse((rate-SE) < 0, 0, (rate-SE)), ymax = rate+SE), 
                position = position_dodge(), width = .2) +
  ggtitle("Bulinus truncatus") +
  ylab("Mean abundance") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1))
p4 <- ggplot(locality_Rad) +
  geom_bar(aes(x = locality, y = rate), col = "black", fill = "lightgrey",
           position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = locality, ymin = ifelse((rate-SE) < 0, 0, (rate-SE)), ymax = rate+SE), 
                position = position_dodge(), width = .2) +
  ggtitle("Radix natalensis") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1)) 

jpeg(file="Figures/count_locality.jpg", width = 170, height = 140, units = "mm", res = 600)
grid.arrange(p3, p2, ncol = 2)
dev.off()

#' FOR SUBMISSION: month abundance figure 
#' 
month_BF  <- data.frame(emmeans(BF_m1,  ~ month, type = "response")) %>% mutate(species = "Bulinus forskalii ")
month_BT  <- data.frame(emmeans(BT_m1,  ~ month, type = "response")) %>% mutate(species = "Bulinus truncatus")

p1 <- ggplot(month_BF) +
  geom_bar(aes(x = month, y = rate), col = "black", fill = "lightgrey",
           position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = site_type, ymin = ifelse((rate-SE) < 0, 0, (rate-SE)), ymax = rate+SE), 
                position = position_dodge(), width = .2) +
  ggtitle("Bulinus forskalii") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank()) +
  ylim(0,13)

p2 <- ggplot(month_BT) +
  geom_bar(aes(x = month, y = rate), col = "black", fill = "lightgrey",
           position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = site_type, ymin = ifelse((rate-SE) < 0, 0, (rate-SE)), ymax = rate+SE), 
                position = position_dodge(), width = .2) +
  ggtitle("Bulinus truncatus") +
  ylim(0,13) +
  ylab("Mean abundance") +
  theme(axis.title.x=element_blank())

##############################
######### interactions water attributes plots

#' Temp_Water:site_type
Temp_Water_site_type_df <-  data.frame(emtrends(Bulinus_m1,  ~ site_type, var = "Temp_Water", type = "response"))
ggplot(Temp_Water_site_type_df) +
  geom_errorbar(aes(x = reorder(site_type, Temp_Water.trend), ymin = Temp_Water.trend-SE, ymax = Temp_Water.trend+SE),
                width = .2) +
  geom_point(aes(x = reorder(site_type, Temp_Water.trend), y = Temp_Water.trend), pch = 23, cex = 4,
             fill = "white") +
  ylab("Estimated slope with Temp_Water") +
  xlab("")

ggplot(Temp_Water_site_type_df) +
  geom_bar(aes(x = reorder(site_type, Temp_Water.trend), y = Temp_Water.trend), col = "black",
           fill = "white", position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = reorder(site_type, Temp_Water.trend), ymin = Temp_Water.trend-SE, ymax = Temp_Water.trend+SE),
                width = .2) +
  ylab("Estimated slope with Temp_Water") +
  xlab("")

#' Cond:site_type
Cond_site_type_df <-  data.frame(emtrends(Bulinus_m1,  ~ site_type, var = "Cond", type = "response"))
ggplot(Cond_site_type_df) +
  geom_errorbar(aes(x = reorder(site_type, Cond.trend), ymin = Cond.trend-SE, ymax = Cond.trend+SE),
                width = .2) +
  geom_point(aes(x = reorder(site_type, Cond.trend), y = Cond.trend), pch = 23, cex = 4,
             fill = "white") +
  ylab("Estimated slope with Cond") +
  xlab("")

#' Bulinus_pos_tot
emmip(Bulinus_m1, ~ Bulinus_pos_tot, cov.reduce = range)

#' locality
locality_df <-  data.frame(emmeans(Bulinus_m1,  ~ locality, type = "response"))
locality_s <- emmeans(Bulinus_m1,  ~ locality)
pairs(locality_s, type = "response", adjust = "tukey")
ggplot(locality_df) +
  geom_bar(aes(x = reorder(locality, rate), y = rate), col = "black",
           fill = "white", position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = reorder(locality, rate), ymin = rate-SE, ymax = rate+SE),
                width = .2) +
  ylab("Estimated marginal mean") +
  xlab("")

#' water_speed_ms
emmip(Bulinus_m1, ~ water_speed_ms, cov.reduce = range)


############prevalence plots

#' Prevalence site_type
prev_site_type_BT <-  data.frame(emmeans(BT_prev_m1,  ~ site_type, type = "response")) %>% mutate(species = "Bulinus truncatus")
prev_site_type_BP <-  data.frame(emmeans(BP_prev_m1,  ~ site_type, type = "response")) %>% mutate(species = "Biomphilaria pfeifferi")
prev_site_type_all <- rbind(prev_site_type_BT, prev_site_type_BP)

ggplot(prev_site_type_BT) +
  geom_bar(aes(x = site_type, y = prob, fill = species), col = "black", 
           position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = site_type, ymin = prob-SE, ymax = prob+SE, group = species), position = position_dodge()) 

#' Prevalence month
prev_month_BT <-  data.frame(emmeans(BT_prev_m1,  ~ month, type = "response")) %>% mutate(species = "Bulinus truncatus")
prev_month_BP <-  data.frame(emmeans(BP_prev_m1,  ~ month, type = "response")) %>% mutate(species = "Biomphilaria pfeifferi")
prev_month_all <- rbind(prev_month_BT, prev_month_BP)

ggplot(prev_month_all) +
  geom_bar(aes(x = month, y = prob, fill = species), col = "black", 
           position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = month, ymin = prob-SE, ymax = prob+SE, group = species), position = position_dodge()) 

#' Prevalence locality
prev_locality_BT <-  data.frame(emmeans(BT_prev_m1,  ~ locality, type = "response")) %>% mutate(species = "Bulinus truncatus")
prev_locality_BP <-  data.frame(emmeans(BP_prev_m1,  ~ locality, type = "response")) %>% mutate(species = "Biomphilaria pfeifferi")
prev_locality_all <- rbind(prev_locality_BT, prev_locality_BP)

tiff(file="Figures/prev_locality.tiff", width = 170, height = 150, units = "mm", res = 200)
ggplot(prev_locality_BT) +
  geom_bar(aes(x = reorder(locality, prob), y = prob), col = "black", fill = "lightgrey",
           position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(x = locality, ymin = prob-SE, ymax = prob+SE), position = position_dodge()) +
  ggtitle("Bulinus truncatus")
dev.off()


#########################################################################################
##############################      Anova tables      ###################################
#########################################################################################

glmmTMB::Anova.glmmTMB(BT_m1)   
summary(BT_m1)

glmmTMB::Anova.glmmTMB(BF_m1)  
summary(BF_m1)

bp_m1

btpos_m2
BFpos_m1
BP_pos_m1_month
BP_pos_m1_site
BP_pos_m1_loc


#########################################################################################
############################      Exporting tables      #################################
#########################################################################################


anova_table_glmmTMB <- glmmTMB::Anova.glmmTMB(BP_pos_m1_loc)
write.csv(anova_table_glmmTMB, "BP_pos_m1_month.csv")

summary_table_glmmTMB <- data.frame(summary(BP_pos_m1_month)$coefficients$cond) %>% 
  mutate(coefs = rownames(.)) %>% 
  dplyr::select(coefs, everything())
write.csv(summary_table_glmmTMB, "BP_pos_m1_month_sum.csv")





