# ============================================================================
# GLM Analysis: Phosphorus Relationships in Diatoms
# Comparing Asterionella and Fragilaria responses to nutrient gradients
# ============================================================================

# Research Questions:
# 1. Do these diatoms show different phosphorus optima?
# 2. How does N:P stoichiometry influence their distribution?
# 3. Does light availability modify phosphorus effects?

library(tidyverse)
library(broom)      # For tidy model outputs
library(patchwork)  # For combining plots
library(DHARMa)     # Check model assumptions

# Load data ----
# (Assumes you've already run the EDA script and have these objects)
# phyto <- read.csv("data/raw/nla2022_phytoplanktoncount_wide.csv")
# env_data <- [your environmental data from EDA script, since it combines several csv files]

# PART 1: PREPARE SPECIES-SPECIFIC DATASETS ----------------------------------

# Extract Asterionella and Fragilaria from phytoplankton data
phyto_v1 <- phyto %>%
  filter(VISIT_NO == 1) %>%
  filter(!is.na(BIOVOLUME))

# Get the two focal species
asterionella <- phyto_v1 %>%
  filter(str_detect(TARGET_TAXON, "ASTERIONELLA")) %>%
  group_by(SITE_ID) %>%
  summarize(
    presence = 1,
    biovolume = sum(BIOVOLUME, na.rm = TRUE),
    density = sum(DENSITY, na.rm = TRUE),
    .groups = "drop"
  )

fragilaria <- phyto_v1 %>%
  filter(str_detect(TARGET_TAXON, "FRAGI")) %>%
  group_by(SITE_ID) %>%
  summarize(
    presence = 1,
    biovolume = sum(BIOVOLUME, na.rm = TRUE),
    density = sum(DENSITY, na.rm = TRUE),
    .groups = "drop"
  )

# How many sites have each species?
nrow(asterionella)
nrow(fragilaria)

# Join with environmental data and create presence/absence
# Include all sites (0 for absence)
modeling_data <- env_data %>%
  select(SITE_ID, 
         # Nutrients
         NTL_RESULT, PTL_RESULT, CHLA_RESULT,
         # Light
         DOC_RESULT, SECCHI_DEPTH, TURB_RESULT,
         # Water chemistry
         PH_RESULT, COND_RESULT,
         # Physical
         ELEVATION, AREA_HA, INDEX_SITE_DEPTH) %>%
  # Add Asterionella
  left_join(asterionella %>% 
              select(SITE_ID, aster_presence = presence, 
                     aster_biovolume = biovolume),
            by = "SITE_ID") %>%
  # Add Fragilaria
  left_join(fragilaria %>% 
              select(SITE_ID, frag_presence = presence,
                     frag_biovolume = biovolume),
            by = "SITE_ID") %>%
  # Convert NAs to 0 for absence
  mutate(
    aster_presence = replace_na(aster_presence, 0),
    frag_presence = replace_na(frag_presence, 0),
    aster_biovolume = replace_na(aster_biovolume, 0),
    frag_biovolume = replace_na(frag_biovolume, 0),
    # Create N:P ratio (molar) - convert from mass
    # N (mg/L) / P (ug/L) -> need to adjust units
    # Molar N:P = (NTL_mg/L / 14) / (PTL_ug/L / 31000)
    NP_ratio = (NTL_RESULT / 14) / (PTL_RESULT / 31000),
    # Log-transform nutrients for modeling (avoid log(0))
    log_PTL = log10(PTL_RESULT + 1),
    log_NTL = log10(NTL_RESULT + 1),
    log_NP = log10(NP_ratio + 1)
  ) %>%
  # Remove sites with missing nutrient data
  filter(!is.na(PTL_RESULT), !is.na(NTL_RESULT))

# Check prevalence
cat("\nPrevalence:\n")
cat("Asterionella:", sum(modeling_data$aster_presence), "/", 
    nrow(modeling_data), "sites =", 
    round(100*mean(modeling_data$aster_presence), 1), "%\n")
cat("Fragilaria:", sum(modeling_data$frag_presence), "/", 
    nrow(modeling_data), "sites =",
    round(100*mean(modeling_data$frag_presence), 1), "%\n")


# PART 2: PRESENCE/ABSENCE MODELS ---------------------------------------------

# Theory: If phosphorus is a key limiting nutrient, we expect to see:
# - Increased probability of occurrence with higher P availability
# - Potential differences between species in P optima
# - Possible non-linear relationships if high P favors competitors

## 2.1 Asterionella presence/absence ----

# Model 1: Phosphorus only
aster_m1 <- glm(aster_presence ~ log_PTL,
                family = binomial(link = "logit"),
                data = modeling_data)

summary(aster_m1)
# weird? Ask Mathew why this result

coef(aster_m1)["log_PTL"] # Negative and signif: aster is less likely at high P
exp(coef(aster_m1)["log_PTL"]) # P increase 10 fold = odds of presence multiplied by 0.27
# That's saying 73% decrease in odds
# Percent decrease in odds:
(1 - exp(coef(aster_m1)["log_PTL"])) * 100

# Pseudo-R²
1 - (aster_m1$deviance / aster_m1$null.deviance)
# P only explains about 5% variation

# Odds ratio of 0.27 means a 10-fold increase in phosphorus reduces the odds of finding Asterionella by 73%
predict(aster_m1, 
        newdata = data.frame(log_PTL = log10(c(10, 100) + 1)),
        type = "response")
# When P increases 10-fold, probability drops from 21% to 7%


# Check model assumptions using DHARMa
aster_m1_resid <- simulateResiduals(aster_m1, n = 1000)
plot(aster_m1_resid)

# Key diagnostic tests
testDispersion(aster_m1_resid)  # Check for over/underdispersion
testZeroInflation(aster_m1_resid)  # Check for zero-inflation


# Model 2: Phosphorus + Stoichiometry
aster_m2 <- glm(aster_presence ~ log_PTL + log_NP,
                family = binomial(link = "logit"),
                data = modeling_data)

summary(aster_m2)

aster_m2_resid <- simulateResiduals(aster_m2, n = 1000)
plot(aster_m2_resid)
testDispersion(aster_m2_resid)
testOutliers(aster_m2_resid)
testZeroInflation(aster_m2_resid)


# MODEL COMPARISON

# Approach 1: AIC comparison
# Lower AIC = better model (balances fit vs. complexity)
AIC(aster_m1, aster_m2)


# Approach 2: Likelihood Ratio Test (for nested models only)
anova(aster_m1, aster_m2, test = "Chisq")

#Interpretation: Significant p-value means the more complex model
# provides a significantly better fit than the simpler model)


## 2.2 Fragilaria presence/absence ----

# Repeat this section for the other species


# PART 3: BIOVOLUME MODELS (WHERE PRESENT) -----------------------------------

# Theory: Among sites where a species occurs, does biovolume (a proxy for 
# population size) increase with P availability? This tests whether P 
# influences not just occurrence but also population density.

# Subset to sites where each species is present
aster_present <- modeling_data %>%
  filter(aster_presence == 1, aster_biovolume > 0)
nrow(aster_present)

frag_present <- modeling_data %>%
  filter(frag_presence == 1, frag_biovolume > 0)
nrow(frag_present)


## 3.1 Asterionella biovolume ----

# Log-transform biovolume (right-skewed distribution)
hist(aster_present$aster_biovolume)
hist(log10(aster_present$aster_biovolume))

aster_present <- aster_present %>%
  mutate(log_biovolume = log10(aster_biovolume + 1))

# Model 1: Phosphorus only
aster_bv_m1 <- glm(log_biovolume ~ log_PTL,
                   family = gaussian,
                   # family = Gamma(link = "log"), # both are not great
                   data = aster_present)

summary(aster_bv_m1)

# Checks
aster_bv_m1_resid <- simulateResiduals(aster_bv_m1, n = 1000)
plot(aster_bv_m1_resid)

testDispersion(aster_bv_m1_resid)  
testOutliers(aster_bv_m1_resid)   
testZeroInflation(aster_bv_m1_resid)

# Model 2: Phosphorus + Stoichiometry  
aster_bv_m2 <- glm(log_biovolume ~ log_PTL + log_NP,
                   family = gaussian,
                   data = aster_present)

summary(aster_bv_m2)

# Checks
aster_bv_m2_resid <- simulateResiduals(aster_bv_m2, n = 1000)
plot(aster_bv_m2_resid)

testDispersion(aster_bv_m2_resid)  
testOutliers(aster_bv_m2_resid)   
testZeroInflation(aster_bv_m2_resid)


## 3.2 Fragilaria biovolume ----

## Your turn, repeat with other species

# PART 4: MODEL SUMMARIES -----------------------------------------------------

# use summary() to see all the models
tidy(aster_m2) %>% mutate(Species = "Asterionella", Model = "Presence_M2")


# PART 5: VISUALIZE RESULTS ---------------------------------------------------

## 5.1 Presence/Absence vs Phosphorus ----

# Create prediction data
pred_data <- data.frame(
  log_PTL = seq(min(modeling_data$log_PTL, na.rm = TRUE),
                max(modeling_data$log_PTL, na.rm = TRUE),
                length.out = 100),
  log_NP = median(modeling_data$log_NP, na.rm = TRUE)
)

# Predictions for Asterionella
pred_data$aster_pred <- predict(aster_m2, newdata = pred_data, 
                                type = "response")

# Back-transform PTL for plotting
pred_data$PTL <- 10^pred_data$log_PTL - 1

# Plot
p_presence <- ggplot() +
  # Asterionella
  geom_point(data = modeling_data, 
             aes(x = PTL_RESULT, y = aster_presence),
             alpha = 0.2, color = "blue", position = position_jitter(height = 0.02)) +
  geom_line(data = pred_data,
            aes(x = PTL, y = aster_pred),
            color = "blue", linewidth = 1) +
  scale_x_log10() +
  labs(title = "Probability of Occurrence vs Total Phosphorus",
       x = "Total Phosphorus (μg/L, log scale)",
       y = "Probability of Presence") +
  theme_minimal()

p_presence


# N:P effect
pred_data_NP <- data.frame(
  log_PTL = median(modeling_data$log_PTL, na.rm = TRUE),  # Hold P constant
  log_NP = seq(min(modeling_data$log_NP, na.rm = TRUE),
               max(modeling_data$log_NP, na.rm = TRUE),
               length.out = 100)
)

# Predictions for Asterionella
pred_data_NP$aster_pred <- predict(aster_m2, newdata = pred_data_NP, 
                                type = "response")

# Back-transform N:P to original scale
# We created: log_NP = log10(NP_ratio + 1)
# So reverse: NP_ratio = 10^log_NP - 1
pred_data_NP$NP_ratio <- 10^pred_data_NP$log_NP - 1

# Plot
np_presence <- ggplot() +
  # Asterionella
  geom_point(data = modeling_data, 
             aes(x = NP_ratio, y = aster_presence),
             alpha = 0.2, color = "blue", position = position_jitter(height = 0.02)) +
  geom_line(data = pred_data_NP,
            aes(x = NP_ratio, y = aster_pred),
            color = "blue", linewidth = 1) +
  scale_x_log10() +
  labs(title = "Asterionella: Probability of Occurrence vs N:P Ratio",
       x = "N:P Molar Ratio (log scale)",
       y = "Probability of Presence") +
  theme_minimal()

np_presence
