# ============================================================================
# NLA 2022 Comprehensive Data Exploration
# Advanced Ecology Lab 1 - Systematic EDA Framework
# ============================================================================

# Load packages ----
library(tidyverse)
library(sf)           
library(corrplot)     
library(GGally)       
library(patchwork)    


# PART 1: LOAD ALL DATA  --------------------------------------------------

# Site information
sites <- read.csv("data/raw/nla22_siteinfo.csv")

# Basins (spatial watersheds)
basins <- read_sf("data/raw/nla2022_basins/NLA2022_basins.shp")

# Landscape data (watershed characteristics)
landscape <- read.csv("data/raw/nla2022_landscape_wide_0.csv")

# Water chemistry
chem <- read.csv("data/raw/nla22_waterchem_wide.csv")

# Secchi (water clarity)
secchi <- read.csv("data/raw/nla22_secchi.csv")

# Physical habitat
habitat <- read.csv("data/raw/nla2022_phab_wide.csv")

# Biological data
zoop <- read.csv("data/raw/nla22_zooplanktoncount_wide.csv")

benthos <- read.csv("data/raw/nla22_benthic_counts.csv")

phyto <- read.csv("data/raw/nla2022_phytoplanktoncount_wide.csv")


# PART 2: UNDERSTAND DATA STRUCTURE ---------------------------------------


# Check dimensions
dim(sites)
dim(basins)
dim(landscape)

# Biological data - how many sites and taxa?
n_distinct(zoop$SITE_ID)
n_distinct(zoop$TARGET_TAXON)

n_distinct(benthos$SITE_ID)
n_distinct(benthos$TARGET_TAXON)

n_distinct(phyto$SITE_ID)
n_distinct(phyto$TARGET_TAXON)

# PART 3: CHECK DATA OVERLAP (SITE MATCHING) ------------------------------


# Get unique sites from each dataset
sites_main <- unique(sites$SITE_ID)
sites_basin <- unique(basins$SITE_ID)
sites_landscape <- unique(landscape$SITE_ID)
sites_chem <- unique(chem$SITE_ID)
sites_zoop <- unique(zoop$SITE_ID)
sites_benthos <- unique(benthos$SITE_ID)
sites_phyto <- unique(phyto$SITE_ID)

# Check overlaps
sum(sites_main %in% sites_basin)
sum(sites_main %in% sites_landscape)
sum(sites_main %in% sites_chem)
sum(sites_main %in% sites_zoop)
sum(sites_main %in% sites_benthos)
sum(sites_main %in% sites_phyto)

# Sites with ALL THREE taxonomic groups
sites_all_bio <- Reduce(intersect, list(sites_zoop, sites_benthos, sites_phyto))
length(sites_all_bio)


# PART 4: EXPLORE SITE CHARACTERISTICS ------------------------------------


# Lake size distribution
summary(sites$AREA_HA)

p_size <- ggplot(sites, aes(x = AREA_HA)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  scale_x_log10() +
  labs(title = "Distribution of Lake Sizes",
       x = "Lake Area (ha, log scale)",
       y = "Count") +
  theme_minimal()

p_size

# Lake origin
table(sites$LAKE_ORGN)

# Regional distribution
table(sites$AG_ECO9_NM)

# Elevation
summary(sites$ELEVATION)


# PART 5: PROCESS WATER CHEMISTRY DATA ------------------------------------


# Key nutrients for lake ecology
# NTL = Total Nitrogen, PTL = Total Phosphorus (primary limiting nutrients)
# CHLA = Chlorophyll-a (phytoplankton biomass proxy)
# Nitrate/Ammonia = Nitrogen forms available to algae
# pH, Conductivity = Water chemistry baseline
# DOC = Dissolved Organic Carbon (carbon availability, light attenuation)
# Turbidity = Water clarity/suspended particles
# ANC = Acid Neutralizing Capacity (buffering capacity)
key_nutrients <- c("NTL_RESULT", "PTL_RESULT", "CHLA_RESULT", 
                   "NITRATE_NITRITE_N_RESULT", "AMMONIA_N_RESULT",
                   "PH_RESULT", "COND_RESULT", "DOC_RESULT",
                   "TURB_RESULT", "ANC_RESULT")

# Create simplified chemistry dataset (first visit only)
chem_simple <- chem %>%
  filter(VISIT_NO == 1 | is.na(VISIT_NO)) %>%
  select(SITE_ID, any_of(key_nutrients)) %>%
  group_by(SITE_ID) %>%
  summarize(across(everything(), ~mean(.x, na.rm = TRUE)))

summary(chem_simple)


# PART 6: PROCESS SECCHI DATA (WATER CLARITY) -----------------------------


# Calculate average Secchi depth
secchi_processed <- secchi %>%
  filter(VISIT_NO == 1 | is.na(VISIT_NO)) %>%
  mutate(SECCHI_DEPTH = (DISAPPEARS + REAPPEARS) / 2) %>%
  select(SITE_ID, SECCHI_DEPTH, CLEAR_TO_BOTTOM, INDEX_SITE_DEPTH)

summary(secchi_processed$SECCHI_DEPTH)
table(secchi_processed$CLEAR_TO_BOTTOM)


# PART 7: CREATE MASTER ENVIRONMENTAL DATASET -----------------------------

env_data <- sites %>%
  select(SITE_ID, UNIQUE_ID, PSTL_CODE, 
         LAT_DD83, LON_DD83, ELEVATION,
         AREA_HA, LAKE_ORGN,
         AG_ECO9, AG_ECO9_NM,
         US_L3NAME,
         URBN_NLA22) %>%
  # Add basin area
  left_join(basins %>% 
              st_drop_geometry() %>% 
              select(SITE_ID, WSAREASQKM = AreaSqKm), 
            by = "SITE_ID") %>%
  # Add water chemistry
  left_join(chem_simple, by = "SITE_ID") %>%
  # Add Secchi
  left_join(secchi_processed, by = "SITE_ID") %>%
  # Add key landscape variables
  # Land use (crop, urban, forest), climate (temp, precip), and nutrient inputs (fertilizer, manure)
  left_join(landscape %>% 
              select(SITE_ID, 
                     PCTCROP2021, PCTURBHI2021, PCTDECID2021,
                     TMEAN9120, PRECIP9120,
                     FERT, MANURE), 
            by = "SITE_ID")

dim(env_data)

# How many sites have complete nutrient data?
complete_sites <- env_data %>%
  filter(complete.cases(NTL_RESULT, PTL_RESULT, CHLA_RESULT))

nrow(complete_sites)


# PART 8: EXPLORE ENVIRONMENTAL GRADIENTS ---------------------------------


# Key environmental variables correlation
env_numeric <- env_data %>%
  select(AREA_HA, WSAREASQKM, ELEVATION,
         NTL_RESULT, PTL_RESULT, CHLA_RESULT,
         PH_RESULT, COND_RESULT, DOC_RESULT,
         SECCHI_DEPTH,
         PCTCROP2021, PCTURBHI2021,
         TMEAN9120, PRECIP9120) %>%
  mutate(AREA_HA_log = log10(AREA_HA + 1),
         WSAREASQKM_log = log10(WSAREASQKM + 1)) %>%
  select(-AREA_HA, -WSAREASQKM)

# Correlation matrix
cor_matrix <- cor(env_numeric, use = "pairwise.complete.obs")

# Plot correlation matrix
corrplot(cor_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45,
         title = "Environmental Variable Correlations",
         mar = c(0,0,2,0))


# Nutrient distributions
p_nutrients <- env_data %>%
  select(NTL_RESULT, PTL_RESULT, CHLA_RESULT) %>%
  pivot_longer(everything(), names_to = "nutrient", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50, fill = "forestgreen", alpha = 0.7) +
  scale_x_log10() +
  facet_wrap(~nutrient, scales = "free") +
  labs(title = "Distribution of Key Nutrients",
       x = "Concentration (log scale",
       y = "Count") +
  theme_minimal()

p_nutrients


#  PART 9: EXPLORE ZOOPLANKTON COMMUNITY ----------------------------------



# Remove records with missing density values (NA = measurement issue, not absence)
zoop_v1 <- zoop %>%
  filter(VISIT_NO == 1) %>%
  filter(!is.na(DENSITY))

# Species richness per site
zoop_richness <- zoop_v1 %>%
  group_by(SITE_ID) %>%
  summarize(richness = n_distinct(TARGET_TAXON),
            total_density = sum(DENSITY, na.rm = TRUE))

summary(zoop_richness$richness)

# Most common taxa
zoop_prevalence <- zoop_v1 %>%
  group_by(TARGET_TAXON) %>%
  summarize(n_sites = n_distinct(SITE_ID),
            mean_density = mean(DENSITY, na.rm = TRUE)) %>%
  arrange(desc(n_sites)) %>%
  head(10)

# Body size patterns (for cladocerans)
if("CLADOCERA_SIZE" %in% names(zoop_v1)) {
  table(zoop_v1$CLADOCERA_SIZE, useNA = "ifany")
}

# Functional groups
if("FFG" %in% names(zoop_v1)) {
  table(zoop_v1$FFG, useNA = "ifany")
}

# Visualization
p_zoo_rich <- ggplot(zoop_richness, aes(x = richness)) +
  geom_histogram(bins = 30, fill = "coral", alpha = 0.7) +
  labs(title = "Zooplankton Species Richness",
       x = "Number of Taxa",
       y = "Number of Sites") +
  theme_minimal()

p_zoo_rich


# PART 10: EXPLORE BENTHIC MACROINVERTEBRATE COMMUNITY --------------------



# Filter to first visit
benthos_v1 <- benthos %>%
  filter(VISIT_NO == 1) %>%
  filter(!is.na(TOTAL))

# Species richness per site
benthos_richness <- benthos_v1 %>%
  group_by(SITE_ID) %>%
  summarize(richness = n_distinct(TARGET_TAXON),
            total_count = sum(TOTAL, na.rm = TRUE))

summary(benthos_richness$richness)


# Most common taxa
benthos_prevalence <- benthos_v1 %>%
  group_by(TARGET_TAXON) %>%
  summarize(n_sites = n_distinct(SITE_ID),
            mean_count = mean(TOTAL, na.rm = TRUE)) %>%
  arrange(desc(n_sites)) %>%
  head(10)

benthos_prevalence

# Taxonomic composition
benthos_orders <- benthos_v1 %>%
  group_by(ORDER) %>%
  summarize(n_taxa = n_distinct(TARGET_TAXON),
            n_sites = n_distinct(SITE_ID)) %>%
  arrange(desc(n_sites)) %>%
  head(10)

benthos_orders

# Visualization
p_bent_rich <- ggplot(benthos_richness, aes(x = richness)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
  labs(title = "Benthic Macroinvertebrate Richness",
       x = "Number of Taxa",
       y = "Number of Sites") +
  theme_minimal()
p_bent_rich


# PART 11: EXPLORE PHYTOPLANKTON COMMUNITY --------------------------------


# Filter to first visit
phyto_v1 <- phyto %>%
  filter(VISIT_NO == 1) %>%
  filter(!is.na(BIOVOLUME))

# Species richness per site
phyto_richness <- phyto_v1 %>%
  group_by(SITE_ID) %>%
  summarize(richness = n_distinct(TARGET_TAXON),
            total_density = sum(DENSITY, na.rm = TRUE),
            total_biovolume = sum(BIOVOLUME, na.rm = TRUE))

summary(phyto_richness$richness)
summary(phyto_richness$total_biovolume)

# Algal group composition
phyto_groups <- phyto_v1 %>%
  group_by(ALGAL_GROUP) %>%
  summarize(n_taxa = n_distinct(TARGET_TAXON),
            n_sites = n_distinct(SITE_ID),
            mean_biovolume = mean(BIOVOLUME, na.rm = TRUE)) %>%
  arrange(desc(n_sites))

# Most common taxa
phyto_prevalence <- phyto_v1 %>%
  group_by(TARGET_TAXON, ALGAL_GROUP) %>%
  summarize(n_sites = n_distinct(SITE_ID),
            mean_biovolume = mean(BIOVOLUME, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(n_sites)) %>%
  head(10)

# Visualizations
p_phyto_rich <- ggplot(phyto_richness, aes(x = richness)) +
  geom_histogram(bins = 30, fill = "darkblue", alpha = 0.7) +
  labs(title = "Phytoplankton Species Richness",
       x = "Number of Taxa",
       y = "Number of Sites") +
  theme_minimal()

p_phyto_rich

p_phyto_groups <- ggplot(phyto_groups, aes(x = reorder(ALGAL_GROUP, n_sites), 
                                            y = n_sites)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(title = "Phytoplankton Prevalence by Algal Group",
       x = "Algal Group",
       y = "Number of Sites") +
  theme_minimal()

p_phyto_groups


# PART 12: COMPARE RICHNESS ACROSS TROPHIC LEVELS -------------------------



# Combine richness data
richness_combined <- zoop_richness %>%
  select(SITE_ID, zoo_richness = richness) %>%
  full_join(benthos_richness %>% 
              select(SITE_ID, benthos_richness = richness),
            by = "SITE_ID") %>%
  full_join(phyto_richness %>% 
              select(SITE_ID, phyto_richness = richness),
            by = "SITE_ID")

# Sites with all three
richness_complete <- richness_combined %>%
  filter(complete.cases(zoo_richness, benthos_richness, phyto_richness))

nrow(richness_complete)

# Correlation between groups
if(nrow(richness_complete) > 0) {
  cor_richness <- cor(richness_complete %>% select(-SITE_ID), 
                      use = "complete.obs")
  round(cor_richness, 3)
}

# Visualize richness relationships
p_rich_compare <- richness_complete %>%
  pivot_longer(cols = c(zoo_richness, benthos_richness, phyto_richness),
               names_to = "group",
               values_to = "richness") %>%
  ggplot(aes(x = group, y = richness, fill = group)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c("coral", "darkgreen", "darkblue")) +
  labs(title = "Species Richness Across Trophic Levels",
       x = "",
       y = "Number of Taxa") +
  theme_minimal() +
  theme(legend.position = "none")

p_rich_compare


# PART 13: LINK COMMUNITIES TO ENVIRONMENT --------------------------------



# Combine richness with environmental data
zoo_env <- zoop_richness %>%
  left_join(env_data, by = "SITE_ID") %>%
  filter(complete.cases(richness, NTL_RESULT, PTL_RESULT))

benthos_env <- benthos_richness %>%
  left_join(env_data, by = "SITE_ID") %>%
  filter(complete.cases(richness, NTL_RESULT, PTL_RESULT))

phyto_env <- phyto_richness %>%
  left_join(env_data, by = "SITE_ID") %>%
  filter(complete.cases(richness, NTL_RESULT, PTL_RESULT))

# Richness vs nutrients

# Zooplankton
p_zoo_nutrients <- zoo_env %>%
  ggplot(aes(x = PTL_RESULT, y = richness)) +
  geom_point(alpha = 0.5, color = "coral") +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_log10() +
  labs(title = "Zooplankton Richness vs Total Phosphorus",
       x = "Total Phosphorus (μg/L, log scale",
       y = "Zooplankton Richness") +
  theme_minimal()

p_zoo_nutrients

# Phytoplankton
p_phyto_nutrients <- phyto_env %>%
  ggplot(aes(x = PTL_RESULT, y = richness)) +
  geom_point(alpha = 0.5, color = "darkblue") +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_log10() +
  labs(title = "Phytoplankton Richness vs Total Phosphorus",
       x = "Total Phosphorus (μg/L, log scale)",
       y = "Phytoplankton Richness") +
  theme_minimal()

p_phyto_nutrients


# PART 14: CREATE SITE × SPECIES MATRICES FOR JSDM ------------------------


# Helper function to create matrix
create_species_matrix <- function(data, response_var = "DENSITY") {
  data %>%
    filter(VISIT_NO == 1) %>%
    filter(!is.na(.data[[response_var]])) %>%  # Remove NAs
    select(SITE_ID, TARGET_TAXON, all_of(response_var)) %>%
    pivot_wider(names_from = TARGET_TAXON,
                values_from = all_of(response_var),
                values_fill = 0,
                values_fn = mean) # In case of duplicates, might change this
}

# Create matrices
zoo_matrix <- create_species_matrix(zoop, "DENSITY")
dim(zoo_matrix)

benthos_matrix <- create_species_matrix(benthos, "TOTAL")
dim(benthos_matrix)

phyto_matrix <- create_species_matrix(phyto, "BIOVOLUME")  
dim(phyto_matrix)

# Check sparsity
calc_zeros <- function(mat) {
  vals <- as.matrix(mat %>% select(-SITE_ID))
  sum(vals == 0) / length(vals) * 100
}


calc_zeros(zoo_matrix)
calc_zeros(benthos_matrix)
calc_zeros(phyto_matrix)

