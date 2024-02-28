# Rscript 'Fire and Summer interact to shape seed dormancy thresholds'
# Zomer et al. 2022

#check working directory
getwd()

#source functions for script
source(file="R/functions_seeds.R")

# load libraries
library(DT)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(rcartocolor)
library(drc)
library(drc)
library(lmtest)
library(sandwich)
library(cowplot)
library(performance)


### Read germination data -----
seeds_data = read.csv("input/germination_data.csv")


# 1. Germination Indices -----

#remove empty seeds from total seeds, and calculate percentage germinated of new total
seeds_data <- seeds_data %>%
  mutate(new_total_seeds = total_seeds - coalesce(empty_seeds, 0),
         germination_percentage = (germinated_seeds / new_total_seeds) * 100)

# Show data in a table
datatable(seeds_data, filter = "top") %>%
  formatRound('germination_percentage', digits = 3)


#Standardize for non-dormant seeds ------
# The lowest heat treatment of 30 °C, equivalent to the mean temperature during the summer, 
# was used to standardize the germination curve. 
# This temperature has been previously shown to not break seed dormancy for our species (Moreira and Pausas, 2012), 
# therefore germination at this heat treatment represents the fraction of non-dormant seeds. 
# For each population, we removed this non-dormant fraction from the final analysis. 


# Calculate mean germination percentage for the '30' heat treatment: represent non-dormant seeds
mean_gp_control <- seeds_data %>%
  filter(heat_treatment == '30') %>%
  group_by(species, population) %>%
  summarise(gp_control = mean(germination_percentage, na.rm = TRUE)) %>%
  ungroup()

# Merge mean germination data back into the original data frame
#Standardize: germination_percentage - % non dormant seeds) / *(100 %  - non dormant seeds%) *100
#*Subtracting the non-dormant seeds percentage from 100% provides the percentage of viable seeds that could potentially germinate under optimal conditions.
#Formula scales the difference in germination percentage by the potential maximum germination percentage, 
#Calculate both percentage gp and decimal form gd

merged_data <- seeds_data %>%
  left_join(mean_gp_control, by = c("species", "population")) %>%
  mutate(gp_standardized = (germination_percentage - gp_control) / (100 - gp_control) * 100,
         gd_standardized = (germination_percentage - gp_control) / (100 - gp_control))


#Identify maximum germination and the corresponding heat treatment for each population -----

# Calculate the mean_standardized germination of the three replicates for each species, population, and heat treatment
mean_standardized <- merged_data %>%
  group_by(species, population, heat_treatment) %>%
  summarise(mean_gp_standardized = mean(gp_standardized))

# Identify the maximum mean_standardized germination for each species and population
#CSA popluation 11 has two equal max germination percentage - choose highest correspondin temperature
max_standardized <- mean_standardized %>%
  group_by(species, population) %>%
  filter(mean_gp_standardized == max(mean_gp_standardized)) %>%
  slice_max(heat_treatment) %>%
  ungroup()


# Rename 'heat_treatment' column in 'max_standardized' to 'max_heat_treatment' and mean_gp to max_gp_mean
max_standardized <- max_standardized %>%
  rename(max_heat_treatment = heat_treatment,
         max_gp_mean = mean_gp_standardized)

# Merge 'max_standardized' back into 'merged_data'
final_data <- merged_data %>%
  left_join(max_standardized, by = c("species", "population"))

#Set gp_standardized and gd_standardized in heat treatment 30 to 0
# subtracting the mean of the control treatment from itself (each replicate) does not provide useful information 
final_data <- final_data %>%
  mutate(gp_standardized = ifelse(heat_treatment == '30', 0, gp_standardized),
         gd_standardized = ifelse(heat_treatment == '30', 0, gd_standardized))

#Save new dataframe as RDS file
saveRDS(final_data, file = "input/final_data.rds")

#Plots ----

#set colour scheme
display_carto_all()
my_palette <- rcartocolor::carto_pal(n = 12, name = "Prism")[c(2, 6)]

# Set order of factor 'species'
max_standardized$species <- factor(max_standardized$species, levels = c("CAL", "CSA"))

# jittered points for original data, set dodge width
dodge <- position_dodge(width = 5)

# box and whisker plot of maximum germination for each population

a<-ggplot(max_standardized, aes(x = species, y = max_gp_mean, color = species, fill = species)) +
  scale_y_continuous(name = "Maximum germination (%)", limits = c(50, 100)) +
  scale_color_manual(values = my_palette, guide = waiver()) +
  scale_fill_manual(values = my_palette, guide = waiver()) +
  scale_x_discrete(name = "", labels = c("CAL" = "C. albidus", "CSA" = "C. salviifolius")) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3) +
  theme(axis.text.x = element_text(size = 14)) +
  theme_classic2(base_size = 14) +
  geom_point(
    position = position_jitter(width = 0.1, seed = 0),
    size = 4, alpha = 0.5,
  )


b <- ggplot(max_standardized, aes(x = species, y = max_heat_treatment, color = species, fill = species)) +
  scale_y_continuous(name = "Heat for maximum germination (°C)", limits = c(70, 120)) +
  scale_color_manual(values = my_palette, guide = waiver()) +
  scale_fill_manual(values = my_palette, guide = waiver()) +
  scale_x_discrete(name = "", labels = c("CAL" = "C. albidus", "CSA" = "C. salviifolius")) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3) +
  theme(axis.text.x = element_text(size = 14)) +
  theme_classic2(base_size = 14) +
  geom_point(position = position_jitter(width = 0.1, seed = 0), size = 4, alpha = 0.5, )


# arrange plots on same page
a + b

# 2. DRC curves ----


# DRC curve for each population (all replicates) to find effective dose 20, 30, 40 , 50, 60 % dormancy release of total dormant seed population

# Four parameters are necessary because the standardized germination % (removing g% of 30) drops below 0 for many populations in treatments 40 - 60
# This reflects mortality of non-dormant seeds with low heat before dormancy release for dormant seeds.


#upload clean data
data <- readRDS(file = "input/final_data.rds")

# Define custom function to select peak and filter data
select_peak_and_filter <- function(data, population, max_heat_treatment) {
  # Print population
  cat("Population:", population, "\n")
  
  # Get unique maximum heat treatment value for the population
  unique_max_heat_treatment <- unique(max_heat_treatment)
  
  # Print the maximum heat treatment value (if it's unique)
  if (length(unique_max_heat_treatment) == 1) {
    cat("Max Treatment:", unique_max_heat_treatment, "\n")
  }
  
  # Calculate mean germination percentage for each heat treatment
  mean_germination <- data %>%
    group_by(species, population, heat_treatment) %>%
    summarise(mean_gp = mean(gp_standardized))
  
  # Plot mean germination percentage
  plot(mean_germination$heat_treatment, mean_germination$mean_gp, 
       type = "o", 
       xlab = "Heat Treatment", 
       ylab = "Mean Germination Percentage",
       main = paste("Mean Germination Percentage vs. Heat Treatment for Population:", population, "\n", "Max Treatment:", unique_max_heat_treatment))
  
  # Customize x-axis to split by 5 degrees
  axis(1, at = seq(0, max(mean_germination$heat_treatment), by = 5))
  
  # Allow user to select peak interactively
  peak <- readline(prompt = "Enter peak (y/n): ")
  
  if (tolower(peak) == "y") {
    peak_value <- as.numeric(readline(prompt = "Enter peak value: "))

# Filter data based on peak heat treatment
filtered_data <- data[data$heat_treatment <= peak_value, ]
return(filtered_data)
  } else {
    # If peak is not selected, use the original max_heat_treatment value for filtering
    max_heat_treatment <- max(max_heat_treatment)
    filtered_data <- data[data$heat_treatment <= max_heat_treatment, ]
    return(filtered_data)
  }
}

# Split data into separate data frames for each population of Cistus salviifolius
csapop_split <- split(data[data$species == "CSA", ], data[data$species == "CSA", ]$population)

# Split data into separate data frames for each population of Cistus albidus
calpop_split <- split(data[data$species == "CAL", ], data[data$species == "CAL", ]$population)

# select peak
#pass through population plots one by one
#approve or manually set peak heat

#csa 3 change to 120 ;csa 10 to 100
csa_populations <- Map(select_peak_and_filter, csapop_split, names(csapop_split), lapply(csapop_split, function(data) data$max_heat_treatment))


#cal 3 change to 95 ;cal 8 to 90; cal 11 to 95; cal 13 to 90; cal 15 to 105
# cal 18 to 90; cal 19 to 85; cal 20 to 90
cal_populations <- Map(select_peak_and_filter, calpop_split, names(calpop_split), lapply(calpop_split, function(data) data$max_heat_treatment))


combined_csa <- bind_rows(csa_populations, .id = "population_id")
combined_cal <- bind_rows(cal_populations, .id = "population_id")

#Save filtered data as RDS file
saveRDS(combined_csa, file = "input/filtered_csa_data.rds")
saveRDS(combined_cal, file = "input/filtered_cal_data.rds")



#upload filtered data-----
csa <- readRDS(file = "input/filtered_csa_data.rds")
cal <- readRDS(file = "input/filtered_cal_data.rds")


#unfortunately there is a bug when mselect function is applied to a list or used in a function, 
#so we have to apply drm and mselect functions manually
#Step 1: models drm()
#Step 2: mselect (model) to select most parsimonous modely type
#Step 3: rerun model with selected model type
# germination decimal format used in order to later use effective dose (ED function) type 'absolute'


### Models----
# (correct model already selected using mselect function - see section below)

model_csa_1 <- drm(gd_standardized ~ heat_treatment, data = csa[csa$population == "1",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_csa_2 <- drm(gd_standardized ~ heat_treatment, data = csa[csa$population == "2",], fct = W2.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_csa_3 <- drm(gd_standardized ~ heat_treatment, data = csa[csa$population == "3",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_csa_4 <- drm(gd_standardized ~ heat_treatment, data = csa[csa$population == "4",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_csa_5 <- drm(gd_standardized ~ heat_treatment, data = csa[csa$population == "5",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_csa_6 <- drm(gd_standardized ~ heat_treatment, data = csa[csa$population == "6",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_csa_7 <- drm(gd_standardized ~ heat_treatment, data = csa[csa$population == "7",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_csa_8 <- drm(gd_standardized ~ heat_treatment, data = csa[csa$population == "8",], fct = W2.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_csa_9 <- drm(gd_standardized ~ heat_treatment, data = csa[csa$population == "9",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_csa_10 <- drm(gd_standardized ~ heat_treatment, data = csa[csa$population == "10",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_csa_11 <- drm(gd_standardized ~ heat_treatment, data = csa[csa$population == "11",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

model_cal_1 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "1",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_2 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "2",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_3 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "3",], fct = W2.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_4 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "4",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_5 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "5",], fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_6 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "6",], fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_7 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "7",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_8 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "8",], fct = W2.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_9 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "9",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_10 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "10",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_11 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "11",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_12 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "12",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_13 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "13",], fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_14 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "14",], fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_15 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "15",], fct = W2.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_16 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "16",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_17 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "17",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_18 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "18",], fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_19 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "19",], fct = W2.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
model_cal_20 <- drm(gd_standardized ~ heat_treatment, data = cal[cal$population == "20",], fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))


### MySelect ------
# select most parsimonous model based on lowest AIC and highest log-likelihood values
# four parameter log logistic curve (LL.4 or LL2.4) or four parameter Weibull curve (W1.4 or W2.4)

mselect(model_csa_1, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_csa_2, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_csa_3, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_csa_4, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_csa_5, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_csa_6, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_csa_7, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_csa_8, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_csa_9, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_csa_10, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_csa_11, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)

mselect(model_cal_1, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_2, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_3, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_4, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_5, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_6, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_7, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_8, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_9, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_10, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_11, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_12, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_13, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_14, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_15, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_16, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_17, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_18, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_19, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)
mselect(model_cal_20, fctList = list(LL.4(), W1.4(), W2.4(), LL2.4()), linreg = TRUE)


# model list
models <- list(
  model_csa_1 = model_csa_1, model_csa_2 = model_csa_2, model_csa_3 = model_csa_3, model_csa_4 = model_csa_4, model_csa_5 = model_csa_5, model_csa_6 = model_csa_6,
  model_csa_7 = model_csa_7, model_csa_8 = model_csa_8, model_csa_9 = model_csa_9, model_csa_10 = model_csa_10, model_csa_11 = model_csa_11,
  model_cal_1 = model_cal_1, model_cal_2 = model_cal_2, model_cal_3 = model_cal_3, model_cal_4 = model_cal_4, model_cal_5 = model_cal_5, model_cal_6 = model_cal_6, model_cal_7 = model_cal_7,
  model_cal_8 = model_cal_8, model_cal_9 = model_cal_9, model_cal_10 = model_cal_10, model_cal_11 = model_cal_11, model_cal_12 = model_cal_12, model_cal_13 = model_cal_13,
  model_cal_14 = model_cal_14, model_cal_15 = model_cal_15, model_cal_16 = model_cal_16, model_cal_17 = model_cal_17, model_cal_18 = model_cal_18, model_cal_19 = model_cal_19,
  model_cal_20 = model_cal_20
)

# Iterate over each model in the list and plot them
for (model_name in names(models)) {
  plot(models[[model_name]], main = model_name)
  # Wait for the plot to be closed
  input <- readline(prompt = "Press Enter to continue...")
}

### Effective Dose ED function ----
# Estimated effective doses (experimental heat) for the 20th – 60th percentiles of dormancy release
# obtained from ‘ED’ function in drc package;
# considering the whole seed population rather than only the germination fraction (type = “absolute”);
# with 95% confidence intervals calculated using the delta method (interval = “delta”).
# ie. ED (model, c(0.2, 0.3, 0.4, 0.5, 0.6),  type = "absolute", interval = "delta")

# function for 'ed'
ed <- function(x) {
  ED(x, c(0.2, 0.3, 0.4, 0.5, 0.6), type = "absolute", interval = "delta")
}

# apply ed function to the list of all models
ed.all <- sapply(models, FUN = ed, simplify = FALSE, USE.NAMES = TRUE)

# convert to dataframe
ed.df <- do.call(rbind.data.frame, ed.all)

# export to csv
#write.csv(ed.df, file = "input/ed.all.csv")

# note: in excel model name and dormancy release level are combined in first column.
# * reformatted csv loaded in next section
# e:1:0.2 = 0.2 or 20% dormancy release
# e:1:0.3 = 30%
# e:1:0.4 = 40%
# e:1:0.5 = 50%
# e:1:0.6 = 60%


# 3. Regression effective dose and climate ----

# reformatted effective dose csv with columns separated into species, pop, dormancy release, estimate
ed.estimate <- read.csv(("input/ed.all.csv"), header = TRUE, sep = ",", na.strings = c("", "NA"))

# site characteristics, climate, fire, and seed mass variables
sites <- read.csv(("input/site_data.csv"), header = TRUE, sep = ",", na.strings = c("", "NA"))

#merge
ed_climate <- merge(sites, ed.estimate, by = c("species", "pop"))

# clean data
ed_climate <- na.omit(ed_climate)
ed_climate$Estimate <- as.numeric(ed_climate$Estimate)

# export to csv
#write.csv(ed_climate, file = "input/ed_climate.csv")


#upload merged data
ed_climate = read.csv("input/ed_climate.csv")

#  colinearity of  variables
vif.mod <- lm(Estimate ~ species + dormancy_release + BIO_05 + BIO_01 + BIO_12 + elevation + seed_mass + burned_area_percent + lat, ed_climate)
check_collinearity(vif.mod)

# Model Selection
# add variables into model in order of contribution to variance
# compare to null model using ANOVA


# Create a null model
null_model <- lm(Estimate ~ 1, data = ed_climate)
AIC_null <- AIC(null_model)
cat("Null model AIC:", AIC_null, "\n")

# List of variables
variables <- c("dormancy_release", "species", "BIO_05", "BIO_01", "elevation", "seed_mass", "BIO_12", "aridity_index", "burned_area_percent", "annual_number_fires", "lat")

# Perform ANOVA for each variable
for (var in variables) {
  # Update the null model with the current variable
  updated_model <- update(null_model, reformulate(var, response = "Estimate"), data = ed_climate)
  
  # Perform ANOVA
  anova_result <- anova(null_model, updated_model)
  
  # Extract p-value
  p_value <- anova_result[2, "Pr(>F)"]
  
  # Display results
  cat("Variable:", var, "\n")
  cat("AIC:", AIC(updated_model), "\n")
  cat("p-value:", format(p_value, scientific = TRUE), "\n")
  if (p_value < 0.001) {
    cat("***\n")
  } else if (p_value < 0.01) {
    cat("**\n")
  } else if (p_value < 0.05) {
    cat("*\n")
  } else {
    cat("ns\n")
  }
  cat("\n")
}


#Iteration 2
# Create a null model
null_model_2 <- lm(Estimate ~ dormancy_release, data = ed_climate)
AIC_null_2 <- AIC(null_model_2)
cat("Null model AIC:", AIC_null_2, "\n")

# List of variables
variables_2 <- c( "species", "BIO_05", "BIO_01", "elevation", "seed_mass", "BIO_12", "aridity_index", "burned_area_percent", "annual_number_fires", "lat")

# Perform ANOVA for each variable
for (var in variables_2) {
  # Update the null model with the current variable
  updated_model_2 <- update(null_model_2, reformulate(c("dormancy_release", var), response = "Estimate"), data = ed_climate)
  
  # Perform ANOVA
  anova_result_2 <- anova(null_model_2, updated_model_2)
  
  # Extract p-value
  p_value <- anova_result_2[2, "Pr(>F)"]
  
  # Display results
  cat("Variable:", var, "\n")
  cat("AIC:", AIC(updated_model_2), "\n")
  cat("p-value:", format(p_value, scientific = TRUE), "\n")
  if (p_value < 0.001) {
    cat("***\n")
  } else if (p_value < 0.01) {
    cat("**\n")
  } else if (p_value < 0.05) {
    cat("*\n")
  } else {
    cat("ns\n")
  }
  cat("\n")
}


#Iteration 3

# Create a null model
null_model_3 <- lm(Estimate ~ dormancy_release + species, data = ed_climate)
AIC_null_3 <- AIC(null_model_3)
cat("Null model AIC:", AIC_null_3, "\n")

# List of variables
variables_3 <- c( "BIO_05", "BIO_01", "elevation", "seed_mass", "BIO_12", "aridity_index", "burned_area_percent", "annual_number_fires", "lat")

# Perform ANOVA for each variable
for (var in variables_3) {
  # Update the null model with the current variable
  updated_model_3 <- update(null_model_3, reformulate(c("dormancy_release", "species",var), response = "Estimate"), data = ed_climate)
  
  # Perform ANOVA
  anova_result_3 <- anova(null_model_3, updated_model_3)
  
  # Extract p-value
  p_value <- anova_result_3[nrow(anova_result_3), "Pr(>F)"]
  
  # Display results
  cat("Variable:", var, "\n")
  cat("AIC:", AIC(updated_model_3), "\n")
  cat("p-value:", format(p_value, scientific = TRUE), "\n")
  if (p_value < 0.001) {
    cat("***\n")
  } else if (p_value < 0.01) {
    cat("**\n")
  } else if (p_value < 0.05) {
    cat("*\n")
  } else {
    cat("ns\n")
  }
  cat("\n")
}

#Iteration 4
# Create a null model
null_model_4 <- lm(Estimate ~ dormancy_release + species +BIO_05, data = ed_climate)
AIC_null_4 <- AIC(null_model_4)
cat("Null model AIC:", AIC_null_4, "\n")

# List of variables
variables_4 <- c( "BIO_01", "elevation", "seed_mass", "BIO_12", "aridity_index", "burned_area_percent", "annual_number_fires", "lat")

# Perform ANOVA for each variable
for (var in variables_4) {
  # Update the null model with the current variable
  updated_model_4 <- update(null_model_4, reformulate(c("dormancy_release", "species","BIO_05",var), response = "Estimate"), data = ed_climate)
  
  # Perform ANOVA
  anova_result_4 <- anova(null_model_4, updated_model_4)
  
  # Extract p-value
  p_value <- anova_result_4[nrow(anova_result_4), "Pr(>F)"]
  
  # Display results
  cat("Variable:", var, "\n")
  cat("AIC:", AIC(updated_model_4), "\n")
  cat("p-value:", format(p_value, scientific = TRUE), "\n")
  if (p_value < 0.001) {
    cat("***\n")
  } else if (p_value < 0.01) {
    cat("**\n")
  } else if (p_value < 0.05) {
    cat("*\n")
  } else {
    cat("ns\n")
  }
  cat("\n")
}




# null model
m0 <- lm(Estimate ~ 1, ed_climate) # null model
AIC(m0)
summary(m0)

# test all variables
mdorm <- update(m0, . ~ . + dormancy_release) # % dormancy release
mspec <- update(m0, . ~ . + species) # species
mmax <- update(m0, . ~ . + BIO_05) # max temp warmest month
mtas <- update(m0, . ~ . + BIO_01) # annual temp
mel <- update(m0, . ~ . + elevation) # elevation
ms <- update(m0, . ~ . + seed_mass) # seed mass
mpr <- update(m0, . ~ . + BIO_12) # annual prec
mai <- update(m0, . ~ . + aridity_index) # aridity index
mperc <- update(m0, . ~ . + burned_area_percent) # percentage burned area in buffer circle
mcount <- update(m0, . ~ . + annual_number_fires)
mlat <- update(m0, . ~ . + lat)

anova(m0, mdorm) # 1.031e-14 ***
anova(m0, mspec) # 4.249e-10 ***
anova(m0, mmax) # 2.096e-07 ***
anova(m0, mtas) # 5.114e-05 ***
anova(m0, mel) #* 0.0001153 ***
anova(m0, ms) # ns
anova(m0, mpr) # ns
anova(m0, mai) # ns
anova(m0, mperc) # ns
anova(m0, mcount) # ns
anova(m0, mlat) # 0.008852**

# add variables into model in order of contribution to variance

m1 <- update(m0, . ~ . + dormancy_release)
AIC(m1) # 931.8578
anova(m0, m1) # 1.031e-14 ***
summary(m1)
anova(m1)

m2 <- update(m1, . ~ . + species)
AIC(m2) # 872.3845
anova(m1, m2) # 9.548e-15 ***
summary(m2)
anova(m2)

m3 <- update(m2, . ~ . + BIO_05) # max temp warmest month
AIC(m3) # 831.2089
anova(m2, m3) # 9.955e-11 ***
summary(m3)
anova(m3)

m4 <- update(m3, . ~ . + BIO_01) # annual temp
AIC(m4) # higher than m3 : 833.1821
anova(m3, m4) # not significant
summary(m4)
anova(m4)

m5 <- update(m3, . ~ . + elevation) # elevation
AIC(m5) # higher than m3 : 833.1272
anova(m3, m5) # not significant
summary(m5)
anova(m5)

m6 <- update(m3, . ~ . + lat) # latitude
AIC(m6) # higher than m3 : 833.1272
anova(m3, m6) # not significant
summary(m6)
anova(m6)


# interactions with three significant variables
m_int1 <- update(m3, . ~ . + BIO_05 * dormancy_release)
AIC(m_int1) # 824.1229
anova(m3, m_int1) # 0.003114 **
summary(m_int1)
anova(m_int1)

m_int2 <- update(m_int1, . ~ . + BIO_05 * species)
AIC(m_int2) # 798.5907
anova(m_int1, m_int2) # 2.936e-07 ***
anova(m_int2)
summary(m_int2)

m_int3 <- update(m_int2, . ~ . + dormancy_release * species) # not significant
AIC(m_int3) # higher than m_int2: 800.0697
anova(m_int2, m_int2) # not significant
summary(m_int3)

# Check other variables again adding on to m_int2
mel_2 <- update(m_int2, . ~ . + elevation) # elevation
ms_2 <- update(m_int2, . ~ . + seed_mass) # seed mass
mpr_2 <- update(m_int2, . ~ . + BIO_12) # annual prec
mtas_2 <- update(m_int2, . ~ . + BIO_01) # annual temp
mai_2 <- update(m_int2, . ~ . + aridity_index) # aridity index
mperc_2 <- update(m_int2, . ~ . + burned_area_percent)
mcount_2 <- update(m_int2, . ~ . + annual_number_fires)
mlat_2 <- update(m_int2, . ~ . + lat)
anova(m_int2, mel_2) # ns
anova(m_int2, ms_2) # ns
anova(m_int2, mpr_2) # ns
anova(m_int2, mtas_2) # ns
anova(m_int2, mai_2) # ns
anova(m_int2, mperc_2) # ns
anova(m_int2, mcount_2) # ns
anova(m_int2, mlat_2) # ns

# Final linear model = m_int2:
lm(Estimate ~ dormancy_release + species + BIO_05 + dormancy_release * BIO_05 + species * BIO_05, ed_climate)

### Predictions plot----

hist(ed_climate$Estimate)
shapiro.test(ed_climate$Estimate) #>0.05, normal distribution

# linear regression with maximum temperature of warmest month (BIO_05), dormancy release %, and species as fixed factors
lm_ed <- lm(Estimate ~ BIO_05 * dormancy_release + BIO_05 * species, ed_climate)

# lm parameters
anova(lm_ed)
summary(lm_ed)

par(mfrow = c(2, 2))
plot(lm_ed) # check residuals
par(mfrow = c(1, 1))

#' predict' r base to calculate predictions and confidence intervals
ed_climate$predlm <- predict(lm_ed)

# remove NA
ed_climate <- na.omit(ed_climate)

predslm.ed <- predict(lm_ed, interval = "confidence")
head(predslm.ed)

datlm <- cbind(ed_climate, predslm.ed)
head(datlm)


# colours
library(rcartocolor)
display_carto_all()

# set colour scheme
my_pal3 <- rcartocolor::carto_pal(n = 12, name = "Prism")[c(6, 2)]

# rename dormancy release levels
new_labels2 <- c("20" = "20", "30" = "30  ", "40" = "40", "50" = "50", "60" = "60")

# set species order
datlm$species <- factor(datlm$species, levels = c("CSA", "CAL"))

# plot predicted linear regression with confidence intervals and raw data (points)
plot <- ggplot(datlm, aes(x = BIO_05 / 10, y = Estimate, color = species, fill = species)) +
  geom_point(alpha = 0.27) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = NULL), alpha = .2) +
  geom_line(aes(y = fit), size = 0.7) +
  scale_colour_manual(values = my_pal3, name = "", labels = c("C. salviifolius", "C. albidus")) +
  scale_fill_manual(values = my_pal3, name = "", labels = c("C. salviifolius", "C. albidus")) +
  theme_classic()

# facet grid to divide by dormancy release %
plot <- plot + facet_grid(~dormancy_release, labeller = labeller(germ_per = new_labels2))

# axes titles
plot <- plot + xlab(expression("Maximum temperature of warmest month (°C)")) +
  ylab(expression(" Lower heat threshold (°C)"))

# format legend
plot <- plot + theme(panel.background = element_rect(fill = NA, color = "black")) +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = "NA")
  ) +
  theme(legend.position = c(0.19, 0.8), legend.direction = "horizontal")

# legend position
plot + theme(strip.text.x = element_text(face = "bold")) + guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE)) +
  ggtitle("Dormancy release (%)") +
  theme(plot.title = element_text(hjust = 0.07, vjust = 0.5, size = 9, face = "bold"))


library(Cairo)
ggsave("Fig4.eps", device = cairo_ps, fallback_resolution = 600, width = 17.3, height = 12, units = "cm")



## plots for powerpoint transparent -----

library(cowplot)

# boxplots
colours <- c("#7E351F", "#C1A09D")

colours2 <- c("#C1A09D", "#7E351F")


germ1$species <- factor(germ1$species, level = c("CAL", "CSA"))

colours4 <- c("#C49B33", "#97BC92")

c5 <- c("steelblue3", "firebrick2")

e <- ggplot(germ1, aes(x = species, y = max_gp_st, fill = species)) +
  scale_fill_manual(values = colours4) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, colour = "white") +
  geom_jitter(color = "white", size = 3, alpha = 0.4, show.legend = FALSE) +
  theme_classic() +
  scale_y_continuous(name = "Maximum germination (%)", limits = c(50, 100)) +
  scale_x_discrete(name = "", labels = c("CAL" = "C. albidus", "CSA" = "C. salviifolius")) +
  theme(text = element_text(size = 24, color = "white")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 20, color = "white")) +
  theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )
e

f <- ggplot(germ1, aes(x = species, y = max_treat_st, fill = species)) +
  scale_fill_manual(values = c6) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, colour = "black", alpha = 1) +
  theme_classic() +
  scale_y_continuous(name = "Heat for maximum germination (°C)", limits = c(70, 120)) +
  scale_x_discrete(name = "", labels = c("CAL" = "C. albidus", "CSA" = "C. salviifolius")) +
  theme(text = element_text(size = 24, color = "black")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 20, color = "black"))

f
o4 <- cowplot::plot_grid(e, f, nrow = 1)
o4
ggsave(
  plot = o4, file = "seedbox12.png",
  type = "cairo-png",
  width = 20, height = 15, units = "cm", dpi = 800
)


geom_jitter(size = 5, alpha = 0.4, show.legend = FALSE) +
  +
    theme(
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", colour = NA)
    )
f

o3 <- cowplot::plot_grid(e, f, nrow = 1)
o3
ggsave(
  plot = o3, file = "seedbox4.png",
  type = "cairo-png", bg = "transparent",
  width = 20, height = 15, units = "cm", dpi = 800
)


# threshold linear regression
datlm$species <- factor(datlm$species, levels = c("CSA", "CSA"))
colours3 <- c("#97BC92", "#C49B33")

c4 <- c("#B4A4A1", "#EDBE60")

c6 <- c("#66B5C0", "#F2B15D")

# plot predicted linear regression with confidence intervals and raw data (points)
plot_t <- ggplot(datlm, aes(x = BIO_05 / 10, y = Estimate, color = species, fill = species)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = NULL), alpha = 0.9, show.legend = FALSE) +
  geom_line(aes(y = fit), size = 0.5, show.legend = FALSE, color = "black") +
  scale_colour_manual(values = c6, name = "", labels = c("C. salviifolius", "C. albidus")) +
  scale_fill_manual(values = c6, name = "", labels = c("C. salviifolius", "C. albidus")) +
  theme_classic()

# facet grid to divide by dormancy release %
plot_t <- plot_t + facet_grid(~dormancy_release, labeller = labeller(germ_per = new_labels2))

# axes titles
plot_t <- plot_t + xlab(expression("Maximum temperature of warmest month (°C)")) +
  ylab(expression(" Lower heat threshold (°C)"))
plot_t

# format legend
plot_t <- plot_t + theme(panel.background = element_rect(fill = NA, color = "black"))

# legend position
pt <- plot_t + theme(strip.text.x = element_text(face = "bold")) + guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE)) +
  ggtitle("Dormancy release (%)") +
  theme(plot.title = element_text(hjust = 0.07, vjust = 0.5, size = 9, face = "bold"))


pt

library(Cairo)
ggsave("seed_reg15.eps", device = cairo_ps, fallback_resolution = 600, width = 17.3, height = 12, units = "cm")



# R core and R package citations -----
library(knitcitations)
knitr::write_bib(c(.packages()), "packages.bib")

citation()
