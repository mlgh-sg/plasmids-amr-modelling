library(brms)    # For Bayesian modeling
library(dplyr)   # For data manipulation
library(readr)   # For data reading
library(ggplot2)
library(here)
library(cmdstanr)
# Read the data
data <- read_csv(here("output","modelling_data.csv"))
# Inspect the data structure
str(data)

# Create Donor-Recipient Pair
data <- data %>%
  mutate(Donor_Recipient_Pair = paste(Conjugation_Donor, Conjugation_Recipient, sep = "_")) %>%
  mutate(Donor_Recipient_Pair = as.factor(Donor_Recipient_Pair)) %>%
  mutate(
    Plasmid = as.factor(Plasmid),
    Conjugation_Donor = as.factor(Conjugation_Donor),
    Conjugation_Recipient = as.factor(Conjugation_Recipient),
    log_Change_Transconjugants = log(Change_Transconjugants + 1),
    log_Donor_Count = log(Donor_Count + 1),
    log_True_Recipient_Count = log(True_Recipient_Count + 1),
    log_Transconjugants_Count = log(Transconjugants_Count + 1)
  )

# Define the model formula with interactions
model_formula <- bf(
  Change_Transconjugants ~ 
    Donor_Count + True_Recipient_Count + Transconjugants_Count + (1 | Plasmid) + (1+Donor_Recipient_Pair)
)

# Fit the model
model <- brm(
  formula = model_formula,
  data = data,
  family = gaussian(),  # Since we're modeling log-transformed data
  prior = c(
    set_prior("normal(0, 10)", class = "b")    # Priors for fixed effects
  ),
  cores = 4,            # Adjust based on your machine
  chains = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "cmdstanr",
  seed = 123            # For reproducibility
)
# Summarize the model
summary(model)
