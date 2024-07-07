# =============================================================================
# TO DO: Data preparation of Plague outbreak dataset
# Author: Andrea Pinke
# =============================================================================

# Load required packages
library(tidyverse)

# Load the data
df_plague <- read.csv(file = "Data/data_plague_2017_madagascar.csv")

# Introduce total infected number as pneumonic + bubonic cases
df_plague %>%
  mutate(total = pneumonic+bubonic) -> df_plague

df_plague_relevant <- select(df_plague, c("date", "total"))

## Preparations for dealing with date column
# convert date to proper format
df_plague_relevant$date <- as.Date(df_plague_relevant$date, format = "%d/%m/%Y")

# as.numeric gives number of days since "1970-01-01"
df_plague_relevant$date <- as.numeric(df_plague_relevant$date)

# rename columns
df_plague_relevant <- rename(df_plague_relevant, c(time = date, I1 = total))

# re-scale time for ode solver to start with day 1
df_plague_relevant$time <- df_plague_relevant$time - 17378

# ahift to 53th day (=start of large epidemic wave) --> same as in paper
df_plague_relevant <- filter(df_plague_relevant, time>52)

# re-scale time for ode solver to start with day 1
df_plague_relevant$time <- df_plague_relevant$time - 52

# save dataset
df_plague_relevant %>%
  write.csv("Data/relevant_data_plague_2017_madagascar.csv", row.names = FALSE)