#load packages
library(ggplot2)

# Set working directory
setwd("C:/Users/lihel/Documents/master/so24/case_studies/CaseStudiesLifeSciences/")

# Load data
df_plague <- read.csv(file = "data_plague_2017_madagascar.csv")
df_plague_relevant <- read.csv(file="data_plague_relevant.csv")


# =============================================================================
#plot 1
#plot of data set (number infected over time)
df_plague_relevant %>%
  ggplot(aes(x = time, y = I1)) +
  geom_point(color="lightseagreen", size=2) +
  labs(title="Visualization of the available data set",
       subtitle="Madagascar plague epedemic 2017",
       x="days since outbreak",
       y="infected individuals") +
  theme_gray(base_size=20)

# =============================================================================
#plot 3
#plot visualizing the fitting methods

