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
  theme_gray(base_size=15) 

# =============================================================================
#plot 3
#plot visualizing the fitting methods

# Data frame to specify the arrows
df_arrows_up <- data.frame(
  time = c(15, 23, 40),     # Specify the x positions for arrows
  I1 = c(64, 50, 24),        # Specify the y positions for arrows
  arrow_length = c(7, 6, 6) # Length of the arrows
)

df_arrows_down <- data.frame(
  time = c(7),     # Specify the x positions for arrows
  I1 = c(40),        # Specify the y positions for arrows
  arrow_length = c(11) # Length of the arrows
)

# Data frame for labels
df_labels <- data.frame(
  time = 35,              # x position for the label
  I1 = 47              # y position for the label
  #,label = "y - hat(y)"  # Text for the label
)

df_plague_relevant %>%
  ggplot(aes(x = time, y = I1)) +
  geom_point(color="lightseagreen", size=2) +
  labs(title="Fit curve to data",
       subtitle="Minimize vertical distance to data points",
       x="days since outbreak",
       y="infected individuals") +
  theme_gray(base_size=15) +
  geom_line(aes(x=time, y=out[,"I1"]), size=1) +
  geom_segment(data = df_arrows_up, aes(x = time, xend = time, y = I1 - arrow_length, yend = I1 + arrow_length),
             arrow = arrow(length = unit(0.4, "cm")), color = "orange", size=1) +
  geom_segment(data = df_arrows_down, aes(x = time, xend = time, yend = I1 - arrow_length, y = I1 + arrow_length),
               arrow = arrow(length = unit(0.4, "cm")), color = "orange", size=1) #+
  # geom_text(data = df_labels, 
  #           aes(x = time, label = label), 
  #           vjust = -1, color = "orange", size = 5, na.rm = TRUE) # Add label 

