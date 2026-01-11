# Generate phenotypic variation data
# Author: Elvira Castillo-Almansa
# E-mail: castilloalmansae@gmail.com

###################

# Load packages
library(stats)
library(vegan)
library(terra)
library(car)
library(caret)

# Load data
dataset <- read.csv2("/file.csv")


# Chi-squarde test
table <- xtabs(formula = ~ morphotype + wei_ggroup, data = dataset)
chi_test <- chisq.test(table, correct = TRUE)
chi_test$stdres
print(chi_test)


# Violin plots 

# Load data with coordinates
df <- read.csv2("/coordenates.csv")

# Load variables: bio1, bio6, bio12 y altitude
temp <- rast("/bio_1.tif")
cold <- rast("/bio_6.tif")
prec <- rast("/bio_12.tif")
altitude <- rast("/elevation_map.tif")

# Extract variable values for each coordenates and add it to the dataframe
elev_values <- extract(altitude, coords_df)
coords_elev <- cbind(coords_df, elev = elev_values[,2])
df$altitude <- coords_elev$elev

temp_values <- extract(temp, coords_df)
coords_temp <- cbind(coords_df, temp = temp_values[,2])
df$temperature <- coords_temp$temp

cold_values <- extract(cold, coords_df)
coords_cold <- cbind(coords_df, cold = cold_values[,2])
df$cold_intensity <- coords_cold$cold

prec_values <- extract(prec, coords_df)
coords_prec <- cbind(coords_df, prec = prec_values[,2])
df$prec <- coords_prec$prec

df <- na.omit(df) # omit NA rows
df$altitude[df$altitude < 0] <- 0 # omit elevation values under zero

ggplot(df, aes(x = morphotyoe, y = altitude, fill = morhotype)) +
  geom_violin(trim = TRUE) +
  scale_fill_brewer(palette = "Set2") +
  xlab("Morphotype") +
  ylab("Elevation (m)") +
  theme_classic() +
  theme(legend.position = "none")


