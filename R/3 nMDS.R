# Open libraries y clean console ####
pacman::p_load(tidyverse, # Data wrangling
               vegan, # nMDS (distance calculation)
               gridExtra, # Multiple plots
               ggrepel) # Plotting text format

rm(list = ls())
shell('cls') # For Windows users
system2('clear') # For Mac users

# Open database ####
origin<- read.csv('Data/Abundance_data.csv',
                  header = T, stringsAsFactors = T)
#origin <- read.csv(here:: here('Data/Diversity_indices.csv'),
#                   header = T, stringsAsFactors = T)

# Set variable "Year" to factor
origin$Year<- as.factor(origin$Year)
as_tibble(origin)

# Set depth zones based on the light processing file ####
Data<- mutate(origin, Zone =
                ifelse(as.numeric(as.character(Initial_depth)) >= 21 &
                         as.numeric(as.character(Final_depth)) >= 21 |
                         as.numeric(as.character(Initial_depth)) >= 21 &
                         Final_depth == '>30' |
                         Initial_depth == '>30' & Final_depth == '>30',
                       'Mesophotic', 'Shallow'))

# Transform NAs to Shallow
Data$Zone[is.na(Data$Zone)]<- 'Shallow'
which(is.na(Data$Zone))

# Convert Zone to factor
Data$Zone<- as.factor(Data$Zone)

# Function to delete columns with abudance 0 ####
delete_sp_cero <- function(df) {
  # Find columns
  columns_to_stay <- apply(df[, 7:ncol(df)], 2, function(col) any(col >= 1))
  # Save data
  columns_to_stay <- c(rep(TRUE, 6), columns_to_stay)
  # Select columns
  df_filtrado <- df[, columns_to_stay]
  
  return(df_filtrado)
}

# Delete columns with abundance 0 ####
Data_filt <- delete_sp_cero(Data)

# Save zone column
last_column <- ncol(Data_filt)

# Arrange columns
data_nMDS <- Data_filt[, c(1:6, last_column, 7:(last_column - 1))]

# Set factors ####
factors<- data_nMDS[, c(1:5, 6:7)]

# Select abundance data ####
abundance<- data_nMDS[, c(8:103)]
rownames(abundance)<- factors[, 6]

which(colSums(abundance) == 0)
which(rowSums(abundance) == 0)

# Filter rows >0
abundance_filtered <- abundance[rowSums(abundance) > 0, ]

which(colSums(abundance_filtered) == 0)
which(rowSums(abundance_filtered) == 0)

# Select factors rows filtered (rows 40 85 91)
rows_to_delete <- c(40, 85, 91)

# Delete rows
factors_filtered <- factors[-rows_to_delete, ]

# Join data ####
full_data <- cbind(factors_filtered, abundance_filtered)
#write_csv(full_data, "Data/nMDS_data.csv")

# Add substrate data
substrate_data <- read.csv('Data/Substrate_data.csv', 
                           sep = ",", header = TRUE, stringsAsFactors = FALSE)
data_combined <- merge(full_data, substrate_data, by = "Video.transect", 
                       all.x = TRUE)

# Substrate matrix
substrate<- data_combined[, c(104:110)]
rownames(substrate)<- data_combined[, 1]

# Fourth root transformation ####
transf<- abundance_filtered^0.25

# Bray-Curtis distance calculation ####
distance<- vegdist(transf, method = 'bray')

# nMDS calculation ####
nMDS.drrh<- metaMDS(comm = distance, trace = F,
                    autotransform = F)
stress<- round(nMDS.drrh$stress, 2)
stress

# Convert NAs
substrate_clean <- substrate
substrate_clean[is.na(substrate_clean)] <- 0

# Luego ejecuta envfit con los datos limpios
substrate_fit <- envfit(nMDS.drrh, substrate_clean, perm = 999)

# Results
print(substrate_fit)

# Extract vectors coordenates 
substrate_scores <- as.data.frame(scores(substrate_fit, display = "vectors"))
substrate_scores$variables <- rownames(substrate_scores)

# Scale factor
scaling_factor <- 7 

# Scale vectors coordenates
substrate_scores$NMDS1 <- substrate_scores$NMDS1 * scaling_factor
substrate_scores$NMDS2 <- substrate_scores$NMDS2 * scaling_factor

# Add data to MDS.coord df
MDS.coord<- data.frame(nMDS.drrh$points)

MDS.coord$Site<- factors_filtered$Site
MDS.coord$Zone<- factors_filtered$Zone
MDS.coord$Video.transect<- factors_filtered$Video.transect
MDS.coord$data_nMDS<- data_nMDS$data_nMDS

# Plots ####
gg <- merge(MDS.coord,
            aggregate(cbind(mean.x = MDS1, mean.y = MDS2) ~ Zone,
                      MDS.coord, mean), by = "Zone")

# Zone Plot
zone_plot <- ggplot() + 
  geom_point(data = gg, aes(MDS1, MDS2, color = Zone), size = 2.5) + 
  stat_ellipse(data = gg, aes(MDS1, MDS2, fill = Zone), level = 0.95, 
               geom = "polygon", alpha = 0.3, color = NA) +
  scale_colour_manual(values = c("#27408B", "#8B0000")) + 
  scale_fill_manual(values = c("#27408B", "#8B0000")) +
  labs(x="nMDS1", y="nMDS2")+ 
  annotate("text", x = 1.25, y = 0.85, label = paste("Stress =", stress),
           size = 10) + 
  xlim(-2, 2) + ylim(-1, 1) + 
  theme_bw(base_size = 25) + 
  theme(panel.grid = element_blank())

# Plot substrate vectors 

zone_substrate_plot <- zone_plot + 
  geom_segment(data = substrate_scores, aes(x = 0, xend = NMDS1, 
                                            y = 0, yend = NMDS2), 
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", 
               size = 1) +
  geom_text_repel(data = substrate_scores, aes(x = NMDS1 * (1.2), 
                                               y = NMDS2 * (1.2), 
                                               label = variables), 
                  size = 7, nudge_x = 0.06, nudge_y = 0.02,
                  segment.color = NA)

ggsave('Figs/Figure 2.tiff', plot = zone_substrate_plot, width = 4350,
       height = 3000, units = 'px', dpi = 320, compression = "lzw")
#ggsave(here:: here('Figs/Figure 2.tiff'), plot = multi, width = 4350,
#       height = 3000, units = 'px', dpi = 320, compression = "lzw")
