# Open libraries y clean console ####
pacman::p_load(tidyverse, # Data wrangling
               vegan, # nMDS (distance calculation)
               ggplot2) # Graphs

rm(list = ls())
# system2("clear")
# shell('cls')

# Open database ####
origin<- read.csv('Data/Abundance data.csv',
                  header = T, stringsAsFactors = T)

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
write_csv(full_data, "Data/abundance_nMDS_data.csv")

# Fourth root transformation ####
transf<- abundance_filtered^0.25

# Bray-Curtis distance calculation ####
distance<- vegdist(transf, method = 'bray')

# nMDS calculation ####
nMDS.drrh<- metaMDS(comm = distance, trace = F,
                    autotransform = F)
stress<- round(nMDS.drrh$stress, 2)
stress

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
zone_plot <- ggplot(gg, aes(MDS1, MDS2, color = Zone))+
  scale_colour_manual(values = c("#27408B", "#8B0000"))+
  geom_point(size = 2)+
  geom_segment(aes(x = mean.x, y = mean.y, xend = MDS1, yend = MDS2))+
  geom_point(aes(x = mean.x, y = mean.y), size = 2.5, color = "black")+
  stat_ellipse(level = 0.95)+
  labs(x = "nMDS1", y = "nMDS2")+
  annotate("text", x = 5, y = -1, label = paste("Stress =", stress),
           size = 15) +
  xlim(-2.5, 7.5) + ylim(-1.1, 1.1) +
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank())

forms <- c("Mesophotic" = 16, "Shallow" = 17)

# Sites plot
Sites_plot <- ggplot(gg, aes(MDS1, MDS2, color = Site, shape = Zone)) +
  scale_colour_manual(values = c("#27408B", "#8B0000", "#989823")) +
  geom_point(size = 4) +
  scale_shape_manual(values = forms) +
  labs(x = "nMDS1", y = "nMDS2") +
  annotate("text", x = 5, y = -1, label = paste("Stress =", stress),
           size = 15) +
  xlim(-2.5, 7.5) + ylim(-1.1, 1.1)+
  labs(shape = "Site") + guides(shape = guide_legend(title = "Zone")) +
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank())

# Basic Zone Plot 
bzone_plot <- ggplot(gg, aes(MDS1, MDS2, color = Zone))+
  scale_colour_manual(values = c("#27408B", "#8B0000"))+
  geom_point(size = 3)+
  labs(x = "nMDS1", y = "nMDS2")+
  annotate("text", x = 5, y = -1, label = paste("Stress =", stress),
           size = 15) +
  xlim(-2.5, 7.5) + ylim(-1.1, 1.1) +
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank())
