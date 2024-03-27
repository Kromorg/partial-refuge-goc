# Open libraries y clean console ####
pacman:: p_load(tidyverse, # Data wrangling
                ncdf4, # NetCDF file
                raster, # Raster objects and extract values
                gridExtra) # Multiple plots

rm(list = ls())
shell('cls')

# Either create a new project to set your working directory more easily
# or use 'here' package.
# install.packages('here')
# library(here)


# NetCDF files ####
# Object of the file that will be used
ncfile<- 'Data/Water transparency.nc'
# ncfile<- here:: here('Data/Water transparency.nc')

# Open NetCDF file
cmes_data <- nc_open(ncfile)

# File information
print(cmes_data)
attributes(cmes_data$var)

# Multiband ####
# Raster object using the NetCDF file
multi_transp <- brick(ncfile, varname = 'KD490')

# Analysis using the raster object ####

# Statistical summary
# Mean and standard deviation estimates
polygon_mean <- calc(multi_transp, fun = mean)
polygon_sd <- calc(multi_transp, fun = sd)

# Plot multiband data
plot(polygon_mean, main = "Average Kd490")
plot(polygon_sd, main = "Standard deviation Kd490")

# Date extraction ####
dates <- substr(names(multi_transp), 2, 11)
dates <- dates %>% parse_date_time("Ymd")


# Extracting values from numeric model ####
sites_coords<- read.csv('Data/ROV coordinates.csv', header = T,
                        stringsAsFactors = T)
#sites_coords<- read.csv(here:: here('Data/ROV coordinates.csv'),
#                                    header = T, stringsAsFactors = T)
coords<- sites_coords[, -c(1:2)]

coordinates(coords)<- ~Long + Lat
crs(coords)<- crs(multi_transp)

# Extract mean values per month
mean_transp <-  raster:: extract(multi_transp, coords,
                        fun = mean, na.rm = F)
head(mean_transp, c(6, 4)) # Preview values

rownames(mean_transp)<- sites_coords$Site # Set site names
rownames(mean_transp)

# Analysis per site ####
# Los Islotes ####

# Create a dataframe with two columns (time and mean)
mean_islotes<- mean_transp[8, ]

islotes_transp <- data.frame(date = dates,
                             mean_transp = c(mean_islotes))

# Interanual mean KDPAR
interanual_islotes <- data.frame(islotes_transp,
           k = 0.0864+ (0.884*islotes_transp$mean_transp)-
             (0.00137*islotes_transp$mean_transp^-1))

# Selecting values from 2021-2022
transp.islotes<- interanual_islotes %>%
  filter(year(date) >= '2021' & year(date) <= '2022')

# Interanual mean and standard deviation of coefficients
# Kd490
round(mean(transp.islotes$mean_transp), 2)
round(sd(transp.islotes$mean_transp), 2)

transp.islotes %>%
  group_by(month(date)) %>%
  summarise(Media = mean(mean_transp),
            Maximo = max(mean_transp),
            Minimo = min(mean_transp))

#KdPAR
round(mean(transp.islotes$k), 2)
round(sd(transp.islotes$k), 2)

# Visualize dataframe
as_tibble(transp.islotes)

# Mesophotic zone at Los Islotes ####
mesophotic.islotes<- transp.islotes %>%
  group_by(month(date)) %>%
  summarise(Media = mean(k),
            Maximo = max(k),
            Minimo = min(k),
            Z10 = 2.3/max(k),
            Z1 = (4.6/((max(k) + min(k))/2)),
            Z0.1 = (6.9/min(k)))
mesophotic.islotes

# Graph
graph.islotes<- ggplot(transp.islotes, aes(x = date)) +
  geom_line(aes(y = mean_transp, color = "deepskyblue4"),
            linewidth = 1, alpha = 0.7, show.legend = F)+ 
  geom_line(aes(y = k, color = "peru"), linewidth = 1,
            show.legend = F)+
  scale_x_datetime(date_breaks = "2 months",
                   date_labels = "%m/%y")+
  labs(subtitle = expression(Los~Islotes),
       x = '',
       y = '')+
  scale_color_manual(labels = c(bquote(K[d490]), bquote(K[dPAR])),
                     values = c("deepskyblue4", "peru"))+
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = seq(0, 0.25, by = 0.05))+
  theme_classic(base_size = 25)


# Punta Lobos ####
# Create a dataframe with two columns(time and mean)
mean_punta_lobos<- mean_transp[15, ]

punta_lobos_transp <- data.frame(date = dates,
                             mean_transp = c(mean_punta_lobos))


# Interanual mean KDPAR
interanual_punta_lobos <- data.frame(punta_lobos_transp,
        k = 0.0864+ (0.884*punta_lobos_transp$mean_transp)-
          (0.00137*punta_lobos_transp$mean_transp^-1))

# Selecting values from 2021-2022
transp.lobos<- interanual_punta_lobos %>%
  filter(year(date) >= '2021' & year(date) <= '2022')

# Interanual mean and standard deviation of coefficients
# Kd490
round(mean(transp.lobos$mean_transp), 2)
round(sd(transp.lobos$mean_transp), 2)

transp.lobos %>%
  group_by(month(date)) %>%
  summarise(Media = mean(mean_transp),
            Maximo = max(mean_transp),
            Minimo = min(mean_transp))

# KdPAR
round(mean(transp.lobos$k), 2)
round(sd(transp.lobos$k), 2)

# Visualize dataframe
as_tibble(transp.lobos)

# Mesophotic zone at Punta Lobos ####
mesophotic.punta.lobos<- transp.lobos %>%
  group_by(month(date)) %>%
  summarise(Media = mean(k),
            Maximo = max(k),
            Minimo = min(k),
            Z10 = 2.3/max(k),
            Z1 = (4.6/((max(k) + min(k))/2)),
            Z0.1 = (6.9/min(k)))
mesophotic.punta.lobos

# Graph
graph.lobos<- ggplot(transp.lobos, aes(x = date)) +
  geom_line(aes(y = mean_transp, color = "deepskyblue4"),
            linewidth = 1, alpha = 0.7)+ 
  geom_line(aes(y = k, color = "peru"), linewidth = 1)+
  scale_x_datetime(date_breaks = "2 months",
                   date_labels = "%m/%y")+
  labs(subtitle = expression(Punta~Lobos),
       x = '',
       y = expression(Attenuation~coefficient~(m^-1)),
       color = '')+
  scale_color_manual(labels = c(bquote(K[d490]), bquote(K[dPAR])),
                     values = c("deepskyblue4", "peru"))+
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = seq(0, 0.25, by = 0.05))+
  theme_classic(base_size = 25)


# El Bajo de Espiritu Santo ####
# Create a dataframe with two columns(time and mean)
mean_ebes<- mean_transp[1, ]

ebes_transp <- data.frame(date = dates,
                          mean_transp = c(mean_ebes))

# Interanual mean KDPAR 
interanual_ebes <- data.frame(ebes_transp,
            k = 0.0864+ (0.884*ebes_transp$mean_transp)-
              (0.00137*ebes_transp$mean_transp^-1))

# Selecting values from 2021-2022
transp.ebes<- interanual_ebes %>%
  filter(year(date) >= '2021' & year(date) <= '2022')

# Interanual mean and standard deviation of coefficients
# Kd490
round(mean(transp.ebes$mean_transp), 2)
round(sd(transp.ebes$mean_transp), 2)

transp.ebes %>%
  group_by(month(date)) %>%
  summarise(Media = mean(mean_transp),
            Maximo = max(mean_transp),
            Minimo = min(mean_transp))

# KdPAR
round(mean(transp.ebes$k), 2)
round(sd(transp.ebes$k), 2)


# Visualize dataframe
as_tibble(transp.ebes)

# Mesophotic zone at El Bajo de Espiritu Santo ####
mesophotic.ebes<- transp.ebes %>%
  group_by(month(date)) %>%
  summarise(Media = mean(k),
            Maximo = max(k),
            Minimo = min(k),
            Z10 = 2.3/max(k),
            Z1 = (4.6/((max(k) + min(k))/2)),
            Z0.1 = (6.9/min(k)))
mesophotic.ebes

# Graph
graph.ebes<- ggplot(transp.ebes, aes(x = date)) +
  geom_line(aes(y = mean_transp, color = "deepskyblue4"),
            linewidth = 1, alpha = 0.7, show.legend = F)+ 
  geom_line(aes(y = k, color = "peru"), linewidth = 1,
            show.legend = F)+
  scale_x_datetime(date_breaks = "2 months",
                   date_labels = "%m/%y")+
  labs(subtitle = expression(EBES),
       x = "Time",
       y = '')+
  scale_color_manual(labels = c(bquote(K[d490]), bquote(K[dPAR])),
                     values = c("deepskyblue4", "peru"))+
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = seq(0, 0.25, by = 0.05))+
  theme_classic(base_size = 25)

grid.arrange(graph.islotes, graph.lobos, graph.ebes,
             ncol = 1, nrow = 3)
