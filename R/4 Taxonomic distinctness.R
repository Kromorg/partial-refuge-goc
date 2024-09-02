# Open libraries y clean console ####
pacman:: p_load(tidyverse, # Data wrangling
                vegan) # Taxonomic distinctness

rm(list = ls())
shell('cls')

# Open databases ####
origin <- read.csv('Data/Abundance data.csv',
                   header = T, stringsAsFactors = T)
#origin <- read.csv(here:: here('Data/Abundance data.csv'),
#                   header = T, stringsAsFactors = T)

as_tibble(origin)

# Set variable "Year" to factor
origin$Year<- as.factor(origin$Year)
as_tibble(origin)

# Data between zones ####
data.between<- origin %>%
  filter(as.numeric(Year) >= '2021')


# Set depth zones based on the light processing file ####
data.between<- mutate(data.between,
                      Zone = ifelse(as.numeric(as.character(Initial_depth)) >= 21 &
                                      as.numeric(as.character(Final_depth)) >= 21 |
                                      as.numeric(as.character(Initial_depth)) >= 21 &
                                      Final_depth == '>30' |
                                      Initial_depth == '>30' & Final_depth == '>30',
                                    'Mesophotic', 'Shallow'))

# Convert NAs to "shallow" zone
data.between$Zone[is.na(data.between$Zone)]<- 'Shallow'
which(is.na(data.between$Zone))

# Change zone data into factor type
data.between$Zone<- factor(data.between$Zone,
                           levels = c('Shallow', 'Mesophotic'),
                           ordered = T)

# Taxonomic classification
taxonomy<- read.csv('Data/Classification.csv', header = TRUE,
                    row.names = 1)

# Indices estimation ####

# Estimation of taxonomic distances 
taxdis<-taxa2dist(taxonomy, varstep = TRUE)
spp_per_transt<- pres_data[, c(4:76)]
rownames(spp_per_transt)<- pres_data[, 3]
mod<- taxondive(spp_per_transt, taxdis)
delta<- mod$D
delta[is.na(mod$D)]<- 0
round(mean(delta),2); round(sd(delta)/sqrt(length(delta)),2)

# Transectos de Espíritu Santo ####
transt_ES<- spp_per_transt[36:66, ]
mod<- taxondive(transt_ES, taxdis)
delta<- mod$D
delta[is.na(delta)]<- 0
round(mean(delta),2); round(sd(delta)/sqrt(length(delta)),2)

# Transectos de Cerralvo ####
transt_CE<- spp_per_transt[1:10, ]
mod<- taxondive(transt_CE, taxdis)
delta<- mod$D
delta[is.na(delta)]<- 0
round(mean(delta),2); round(sd(delta)/sqrt(length(delta)),2)

# Transectos de San Benedicto
transt_SB<- spp_per_transt[67:75, ]
mod<- taxondive(transt_SB, taxdis)
delta<- mod$D
delta[is.na(delta)]<- 0
round(mean(delta),2); round(sd(delta)/sqrt(length(delta)),2)

# Transectos de Socorro ####
transt_SOC<- spp_per_transt[76:78, ]
mod<- taxondive(transt_SOC, taxdis)
delta<- mod$D
delta[is.na(delta)]<- 0
round(mean(delta),2); round(sd(delta)/sqrt(length(delta)),2)

# Transectos de Clarión ####
transt_CL<- spp_per_transt[11:35, ]
mod<- taxondive(transt_CL, taxdis)
delta<- mod$D
delta[is.na(delta)]<- 0
round(mean(delta),2); round(sd(delta)/sqrt(length(delta)),2)

# Transectos islas continentales
transt_ic<- spp_per_transt[c(1:10, 36:66), ]
mod<- taxondive(transt_ic, taxdis)
delta<- mod$D
delta[is.na(delta)]<- 0
round(mean(delta),2); round(sd(delta)/sqrt(length(delta)),2)

# Transectos islas oceánicas
transt_io<- spp_per_transt[c(11:35, 67:78), ]
mod<- taxondive(transt_io, taxdis)
delta<- mod$D
delta[is.na(delta)]<- 0
round(mean(delta),2); round(sd(delta)/sqrt(length(delta)),2)
