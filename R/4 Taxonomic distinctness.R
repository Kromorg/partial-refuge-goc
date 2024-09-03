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
origin$Year <- as.factor(origin$Year)
as_tibble(origin)

# Data between zones ####
data.between <- origin %>%
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
data.between$Zone <- factor(data.between$Zone,
                           levels = c('Shallow', 'Mesophotic'),
                           ordered = T)

# Taxonomic classification
taxonomy <- read.csv('Data/Classification.csv', header = TRUE,
                    row.names = 1)


# Indices estimation ####
# Select data from shallow and mesophotic reefs ####
# Shallow reefs
data.shallow <- data.between %>% filter(Zone == 'Shallow')%>% 
  select(Abudefduf_troschelii:Zapteryx_exasperata) %>% droplevels()
rownames(data.shallow)<-  data.between[data.between$Zone == 'Shallow', 6]

# Looking for absent species
which(colSums(data.shallow) == 0)
data.shallow <- data.shallow[, -c(5, 11, 17, 37:38, 41, 48, 50, 57,
                                  65:67, 71, 75, 83, 85, 88:89, 98:99,
                                  101, 103)]

# Mesophotic reefs
data.meso <- data.between %>% filter(Zone == 'Mesophotic')%>% 
  select(Abudefduf_troschelii:Zapteryx_exasperata) %>% droplevels()
rownames(data.meso)<- data.between[data.between$Zone == 'Mesophotic', 6]

# Looking for absent species
which(colSums(data.meso) == 0)
data.meso <- data.meso[, -c(1, 3, 5, 16:17, 19:20, 25, 30:31, 33:35,
                            43, 45, 48, 50, 53:55, 57, 63:64, 67:68,
                            74, 77:79, 80:82, 89, 91, 93, 102)]

# Estimation of taxonomic distances
# Shallow reefs
taxdis.shallow <- taxonomy[rownames(taxonomy) %in%
                             colnames(data.shallow), ] %>% 
  taxa2dist(., varstep = TRUE)

# Mesophotic reefs
taxdis.meso <- taxonomy[rownames(taxonomy) %in%
                          colnames(data.meso), ] %>% 
  taxa2dist(., varstep = TRUE)

# -------------------------------------------------------------- ####
# "varstep" argument determines whether the path length between two
# randomly chosen species levels in the taxonomic classification will
# be separated equally or proportionally as taxon richness decreases
# at each step until the linking division.
# -------------------------------------------------------------- ####

# Estimation of taxonomic distinctness

# Shallow reefs
tax.shallow <- taxondive(data.shallow, taxdis.shallow)
delta.shallow<- tax.shallow$Dstar
round(mean(delta.shallow),2); round(sd(delta.shallow)/sqrt(length(delta.shallow)),2)

# Mesophotic reefs
tax.meso <- taxondive(data.meso, taxdis.meso)
delta.meso<- tax.meso$Dstar
delta.meso[is.na(tax.meso$Dstar)]<- 0
round(mean(delta.meso),2); round(sd(delta.meso)/sqrt(length(delta.meso)),2)
