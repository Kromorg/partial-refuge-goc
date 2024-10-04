# Open libraries y clean console ####
pacman:: p_load(tidyverse, # Data wrangling
                mFD) # Functional diversity

rm(list = ls())
shell('cls')

# Open database ####
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


# Functional diversity analysis ####
# Subset abundance data at each strata
abundance.shallow<- data.between %>% filter(Zone == 'Shallow') %>% 
  select(Abudefduf_troschelii:Zapteryx_exasperata) %>% droplevels()
rownames(abundance.shallow)<- data.between[data.between$Zone == 'Shallow', 6]

abundance.mesophotic<- data.between %>% filter(Zone == 'Mesophotic') %>% 
  select(Abudefduf_troschelii:Zapteryx_exasperata) %>% droplevels()
rownames(abundance.mesophotic)<- data.between[
  data.between$Zone == 'Mesophotic', 6]

# Load database of biological traits ####
traits<- read.csv('Data/Traits.csv',
                  header = T, stringsAsFactors = T,
                  row.names = 1)
#traits<- read.csv(here:: here('Data/Traits.csv'),
#                  header = T, stringsAsFactors = T,
#                  row.names = 1)

as_tibble(traits)

# Set ordinal traits as "ordinal" data in R language
traits$Size<- factor(traits$Size,
                     levels = c('1', '2', '3', '4', '5', '6'),
                     ordered = T)
traits$Mobility<- factor(traits$Mobility,
                         levels = c('1', '2', '3', '4'),
                         ordered = T)
traits$Gregariousness<- factor(traits$Gregariousness,
                               levels = c('1', '2', '3', '4'),
                               ordered = T)
traits$Position<- factor(traits$Position,
                         levels = c('1', '2', '3'),
                         ordered = T)
as_tibble(traits)


# Shallow zone ####
# Determine number of unique atributes per video-transect
FE <- matrix(nrow = nrow(abundance.shallow), ncol = 2)
a <- for (m in 1:nrow(abundance.shallow)){
  FE[m, ]<- dim(unique(traits[
    rownames(traits) %in%
      colnames(abundance.shallow[m, apply(abundance.shallow[m, ], 2, sum)> 0]),]))
}
head(FE)

# Consider those in which the number of species recorded was higher than
# the traits used (S> 6)
# If those records aren't deleted, the following error will eventually
# appear:

# 1) Number of species/FEs must be higher than the number of
# axes needed to compute the convex hull.

# 2) Assemblages that have no species.

nbFE <- FE[, 1] # Add the number of traits used (6)
sum(nbFE < 6) # 3 video-transects (S< 6)

sub.abundance.shallow<- cbind(abundance.shallow, nbFE) # Add number of unique FE
which(sub.abundance.shallow$nbFE < 6) # Rows with unique FE < 6

sub.abundance.shallow <- sub.abundance.shallow[!(sub.abundance.shallow$nbFE < 6), ] # Deleting the rows
sum(sub.abundance.shallow$nbFE < 6) # Number of video-transects S< 6

sub.abundance.shallow<- sub.abundance.shallow[, -104] # Delete row "nbFE"

which(colSums(sub.abundance.shallow) == 0) # Which species are absent due to deleting
sub.abundance.shallow<- sub.abundance.shallow[,-c(5, 11, 17, 37:38, 41,
                           48, 50, 57, 65:67, 71, 75, 83, 85, 88:89, 98,
                           99, 101, 103)] # Delete those species
which(colSums(sub.abundance.shallow) == 0) # No longer absent species


sub.traits.shallow<- traits[-c(5, 11, 17, 37:38, 41, 48, 50, 57, 65:67,
                               71, 75, 83, 85, 88:89, 98, 99, 101, 103), ]


# Load database with traits information ####
types<- read.csv('Data/Traits info.csv',
                 header = T, stringsAsFactors = T)
#types<- read.csv(here:: here('Data/Traits info.csv'),
#                 header = T, stringsAsFactors = T)

sp.tr.summary(tr_cat = types,
              sp_tr = sub.traits.shallow)

sp.to.fe(sub.traits.shallow, types)


# Mesophotic zone ####
# Determine number of unique atributes per video-transect
FE <- matrix(nrow = nrow(abundance.mesophotic), ncol = 2)
a <- for (m in 1:nrow(abundance.mesophotic)){
  FE[m, ]<- dim(unique(traits[
    rownames(traits) %in%
      colnames(abundance.mesophotic[m, apply(abundance.mesophotic[m, ], 2, sum)> 0]),]))
}
head(FE)

# Consider those in which the number of species recorded was
# higher than the traits used (S> 6)

nbFE <- FE[, 1] # Add the number of traits used (6)
sum(nbFE < 6) # 49 video-transects (S< 6)

sub.abundance.mesophotic<- cbind(abundance.mesophotic, nbFE) # Add number of unique FE
which(sub.abundance.mesophotic$nbFE < 6) # Rows with unique FE < 6

sub.abundance.mesophotic <- sub.abundance.mesophotic[!(sub.abundance.mesophotic$nbFE < 6), ] # Deleting the rows
sum(sub.abundance.mesophotic$nbFE < 6) # Number of video-transects S< 6

sub.abundance.mesophotic<- sub.abundance.mesophotic[, -104] # Delete row "nbFE"

which(colSums(sub.abundance.mesophotic) == 0) # Which species are absent due to deleting
sub.abundance.mesophotic<- sub.abundance.mesophotic[,-c(1, 3, 5, 16:17,
                             19:20, 25, 28, 30:31, 33:35, 41, 43:45, 48,
                             50, 53:55, 57, 63:64, 67:68, 74, 77:79,
                             80:83, 85, 89, 90:91, 93, 102)] # Delete those species
which(colSums(sub.abundance.mesophotic) == 0) # No longer absent species


sub.traits.mesophotic<- traits[-c(1, 3, 5, 16:17, 19:20, 25, 28, 30:31,
                                  33:35, 41, 43:45, 48, 50, 53:55, 57,
                                  63:64, 67:68, 74, 77:79, 80:83, 85,
                                  89, 90:91, 93, 102), ]

sp.tr.summary(tr_cat = types,
              sp_tr = sub.traits.mesophotic)

sp.to.fe(sub.traits.mesophotic, types)
