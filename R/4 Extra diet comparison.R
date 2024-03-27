# Open libraries y clean console ####
pacman:: p_load(tidyverse, # Data wrangling
                ggplot2, # Graphs
                gridExtra, # Multiple plots
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
# Subset abundance data
abundance<- data.between[, c(7:109)]
rownames(abundance)<- data.between[, 6]

# Abundance summary
asb.sp.summary(as.matrix(abundance))

# Warning due to video-transects with no records (S= 0)

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

# Determine number of unique atributes per video-transect
FE <- matrix(nrow = nrow(abundance), ncol = 2)
a <- for (m in 1:nrow(abundance)){
  FE[m, ]<- dim(unique(traits[
    rownames(traits) %in%
      colnames(abundance[m, apply(abundance[m, ], 2, sum)> 0]),]))
}
head(FE)

# Consider those in which the number of species recorded was
# higher than the traits used (S> 6)
# If those records aren't deleted, the following error will 
# eventually appear:

# 1) Number of species/FEs must be higher than the number of
# axes needed to compute the convex hull.

# 2) Assemblages that have no species.

nbFE <- FE[, 1] # Add the number of traits used (6)
sum(nbFE < 6) # 52 video-transects (S< 6)

sub.abundance<- cbind(abundance, nbFE) # Add number of unique FE
which(sub.abundance$nbFE < 6) # Rows with unique FE < 6

sub.abundance <- sub.abundance[!(sub.abundance$nbFE < 6), ] # Deleting the rows
sum(sub.abundance$nbFE < 6) # Number of video-transects S< 6

sub.abundance<- sub.abundance[, -104] # Delete row "nbFE"

which(colSums(sub.abundance) == 0) # Which species are absent due to deleting
sub.abundance<- sub.abundance[,-c(5, 17, 41, 48, 50, 57, 67, 83,
                                  85, 89)] # Delete those species
which(colSums(sub.abundance) == 0) # No longer absent species


# We also delete the species Hypanus dipterurus, 
# Seriola lalandi, and Scorpaena guttata in the "trait"
# database. The rest are already gone
which(rownames(traits) == c('Hypanus_dipterurus'))
which(rownames(traits) == c('Scorpaena_guttata'))
which(rownames(traits) == c('Seriola_lalandi'))
sub.traits<- traits[-c(39, 77, 79), ]


# Load database with traits information ####
types<- read.csv('Data/Traits info.csv',
                 header = T, stringsAsFactors = T)
#types<- read.csv(here:: here('Data/Traits info.csv'),
#                 header = T, stringsAsFactors = T)

sp.tr.summary(tr_cat = types,
              sp_tr = sub.traits)


# Gower distance estimation ####
sub.gower<- funct.dist(sp_tr = sub.traits,
                       tr_cat = types,
                       metric = 'gower')

# Warning: distance between speces is equal to 0
# due to similar FE

# Functional space quality estimation ####
qual<- quality.fspaces(sp_dist = sub.gower,
                       maxdim_pcoa = 16, # Maximum
                       deviation_weighting = 'absolute')

# Convex hull ####
# Subset site and zone information
transects<- data.between[, c(1, 110)]
rownames(transects)<- data.between[, 6]
zone<- c('Shallow', 'Mesophotic')

presence.conditions <- lapply(zone, function(x) {
  
  vid.trans <- rownames(transects[transects$Zone == x,])
  
  colSums(sub.abundance[rownames(sub.abundance) %in% vid.trans, ])
  
})#eo lapply

presence.conditions <- do.call(rbind, presence.conditions)

rownames(presence.conditions) <- zone


# Diet comparison ####
# Species coordinates
fd.coord<- qual$details_fspaces$sp_pc_coord %>%
  .[, c('PC1', 'PC2', 'PC3', 'PC4')]

trait.space<- data.frame(fd.coord, sub.traits) %>%
  data.frame(., t(presence.conditions))
trait.space<- trait.space[, -c(5:9)]

# Shallow reefs
sh.herb.detr <- trait.space[trait.space$Diet == 'HD' &
                              trait.space$Shallow >= 1, ]
sh.sess.inv <- trait.space[trait.space$Diet == 'IS' &
                             trait.space$Shallow >= 1, ]
sh.mob.inv <- trait.space[trait.space$Diet == 'IM' &
                            trait.space$Shallow >= 1, ]
sh.plank <- trait.space[trait.space$Diet == 'Pk' &
                          trait.space$Shallow >= 1, ]
sh.pisc <- trait.space[trait.space$Diet == 'FC' &
                         trait.space$Shallow >= 1, ]
sh.omni <- trait.space[trait.space$Diet == 'OM' &
                         trait.space$Shallow >= 1, ]

# Mesophotic reefs
mes.herb.detr <- trait.space[trait.space$Diet == 'HD' &
                               trait.space$Mesophotic >= 1, ]
mes.sess.inv <- trait.space[trait.space$Diet == 'IS' &
                              trait.space$Mesophotic >= 1, ]
mes.mob.inv <- trait.space[trait.space$Diet == 'IM' &
                             trait.space$Mesophotic >= 1, ]
mes.plank <- trait.space[trait.space$Diet == 'Pk' &
                           trait.space$Mesophotic >= 1, ]
mes.pisc <- trait.space[trait.space$Diet == 'FC' &
                          trait.space$Mesophotic >= 1, ]
mes.omni <- trait.space[trait.space$Diet == 'OM' &
                          trait.space$Mesophotic >= 1, ]


# Convex hull based on the axes 1 and 2 ####
find_hull <- function(x)x[chull(x$PC1,x$PC2),]
hulls_all <- find_hull(trait.space)

# Shallow reefs
hull.sh.hd <- find_hull(sh.herb.detr)
hull.sh.is <- find_hull(sh.sess.inv)
hull.sh.im <- find_hull(sh.mob.inv)
hull.sh.pk <- find_hull(sh.plank)
hull.sh.fc <- find_hull(sh.pisc)
hull.sh.om <- find_hull(sh.omni)
shallow.hulls<- rbind(hull.sh.hd, hull.sh.is, hull.sh.im,
                      hull.sh.pk, hull.sh.fc, hull.sh.om)
rm(sh.herb.detr, sh.sess.inv, sh.mob.inv, sh.plank, sh.pisc,
   sh.omni, hull.sh.hd, hull.sh.is, hull.sh.im, hull.sh.pk,
   hull.sh.fc, hull.sh.om)

# Mesophotic reefs
hull.mes.hd <- find_hull(mes.herb.detr)
hull.mes.is <- find_hull(mes.sess.inv)
hull.mes.im <- find_hull(mes.mob.inv)
hull.mes.pk <- find_hull(mes.plank)
hull.mes.fc <- find_hull(mes.pisc)
hull.mes.om <- find_hull(mes.omni)
meso.hulls<- rbind(hull.mes.hd, hull.mes.is, hull.mes.im,
                   hull.mes.pk, hull.mes.fc, hull.mes.om)
rm(mes.herb.detr, mes.sess.inv, mes.mob.inv, mes.plank, mes.pisc,
   mes.omni, hull.mes.hd, hull.mes.is, hull.mes.im, hull.mes.pk,
   hull.mes.fc, hull.mes.om)

# Graphics ####
# Shallow reefs
sh.diets<- ggplot(trait.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = hulls_all, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = 'black')+
  geom_point(data = trait.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = 'grey')+
  geom_polygon(data = shallow.hulls, aes (x = PC1, y = PC2,
                                          fill = Diet, colour = Diet),
               linetype = 1, linewidth = 1, alpha = 0.3)+
  scale_color_manual(values = c('HD' = 'greenyellow',
                                'IS' = 'skyblue2',
                                'IM' = 'gold',
                                'Pk' = 'mediumorchid',
                                'FC' = 'coral2',
                                'OM' = 'chocolate1'))+
  scale_fill_manual(values = c('HD' = 'greenyellow',
                               'IS' = 'skyblue2',
                               'IM' = 'gold',
                               'Pk' = 'mediumorchid',
                               'FC' = 'coral2',
                               'OM' = 'chocolate1'))+
  labs(subtitle = expression(Shallow~diet~traits),
       x = NULL,
       y = expression(PCoA~2~('22.12%')))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 15))

# Mesophotic reefs
mes.diets<- ggplot(trait.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = hulls_all, aes (x = PC1, y = PC2),
               linetype = 1, size = 1, alpha = 0, col = "black")+
  geom_point(data = trait.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = meso.hulls, aes (x = PC1, y = PC2,
                                       fill = Diet, colour = Diet),
               linetype = 1, linewidth = 1, alpha = 0.3)+
  scale_color_manual(values = c('HD' = 'springgreen4',
                                'IS' = 'steelblue4',
                                'IM' = 'goldenrod3',
                                'Pk' = 'lightseagreen',
                                'FC' = 'firebrick',
                                'OM' = 'darkorange3'))+
  scale_fill_manual(values = c('HD' = 'springgreen4',
                               'IS' = 'steelblue4',
                               'IM' = 'goldenrod3',
                               'Pk' = 'lightseagreen',
                               'FC' = 'firebrick',
                               'OM' = 'darkorange3'))+
  labs(subtitle = expression(Mesophotic~diet~traits),
       x = expression(PCoA~1~('29.45%')),
       y = expression(PCoA~2~('22.12%')))+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size = 15))

# All in one graph
x11()
grid.arrange(sh.diets, mes.diets, nrow = 2)
