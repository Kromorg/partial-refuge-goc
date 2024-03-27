# Open libraries y clean console ####
pacman:: p_load(tidyverse, # Data wrangling
                ggplot2, # Graphs
                gridExtra, # Multiple plots
                mFD, elbow, # Functional diversity
                geometry, # Compute convex hull
                tripack) # Vertices triangulation (plot)

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

# Elbow method
var.expl<- qual$details_fspaces$pc_eigenvalues$Eigenvalues/
  sum(qual$details_fspaces$pc_eigenvalues$Eigenvalues)

round(sum(var.expl[1:5]), 2) # Quality using five axes
round(sum(var.expl[1:4]), 2) # Quality using four axes

inflection_points<- data.frame(
  'Dimensions' = 1:length(var.expl),
  'Explained variance' = cumsum(var.expl))
perf_PNZMAES<- elbow(
  data = inflection_points[, c('Dimensions',
                               'Explained.variance')])

# Functional indices estimation ####
# Species coordinates
coord<- qual$details_fspaces$sp_pc_coord

# Indices
indices<- alpha.fd.multidim(
  sp_faxes_coord = coord[,c('PC1', 'PC2', 'PC3', 'PC4')],
  asb_sp_w = as.matrix(sub.abundance),
  ind_vect = c('fric', 'fori', 'fdiv'),
  scaling = T, check_input = T, details_returned = T)

funct_diversity<- indices$functional_diversity_indices


# Statystical analysis ####
# Add factors to assess considering video-transects with less
# than six unique FE
Year<- data.between[data.between$Video.transect %in%
             rownames(funct_diversity), 2]
Zone<- data.between[data.between$Video.transect %in%
              rownames(funct_diversity), 110]
Season<-  data.between[data.between$Video.transect %in%
                    rownames(funct_diversity), 3]
Site<- data.between[data.between$Video.transect %in%
               rownames(funct_diversity), 1]

# Create an object including the indices and factors 
funct_diversity<- cbind(Year, Zone, Season, Site,
                        funct_diversity)

glimpse(funct_diversity)

# Summary per zone
indices.zones<- funct_diversity %>%
  group_by(Zone) %>%
  summarise(MeanFRic = round(mean(fric), 2),
            SD.FRic = round(sd(fric), 2),
            MeanFOri = round(mean(fori), 2),
            SD.FOri = round(sd(fori), 2),
            MeanFDiv = round(mean(fdiv), 2),
            SD.FDiv = round(sd(fdiv), 2))
as_tibble(indices.zones)

# Summary per zone and site
indices.sites<- funct_diversity %>%
  group_by(Zone, Site) %>%
  summarise(MeanFRic = round(mean(fric), 2),
            SD.FRic = round(sd(fric), 2),
            MeanFOri = round(mean(fori), 2),
            SD.FOri = round(sd(fori), 2),
            MeanFDiv = round(mean(fdiv), 2),
            SD.FDiv = round(sd(fdiv), 2))
as_tibble(indices.sites)

# Functional richness ####
# Validation of assumptions
# Normality
shapiro.test(funct_diversity$fric[funct_diversity$Zone == 'Shallow'])
shapiro.test(funct_diversity$fric[funct_diversity$Zone == 'Mesophotic'])

# Homoscedasticity
var(funct_diversity$fric[funct_diversity$Zone == 'Shallow'])/
  var(funct_diversity$fric[funct_diversity$Zone == 'Mesophotic'])

var.test(funct_diversity$fric[funct_diversity$Zone == 'Shallow'],
         funct_diversity$fric[funct_diversity$Zone == 'Mesophotic'],
         ratio = 1, alternative = 'g') # p< 0.001

# Statistical test
wilcox.test(funct_diversity$fric[funct_diversity$Zone == 'Shallow'],
            funct_diversity$fric[funct_diversity$Zone == 'Mesophotic'],
            paired = F, alternative = 'g', correct = T) # p< 0.001


# Functional originality ####
# Validation of assumptions
# Normality
shapiro.test(funct_diversity$fori[funct_diversity$Zone == 'Shallow'])
shapiro.test(funct_diversity$fori[funct_diversity$Zone == 'Mesophotic'])

# Homoscedasticity
var(funct_diversity$fori[funct_diversity$Zone == 'Shallow'])/
  var(funct_diversity$fori[funct_diversity$Zone == 'Mesophotic'])

var.test(funct_diversity$fori[funct_diversity$Zone == 'Shallow'],
         funct_diversity$fori[funct_diversity$Zone == 'Mesophotic'],
         ratio = 1, alternative = 'l') # p> 0.05

# Statistical test
t.test(funct_diversity$fori[funct_diversity$Zone == 'Shallow'],
       funct_diversity$fori[funct_diversity$Zone == 'Mesophotic'],
       paired = F, alternative = 'g', correct = T)


# Functional divergence ####
# Validation of assumptions
# Normality
shapiro.test(funct_diversity$fdiv[funct_diversity$Zone == 'Shallow'])
shapiro.test(funct_diversity$fdiv[funct_diversity$Zone == 'Mesophotic'])

# Homoscedasticity
var(funct_diversity$fdiv[funct_diversity$Zone == 'Shallow'])/
  var(funct_diversity$fdiv[funct_diversity$Zone == 'Mesophotic'])

var.test(funct_diversity$fdiv[funct_diversity$Zone == 'Shallow'],
         funct_diversity$fdiv[funct_diversity$Zone == 'Mesophotic'],
         ratio = 1, alternative = 'l') # p< 0.001

# Statistical test
wilcox.test(funct_diversity$fdiv[funct_diversity$Zone == 'Shallow'],
            funct_diversity$fdiv[funct_diversity$Zone == 'Mesophotic'],
            paired = F, alternative = 'g', correct = T) # p> 0.05

# Data visualization ####
# Functional richness graph
fric<- ggplot(funct_diversity,
              aes(x = Site, y = fric, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  ylim(0, 0.40)+ labs(y = 'FRic')+
  annotate('text', x = 1, y = 0.4, label = 'A)*',
           size = 15)+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30))

# Functional originality graph
fori<- ggplot(funct_diversity,
              aes(x = Site, y = fori, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  ylim(0, 0.40)+ labs(y = 'FOri')+
  stat_summary(fun = mean, geom = "point", col = "black",
               size = 9, fill = 'black', shape = 23,
               aes(group = Zone))+
  annotate('text', x = 1, y = 0.4, label = 'B)*',
           size = 15)+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30))

# Functional divergence graph
fdiv<- ggplot(funct_diversity,
              aes(x = Site, y = fdiv, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  ylim(0.60, 1)+ labs(y = 'FDiv')+
  annotate('text', x = 1, y = 1, label = 'C)',
           size = 15)+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30))

multi <- grid.arrange(fric, fori, fdiv,
             ncol = 3, nrow = 1)

ggsave('Figs/Figure 3.tiff', plot = multi,
       width = 6560, height = 3440, units = 'px', dpi = 320)
#ggsave(here:: here('Figs/Functional diversity between zones.tiff'),
#       plot = multi, width = 6560, height = 3440, units = 'px',
#       dpi = 320)


# Remove objects
rm(FE, nbFE, perf_PNZMAES, inflection_points, var.expl,
   Year, Zone, Season, Site, funct_diversity,
   indices.zones, indices.sites, fric, fori, fdiv)


# Convex hull ####
# Species coordinates
fd.coord<- coord[, c('PC1', 'PC2', 'PC3', 'PC4')]
# Load Species and FE data
spp_fes<- data.frame(FE = str_c(sub.traits$Size, '',
                                sub.traits$Mobility, '',
                                sub.traits$Activity, '',
                                sub.traits$Gregariousness, '',
                                sub.traits$Position, '',
                                sub.traits$Diet))
rownames(spp_fes)<- rownames(sub.traits)
as_tibble(spp_fes)

# Avoiding repeated FE and changing row names from
# fd.coord object
rownames(fd.coord) = spp_fes[, ]
fd.coord<- fd.coord[!duplicated(rownames(fd.coord)), ]

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


FSpace <- lapply(zone, function (x) {
  
  extant <- colnames(presence.conditions)[which(
    presence.conditions[x, ] > 0)]
  
  fes_cond <- spp_fes[rownames(spp_fes) %in% extant, ]
  
  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
  
  ch <- convhulln(m, options = "FA")
  
  chg <- convhulln(fd.coord, options = "FA")
  
  c(length(extant), length(extant)/93*100, dim(m)[1], 
    dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
  # Accounting for 93 species in total
})

names(FSpace) <- zone


# FSpace contains the number of species(NbSp) and FEs (NbFEs),
# relative percentages (NbSpP, NbFEsP), and the volume among
# shallow and mesophotic zones

FSpace <- do.call(rbind, FSpace)
colnames(FSpace) <- c("NbSp", "NbSpP", "NbFEs",
                      "NbFEsP", "Vol")


# Plot covex hull ####
cols <- c('darkred', 'royalblue4')
names(cols) <- zone

# Species and functional diversity changes among zones
# All volumes in distinct plots

labels_fig_cv <- c("No. spp.", "No. FE", "Vol.")

tiff(filename = "Figs/Figure 2.tiff",
     height = 7500, width = 5700,
     antialias = 'cleartype', res = 320, pointsize = 36)
#tiff(filename = here:: here("Figs/Figure 2.tiff"),
#     height = 7500, width = 5700,
#     antialias = 'cleartype', res = 320, pointsize = 36)

# x11()
par(mfrow = c(3, 2))

for (i in zone) {
  
  midpoints <-  barplot(FSpace[i, c(2, 4, 5)],
                        ylab= "Relative richness (%)",
                        names.arg = labels_fig_cv,
                        col = cols[i], ylim = c(0, 105),
                        main = paste0(i), 
                        col.main = 'black', cex.main = 1.5,
                        cex.names = 1.3, cex.lab = 1.3) 
  
  text(midpoints[c(1, 2), ],
       FSpace[i, c(2,4)] + 8,
       labels = round(FSpace[i, c(1,3)]), 
       col= 'black',
       cex = 1.3,
       font = 2)
}

for (i in zone) {
  
  extant <- colnames(presence.conditions)[which(
    presence.conditions[i,] > 0)]
  
  fes_cond <- spp_fes[rownames(spp_fes) %in% extant, ]
  
  m <- fd.coord[rownames(fd.coord) %in% fes_cond, ]
  
  tr_all<- tri.mesh(fd.coord[, 1], fd.coord[, 2],
                    duplicate = 'remove')
  ch_all<- convex.hull(tr_all)
  tr <-tri.mesh(m[, 1], m[, 2], duplicate = 'remove')
  ch <- convex.hull(tr)
  
  plot(fd.coord[, 1], fd.coord[, 2],
       xlab = "PCoA 1 (29%)", ylab = "PCoA 2 (22%)",
       type="n", cex.lab = 1.3, cex.axis = 1.1)
  polygon(ch_all, col = 'white', border = 'black', lwd = 2)
  polygon(ch, col = 'gainsboro', border = 'black')
  points(m[, 1:2], pch = 16, col = cols[i])
  
}

for (i in zone) {
  
  extant <- colnames(presence.conditions)[which(
    presence.conditions[i,] > 0)]
  
  fes_cond <- spp_fes[rownames(spp_fes) %in% extant, ]
  
  m <- fd.coord[rownames(fd.coord) %in% fes_cond, ]
  
  tr_all<- tri.mesh(fd.coord[, 3], fd.coord[, 4],
                    duplicate = 'remove')
  ch_all<- convex.hull(tr_all)
  tr <-tri.mesh(m[, 3], m[, 4], duplicate = 'remove')
  ch <- convex.hull(tr)
  
  plot(fd.coord[, 3], fd.coord[, 4],
       xlab = 'PCoA 3 (13%)', ylab = 'PCoA 4 (11%)',
       type = 'n', cex.lab = 1.3, cex.axis = 1.1)
  polygon(ch_all, col = 'white', border = 'black', lwd = 2)
  polygon(ch, col = 'gainsboro', border = 'black')
  points(m[, 3:4], pch = 16, col = cols[i])
}

dev.off()

# Functional divergence comparison graph ####
fdiv.plot<- alpha.fd.multidim(
  sp_faxes_coord = coord[, c('PC1', 'PC2', 'PC3', 'PC4')],
  asb_sp_w = presence.conditions, ind_vect = 'fdiv',
  scaling = T, check_input = T, details_returned = T) %>%
  alpha.multidim.plot(output_alpha_fd_multidim = .,
    plot_asb_nm = c('Shallow', 'Mesophotic'),
    ind_nm = 'fdiv', faxes = NULL, faxes_nm = NULL,
    color_bg = '#FFFFFF',
    # Convex parameters
    color_ch = c(pool = '#000000', asb1 = '#8B0000',
                 asb2 =  '#27408B'),
    fill_ch = c(pool = '#DCDCDC', asb1 = '#8B0000',
                asb2 =  '#27408B'),
    alpha_ch = c(pool = 0.8, asb1 = 1, asb2 = 1),
    # Vertices parameters
    color_vert = c(pool = '#000000', asb1 = '#8B0000',
                   asb2 =  '#27408B'),
    fill_vert = c(pool = NA, asb1 = '#8B0000',
                  asb2 =  '#27408B'))
fdiv.plot

ggsave('Figs/Appendix S3.tiff', width = 4280, height = 3536,
       units = 'px',  dpi = 320)
#ggsave(here:: here('Figs/FEs distribution.tiff'), width = 4280, height = 3536,
#       units = 'px',  dpi = 320)


# Remove objects
rm(coord, spp_fes, transects, zone, FSpace, cols, labels_fig_cv,
   m, midpoints, ch, ch_all, tr, tr_all, a, extant, fes_cond, i,
   fd.coord)

# Hill numbers ####

# According to Loiseau et al. (2022): it is possible to compute
# taxonomic alpha and beta with beta.fd.hill but functional
# distance must be higher than 0 to avoid transformation of
# species in FE (reducing the real value of alpha)

# Taxonomic diversity ####

dist.taxo <- sub.gower
dist.taxo[dist.taxo == 0]<- 0.0001

# Species richness
hill.tax.q0<- alpha.fd.hill(asb_sp_w = as.matrix(sub.abundance),
                            sp_dist = dist.taxo,
                            q = 0,
                            tau = 'min')$asb_FD_Hill

# Entropy (Shannon index)
hill.tax.q1<- alpha.fd.hill(asb_sp_w = as.matrix(sub.abundance),
                            sp_dist = dist.taxo,
                            q = 1,
                            tau = 'min')$asb_FD_Hill

# Functional diversity ####

# Effective number of functional groups
hill.fd.q0<- alpha.fd.hill(asb_sp_w = as.matrix(sub.abundance),
              sp_dist = sub.gower,
              q = 0,
              tau = 'mean')$asb_FD_Hill

# Entropy
hill.fd.q1<- alpha.fd.hill(asb_sp_w = as.matrix(sub.abundance),
              sp_dist = sub.gower,
              q = 1,
              tau = 'mean')$asb_FD_Hill

# Merge ---- Fallo por desigualdad de filas
hill.div <- data.frame(tax.rich = hill.tax.q0[, 1],
                       tax.entro = hill.tax.q1[, 1],
                       fd.rich = hill.fd.q0[, 1],
                       fd.entro = hill.fd.q1[, 1])
hill.div


# Statystical analysis ####
Zone<- data.between[data.between$Video.transect %in%
                      rownames(hill.div), 110]
Site<- data.between[data.between$Video.transect %in%
                      rownames(hill.div), 1]

# Create an object including the indices and factors 
hill.div<- cbind(Zone, Site, hill.div)

glimpse(hill.div)

# Summary per zone
hill.zones<- hill.div %>%
  group_by(Zone) %>%
  summarise(Taxon_richness = round(mean(tax.rich), 2),
            TaxRich_SD = round(sd(tax.rich), 2),
            Taxon_entropy = round(mean(tax.entro), 2),
            TaxEntro_SD = round(sd(tax.entro), 2),
            Funct_richness = round(mean(fd.rich), 2),
            FuncRich_SD = round(sd(fd.rich), 2),
            Funct_entropy = round(mean(fd.entro), 2),
            FuncEntro_SD = round(sd(fd.entro), 2))

# Summary per site
hill.sites<- hill.div %>%
  group_by(Zone, Site) %>%
  summarise(Taxon_richness = round(mean(tax.rich), 2),
            TaxRich_SD = round(sd(tax.rich), 2),
            Taxon_entropy = round(mean(tax.entro), 2),
            TaxEntro_SD = round(sd(tax.entro), 2),
            Funct_richness = round(mean(fd.rich), 2),
            FuncRich_SD = round(sd(fd.rich), 2),
            Funct_entropy = round(mean(fd.entro), 2),
            FuncEntro_SD = round(sd(fd.entro), 2))

# Taxonomic diversity ####

# Richness
# Validation of assumptions
# Normality
shapiro.test(hill.div$tax.rich[hill.div$Zone == 'Shallow'])
shapiro.test(hill.div$tax.rich[hill.div$Zone == 'Mesophotic'])

# Homoscedasticity
var(hill.div$tax.rich[hill.div$Zone == 'Shallow'])/
  var(hill.div$tax.rich[hill.div$Zone == 'Mesophotic'])

var.test(hill.div$tax.rich[hill.div$Zone == 'Shallow'],
         hill.div$tax.rich[hill.div$Zone == 'Mesophotic'],
         ratio = 1, alternative = 'g') # p< 0.001

# Statistical test
wilcox.test(hill.div$tax.rich[hill.div$Zone == 'Shallow'],
            hill.div$tax.rich[hill.div$Zone == 'Mesophotic'],
            paired = F, alternative = 'g', correct = T) # p< 0.001

# Entropy
# Validation of assumptions
# Normality
shapiro.test(hill.div$tax.entro[hill.div$Zone == 'Shallow'])
shapiro.test(hill.div$tax.entro[hill.div$Zone == 'Mesophotic'])

# Homoscedasticity
var(hill.div$tax.entro[hill.div$Zone == 'Shallow'])/
  var(hill.div$tax.entro[hill.div$Zone == 'Mesophotic'])

var.test(hill.div$tax.entro[hill.div$Zone == 'Shallow'],
         hill.div$tax.entro[hill.div$Zone == 'Mesophotic'],
         ratio = 1, alternative = 'g') # p< 0.005

# Statistical test
wilcox.test(hill.div$tax.entro[hill.div$Zone == 'Shallow'],
            hill.div$tax.entro[hill.div$Zone == 'Mesophotic'],
            paired = F, alternative = 'g', correct = T) # p< 0.001

# Functional diversity ####
# Richness
# Validation of assumptions
# Normality
shapiro.test(hill.div$fd.rich[hill.div$Zone == 'Shallow'])
shapiro.test(hill.div$fd.rich[hill.div$Zone == 'Mesophotic'])

# Homoscedasticity
var(hill.div$fd.rich[hill.div$Zone == 'Shallow'])/
  var(hill.div$fd.rich[hill.div$Zone == 'Mesophotic'])

var.test(hill.div$fd.rich[hill.div$Zone == 'Shallow'],
         hill.div$fd.rich[hill.div$Zone == 'Mesophotic'],
         ratio = 1, alternative = 'l') # p> 0.05

# Statistical test
t.test(hill.div$fd.rich[hill.div$Zone == 'Shallow'],
       hill.div$fd.rich[hill.div$Zone == 'Mesophotic'],
       paired = F, alternative = 'g', correct = T) # p> 0.05

# Entropy
# Validation of assumptions
# Normality
shapiro.test(hill.div$fd.entro[hill.div$Zone == 'Shallow'])
shapiro.test(hill.div$fd.entro[hill.div$Zone == 'Mesophotic'])

# Homoscedasticity
var(hill.div$fd.entro[hill.div$Zone == 'Shallow'])/
  var(hill.div$fd.entro[hill.div$Zone == 'Mesophotic'])

var.test(hill.div$fd.entro[hill.div$Zone == 'Shallow'],
         hill.div$fd.entro[hill.div$Zone == 'Mesophotic'],
         ratio = 1, alternative = 'g') # p> 0.05

# Statistical test
t.test(hill.div$fd.entro[hill.div$Zone == 'Shallow'],
       hill.div$fd.entro[hill.div$Zone == 'Mesophotic'],
       paired = F, alternative = 'g', correct = T) # p> 0.05

# Data visualization ####
# Taxonomic diversity graph
rich<- ggplot(hill.div,
              aes(x = Site, y = tax.rich, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  scale_y_continuous(limits = c(5, 30), breaks = seq(5, 30, 5))+
  labs(y = 'Richness')+
  labs(title = expression(Taxonomic~diversity))+
  annotate('text', x = 1, y = 30, label = 'A)*',
           size = 15)+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Taxonomic entropy graph
entro<- ggplot(hill.div,
              aes(x = Site, y = tax.entro, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  scale_y_continuous(limits = c(1, 13), breaks = seq(1, 13, 3))+
  labs(y = 'Entropy')+
  annotate('text', x = 1, y = 13, label = 'C)*',
           size = 15)+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30))

# Functional diversity
fd.rich<- ggplot(hill.div,
                 aes(x = Site, y = fd.rich, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  scale_y_continuous(limits = c(1, 6), breaks = seq(1, 6, 1))+
  labs(title = expression(Functional~diversity))+
  stat_summary(fun = mean, geom = "point", col = "black",
               size = 9, fill = 'black', shape = 23,
               aes(group = Zone))+
  annotate('text', x = 1, y = 6, label = 'B)',
           size = 15)+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

fd.entro<- ggplot(hill.div,
                 aes(x = Site, y = fd.entro, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  ylim(1, 5)+
  stat_summary(fun = mean, geom = "point", col = "black",
               size = 9, fill = 'black', shape = 23,
               aes(group = Zone))+
  annotate('text', x = 1, y = 5, label = 'D)',
           size = 15)+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30),
        axis.title.y = element_blank())

hill <- grid.arrange(rich, fd.rich, entro, fd.entro,
             ncol = 2, nrow = 2)

ggsave('Figs/Figure 4.tiff', plot = hill,
       width = 5500, height = 4700, units = 'px', dpi = 320)
#ggsave(here:: here('Figs/Figure 4.tiff'), plot = hill,
#       width = 5500, height = 4700, units = 'px', dpi = 320)

# Remove objects
rm(entro, fd.entro, fd.rich, fdiv.plot, hill.div, hill.fd.q0,
   hill.fd.q1, hill.sites, hill.tax.q0, hill.tax.q1, hill.zones,
   dist.taxo, Site, sub.gower, Zone)

