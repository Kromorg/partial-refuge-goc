# Open libraries y clean console ####
pacman:: p_load(tidyverse, # Data wrangling
                mvabund, lme4, # Multivariate analysis and GLMMs
                labdsv, # Convert wide to long format
                vegan, # Taxonomic distinctness
                mFD, elbow, # Functional diversity
                geometry, # Compute convex hull
                tripack, # Vertices triangulation (plot)
                ggResidpanel) # Model assumption graphs

rm(list = ls())
shell('cls') # For Windows users
system2('clear') # For Mac users


# Open database ####
origin <- read.csv('Data/Abundance data.csv',
                   header = T, stringsAsFactors = T)
#origin <- read.csv(here:: here('Data/Abundance data.csv'),
#                   header = T, stringsAsFactors = T)

as_tibble(origin)

# Set variable "Year" to factor
origin$Year<- as.factor(origin$Year)
as_tibble(origin)


# Set depth zones based on the light processing file ####
data.between<- mutate(origin,
    Zone = ifelse(as.numeric(as.character(Initial_depth)) >= 21 &
              as.numeric(as.character(Final_depth)) >= 21 |
              as.numeric(as.character(Initial_depth)) >= 21 &
              Final_depth == '>30' |
              Initial_depth == '>30' & Final_depth == '>30',
              'Mesophotic', 'Shallow')) %>% droplevels()

# Convert NAs to "shallow" zone
data.between$Zone[is.na(data.between$Zone)]<- 'Shallow'
which(is.na(data.between$Zone))

# Change zone data into factor type
data.between$Zone<- factor(data.between$Zone,
                           levels = c('Shallow', 'Mesophotic'),
                           ordered = T)


# Multivariate analysis ####
# Creating objects to run the functions ####

# Subset abundance data
abundance <- data.between[, c(7:109)]
rownames(abundance) <- data.between[, 6]

# Testing if environmental variables have an effect on the abundance
# of fish species
mv.data <- data.between %>% select(Zone, Year, Season, Site) %>%
  droplevels() %>% as.list()
mv.data$Fish <- abundance %>%  as.matrix()

# Multivariate model
species.mod <- manyglm(Fish ~ Zone*Year*Season*Site,
                       data = mv.data, family = 'negative.binomial')

# Model visual assumptions
plot(species.mod)

# Model results

# ------------------------------------------------------- Estimation
# Results and summary having been saved as a .RDS file because
# it taks a lot of computation time. If you decide to rerun these 
# line codes, be patient (> 3 hrs) and continue this section.
# Otherwise, skip to "Results" below.
results.mod <- anova.manyglm(species.mod, p.uni = 'none')
saveRDS(results.mod, file = 'Data/Multivariate_results.RDS')

summary.mod <- summary.manyglm(species.mod)
saveRDS(summary.mod, file = 'Data/Multivariate_summary.RDS')

# ------------------------------------------------------- Results
results.mod <- readRDS('Data/Multivariate_results.RDS')
results.mod

summary.mod <- readRDS('Data/Multivariate_summary.RDS')
summary.mod

# Create base object to save diversity indices ####
diversity.metrics <- data.between %>% 
  select(Year, Season, Site, Zone, Video.transect)

long.abund <- dematrify(abundance) %>% 
  rename(Video.transect = sample, Species = species,
         Abundance = abundance) %>% group_by(Video.transect) %>% 
  reframe(Richness = length(unique(Species)),
          Abundance = sum(Abundance))

# Merging data
diversity.metrics <- diversity.metrics %>%
  left_join(long.abund, by = 'Video.transect')
diversity.metrics[is.na(diversity.metrics)]<- 0

# Remove objects
rm(mv.data, species.mod, long.abund, results.mod, summary.mod)


# Taxonomic distinctness ####
# Taxonomic classification
taxonomy <- read.csv('Data/Classification.csv', header = TRUE,
                     row.names = 1)
#taxonomy <- read.csv(here:: here('Data/Classification.csv'),
#header = TRUE, row.names = 1)

# Select data from shallow and mesophotic reefs
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
delta.shallow<- tax.shallow$Dstar %>% as.data.frame()

# Mesophotic reefs
tax.meso <- taxondive(data.meso, taxdis.meso)
delta.meso <- tax.meso$Dstar %>% as.data.frame()
delta.meso[is.na(delta.meso)] <- 0.0001 # To avoid "non-positive" values in glmm

# Merging data
delta.data <- rbind(delta.shallow, delta.meso) %>%
  rownames_to_column(., 'Video.transect') %>%
  rename(Distinctness = ".")
diversity.metrics <- diversity.metrics %>%
  left_join(delta.data, by = 'Video.transect')

# Remove objects
rm(taxonomy, data.shallow, taxdis.shallow, data.meso, taxdis.meso,
   delta.shallow, delta.meso, delta.data, tax.shallow, tax.meso)


# Functional diversity analysis ####
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

funct.diversity<- indices$functional_diversity_indices

# Remove objects
rm(FE, nbFE, perf_PNZMAES, inflection_points, var.expl, indices)


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

tiff(filename = "Figs/Figure 2.tiff", height = 7500, width = 5700,
     compression = "lzw", antialias = 'cleartype', res = 320,
     pointsize = 36)
#tiff(filename = here:: here("Figs/Figure 2.tiff"), height = 7500,
#     width = 5700, compression = "lzw", antialias = 'cleartype',
#     res = 320, pointsize = 36)

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
       compression = "lzw", units = 'px',  dpi = 320)
#ggsave(here:: here('Figs/FEs distribution.tiff'), width = 4280,
#        height = 3536, compression = "lzw", units = 'px',  dpi = 320)


# Remove objects
rm(coord, spp_fes, transects, zone, FSpace, cols, labels_fig_cv,
   m, midpoints, ch, ch_all, tr, tr_all, a, extant, fes_cond, i,
   fd.coord, fdiv.plot)

# Hill numbers ####

# According to Loiseau et al. (2022): it is possible to compute
# taxonomic alpha and beta with beta.fd.hill but functional
# distance must be higher than 0 to avoid transformation of
# species in FE (reducing the real value of alpha)

# Taxonomic diversity ####

dist.taxo <- sub.gower
dist.taxo[dist.taxo == 0]<- 0.0001


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

# Merging data
distance.div <- data.frame(funct.diversity[2:4], hill.tax.q1,
                       hill.fd.q0, hill.fd.q1) %>%
  rownames_to_column(., 'Video.transect') %>%
  rename(FRic = fric, FOri = fori, FDiv = fdiv,
         Tax.Entro = FD_q1, FD.Rich = FD_q0, FD.Entro = FD_q1.1)

diversity.metrics <- diversity.metrics %>%
  left_join(distance.div, by = 'Video.transect')

# Particularly for these data, we keep NAs because the value is not 0.

# Saving diversity indices in a new file
write.csv(diversity.metrics, 'Data/Diversity_indices.csv',
          row.names = F)

# Summary functional diversity per zone
func.indices.zones<- diversity.metrics %>%
  group_by(Zone) %>% na.omit() %>% 
  summarise(MeanFRic = round(mean(FRic), 2),
            SD.FRic = round(sd(FRic), 2),
            MeanFOri = round(mean(FOri), 2),
            SD.FOri = round(sd(FOri), 2),
            MeanFDiv = round(mean(FDiv), 2),
            SD.FDiv = round(sd(FDiv), 2))
as_tibble(func.indices.zones)

# Summary functional diversity per zone and site
func.indices.site <- diversity.metrics %>%
  group_by(Zone, Site) %>% na.omit() %>% 
  summarise(MeanFRic = round(mean(FRic), 2),
            SD.FRic = round(sd(FRic), 2),
            MeanFOri = round(mean(FOri), 2),
            SD.FOri = round(sd(FOri), 2),
            MeanFDiv = round(mean(FDiv), 2),
            SD.FDiv = round(sd(FDiv), 2))
as_tibble(func.indices.site)

# Summary Hill numbers per zone
hill.indices.zones<- diversity.metrics %>%
  group_by(Zone) %>% na.omit() %>% 
  summarise(Taxq0 = round(mean(Richness), 2),
            Taxq0.SD = round(sd(Richness), 2),
            Taxq1 = round(mean(Tax.Entro), 2),
            Taxq1.SD = round(sd(Tax.Entro), 2),
            Functq0 = round(mean(FD.Rich), 2),
            Functq0.SD = round(sd(FD.Rich), 2),
            Functq1 = round(mean(FD.Entro), 2),
            Functq1.SD = round(sd(FD.Entro), 2))
as_tibble(hill.indices.zones)

# Summary Hill numbers per site
hill.indices.sites<- diversity.metrics %>%
  group_by(Zone, Site) %>% na.omit() %>% 
  summarise(Taxq0 = round(mean(Richness), 2),
            Taxq0.SD = round(sd(Richness), 2),
            Taxq1 = round(mean(Tax.Entro), 2),
            Taxq1.SD = round(sd(Tax.Entro), 2),
            Functq0 = round(mean(FD.Rich), 2),
            Functq0.SD = round(sd(FD.Rich), 2),
            Functq1 = round(mean(FD.Entro), 2),
            Functq1.SD = round(sd(FD.Entro), 2))
as_tibble(hill.indices.sites)

# Remove objects
rm(origin, abundance, data.between, presence.conditions, qual,
   sub.abundance, sub.traits, traits, types, funct.diversity, dist.taxo,
   sub.gower, hill.tax.q1, hill.fd.q0, hill.fd.q1, distance.div,
   func.indices.zones, func.indices.site, hill.indices.zones,
   hill.indices.sites)


# Generalized Linear Mixed Models ####
# Richness/Taxonomic richness (q= 0)
glmm.rich.mod <- diversity.metrics %>% 
  glmer(Richness ~ Zone + (1|Year) + (1|Season) + (1|Site),
        family = 'poisson', data = .)

summary(glmm.rich.mod) # Significant effect
resid_panel(glmm.rich.mod)

# Taxonomic Distinctness
glmm.distinc.mod <- diversity.metrics %>% 
  glmer(Distinctness ~ Zone + (1|Year) + (1|Season) + (1|Site),
        family = 'Gamma', data = .)

summary(glmm.distinc.mod) # Non statistical effect
resid_panel(glmm.distinc.mod)

# Functional richness (FRic)
glmm.fric.mod <- diversity.metrics %>% na.omit() %>% 
  glmer(FRic ~ Zone + (1|Year) + (1|Season) + (1|Site),
        family = 'Gamma', data = .)

summary(glmm.fric.mod) # Significant effect
resid_panel(glmm.fric.mod)

# Functional originality (FOri)
glmm.fori.mod <- diversity.metrics %>% na.omit() %>% 
  glmer(FOri ~ Zone + (1|Year) + (1|Season) + (1|Site),
        family = 'Gamma', data = .)

summary(glmm.fori.mod) # Non statistical effect
resid_panel(glmm.fori.mod)

# Functional divergence (FDiv)
glmm.fdiv.mod <- diversity.metrics %>% na.omit() %>% 
  glmer(FDiv ~ Zone + (1|Year) + (1|Season) + (1|Site),
        family = 'Gamma', data = .)

summary(glmm.fdiv.mod) # Non statistical effect
resid_panel(glmm.fdiv.mod)

# Taxonomic entropy (q= 1)
glmm.taxq1.mod <- diversity.metrics %>% na.omit() %>% 
  glmer(Tax.Entro ~ Zone + (1|Year) + (1|Season) + (1|Site),
        family = 'Gamma', data = .)

summary(glmm.taxq1.mod) # Significant effect
resid_panel(glmm.taxq1.mod)

# Functional richness (q= 0)
glmm.funcq0.mod <- diversity.metrics %>% na.omit() %>% 
  glmer(FD.Rich ~ Zone + (1|Year) + (1|Season) + (1|Site),
        family = 'Gamma', data = .)

summary(glmm.funcq0.mod) # Non statistical effect
resid_panel(glmm.funcq0.mod)

# Functional entropy (q= 1)
glmm.funcq1.mod <- diversity.metrics %>% na.omit() %>% 
  glmer(FD.Entro ~ Zone + (1|Year) + (1|Season) + (1|Site),
        family = 'Gamma', data = .,
        control = glmerControl(optimizer = 'bobyqa'))

summary(glmm.funcq1.mod) # Non statistical effect
resid_panel(glmm.funcq1.mod)
