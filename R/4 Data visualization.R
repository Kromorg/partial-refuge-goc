# Open libraries y clean console ####
pacman::p_load(tidyverse,  # Data wrangling
               labdsv, # Convert wide to long format
               ggVennDiagram, # Venn diagram
               gridExtra) # Multiple plots

rm(list = ls())
shell('cls') # For Windows users
system2('clear') # For Mac users

# Venn diagram ####
origin <- read.csv('Data/Abundance_data.csv',
                   header = T, stringsAsFactors = T)
#origin <- read.csv(here:: here('Data/Abundance_data.csv'),
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
# Subset abundance data
abundance <- data.between %>%
  select(Abudefduf_troschelii:Zapteryx_exasperata)
rownames(abundance) <- data.between[, 6]

# Long format
long.abund <- dematrify(abundance) %>% 
  rename(Video.transect = sample, Species = species,
         Abundance = abundance)

# Add zone data
long.abund$Zone <- for (n in 1:length(long.abund$Video.transect)) {
  actual <- long.abund$Video.transect[n]
  extant <- which(data.between$Video.transect == actual)
  if(length(extant) > 0){
    value <- data.between$Zone[extant]
    long.abund[n, 4] <- value
  } else {
    long.abund[n, 4] <-  NA
  }
}

# Add site data
long.abund$Site <- for (n in 1:length(long.abund$Video.transect)) {
  actual <- long.abund$Video.transect[n]
  extant <- which(data.between$Video.transect == actual)
  if(length(extant) > 0){
    value <- data.between$Site[extant]
    long.abund[n, 5] <- value
  } else {
    long.abund[n, 5] <-  NA
  }
}

long.abund <- long.abund %>% rename(Zone = 'V4', Site = 'V5') %>% 
  select(Zone, Site, Species)

venn.zone <- list(Shallow = unique(long.abund[long.abund$Zone == 'Shallow', 3]),
                   Mesophotic = unique(long.abund[long.abund$Zone == 'Mesophotic', 3]))

venn.data <- venn.zone %>% Venn() %>% process_data()
venn.plot.zone <- ggplot()+
  geom_polygon(aes(X, Y, fill = id, group = id), 
               data = venn_regionedge(venn.data), 
               show.legend = F)+
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(venn.data), 
            show.legend = F)+
  geom_text(aes(X, Y, label = name), 
            data = venn_setlabel(venn.data), size = 15)+
  geom_label(aes(X, Y, label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), 
             data = venn_regionlabel(venn.data), size = 15) +
  scale_x_continuous(limits = c(-10, 5))+
  scale_fill_manual(values = c('darkred', 'darkslateblue', 'royalblue4'))+
  coord_equal()+
  theme_void()

venn.site <- list(Bajo_Shallow = unique(long.abund[long.abund$Zone == 'Shallow' &
                                                     long.abund$Site == 'El Bajo', 3]),
                  Lobos_Shallow = unique(long.abund[long.abund$Zone == 'Shallow'&
                                                      long.abund$Site == 'Punta Lobos', 3]),
                  Islotes_Shallow = unique(long.abund[long.abund$Zone == 'Shallow'&
                                                      long.abund$Site == 'Los Islotes', 3]),
                  Bajo_Meso = unique(long.abund[long.abund$Zone == 'Mesophotic'&
                                                   long.abund$Site == 'El Bajo', 3]),
                  Lobos_Meso = unique(long.abund[long.abund$Zone == 'Mesophotic'&
                                                  long.abund$Site == 'Punta Lobos', 3]),
                  Islotes_Meso = unique(long.abund[long.abund$Zone == 'Mesophotic'&
                                                  long.abund$Site == 'Los Islotes', 3]))

venn.data <- venn.site %>% Venn() %>% process_data()
venn.plot.site <- ggplot()+
  geom_polygon(aes(X, Y, fill = id, group = id), 
               data = venn_regionedge(venn.data), 
               show.legend = F)+
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(venn.data), 
            show.legend = F)+
  geom_text(aes(X, Y, label = name), 
            data = venn_setlabel(venn.data), size = 5)+
  geom_label(aes(X, Y, label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), 
             data = venn_regionlabel(venn.data), size = 5) +
  scale_x_continuous(limits = c(-125, 1250))+
  coord_equal()+
  theme_void()

multi <- grid.arrange(venn.plot.zone, venn.plot.site,
                      ncol = 1, nrow = 2)

ggsave('Figs/Appendix S5.tiff', plot = multi, width = 5000,
       height = 7000, units = 'px', dpi = 320, compression = "lzw")
#ggsave('Figs/Appendix S5.tiff', plot = multi, width = 1500,
#       height = 3000, units = 'px', dpi = 320, compression = "lzw")

# Open database ####
diversity.metrics <- read.csv('Data/Diversity_indices.csv',
                   header = T, stringsAsFactors = T)
#origin <- read.csv(here:: here('Data/Diversity_indices.csv'),
#                   header = T, stringsAsFactors = T)

# Change zone data into factor type
diversity.metrics$Zone<- factor(diversity.metrics$Zone,
                           levels = c('Shallow', 'Mesophotic'),
                           ordered = T)

# Functional richness graph
fric<- ggplot(diversity.metrics,
              aes(x = Site, y = FRic, colour = Zone))+
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
fori<- ggplot(diversity.metrics,
              aes(x = Site, y = FOri, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  ylim(0, 0.40)+ labs(y = 'FOri')+
  annotate('text', x = 1, y = 0.4, label = 'B)',
           size = 15)+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30))

# Functional divergence graph
fdiv<- ggplot(diversity.metrics,
              aes(x = Site, y = FDiv, colour = Zone))+
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

multi <- grid.arrange(fric, fori, fdiv, ncol = 3, nrow = 1)

ggsave('Figs/Figure 4.tiff', plot = multi, width = 6560, height = 3440,
       units = 'px', dpi = 320, compression = "lzw")
#ggsave(here:: here('Figs/Figure 3.tiff'), plot = multi, width = 6560,
#       height = 3440, units = 'px', dpi = 320, compression = "lzw")


# Taxonomic diversity graph
rich <- ggplot(diversity.metrics,
              aes(x = Site, y = Richness, colour = Zone))+
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
entro <- ggplot(diversity.metrics,
               aes(x = Site, y = Tax.Entro, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  scale_y_continuous(limits = c(1, 13), breaks = seq(1, 13, 3))+
  labs(y = 'Entropy')+
  annotate('text', x = 1, y = 13, label = 'B)*',
           size = 15)+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30))

# Taxonomic distinctness graph
delta <- ggplot(diversity.metrics,
               aes(x = Site, y = Distinctness, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
  labs(title = expression(Phylogenetic~diversity))+
  labs(y = 'Taxonomic distinctness')+
  annotate('text', x = 1, y = 100, label = 'E)',
           size = 15)+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30))

# Taxonomic distinctness funnel plot
funnel <- ggplot(diversity.metrics,
               aes(x = Richness, y = Distinctness, colour = Zone))+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
  scale_x_continuous(limits = c(0, 27), breaks = seq(0, 27, 5))+
  geom_hline(yintercept = diversity.metrics$EDstar,
            color = 'black',
            linetype = 'dashed')+
  geom_ribbon(aes(x = Richness, colour = NULL,
                  ymin = EDstar - 2*sd.Dplus,
                  ymax = EDstar + 2*sd.Dplus),
              alpha = 0.2, show.legend = F)+
  labs(x = 'Species richness', y = 'Taxonomic distinctness')+
  annotate('text', x = 5, y = 100, label = 'F)',
           size = 15)+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30))

# Functional diversity
fd.rich<- ggplot(diversity.metrics,
                 aes(x = Site, y = FD.Rich, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  scale_y_continuous(limits = c(1, 6), breaks = seq(1, 6, 1))+
  labs(title = expression(Functional~diversity))+
  annotate('text', x = 1, y = 6, label = 'C)',
           size = 15)+
  labs(y = 'Richness')+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

fd.entro<- ggplot(diversity.metrics,
                  aes(x = Site, y = FD.Entro, colour = Zone))+
  geom_boxplot(show.legend = F, size = 1.5)+
  geom_point(size = 5.5, alpha = 0.4, position = 'jitter',
             show.legend = F)+
  scale_colour_manual(values = c('darkred', 'royalblue4'))+
  ylim(1, 5)+
  annotate('text', x = 1, y = 5, label = 'D)',
           size = 15)+
  labs(y = 'Entropy')+
  theme_classic(base_size = 25)+
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 30))

hill <- grid.arrange(rich, fd.rich, delta, entro, fd.entro, funnel,
                     ncol = 3, nrow = 2)

ggsave('Figs/Figure 5.tiff', plot = hill, width = 7000, height = 4700,
       compression = "lzw", units = 'px', dpi = 320)
#ggsave(here:: here('Figs/Figure 4.tiff'), plot = hill, width = 5500,
#       height = 4700, compression = "lzw", units = 'px', dpi = 320)
