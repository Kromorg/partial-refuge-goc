# Open libraries y clean console ####
pacman::p_load(tidyverse,  # Data wrangling
               gridExtra) # Multiple plots

rm(list = ls())
shell('cls') # For Windows users
system2('clear') # For Mac users

# Open database ####
diversity.metrics <- read.csv('Data/Diversity_indices.csv',
                   header = T, stringsAsFactors = T)
#origin <- read.csv(here:: here('Data/Diversity_indices.csv'),
#                   header = T, stringsAsFactors = T)

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
  stat_summary(fun = mean, geom = "point", col = "black",
               size = 9, fill = 'black', shape = 23,
               aes(group = Zone))+
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

ggsave('Figs/Figure 3.tiff', plot = multi, width = 6560, height = 3440,
       units = 'px', dpi = 320, compression = "lzw")
#ggsave(here:: here('Figs/Figure 3.tiff'), plot = multi, width = 6560,
#       height = 3440, units = 'px', dpi = 320, compression = "lzw")


# Taxonomic diversity graph
rich<- ggplot(diversity.metrics,
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
entro<- ggplot(diversity.metrics,
               aes(x = Site, y = Tax.Entro, colour = Zone))+
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
fd.rich<- ggplot(diversity.metrics,
                 aes(x = Site, y = FD.Rich, colour = Zone))+
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

fd.entro<- ggplot(diversity.metrics,
                  aes(x = Site, y = FD.Entro, colour = Zone))+
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

ggsave('Figs/Figure 4.tiff', plot = hill, width = 5500, height = 4700,
       compression = "lzw", units = 'px', dpi = 320)
#ggsave(here:: here('Figs/Figure 4.tiff'), plot = hill, width = 5500,
#       height = 4700, compression = "lzw", units = 'px', dpi = 320)
