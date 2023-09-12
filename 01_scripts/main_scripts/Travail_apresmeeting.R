#### Importation des packages ----
library(vegan)
library(dplyr)
library(tidyverse)
library(ade4)
library(adespatial)
library(adegraphics)
library(sf)
library(tmap)
library(spdep)
library(ggplot2)

#### Importation des données ----
# Arbres en fonction des placettes
spe <- read.csv("data_output/spe.csv", row.names=1) # rn pour rownames (garde les grid_id)
head(spe)

# Placettes forestières et données environnementales
grid_topo <- read_sf("data_output/grid_topo.geojson") # on lit le fichier

grid.env2 <- read.csv("data_output/grid.env2.csv", row.names=1) # rn pour rownames (garde les grid_id)
head(grid.env2)

#### Transformation des données environnementales ----
# Scatter plots for all pairs of environmental variables (NEwR, chap.2) 
source('functions/panelutils.R')
# Bivariate plots with histograms on the diagonal and smooth fitted curves
dev.new(
  title = "Bivariate descriptor plots - env",
  width = 10,
  height = 10,
  noRStudioGD = TRUE)
pairs(grid.env2, 
      panel = panel.smooth, 
      diag.panel = panel.hist,
      main = "Bivariate Plots with Histograms and Smooth Curves")

# Standardization of all environmental variables (NEwR, chap.2)
# Center and scale = standardize the variables (z-scores)
grid.env2 <- as.data.frame(scale(grid.env2)) # means = 0 // standard deviations = 1

#### Exploration de la table de données d'arbres - après correction ----
# From NEwR: # Exploration of a data frame using basic R functions (Chap.2)
spe[1:5, 1:10]            # Display only 5 lines and 10 columns
head(spe)                 # Display only the first 6 lines
tail(spe)                 # Display only the last 6 rows
nrow(spe)                 # Number of rows (sites)
ncol(spe)                 # Number of columns (species)
dim(spe)                  # Dimensions of the data frame (rows, columns)
colnames(spe)             # Column labels (descriptors = species)
summary(spe)              # Descriptive statistics for columns

# Création de carte choroplètes montrant des statistiques d'abondance / diversité ----
# Ajout de statistiques dans un tableau
spe_stat <- spe %>% 
  mutate(total_n = rowSums(.)) %>% # nombre d'individus
  mutate(n_sp = apply(spe > 0, 1, sum)) %>% # nombre d'espèces
  mutate(grid_id = row.names(spe)) # grid_id
head(spe_stat)

# On joint les statistiques avec les placettes pour pouvoir les visualiser
grid_arbres <- left_join(grid_topo, spe_stat)

# Création des cartes choroplètes
carte_total_n <- tm_shape(grid_arbres) +
  tm_fill(col = "total_n",
          title = "Nombre d'individus") +
  tm_borders() 

carte_n_sp <- tm_shape(grid_arbres) +
  tm_fill(col = "n_sp",
          title = "Nombre d'espèces") +
  tm_borders() 

nb_sp_nb_ind <- tmap_arrange(carte_total_n, carte_n_sp, ncol = 2)
tmap_save(tm = nb_sp_nb_ind, # objet tmap qu'on veut enregistrer
          filename = "figures/nb_sp_nb_ind.png", # nom et chemin du fichier qu'on enregistre, incluant l'extension
          height = 4, # hauteur de la carte
          width = 8, # largeur de la carte
          units = "in", # unités utilisées pour la hauteur et la largeur
          dpi = 300) # par défault dpi = 300
# Distribution des abondances ----
# From NEwR: Overall distribution of abundances (dominance codes) (Chap.2)

# Minimum and maximum of abundance values in the whole data set
range(spe)
# Minimum and maximum value for each species
apply(spe, 2, range)
# Count the cases for each abundance class
(ab <- table(unlist(spe)))
ab_nozero <- ab[2:30]

# Barplot of the distribution, all species confounded
png("figures/classes_abondance_avec0.png", height = 4, width = 5, 
    units = "in", # unités
    res = 300) # résolution
barplot(ab,
        main = "Distribution des classes d'abondance",
        las = 1,
        xlab = "Classe d'abondance",
        ylab = "Fréquence",
        col = gray(5 : 0 / 5))
dev.off()

# La classe d'abondance comptant les zéros a été enlevée afin de permettre d'avoir une meilleure idée de la répartition des classes d'abondance.
# Il y avait 6075 zéros dans la table de données.
png("figures/classes_abondance.png", height = 4, width = 5, 
    units = "in", # unités
    res = 300) # résolution
barplot(ab_nozero,
        main = "Distribution des classes d'abondance",
        las = 1,
        xlab = "Classe d'abondance",
        ylab = "Fréquence",
        ylim = c(0,600),
        col = gray(5 : 0 / 5))
dev.off()

# Number of absences
sum(spe == 0)
# Proportion of zeros in the community data set
sum(spe == 0) / (nrow(spe) * ncol(spe))

# From NEwR: Compare species: number of occurrences (Chap.2) ----
# Compute the number of sites where each species is present
spe.pres <- apply(spe > 0, 2, sum) # To sum by columns, the second argument of apply(), MARGIN, is set to 2
sort(spe.pres) # Sort the results in increasing order

# Compute percentage frequencies
spe.relf <- 100 * spe.pres/nrow(spe)
round(sort(spe.relf), 1) # Round the sorted output to 1 digit

# Plot the histograms
# Enregistrement de la figure
png("figures/hist_freq.png", height = 4, width = 8, 
    units = "in", # unités
    res = 300) # résolution

# Divide the window horizontally
par(mfrow = c(1,2))
hist(spe.pres, 
     main = "Occurrences des espèces", 
     right = FALSE, 
     las = 1, 
     xlab = "Nombre d'occurences", 
     ylab = "Nombre d'espèces", 
     breaks = seq(0, 500, by = 25),
     col = "bisque")

hist(spe.relf, 
     main = "Fréquences relatives des espèces", 
     right = FALSE, 
     las = 1,
     xlab = "Fréquence des occurrences (%)", 
     ylab = "Nombre d'espèces",
     breaks = seq(0, 100, by = 10),
     col = "bisque")

dev.off()

#### Transformation des données d'arbres ----
# Q-Mode: Transformations pour des données quantitatives d'espèces ----
# Nous voulons conserver les distances entre les éléments pour les analyses statistiques qui suivront.
# Les distances de Hellinger et de chorde sont appropriées pour l’ordination des données de communautés
# parce qu’elles possèdent 9 propriétés importantes pour l’étude de la diversité bêta.

source('functions/box.cox.chord.R') # provient de Appendix 4 (Legendre & Borcard, 2018)

# Sans transformation
spe.d <- dist(spe)
is.euclid(spe.d)

# Transformation en profils d'abondances relatives
spe.total = decostand(spe,"total") ## package vegan
spe.dt <- dist(spe.total)
is.euclid(spe.dt) # True ## package ade4

# Transformation de khi-carré
spe.chisq = decostand(spe,"chi.sq")
spe.dchisq <- dist(spe.chisq)
is.euclid(spe.dchisq) # True

# Chord distance matrix
spe.norm <- decostand(spe, "nor")
spe.dc <- dist(spe.norm)
is.euclid(spe.dc) # True

# Hellinger distance matrix (Box-cox chord : bc.exp = 0.5)
spe.hel <- decostand(spe, "hel") # transfo Hellinger
spe.dh <- dist(spe.hel) # Matrice de distance Hellinger
is.euclid(spe.dh) # True

# Log-chord distance matrix (Box-cox chord : bc.exp = 0)
spe.ln <- log1p(spe)
spe.ln.norm <- decostand(spe.ln, "nor")
spe.logchord <- dist(spe.ln.norm)
is.euclid(spe.logchord) # True

# Double square-root (fourth-root) (Box-cox chord : bc.exp = 0.25)
spe.bc <- box.cox.chord(spe, bc.exp = 0.25)
spe.dbc <- dist(spe.bc)
is.euclid(spe.dbc) # True

# Boxplots of transformed abundances of a common species ----
# (BEPA, species #3)
# Enregistrement de la figure
png("figures/transfo.png", height = 4, width = 8, 
    units = "in", # unités
    res = 300) # résolution
par(mfrow = c(1,2))
lablist1<-as.vector(c("Données brutes", "sqrt", "log"))
boxplot(spe$BEPA,
        sqrt(spe$BEPA), 
        log1p(spe$BEPA),
        las = 1, 
        main = "Transformations simples",
        col = "bisque",
        pars  =  list(xaxt = "n"))
axis(1, at=seq(1,3, by=1), labels = FALSE)
text(seq(1,3, by=1), par("usr")[1]-3.25, labels = lablist1, srt = 35, pos = 1, xpd = TRUE)

lablist2<-as.vector(c("Hellinger", "Chord", "log-Chord", "Khi-carré", "Total", "Double sqrt" ))
boxplot(spe.hel$BEPA, 
        spe.norm$BEPA, 
        spe.ln.norm$BEPA,
        spe.chisq$BEPA,
        spe.total$BEPA,
        spe.bc$BEPA,
        las = 1, 
        main = "Transformations
        (distances euclidiennes)", 
        col = "lightblue",
        pars  =  list(xaxt = "n"))
axis(1, at=seq(1,6, by=1), labels = FALSE)
text(seq(1,6, by=1), par("usr")[1]-0.5, labels = lablist2, srt = 35, pos = 1, xpd = TRUE)
dev.off()

#### Ordination canonique - Choix de la transformation des abondances ----
# Sans transformation
spe.brute.rda <- rda(spe, grid.env2)
brute.R2aj <- RsquareAdj(spe.brute.rda)$adj.r.squared

# Chorde
spe.norm.rda <- rda(spe.norm, grid.env2)
chorde.R2aj <- RsquareAdj(spe.norm.rda)$adj.r.squared 

# Hellinger
spe.hel.rda <- rda(spe.hel, grid.env2)
hel.R2aj <- RsquareAdj(spe.hel.rda)$adj.r.squared 

# Khi-carré
spe.chisq.rda <- rda(spe.chisq, grid.env2)
chisq.R2aj <- RsquareAdj(spe.chisq.rda)$adj.r.squared 

# Total
spe.total.rda <- rda(spe.total, grid.env2)
tot.R2aj <- RsquareAdj(spe.total.rda)$adj.r.squared 

# log-chorde
spe.ln.norm.rda <- rda(spe.ln.norm, grid.env2)
ln.norm.R2aj <- RsquareAdj(spe.ln.norm.rda)$adj.r.squared

# double square root (Box-cox chord, bc.exp = 0.25)
spe.bc.rda <- rda(spe.bc, grid.env2)
boxcox.R2aj <- RsquareAdj(spe.bc.rda)$adj.r.squared 

transfo.R2aj <- data.frame("Transformation" = c("Aucune", "Chorde", "Hellinger",
                                        "Khi-carré", "Profils d'abondances relatives",
                                        "log-Chorde", "Double racine carrée"),
                   "R2ajusté" = c(brute.R2aj, chorde.R2aj, hel.R2aj, chisq.R2aj,
                                   tot.R2aj, ln.norm.R2aj, boxcox.R2aj))


write.csv(transfo.R2aj,"data_output/transfo.R2aj.csv", row.names = FALSE)

# Les transformations de Hellinger et log-Chorde donnent toutes les deux de bons R2.
# Nous allons donc tenter une analyse de groupement avec contrainte de contiguité spatiale
# pour chacune des transformations et voir laquelle produit le résultat le plus logique d'un point de vue écologique.

## Choix de la transformation:
# les noms sont changés par des noms plus généralistes pour limiter
# les modifications dans le cas où la transformation serait à changer
spe.transfo <- spe.ln.norm
spe.dist <- spe.logchord

#### Groupement avec contrainte de contiguité spatiale ----
# Liste de liens ----
# Source du script (liste de liens et retrait des liens trop longs) : https://www.jstatsoft.org/article/view/v103i07

# Matrice Y = spe
# Matrice xy = site coordinates = placettes
# Matrice E = matrix of connecting edges
placettes <- grid_coord[1:2]

# Obtenir une liste de liens
placettes.edge <- placettes %>% # package spdep
  tri2nb %>%
  nb2listw(style = "B") %>%
  listw2sn

# Calcul de la longueur des liens pour enlever ceux trop longs (longueur normale: 20 m)
names(placettes.edge)[3L] <- "distance"
placettes.edge$distance <- placettes %>%
  dist %>%
  as.matrix %>%
  .[placettes.edge[,1L:2L] %>% as.matrix]

# Les liens en diagonale semblent être 28.28 m au lieu de 20 m.
# C'est probablement une erreur de la projection.
# On accepte donc tous les liens sous 30 m.
library("magrittr")
placettes.edge %<>% .[.$distance <= 30,] # package magrittr pour le pipe %<>%

# Groupement (log-Chord) ----
# Source du script : https://search.r-project.org/CRAN/refmans/adespatial/html/constr.hclust.html
# Clustering with a contiguity constraint described by a list of links:
grpWD2cst_constr_hclust <-
  constr.hclust(
    d = spe.dist, # Response dissimilarity matrix (log-Chord distance)
    method = "ward.D2", # Clustering method
    links = placettes.edge, # File of link edges (constraint)
    coords = placettes) # File of geographic coordinates

# Warning message:
# In matrix(links, nl, 2L) :
#   data length differs from size of matrix: [6714 != 2238 x 2]

## To visualize using hclust's plotting method:
#stats:::plot.hclust(grpWD2cst_constr_hclust, hang=-1)

# Generic functions from hclust can be used, for instance to obtain a list of members of each cluster:
spe.3constr.g <- cutree(grpWD2cst_constr_hclust, k=3) # pour 3 communautés
spe.4constr.g <- cutree(grpWD2cst_constr_hclust, k=4) # pour 4 communautés
spe.5constr.g <- cutree(grpWD2cst_constr_hclust, k=5) # pour 5 communautés

# Visualisation des groupements sur une carte de liens avec contrainte de contiguité spatiale ----
# Plot the results on a map with k=3 clusters:
png("figures/carte_liens3_logchord.png",
    height = 4, width = 3, # taille
    units = "in", # unités
    res = 300) # résolution
plot(grpWD2cst_constr_hclust, k=3, links=TRUE, main = "Trois groupements", las=1, xlab="x (m)",
     ylab="y (m)", cex=0.6, lwd=1, cex.axis = 0.5, cex.main = 0.75, cex.lab = 0.75)
dev.off()

# Now with k=4 clusters:
png("figures/carte_liens4_logchord.png",
    height = 4, width = 3, # taille
    units = "in", # unités
    res = 300) # résolution
plot(grpWD2cst_constr_hclust, k=4, links=TRUE, main = "Quatre groupements", las=1, xlab="x (m)",
     ylab="", yaxt="n", cex=0.6, lwd=1, cex.axis = 0.5, cex.main = 0.75, cex.lab = 0.75)
dev.off()

# Now with k=5 clusters:
png("figures/carte_liens5_logchord.png",
    height = 4, width = 3, # taille
    units = "in", # unités
    res = 300) # résolution
plot(grpWD2cst_constr_hclust, k=5, links=TRUE, main = "Cinq groupements", las=1, xlab="x (m)",
     ylab="", yaxt="n",
     #ylab="y (m)",
     cex=0.6, lwd=1, cex.axis = 0.5, cex.main = 0.75, cex.lab = 0.75)
dev.off()

#### Définition des communautés écoforestières ----
# Fréquences relatives des espèces pour nommer les communautés ----
# Tableau des superficies de canopées par espèce par placette
surfs <- read_csv("data_output/surfs.csv")
head(surfs)

surfs$grid_id <- as.character(surfs$grid_id) # On transforme que les grid_id en character

# Création d'une table de type "tibble" - 3 communautés
comm3_tbl <- tibble(grid_id = names(spe.3constr.g), # on crée une variable grid_id avec les noms d'étiquettes
                    comm3 = spe.3constr.g) # on nomme la variable contenant les groupes "comm" pour communauté
comm3_tbl # on examine le résultat

# Création d'une table de type "tibble" - 4 communautés
comm4_tbl <- tibble(grid_id = names(spe.4constr.g), # on crée une variable grid_id avec les noms d'étiquettes
                    comm4 = spe.4constr.g) # on nomme la variable contenant les groupes "comm" pour communauté
comm4_tbl # on examine le résultat

# Création d'une table de type "tibble" - 5 communautés
comm5_tbl <- tibble(grid_id = names(spe.5constr.g), # on crée une variable grid_id avec les noms d'étiquettes
                    comm5 = spe.5constr.g) # on nomme la variable contenant les groupes "comm" pour communauté
comm5_tbl # on examine le résultat

# Jointure des superficies de canopées par espèce par placette avec leurs communautés respectives
surfs_comm <- surfs %>% # on part du tableau de superficies pour espèces sélectionnées
  left_join(comm3_tbl) %>% # on joint les deux tableaux par la variable commune
  left_join(comm4_tbl) %>%
  left_join(comm5_tbl)
surfs_comm

# Calcul des superficies relatives moyennes par espèce par communauté - 3 communautés
surf_moy_comm3 <- surfs_comm %>% 
  #select(grid_id, comm3, Label, sum_surf_m2_ha) %>%  # sélectionner seulement variables pertinentes
  ungroup() %>% # pour enlever les groupes précédents au cas où
  group_by(comm3, Label) %>% 
  summarise(surf_moy = mean(sum_surf_m2_ha)) %>% 
  ungroup() %>% 
  group_by(comm3) %>% 
  mutate(surf_rel_moy = surf_moy / sum(surf_moy))

# Calcul des superficies relatives moyennes par espèce par communauté - 4 communautés
surf_moy_comm4 <- surfs_comm %>% 
  #select(grid_id, comm3, Label, sum_surf_m2_ha) %>%  # sélectionner seulement variables pertinentes
  ungroup() %>% # pour enlever les groupes précédents au cas où
  group_by(comm4, Label) %>% 
  summarise(surf_moy = mean(sum_surf_m2_ha)) %>% 
  ungroup() %>% 
  group_by(comm4) %>% 
  mutate(surf_rel_moy = surf_moy / sum(surf_moy))

# Calcul des superficies relatives moyennes par espèce par communauté - 5 communautés
surf_moy_comm5 <- surfs_comm %>% 
  #select(grid_id, comm5, Label, sum_surf_m2_ha) %>%  # sélectionner seulement variables pertinentes
  ungroup() %>% # pour enlever les groupes précédents au cas où
  group_by(comm5, Label) %>% 
  summarise(surf_moy = mean(sum_surf_m2_ha)) %>% 
  ungroup() %>% 
  group_by(comm5) %>% 
  mutate(surf_rel_moy = surf_moy / sum(surf_moy))

# Définition des noms des communautés ----
# Pour 3 communautés (avec log-Chord) :
comm3_tbl <- comm3_tbl %>% # on peut stocker les noms dans notre objet contenant l'appartenance de chaque station
  mutate(comm3_nom = recode_factor(comm3, # on copie le facteur comm, en renommant les niveaux
                                   `1` = 'Bétulaie à peuplier et érable rouge',
                                   `2` = 'Lariçaie à bouleau et pin blanc',
                                   `3` = 'Bétulaie à peuplier et pin blanc'))

comm3_tbl # on examine le résultat

# Pour 4 communautés (avec log-Chord) :
comm4_tbl <- comm4_tbl %>% # on peut stocker les noms dans notre objet contenant l'appartenance de chaque station
  mutate(comm4_nom = recode_factor(comm4, # on copie le facteur comm, en renommant les niveaux
                                   `1` = 'Bétulaie à peuplier et érable rouge',
                                   `2` = 'Lariçaie à bouleau et pin blanc',
                                   `3` = 'Bétulaie à peuplier et pin blanc',
                                   `4` = 'Bétulaie à peuplier et érable rouge (2)'))

comm4_tbl # on examine le résultat

# Pour 5 communautés (avec log-Chord) :
comm5_tbl <- comm5_tbl %>% # on peut stocker les noms dans notre objet contenant l'appartenance de chaque station
  mutate(comm5_nom = recode_factor(comm5, # on copie le facteur comm, en renommant les niveaux
                                   `1` = 'Bétulaie à peuplier et érable rouge',
                                   `2` = 'Lariçaie à bouleau blanc et pin blanc',
                                   `3` = 'Pinède à peuplier',
                                   `4` = 'Bétulaie à peuplier et érable rouge (2)',
                                   `5` = 'Peupleraie à érable rouge'))

comm5_tbl # on examine le résultat

# On enregistre le fichier avec les noms des communautés
write_csv(comm3_tbl, "data_output/comm3.csv")
write_csv(comm4_tbl, "data_output/comm4.csv")
write_csv(comm5_tbl, "data_output/comm5.csv")

# Figures des fréquences relatives (après la détermination des noms) ----
# Pour 3 communautés (avec log-Chord) :
surf_moy_comm3 <- surf_moy_comm3 %>% # on peut stocker les noms dans notre objet contenant l'appartenance de chaque station
  mutate(comm3_nom = recode_factor(comm3, # on copie le facteur comm, en renommant les niveaux
                                   `1` = 'Bétulaie à peuplier et érable rouge',
                                   `2` = 'Lariçaie à bouleau blanc et pin blanc',
                                   `3` = 'Bétulaie à peuplier et pin blanc'))
surf_moy_comm3 # on examine le résultat

# Pour 4 communautés (avec log-Chord) :
surf_moy_comm4 <- surf_moy_comm4 %>% # on peut stocker les noms dans notre objet contenant l'appartenance de chaque station
  mutate(comm4_nom = recode_factor(comm4, # on copie le facteur comm, en renommant les niveaux
                                   `1` = 'Bétulaie à peuplier et érable rouge',
                                   `2` = 'Lariçaie à bouleau et pin blanc',
                                   `3` = 'Bétulaie à peuplier et pin blanc',
                                   `4` = 'Bétulaie à peuplier et érable rouge (2)'))
surf_moy_comm4 # on examine le résultat

# Pour 5 communautés (avec log-Chord) :
surf_moy_comm5 <- surf_moy_comm5 %>% # on peut stocker les noms dans notre objet contenant l'appartenance de chaque station
  mutate(comm5_nom = recode_factor(comm5, # on copie le facteur comm, en renommant les niveaux
                                   `1` = 'Bétulaie à peuplier et érable rouge',
                                   `2` = 'Lariçaie à bouleau et pin blanc',
                                   `3` = 'Pinède à peuplier',
                                   `4` = 'Bétulaie à peuplier et érable rouge (2)',
                                   `5` = 'Peupleraie à érable rouge'))
surf_moy_comm5 # on examine le résultat

# Figure des fréquences relatives - 3 communautés
graph_comm3 <- ggplot(surf_moy_comm3,
                      aes(x = Label,
                      y = surf_rel_moy)) +
  geom_bar(stat = 'identity') +
  xlab("Code d'espèce") +
  ylab("Superficie relative") +
  coord_flip() +
  facet_wrap(~comm3_nom, ncol = 3)
graph_comm3

# Enregistrement de la figure - 3 communautés
ggsave(graph_comm3, file = 'figures/graph_comm3.png', # nom du fichier avec extension
       width = 8, # largeur en pouces
       height = 3, # hauteur en pouces
       dpi = 300) # résolution en pixels par pouce

# Figure des fréquences relatives - 4 communautés
graph_comm4 <- ggplot(surf_moy_comm4,
                      aes(x = Label,
                      y = surf_rel_moy)) +
  geom_bar(stat = 'identity') +
  xlab("Code d'espèce") +
  ylab("Superficie relative") +
  coord_flip() +
  facet_wrap(~comm4_nom, ncol = 2)
graph_comm4

# Enregistrement de la figure - 4 communautés
ggsave(graph_comm4, file = 'figures/graph_comm4.png', # nom du fichier avec extension
       width = 6, # largeur en pouces
       height = 6, # hauteur en pouces
       dpi = 300) # résolution en pixels par pouce

# Figure des fréquences relatives - 5 communautés
graph_comm5 <- ggplot(surf_moy_comm5,
                           aes(x = Label,
                               y = surf_rel_moy)) +
  geom_bar(stat = 'identity') +
  xlab("Code d'espèce") +
  ylab("Superficie relative") +
  coord_flip() +
  facet_wrap(~comm5_nom, ncol = 3)
graph_comm5

# Enregistrement de la figure - 5 communautés
ggsave(graph_comm5, file = 'figures/graph_comm5.png', # nom du fichier avec extension
       width = 8, # largeur en pouces
       height = 6, # hauteur en pouces
       dpi = 300) # résolution en pixels par pouce

# Cartographie des communautés avec noms ----
# Jointure des communautés avec le fichier des placettes géoréférencées avec statistiques
grid_arbres <- grid_arbres %>% # on part de l'objet grid_arbres
  left_join(comm3_tbl, # on joint le tableau des 3 communautés
            by = "grid_id") %>% # on lie les tableaux par les noms de placettes
  left_join(comm4_tbl, # on joint le tableau des 4 communautés
            by = "grid_id") %>% # on lie les tableaux par les noms de placettes
  left_join(comm5_tbl, # on joint le tableau des 5 communautés
            by = "grid_id") # on lie les tableaux par les noms de placettes
grid_arbres # on vérifie que les nouvelles colonnes sont présentes

plot(grid_arbres["comm3_nom"]) # on visualise le résultat rapido
plot(grid_arbres["comm4_nom"]) # on visualise le résultat rapido
plot(grid_arbres["comm5_nom"]) # on visualise le résultat rapido

dev.off()

# Importation d'une orthomosaïque
require(stars)
ortho <- read_stars('data_input/2021-09-28-sbl-cloutier-z1-P4RTK-MS.tif', 
                    proxy = TRUE,
                    NA_value = 0) %>% # car les valeurs manquantes sont 0 dans ce raster
  st_transform(2950) # on transforme au code 2950
# plot(ortho, rgb = 1:3) # pour une image couleur

st_crs(ortho)$epsg 

# Carte statique des communautés (3)
tmap_mode("plot")

carte_3comm <- tm_shape(ortho, # on charge notre orthomosaïque
                        bbox = grid_arbres) + # on définit l'étendue à partir des communautés
  tm_rgba() + # pour image RGB avec transparence (a)
  tm_shape(grid_arbres) + # on charge les points
  tm_polygons(col = "comm3_nom",
              alpha = 0.5,
              title = "Communautés") +
  tm_layout(legend.position = c('right', 'bottom'),
            inner.margins = c(0.05, .05, .05, .05)) + # on augmente les marges
  tm_compass(type = "8star", position = c("left", "top")) + # rose des vents
  tm_scale_bar() + # l'échelle
  tm_graticules(lwd = 0.3,
                alpha = 0.6)
#carte_3comm # on visualise le résultat

tmap_save(tm = carte_3comm,
          filename = 'figures/carte_3comm.png', # on nomme le fichier # l'extension définit le format
          dpi = 600, # résolution en pixels par pouce
          width = 8) # largeur de l'image en pouces

# Carte statique des communautés (4)
carte_4comm <- tm_shape(ortho, # on charge notre orthomosaïque
                        bbox = grid_arbres) + # on définit l'étendue à partir des communautés
  tm_rgba() + # pour image RGB avec transparence (a)
  tm_shape(grid_arbres) + # on charge les points
  tm_polygons(col = "comm4_nom",
              alpha = 0.5,
              title = "Communautés") +
  tm_layout(legend.position = c('right', 'bottom'),
            inner.margins = c(0.05, .05, .05, .05)) + # on augmente les marges
  tm_compass(type = "8star", position = c("left", "top")) + # rose des vents
  tm_scale_bar() + # l'échelle
  tm_graticules(lwd = 0.3,
                alpha = 0.6)
#carte_4comm # on visualise le résultat

tmap_save(tm = carte_4comm,
          filename = 'figures/carte_4comm.png', # on nomme le fichier # l'extension définit le format
          dpi = 600, # résolution en pixels par pouce
          width = 8) # largeur de l'image en pouces

# Carte statique des communautés (5)
carte_5comm <- tm_shape(ortho, # on charge notre orthomosaïque
                        bbox = grid_arbres) + # on définit l'étendue à partir des communautés
  tm_rgba() + # pour image RGB avec transparence (a)
  tm_shape(grid_arbres) + # on charge les points
  tm_polygons(col = "comm5_nom",
              alpha = 0.5,
              title = "Communautés") +
  tm_layout(legend.position = c('right', 'bottom'),
            inner.margins = c(0.05, .05, .05, .05)) + # on augmente les marges
  tm_compass(type = "8star", position = c("left", "top")) + # rose des vents
  tm_scale_bar() + # l'échelle
  tm_graticules(lwd = 0.3,
                alpha = 0.6)
#carte_5comm # on visualise le résultat

tmap_save(tm = carte_5comm,
          filename = 'figures/carte_5comm.png', # on nomme le fichier # l'extension définit le format
          dpi = 600, # résolution en pixels par pouce
          width = 8) # largeur de l'image en pouces

# EXTRA: Carte interactive - Pour comparaison visuelle entre 3, 4 et 5 communautés ----

## !!!! Ajout MNT et variables env. dans la carte interactive!!!!

# tmap_mode("view")
# 
# carteinteractive <- tm_shape(ortho, # on charge notre orthomosaïque
#                              bbox = grid_arbres,
#                              name = "Orthomosaïque") + # on définit l'étendue à partir des communautés
#   tm_rgba() + # pour image RGB avec transparence (a)
#   tm_shape(grid_arbres,
#            name = "Trois communautés") + # on charge les communautés observées (3)
#   tm_polygons(col = 'comm3_nom',
#               alpha = 0.5,
#               title = "Trois communautés") +
#   tm_shape(grid_arbres,
#            name = "Quatre communautés") + # on charge les communautés observées (4)
#   tm_polygons(col = 'comm4_nom',
#               alpha = 0.5,
#               title = "Quatre communautés") +
#   tm_shape(grid_arbres,
#            name = "Cinq communautés") + # on charge les communautés attendues (5)
#   tm_polygons(col = 'comm5_nom',
#               alpha = 0.5,
#               title = 'Cinq communautés')
# 
# tmap_save(tm = carteinteractive,
#           filename = 'figures/carteinteractive.html', # on nomme le fichier # l'extension définit le format
#           dpi = 600, # résolution en pixels par pouce
#           width = 8) # largeur de l'image en pouces
# tmap_mode("plot")

#### Ordination simple - Analyse en composantes principales ----
# PCA avec rda() de VEGAN ----
spe.acp <- rda(spe.transfo, scale = FALSE) # scale = false parce que les données sont déjà transformées
summary(spe.acp)

(pourc_var1=(spe.acp$CA$eig[1])/(sum(spe.acp$CA$eig))) # 0.2698085  -> 26.98%
(pourc_var2=(spe.acp$CA$eig[2])/(sum(spe.acp$CA$eig))) # 0.1373318  -> 13.73%

pourc_var1+pourc_var2 # 0.4071402  

# Les deux premiers axes expliquent 40.71% de la variance.

spe.good.si <- goodness(spe.acp, display="sites", model="CA") # voir TP en R, p.9

# Graphiques avec ggplot tiré du TP4 - Stage en écologie végétale ----
acp1_summ <- summary(spe.acp, scaling = 1) # on crée un objet contenant le *résumé* de l'ACP, cadrage type 1
acp1_sites <- acp1_summ$sites[, 1:2] # on sauvegarde les coordonnées des stations le long des deux premiers axes
head(acp1_sites) # on examine le résultat

acp1_sp <- acp1_summ$species[, 1:2] # on sauvegarde les coordonnées des vecteurs-espèces le long des deux premiers axes
acp1_sp # on examine le résultat

sites <- as_tibble(acp1_sites) %>% # on convertit en tibble
  mutate(grid_id = rownames(acp1_sites)) %>% # on ajoute une variable station
  left_join(comm3_tbl) # on ajoute les communautés

sites # on examine le résultat

especes <- as_tibble(acp1_sp) %>% # on convertit en tibble
  mutate(codesp = rownames(acp1_sp)) # on ajoute les codes d'espèces

especes # on examine le résultat


library(svglite)
library(spocc)
svglite("figures/ACP_3comm.svg")
# Avec 3 communautés : 
acp_graph_sites <- ggplot(sites) + # on utilise le tableau sites
  geom_text(aes(x = PC1,
                y = PC2,
                label = grid_id, # on utilise les étiquettes de station
                colour = comm3_nom)) +
  xlab('ACP 1 (26.98%)') + # on indique le % de variance pour l'axe 1 et 2
  ylab('ACP 2 (13.73%)') +
  geom_hline(yintercept = 0, linetype = 'dotted') + # on ajoute des lignes à 0,0
  geom_vline(xintercept = 0, linetype = 'dotted') +
  theme_bw() + # thème plus minimaliste
  theme(panel.grid.major = element_blank(), # on enlève les lignes majeures et mineures
        panel.grid.minor = element_blank(),
        legend.position = 'top') + # on bouge la légende au dessus
  scale_colour_discrete(name = NULL) # on enlève le nom de la légende

acp_graph_sites # voyons le résultat

# Ajouter les espèces
exp_sp <- 9 # facteur d'expansion pour les coordonnées des vecteurs-espèces (paramètre esthétique)

acp_graph_sites_sp <- acp_graph_sites + # on utilise le graphique des sites
  geom_text(data = especes, aes(x = PC1 / exp_sp, # x, divisé par facteur expansion
                                y = PC2 / exp_sp, # y , divisé par facteur expansion
                                label = codesp)) + # on ajoute les étiquettes d'espèces
  geom_segment(data = especes, aes(x = 0, y = 0, # on ajoute des lignes pour représenter vecteurs
                                   xend = PC1 / exp_sp, # même chose que pour texte
                                   yend = PC2 / exp_sp))

acp_graph_sites_sp # voyons le résultat

dev.off()

# Enregistrement de la figure ACP
ggsave(acp_graph_sites_sp, file = 'figures/ACP_3comm.png', # nom du fichier avec extension
       width = 8, # largeur en pouces
       height = 6, # hauteur en pouces
       dpi = 300) # résolution en pixels par pouce

#### Ordination canonique - sélection des variables significatives en RDA ----
source('functions/triplot.rda.R')
## Selection of explanatory variables in RDA (Practicals using the R language p.31)
dim(spe.transfo) # Response variables
dim(grid.env2) # Explanatory variables

# Sélection des variables avec forward.sel
var.sel = forward.sel(spe.transfo, grid.env2)

write.csv(var.sel,"data_output/sel.var.csv", row.names = FALSE)

# Run rda on selected variables
(spe.rda.sel <- rda(spe.transfo ~ ., grid.env2[, var.sel$order[1:5]]))

summary(spe.rda.sel)

grid.env2.sel <- grid.env2[, var.sel$order[1:5]] # On crée un sous-jeu de données pour les variables significatives

# Canonical coefficients from the rda object
coef(spe.rda.sel)

# Unadjusted R^2 retrieved from the rda object
(R2 <- RsquareAdj(spe.rda.sel)$r.squared) # 0.1500596

# Adjusted R^2 retrieved from the rda object
(R2adj <- RsquareAdj(spe.rda.sel)$adj.r.squared) # 0.139328

# Global test of the RDA result
anova(spe.rda.sel, permutations = how(nperm = 999)) # significatif

# Tests of all canonical axes
anova(spe.rda.sel, by = "axis", permutations = how(nperm = 999)) 

# Pourcentage de la variance expliqué par les premiers axes
(pourc_var1_rda = R2adj*0.68381) # 0.09527391 ou 9.53%
(pourc_var2_rda = R2adj*0.180955) # 0.02521211 ou 2.52%
  
(pourc_var1_rda+pourc_var2_rda) # 0.120486 ou 12.05%

# Les deux premiers axes expliquent 12.05% de la variance.

## Apply Kaiser-Guttman criterion to residual axes (NEwR: p.221)
#spe.rda.sel$CA$eig[spe.rda.sel$CA$eig > mean(spe.rda.sel$CA$eig)]

## Triplots of the rda results (lc scores) ----
## Site scores as linear combinations of the environmental variables

source('functions/triplot.rda.R')
# Triplots with function triplot.rda(), scalings 1 and 2, lc scores

# Figure : "RDA plot with triplot.rda"
png("figures/RDA_triplot.rda.lc.png",
    height = 8, width = 16, # taille
    units = "in", # unités
    res = 300) # résolution

par(mfrow = c(1, 2))
triplot.rda(spe.rda.sel, 
            site.sc = "lc", 
            scaling = 1, 
            cex.char1 = 1,
            cex.char2 = 1,
            pos.env = 3, 
            pos.centr = 1, 
            mult.arrow = 1.1, 
            mar.percent = 0.05,
            label.sites = FALSE,
)
text(-0.68, 0.58, "a", cex = 2)
triplot.rda(spe.rda.sel, 
            site.sc = "lc", 
            scaling = 2,
            cex.char1 = 1,
            cex.char2 = 1,
            pos.env = 3, 
            pos.centr = 1, 
            mult.arrow = 1.1, 
            mar.percent = 0.05,
            label.sites = FALSE,
)
text(-4.9, 3.7, "b", cex = 2)
dev.off()

library(svglite)
library(spocc)
svglite("figures/RDA_triplot.rda.lc.svg",
        height = 8, width = 16)
par(mfrow = c(1, 2))
triplot.rda(spe.rda.sel, 
            site.sc = "lc", 
            scaling = 1, 
            cex.char1 = 1,
            cex.char2 = 1,
            pos.env = 3, 
            pos.centr = 1, 
            mult.arrow = 1.1, 
            mar.percent = 0.1,
            label.sites = FALSE,
)
text(-0.68, 0.58, "a", cex = 2)
triplot.rda(spe.rda.sel, 
            site.sc = "lc", 
            scaling = 2,
            cex.char1 = 1,
            cex.char2 = 1,
            pos.env = 3, 
            pos.centr = 1, 
            mult.arrow = 1.1, 
            mar.percent = 0.1,
            label.sites = FALSE,
)
text(-4.9, 3.7, "b", cex = 2)
dev.off()

# Triplots of the rda results (wa scores) ----
## Site scores as weighted averages (vegan's default)

# Triplots with function triplot.rda(), scalings 1 and 2, wa scores

# Figure : "RDA plot with triplot.rda"
png("figures/RDA_triplot.rda.wa.png",
    height = 8, width = 16, # taille
    units = "in", # unités
    res = 300) # résolution

par(mfrow = c(1, 2))
triplot.rda(spe.rda.sel, 
            site.sc = "wa", 
            scaling = 1, 
            cex.char1 = 1,
            cex.char2 = 1, 
            pos.env = 3, 
            pos.centr = 1, 
            mult.arrow = 1.1, 
            mar.percent = 0.05,
            label.sites = FALSE,
)
text(-0.97, 0.8, "a", cex = 2)
triplot.rda(spe.rda.sel, 
            site.sc = "wa", 
            scaling = 2, 
            cex.char1 = 1,
            cex.char2 = 1, 
            pos.env = 3, 
            pos.centr = 1, 
            mult.arrow = 1.1, 
            mar.percent = 0.05,
            label.sites = FALSE,
)
text(-7.92, 5.4, "b", cex = 2)
dev.off()

#### Espèces communes ----
sp_communes <- surfs %>% # on part du tableau de surface terrières par espèce/station
  ungroup() %>% # on défait les groupes définis précédemment
  group_by(Label) %>% # on groupe par espèce
  count(.) %>% # On dénombre les stations où se retrouve chaque espèce
  # Le . est pour passer l'objet temporaire au complet
  mutate(freq_rel =  round(n / 402, digits = 2)) %>% # on arrondit à la hausse
  filter(freq_rel >= 0.05) # on filtre les espèces présentes dans au moins 5% des stations
sp_communes

sp_rares <- c('FRNI', 'POBA', 'POTR', 'QURU')

#### Préparation des données pour les analyses dbMEM, RLQ et fourth corner ----
# Construction d'un fichier grid_coord avec des coordonnées en latitude, longitude
grid_lonlat <- grid %>% # on lit le fichier
  st_transform(4617) # on transforme en NAD83(CSRS)
rownames(grid_lonlat) <- grid_lonlat$grid_id # on utilise les grid_id comme noms de lignes

grid_coord_lonlat = st_centroid(grid_lonlat) %>% # On extrait les centroïdes des placettes
  st_coordinates() # On extrait uniquement la matrice des coordonnées (pas d'objet sf)

grid_coord_lonlat <-  data.frame(X = grid_coord_lonlat[,1],
                                 Y = grid_coord_lonlat[,2],
                                 grid_id = grid_id)
rownames(grid_coord_lonlat) <- grid_coord_lonlat$grid_id # on utilise les grid_id comme noms de lignes

# Transformation des coordonnées latlon en coordonnées cartésiennes
library(SoDA)
xygrid = geoXY(grid_coord_lonlat[,2], grid_coord_lonlat[,1]) # package SoDA

# Visualisation des données
head(spe.transfo) # données transformées des abondances d'espèces
head(xygrid) # coordonnées cartésiennes des placettes
head(grid.env2.sel) # données environnementales et topographiques des placettes

#### RLQ and fourth-corner analysis (Chap.6 - NEwR) ----
# Load data ----
# On enlève les espèces rares
spe.transfo.sel <- spe.transfo[c("ABBA","ACRU","BEPA","ACSA","THOC","FAGR","ACPE","Picea","TSCA","POGR","BEAL","PIST","LALA")]
head(spe.transfo.sel)

spe.traits <- read.csv("traits_output/traits_short.csv", row.names=1) # rn pour rownames (garde les grid_id)

spe.traits.sel <- spe.traits[!(row.names(spe.traits) %in% sp_rares),] # On enlève les espèces rares
head(spe.traits.sel) # données des traits des espèces

# Standardization of all traits variables
# Center and scale = standardize the variables (z-scores)
spe.traits.sel <- as.data.frame(scale(spe.traits.sel)) # means = 0 // standard deviations = 1

# Nos 3 tables de départ sont donc:
# spe.transfo.sel (L), grid.env2.sel (R), spe.traits.sel (matrix Q)

# Tri dans les traits ----
# Corrélation des traits foliaires
traits.cor = cor(spe.traits.sel, method = c("spearman"))

library("Hmisc")
library("corrplot")

png("figures/corrplot_traits.png",
    height = 6, width = 6, # taille
    units = "in", # unités
    res = 300) # résolution
corrplot(traits.cor, 
         method = "color", 
         type = 'lower', 
         insig='blank',
         addCoef.col ='black',
         order="hclust", 
         diag=TRUE)
dev.off()

# Scatter plots for all pairs of traits variables (NEwR, chap.2)

# Bivariate plots with histograms on the diagonal and smooth 
# fitted curves
dev.new(
  title = "Bivariate descriptor plots - traits",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
pairs(spe.traits.sel, 
      panel = panel.smooth, 
      diag.panel = panel.hist,
      main = "Bivariate Plots with Histograms and Smooth Curves"
)

# Simple transformation of a trait variable : C:N ==============

range(spe.traits.sel$C.N)
# Log-transformation of the C:N variable (y = ln(x))
# Compare histograms and boxplots of the raw and transformed values
dev.new(title = "Transformation and standardization of variable C:N", 
        noRStudioGD = TRUE)
par(mfrow = c(2, 2))
hist(spe.traits.sel$C.N, 
     col = "bisque", 
     right = FALSE
)
hist(log(spe.traits.sel$C.N), 
     col = "light green", 
     right = FALSE, 
     main = "Histogram of ln(spe.traits.sel$C.N)"
)
boxplot(spe.traits.sel$C.N, 
        col = "bisque", 
        main = "Boxplot of spe.traits.sel$C.N", 
        ylab = "spe.traits.sel$C.N"
)
boxplot(log(spe.traits.sel$C.N), 
        col = "light green", 
        main = "Boxplot of ln(spe.traits.sel$C.N)",
        ylab = "log(spe.traits.sel$C.N)"
)

# Simple transformation of a trait variable : Chla:Chlb ==============

range(spe.traits.sel$Chla.Chlb)
# Log-transformation of the Chla.Chlb variable (y = ln(x))
# Compare histograms and boxplots of the raw and transformed values
dev.new(title = "Transformation and standardization of variable Chla.Chlb", 
        noRStudioGD = TRUE)
par(mfrow = c(2, 2))
hist(spe.traits.sel$Chla.Chlb, 
     col = "bisque", 
     right = FALSE
)
hist(log(spe.traits.sel$Chla.Chlb), 
     col = "light green", 
     right = FALSE, 
     main = "Histogram of ln(spe.traits.sel$Chla.Chlb)"
)
boxplot(spe.traits.sel$Chla.Chlb, 
        col = "bisque", 
        main = "Boxplot of spe.traits.sel$Chla.Chlb", 
        ylab = "spe.traits.sel$Chla.Chlb"
)
boxplot(log(spe.traits.sel$Chla.Chlb), 
        col = "light green", 
        main = "Boxplot of ln(spe.traits.sel$Chla...Chlb)",
        ylab = "log(spe.traits.sel$Chla...Chlb)"
)

# Simple transformation of a trait variable : Carot:Chltot ==============

range(spe.traits.sel$Carot.Chltot)
# Log-transformation of the Carot.Chltot variable (y = ln(x))
# Compare histograms and boxplots of the raw and transformed values
dev.new(title = "Transformation and standardization of variable Carot.Chltot", 
        noRStudioGD = TRUE)
par(mfrow = c(2, 2))
hist(spe.traits.sel$Carot.Chltot, 
     col = "bisque", 
     right = FALSE
)
hist(log(spe.traits.sel$Carot.Chltot), 
     col = "light green", 
     right = FALSE, 
     main = "Histogram of ln(spe.traits.sel$Carot.Chltot)"
)
boxplot(spe.traits.sel$Carot.Chltot, 
        col = "bisque", 
        main = "Boxplot of spe.traits.sel$Carot.Chltot", 
        ylab = "spe.traits.sel$Carot.Chl.tot"
)
boxplot(log(spe.traits.sel$Carot.Chltot), 
        col = "light green", 
        main = "Boxplot of ln(spe.traits.sel$Carot.Chltot)",
        ylab = "log(spe.traits.sel$Carot.Chltot)"
)
# Preliminary analyses: CA, PCA and PCA ----
afcL <- dudi.coa(spe.transfo.sel, scannf = FALSE)
acpR <- dudi.pca(grid.env2.sel, 
                 row.w = afcL$lw,
                 scannf = FALSE)
acpQ <- dudi.pca(spe.traits.sel,
                 row.w = afcL$cw,
                 scannf = FALSE)

# RLQ analysis ----
rlq.sbl <- rlq(       # package ade4
  dudiR = acpR, 
  dudiL = afcL, 
  dudiQ = acpQ,
  scannf = FALSE)

# plot(rlq.sbl) # trop chargé, pas pertinent

# Traits by environment crossed table
rlq.sbl$tab

# Since the plots are crowded, one can plot them one by one 
# in large graphical windows.

# Figure RLQ - site (L) scores
png("figures/RLQ_L_site_scores.png",
    height = 4, width = 4, # taille
    units = "in", # unités
    res = 300) # résolution
s.label(rlq.sbl$lR,             # package adegraphics
        plabels.boxes.draw = FALSE, 
        ppoints.alpha = 0,
        psub.text = "a",
        psub.cex = 2, 
        psub.position = "topleft"
)
dev.off()

# Figure RLQ - species (Q) abundances
png("figures/RLQ_Q_species_scores.png",
    height = 4, width = 4, # taille
    units = "in", # unités
    res = 300) # résolution
s.label(rlq.sbl$lQ,
        plabels.boxes.draw = FALSE, 
        ppoints.alpha = 0,
        psub.text = "b",
        psub.cex = 2, 
        psub.position = "topleft"
)
dev.off()

# Figure RLQ - environmental variables
png("figures/RLQ_environmental_var.png",
    height = 4, width = 4, # taille
    units = "in", # unités
    res = 300) # résolution
s.arrow(rlq.sbl$l1,               # package adegraphics
        psub.text = "c",
        psub.cex = 2, 
        psub.position = "topleft"
)
dev.off()

# Figure RLQ - species traits
png("figures/RLQ_species_traits.png",
    height = 4, width = 4, # taille
    units = "in", # unités
    res = 300) # résolution
s.arrow(rlq.sbl$c1,
        psub.text = "d",
        psub.cex = 2, 
        psub.position = "topleft"
)
dev.off()

# Global test
# model 6: after ter Braak et al. 2012
randtest(rlq.sbl, nrepet = 999, modeltype = 6) # package ade4      # voir ce que ça ferait juste avec modele 2

# H1 = les traits et l'environnement sont reliés.
# H0 = les traits et l'environnement ne sont pas reliés.

# Au niveau de significativité a = 0.05, si p-value <= 0.05, l'hypothèse nulle est rejetée.

# Le premier test résulte en une p-value de 0.001. L'hypothèse nulle est donc rejetée.
# Cela signifie que les liens L - Q sont significatifs.

# Par contre, le deuxième test résulte en une p-value de 0.843. L'hypothèse nulle est donc acceptée.
# Cela signifie que les liens L - R ne sont pas significatifs.

# Article cité:
# ter Braak, C., Cormont, A., Dray, S.: Improved testing of species traits–environment relationships in the fourth corner problem. Ecology. 93, 1525–1526 (2012)

# « Dray and Legendre (2008) identified that under H1 both the
# links L-Q and L-R must exist, and that H0: X = 0 holds if
# at least one of the links is absent. »

# « The overall null hypothesis is thus rejected when both nulls
# (L != Q and L !=R) are rejected. »

# Ainsi, l'hypothèse nulle globale est acceptée, car les deux hypothèses nulles 
# indépendantes doivent être rejetées afin de prouver une significativité globale.
# Cela signifie donc que les traits et l'environnement ne sont pas reliés.

# Fourth-corner analysis (takes time with 49999 permutations!) ----
fourth.sbl <- fourthcorner(     # package ade4
  tabR = grid.env2.sel, 
  tabL = spe.transfo.sel, 
  tabQ = spe.traits.sel,
  modeltype = 6,
  p.adjust.method.G = "none", 
  p.adjust.method.D = "none", 
  nrepet = 999) # c'était 49999 avant, mais pas assez rapide

# Correction for multiple testing, here using FDR
fourth.sbl.adj <- p.adjust.4thcorner(     # package ade4
  fourth.sbl,
  p.adjust.method.G = "fdr", 
  p.adjust.method.D = "fdr", 
  p.adjust.D = "global") 

# Plot significant associations
png("figures/fourth-corner_analysis.png",
    height = 4, width = 4, # taille
    units = "in", # unités
    res = 300) # résolution
plot(fourth.sbl.adj, alpha = 0.05, stat = "D2")
dev.off()

# Il n'y a aucune association significative.

# Biplot combining RLQ and fourth-corner results
png("figures/RLQ_fourth-corner.png",
    height = 4, width = 4, # taille
    units = "in", # unités
    res = 300) # résolution
plot(fourth.sbl.adj, 
     x.rlq = rlq.sbl, 
     alpha = 0.05, 
     stat = "D2", 
     type = "biplot")
dev.off()

# Test RLQ and Fourth-corner associations avec modeltype = 2 ----
# Global test # model 2
randtest(rlq.sbl, nrepet = 999, modeltype = 2) # package ade4 

# Fourth-corner analysis (takes time with 49999 permutations!) ----
fourth.sbl.mod2 <- fourthcorner(     # package ade4
  tabR = grid.env2.sel, 
  tabL = spe.transfo.sel, 
  tabQ = spe.traits.sel,
  modeltype = 2,
  p.adjust.method.G = "none", 
  p.adjust.method.D = "none", 
  nrepet = 999) # c'était 49999 avant, mais pas assez rapide

# Correction for multiple testing, here using FDR
fourth.sbl.mod2.adj <- p.adjust.4thcorner(     # package ade4
  fourth.sbl.mod2,
  p.adjust.method.G = "fdr", 
  p.adjust.method.D = "fdr", 
  p.adjust.D = "global") 

# Plot significant associations
png("figures/fourth-corner_analysis_mod2.png",
    height = 4, width = 4, # taille
    units = "in", # unités
    res = 300) # résolution
plot(fourth.sbl.mod2.adj, alpha = 0.05, stat = "D2")
dev.off()

# Biplot combining RLQ and fourth-corner results
png("figures/RLQ_fourth-corner_mod2.png",
    height = 4, width = 4, # taille
    units = "in", # unités
    res = 300) # résolution
plot(fourth.sbl.mod2.adj, 
     x.rlq = rlq.sbl, 
     alpha = 0.05, 
     stat = "D2", 
     type = "biplot")
dev.off()

#### dbMEM (Chap.7 - NEwR) ----
# Load functions ----
source("functions/sr.value.R")

# Compute the distance-based Moran maps ----
# Creation of the dbMEM eigenfunctions with positive Moran's I 
xygrid.dbmem.tmp <- dbmem(xygrid)
xygrid.dbmem <- as.data.frame(xygrid.dbmem.tmp)
# Count the eigenvalues
length(attributes(xygrid.dbmem.tmp)$values)

# Plot some dbMEM variables using s.value {adegraphics}
png("figures/dbMEM_variables_grid.png",
    height = 6, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
somedbmem <- c(20, 25, 50, 80, 140, 198)
s.value(xygrid, xygrid.dbmem[ ,somedbmem], 
        method = "color", 
        symbol = "circle", 
        ppoints.cex = 0.4)
dev.off()

# Shades of grey:
# values in each eigenvector, from white (largest negative value)
# to black (largest positive value)

## Detrending ----
# Is there a significant linear spatial trend in the data?
anova(rda(spe.transfo, xygrid))	# Result: significant trend

spe.trend.rda <- rda(spe.transfo, xygrid)
(spe.trend.R2a <- RsquareAdj(spe.trend.rda)$adj.r.squared) # 0.1234812

# Computation of linearly detrended spe data - model residuals
spe.transfo.det <- resid(lm(as.matrix(spe.transfo) ~ as.matrix(xygrid)))

## Mantel correlogram ## Ne pas utiliser: Legendre explique que ce n'est pas fiable dans le cours 12 ----
# spe.transfo.D1 <- dist(spe.transfo.det)
# (spe.correlog <- 
#     mantel.correlog(spe.transfo.D1, 
#                     XY = xygrid, 
#                     nperm = 999))
# summary(spe.correlog)
# 
# # Number of classes
# spe.correlog$n.class # or: mite.correlog[2]
# # Break points
# spe.correlog$break.pts # or: mite.correlog[3]
# 
# # Figure : Mantel correlogram of spe data
# png("figures/correlogram.png",
#     height = 5, width = 8, # taille
#     units = "in", # unités
#     res = 300) # résolution
# plot(spe.correlog)
# dev.off()

## Step 1. Construct the matrix of dbMEM variables ----
spe.dbmem.tmp <- dbmem(xygrid, silent = FALSE) # package adespatial
spe.dbmem <- as.data.frame(spe.dbmem.tmp)
# Truncation distance used above:
(thr <- give.thresh(dist(xygrid))) # give.thresh: package adespatial; dist: package stats

# Display and count the eigenvalues
attributes(spe.dbmem.tmp)$values
length(attributes(spe.dbmem.tmp)$values)
# Argument silent = FALSE allows the function to display the truncation level.

## Step 2. Run the global dbMEM analysis on the *detrended* data ----
(spe.dbmem.rda <- rda(spe.transfo.det ~ ., spe.dbmem))
anova(spe.dbmem.rda)

## Step 3. Selection of significant MEM eigenfunctions ----
## Since the R-square is significant, compute the adjusted
## R2 and run a forward selection of the dbmem variables

# Adjusted R-square of the dbMEM model
(spe.R2a <- RsquareAdj(spe.dbmem.rda)$adj.r.squared) # package vegan
# 0.5099192

# Selection of significant MEM eigenfunctions
(spe.dbmem.fwd <- forward.sel(spe.transfo.det, as.matrix(spe.dbmem), # package adespatial
                              adjR2thresh = spe.R2a))

# Number of signif. dbMEM
(nb.sig.dbmem <- nrow(spe.dbmem.fwd)) # 87 - 97

# Identity of the significant dbMEM in increasing order
(dbmem.sign <- sort(spe.dbmem.fwd[ ,2]))

# Write the significant dbMEM to a new object
dbmem.red <- spe.dbmem[ ,c(dbmem.sign)]

## Step 4. New dbMEM analysis with 93 significant dbMEM variables ----
(spe.dbmem.rda2 <- rda(spe.transfo.det ~ ., data = dbmem.red))
anova(spe.dbmem.rda2)

# Adjusted R-square of the dbMEM model after forward selection: R2adj = 0.4738661 --> 47.39%
(spe.fwd.R2a <- RsquareAdj(spe.dbmem.rda2)$adj.r.squared)

# Ordination de la variation spatiale sur les variables environnementales
(spe.dbmemgeneral.env = forward.sel(dbmem.red, grid.env2.sel))
spe.dbmemgeneral.env[,1] # Extraire les variables significatives

# Test of significance of individual axes
(axes.test <- anova(spe.dbmem.rda2, by = "axis"))

# Number of significant axes
(nb.ax <- length(which(axes.test[ ,ncol(axes.test)] <=  0.05))) # 7

# Proportion de la variance expliquée par les 7 axes significatifs
(pourc_var1=(spe.dbmem.rda2$CCA$eig[1])/(sum(spe.dbmem.rda2$CCA$eig))) # 0.2515594
(pourc_var2=(spe.dbmem.rda2$CCA$eig[2])/(sum(spe.dbmem.rda2$CCA$eig))) # 0.155471   
(pourc_var3=(spe.dbmem.rda2$CCA$eig[3])/(sum(spe.dbmem.rda2$CCA$eig))) # 0.1404764 
(pourc_var4=(spe.dbmem.rda2$CCA$eig[4])/(sum(spe.dbmem.rda2$CCA$eig))) # 0.1106959 
(pourc_var5=(spe.dbmem.rda2$CCA$eig[5])/(sum(spe.dbmem.rda2$CCA$eig))) # 0.09790925  
(pourc_var6=(spe.dbmem.rda2$CCA$eig[6])/(sum(spe.dbmem.rda2$CCA$eig))) # 0.07137635 
(pourc_var7=(spe.dbmem.rda2$CCA$eig[7])/(sum(spe.dbmem.rda2$CCA$eig))) # 0.05776293

(pourc_var1a7 <- pourc_var1+pourc_var2+pourc_var3+pourc_var4+pourc_var5+pourc_var6+pourc_var7) # 0.8852512 
(spe.fwd.R2a*pourc_var1a7) # 0.4194906  -> 41.95%

# Les sept axes significatifs expliquent (47.39*0.88525) = 41.95% de la variance totale.

## Step 5. Plot the significant canonical axes ----
spe.rda2.axes <- 
  scores(spe.dbmem.rda2, 
         choices = c(1:nb.ax), 
         display = "lc", 
         scaling = 1)

# Figure : dbMEM analysis of spe data
png("figures/dbMEM_analysis.png",
    height = 8, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
#par(mfrow = c(3,(nb.ax)/2))
par(mfrow = c(3,3))
for(i in 1:nb.ax){
  sr.value(xygrid, spe.rda2.axes[ ,i], 
           sub = paste("RDA",i), 
           csub = 2)
}
dev.off()

# « The dbMEM analysis of the detrended spe data explained 46.12% of the variance
# (see spe.fwd.R2a, the adjusted R2 obtained with 97 dbMEM variables retained by forward selection).
# Six canonical axes, explaining 46.12% X 0.88415 = 40.79 % of the total variance, are significant;
# their fitted site scores have been plotted on a map of the sites.
# These six plots (dbMEM_analysis.png) represent the spatially structured variation
# of the detrended spe data. »

# Interprétation de la variation spatiale ----
# Ordination des axes canoniques significatifs sur les variables environnementales
# Axe 1 :
spe.rda2.axis1.env = forward.sel(spe.rda2.axes[ ,1], grid.env2.sel)
spe.rda2.axis1.env[,1] # Extraire les variables significatives
# "elev_mnt"      "eastness_exp" sont significatifs.

# Axe 2 :
spe.rda2.axis2.env = forward.sel(spe.rda2.axes[ ,2], grid.env2.sel)
spe.rda2.axis2.env[,1] # Extraire les variables significatives
# "twi" "northness_exp" "eastness_exp" sont significatifs.

# Axe 3 :
spe.rda2.axis3.env = forward.sel(spe.rda2.axes[ ,3], grid.env2.sel)
spe.rda2.axis3.env[,1] # Extraire les variables significatives
# "northness_exp""twi" sont significatifs.

# Axe 4 :
spe.rda2.axis4.env = forward.sel(spe.rda2.axes[ ,4], grid.env2.sel)
spe.rda2.axis4.env[,1] # Extraire les variables significatives
# "twi" "northness_exp" sont significatifs.

# Axe 5 :
spe.rda2.axis5.env = forward.sel(spe.rda2.axes[ ,5], grid.env2.sel)
spe.rda2.axis5.env[,1] # Extraire les variables significatives
# Aucune variable significative.

# Axe 6 :
spe.rda2.axis6.env = forward.sel(spe.rda2.axes[ ,6], grid.env2.sel)
spe.rda2.axis6.env[,1] # Extraire les variables significatives
# "northness_exp" est significatif.

# Axe 7 :
spe.rda2.axis7.env = forward.sel(spe.rda2.axes[ ,7], grid.env2.sel)
spe.rda2.axis7.env[,1] # Extraire les variables significatives
# "elev_mnt"  "northness_exp" "slope"     "twi"  sont significatifs.

# Séparation des dbMEM en différentes échelles spatiales ----
# Scalogram of the variance explained by all dbMEM eigenfunctions,
# computed with our homemade function scalog()
source('functions/scalog.R')
png("figures/scalogram.png",
    height = 6, width = 10, # taille
    units = "in", # unités
    res = 300) # résolution
scalog(spe.dbmem.rda)
dev.off()

# Pour ouvrir la figure dans inkscape
library(svglite)
library(spocc)
svglite("figures/scalogram_noir.svg")
scalog(spe.dbmem.rda)
dev.off()

# NEwR p.325
# Broad-scale : 1-14
# Medium-scale : 15-63
# Fine-scale : 68-115
# Very Fine-scale : 121-200

dbmem.sign
dbmem.red
spe.dbmem

## dbMEM analysis of the spe data - broad scale (1-14) ----
broad.scale <- c(1:14)

(spe.dbmem.broad <- 
   rda(spe.transfo.det ~ ., data = spe.dbmem[ ,broad.scale]))
anova(spe.dbmem.broad)

(spe.dbmem.broad.R2a <- RsquareAdj(spe.dbmem.broad)$adj.r.squared)

# Interpreting the broad-scaled spatial variation (general) 
# Ordination de la variation spatiale sur les variables environnementales
(spe.dbmembroad.env = forward.sel(spe.dbmem[ ,broad.scale], grid.env2.sel))
spe.dbmembroad.env[,1] # Extraire les variables significatives
# Toutes les variables

# Interpreting the broad-scaled spatial variation (axes) ----
# Test des axes
(axes.broad <- anova(spe.dbmem.broad, by = "axis"))
# Number of significant axes
(nb.ax.broad <- 
    length(which(axes.broad[ , ncol(axes.broad)] <=  0.05))) # 6

# Plot of the six significant canonical axes
spe.dbmembroad.axes <- 
  scores(spe.dbmem.broad, 
         choices = 1:6, 
         display = "lc", 
         scaling = 1)

# Figure : dbMEM analysis of spe data - broad scale
png("figures/dbMEM_broadscale.png",
    height = 8, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
par(mfrow = c(2, 3))
sr.value(xygrid, spe.dbmembroad.axes[ ,1])
sr.value(xygrid, spe.dbmembroad.axes[ ,2])
sr.value(xygrid, spe.dbmembroad.axes[ ,3])
sr.value(xygrid, spe.dbmembroad.axes[ ,4])
sr.value(xygrid, spe.dbmembroad.axes[ ,5])
sr.value(xygrid, spe.dbmembroad.axes[ ,6])
dev.off()

####
  
# Ordination des axes canoniques significatifs sur les variables environnementales
# Axe 1 :
spe.dbmembroad.ax1.env = forward.sel(spe.dbmembroad.axes[ ,1], grid.env2.sel)
spe.dbmembroad.ax1.env[,1] # Extraire les variables significatives
# "elev_mnt" "eastness_exp" sont significatifs.

# Axe 2 :
spe.dbmembroad.ax2.env = forward.sel(spe.dbmembroad.axes[ ,2], grid.env2.sel)
spe.dbmembroad.ax2.env[,1] # Extraire les variables significatives
# "twi" "northness_exp" "eastness_exp" sont significatifs.

# Axe 3 :
spe.dbmembroad.ax3.env = forward.sel(spe.dbmembroad.axes[ ,3], grid.env2.sel)
spe.dbmembroad.ax3.env[,1] # Extraire les variables significatives
# "northness_exp" "slope"  "elev_mnt" sont significatifs.

# Axe 4 :
spe.dbmembroad.ax4.env = forward.sel(spe.dbmembroad.axes[ ,4], grid.env2.sel)
spe.dbmembroad.ax4.env[,1] # Extraire les variables significatives
# "northness_exp" "twi" "slope"  sont significatifs.

# Axe 5 :
spe.dbmembroad.ax5.env = forward.sel(spe.dbmembroad.axes[ ,5], grid.env2.sel)
spe.dbmembroad.ax5.env[,1] # Extraire les variables significatives
# "eastness_exp""northness_exp" est significatif.

# Axe 6 :
spe.dbmembroad.ax6.env = forward.sel(spe.dbmembroad.axes[ ,6], grid.env2.sel)
spe.dbmembroad.ax6.env[,1] # Extraire les variables significatives
# "elev_mnt"  "northness_exp" est significatif.

## dbMEM analysis of the spe data - medium scale (15-63) ----
med.scale <- c(15,16,17,18,19,20,21,22,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,54,55,57,58,60,61,62,63)

(spe.dbmem.med <- 
   rda(spe.transfo.det ~ ., data = spe.dbmem[ ,med.scale]))
anova(spe.dbmem.med)

(spe.dbmem.med.R2a <- RsquareAdj(spe.dbmem.med)$adj.r.squared)

# Interpreting the medium-scaled spatial variation (general) 
# Ordination de la variation spatiale sur les variables environnementales
(spe.dbmemmed.env = forward.sel(spe.dbmem[ ,med.scale], grid.env2.sel))
spe.dbmemmed.env[,1] # Extraire les variables significatives
# Toutes les variables

# Interpreting the medium-scaled spatial variation (axes) ----
# Test des axes
(axes.med <- anova(spe.dbmem.med, by = "axis"))
# Number of significant axes
(nb.ax.med <- 
    length(which(axes.med[ , ncol(axes.med)] <=  0.05))) # 3

# Plot of the six significant canonical axes
spe.dbmemmed.axes <- 
  scores(spe.dbmem.med, 
         choices = 1:3, 
         display = "lc", 
         scaling = 1)

# Figure : dbMEM analysis of spe data - medium scale
png("figures/dbMEM_mediumscale.png",
    height = 4, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
par(mfrow = c(1, 3))
sr.value(xygrid, spe.dbmemmed.axes[ ,1])
sr.value(xygrid, spe.dbmemmed.axes[ ,2])
sr.value(xygrid, spe.dbmemmed.axes[ ,3])
dev.off()

####

# Ordination des axes canoniques significatifs sur les variables environnementales
# Axe 1 :
spe.dbmemmed.ax1.env = forward.sel(spe.dbmemmed.axes[ ,1], grid.env2.sel)
spe.dbmemmed.ax1.env[,1] # Extraire les variables significatives
# "northness_exp" "twi" sont significatifs.

# Axe 2 :
spe.dbmemmed.ax2.env = forward.sel(spe.dbmemmed.axes[ ,2], grid.env2.sel)
spe.dbmemmed.ax2.env[,1] # Extraire les variables significatives
# "eastness_exp" "elev_mnt"  "twi"  sont significatifs.

# Axe 3 :
spe.dbmemmed.ax3.env = forward.sel(spe.dbmemmed.axes[ ,3], grid.env2.sel)
spe.dbmemmed.ax3.env[,1] # Extraire les variables significatives
# "twi" sont significatifs.

## dbMEM analysis of the spe data - fine scale (68-115) ----
fine.scale <- c(68,69,72,73,75,77,78,79,80,84,85,86,87,90,93,96,99,100,101,104,105,106,107,109,111,114,115)

(spe.dbmem.fine <- 
    rda(spe.transfo.det ~ ., data = spe.dbmem[ ,fine.scale]))
anova(spe.dbmem.fine)

(spe.dbmem.fine.R2a <- RsquareAdj(spe.dbmem.fine)$adj.r.squared)

# Interpreting the fine-scaled spatial variation (general) 
# Ordination de la variation spatiale sur les variables environnementales
(spe.dbmemfine.env = forward.sel(spe.dbmem[ ,fine.scale], grid.env2.sel))
spe.dbmemfine.env[,1] # Extraire les variables significatives
# Aucune variable

# Interpreting the fine-scaled spatial variation (axes) ----
# Test des axes
(axes.fine <- anova(spe.dbmem.fine, by = "axis"))
# Number of significant axes
(nb.ax.fine <- 
    length(which(axes.fine[ , ncol(axes.fine)] <=  0.05))) # 1

# Plot of the six significant canonical axes
spe.dbmemfine.axes <- 
  scores(spe.dbmem.fine, 
         choices = 1, 
         display = "lc", 
         scaling = 1)

# Figure : dbMEM analysis of spe data - fine scale
png("figures/dbMEM_finescale.png",
    height = 4, width = 4, # taille
    units = "in", # unités
    res = 300) # résolution
par(mfrow = c(1, 1))
sr.value(xygrid, spe.dbmemfine.axes[ ,1])
dev.off()

####

# Ordination des axes canoniques significatifs sur les variables environnementales
# Axe 1 :
spe.dbmemfine.ax1.env = forward.sel(spe.dbmemfine.axes[ ,1], grid.env2.sel)
spe.dbmemfine.ax1.env[,1] # Extraire les variables significatives
# "twi" est significatif.

## dbMEM analysis of the spe data - very fine scale (123-193) ----
veryfine.scale <- c(123,124,128,129,134,137,145,147,151,159,165,183,189,193)

(spe.dbmem.veryfine <- 
    rda(spe.transfo.det ~ ., data = spe.dbmem[ ,veryfine.scale]))
anova(spe.dbmem.veryfine)

(spe.dbmem.veryfine.R2a <- RsquareAdj(spe.dbmem.veryfine)$adj.r.squared)

# Interpreting the veryfine-scaled spatial variation (general) 
# Ordination de la variation spatiale sur les variables environnementales
(spe.dbmemveryfine.env = forward.sel(spe.dbmem[ ,veryfine.scale], grid.env2.sel))
spe.dbmemveryfine.env[,1] # Extraire les variables significatives
# Aucune variable

# Interpreting the veryfine-scaled spatial variation (axes) ----
# Test des axes
(axes.veryfine <- anova(spe.dbmem.veryfine, by = "axis"))
# Number of significant axes
(nb.ax.veryfine <- 
    length(which(axes.veryfine[ , ncol(axes.veryfine)] <=  0.05))) # 1

# Plot of the six significant canonical axes
spe.dbmemveryfine.axes <- 
  scores(spe.dbmem.veryfine, 
         choices = 1, 
         display = "lc", 
         scaling = 1)

# Figure : dbMEM analysis of spe data - very fine scale
png("figures/dbMEM_veryfinescale.png",
    height = 4, width = 4, # taille
    units = "in", # unités
    res = 300) # résolution
par(mfrow = c(1, 1))
sr.value(xygrid, spe.dbmemveryfine.axes[ ,1])
dev.off()

####

# Ordination des axes canoniques significatifs sur les variables environnementales
# Axe 1 :
spe.dbmemveryfine.ax1.env = forward.sel(spe.dbmemveryfine.axes[ ,1], grid.env2.sel)
spe.dbmemveryfine.ax1.env[,1] # Extraire les variables significatives
# Il n'y a aucune variable significative.

#### dbMEM variation partitioning (Spe, trend, environment) ----
# 1-4. Prepare data for variation partitioning ----
# 1. Test trend
spe.XY.rda <- rda(spe.transfo, xygrid)
anova(spe.XY.rda) # Significative tend

# 2. Test and forward selection of the environmental variables
head(grid.env2.sel) # Already done

# 3. Test and forward selection of the dbMEM variables
head(dbmem.red) # Already done

# 4. Arbitrarily split the significant dbMEM into 2 differents scales
dbmem.broad <- dbmem.red[ ,1:53] # Broad arbitrary choice: broad + medium
dbmem.fine <- dbmem.red[ ,54:91] # Fine arbitrary choice: fine + very fine

# 5. dbMEM variation partitioning (broad and fine arbitrary scales) ----
(spe.varpart <- 
   varpart(spe.transfo, grid.env2.sel, xygrid, dbmem.broad, dbmem.fine))

# Figure : Spe - environment - dbMEM variation partitioning
png("figures/dbMEM_varpart.png",
    height = 6, width = 12, # taille
    units = "in", # unités
    res = 300) # résolution
# Show the symbols of the fractions and plot their values
par(mfrow = c(1,2))
showvarparts(4, bg = c("red", "blue", "yellow", "green"))
plot(spe.varpart, 
     digits = 2, 
     bg = c("red", "blue", "yellow", "green")
) 
dev.off()

# Pour ouvrir la figure dans inkscape
library(svglite)
library(spocc)
svglite("figures/dbMEM_varpart.svg")
par(mfrow = c(1,2))
showvarparts(4, bg = c("red", "blue", "yellow", "green"))
plot(spe.varpart, 
     digits = 2, 
     bg = c("red", "blue", "yellow", "green")
) 
dev.off()

# Tests of the unique fractions [a], [b], [c] and [d] ----
# Fraction [a], pure environmental
anova(rda(spe.transfo, grid.env2.sel, cbind(xygrid, dbmem.broad, dbmem.fine))) # Significatif

# Fraction [b], pure trend
anova(rda(spe.transfo, xygrid, cbind(grid.env2.sel, dbmem.broad, dbmem.fine))) # Significatif

# Fraction [c], pure broad scale spatial
anova(rda(spe.transfo, dbmem.broad, cbind(grid.env2.sel, xygrid, dbmem.fine))) # Significatif

# Fraction [d], pure fine scale spatial
anova(rda(spe.transfo, dbmem.fine, cbind(grid.env2.sel, xygrid, dbmem.broad))) # Significatif

# Proportion de la variance expliquée par des variables environnementales ou spatiales ----
Rajust <- spe.varpart[["part"]][["fract"]][["Adj.R.square"]]
All.fractions <- length(Rajust)

Rajust[All.fractions]

# 54% de la variance est expliquée par les variables environnementales et spatiales.