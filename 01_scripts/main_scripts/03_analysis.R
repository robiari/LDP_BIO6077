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
spe <- read.csv("02_outdata/spe.csv", row.names=1) # rn pour rownames (garde les grid_id)
head(spe)

# Placettes forestières et données environnementales
grid_topo <- read_sf("02_outdata/grid_topo.geojson") # on lit le fichier

#### Édition du fichier grid_topo pour le faire fonctionner avec les analyses ----
grid_env <- grid_topo %>%
  st_drop_geometry() %>%
  as.data.frame()
rownames(grid_env) <- grid_env$grid_id # on utilise la colonne grid_id comme rownames

ncol(grid_env)
grid_env <- grid_env[4:ncol(grid_env)]

## Corrélation des variables environnementales ----
grid_env.cor = cor(grid_env, method = c("spearman"))

library("Hmisc")
library("corrplot")

png("03_figs/corrplot_env.png",
    height = 6, width = 6, # taille
    units = "in", # unités
    res = 300) # résolution
corrplot(grid_env.cor, 
         method = "color", 
         type = 'lower', 
         insig='blank',
         addCoef.col ='black',
         order="hclust", 
         diag=TRUE)
dev.off()

# Nous devons choisir entre eastness/northness ou eastness exposure/northness exposure.
# Nous choisissons avec exposure, car cela fait plus de sens d'un point de vue écologique.

grid_env <- grid_env %>%
  select(-c(pentes, eastness, northness))

head(grid_env)

#write.csv(grid_env,"02_outdata/grid_env.csv", row.names = TRUE)

# Statistiques sur les variables topographiques ----
env_stat <- data.frame('Moyenne' = apply(grid_env, 2, mean, na.rm = TRUE),
                       'Écart-type' = apply(grid_env, 2, sd, na.rm = TRUE),
                       'Étendue' = apply(grid_env, 2, max, na.rm = TRUE)-apply(grid_env, 2, min, na.rm = TRUE),
                       'Cov' = (apply(grid_env, 2, sd, na.rm = TRUE)/apply(grid_env, 2, mean, na.rm = TRUE))*100) %>%
  t() %>% # transposer le tableau
  as.data.frame()

#### Transformation des données environnementales ----
# Scatter plots for all pairs of environmental variables (NEwR, chap.2) 
source('01_scripts/r_functions/panelutils.R')
# Bivariate plots with histograms on the diagonal and smooth fitted curves
dev.new(
  title = "Bivariate descriptor plots - env",
  width = 10,
  height = 10,
  noRStudioGD = TRUE)
pairs(grid_env, 
      panel = panel.smooth, 
      diag.panel = panel.hist,
      main = "Bivariate Plots with Histograms and Smooth Curves")

# Standardization of all environmental variables (NEwR, chap.2)
# Center and scale = standardize the variables (z-scores)
grid_env <- as.data.frame(scale(grid_env)) # means = 0 // standard deviations = 1

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
          filename = "03_figs/nb_sp_nb_ind.png", # nom et chemin du fichier qu'on enregistre, incluant l'extension
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
png("03_figs/classes_abondance_avec0.png", height = 4, width = 5, 
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
png("03_figs/classes_abondance.png", height = 4, width = 5, 
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
png("03_figs/hist_freq.png", height = 4, width = 8, 
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

source('01_scripts/r_functions/box.cox.chord.R') # provient de Appendix 4 (Legendre & Borcard, 2018)

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
png("03_figs/transfo.png", height = 4, width = 8, 
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
spe.brute.rda <- rda(spe, grid_env)
brute.R2aj <- RsquareAdj(spe.brute.rda)$adj.r.squared

# Chorde
spe.norm.rda <- rda(spe.norm, grid_env)
chorde.R2aj <- RsquareAdj(spe.norm.rda)$adj.r.squared 

# Hellinger
spe.hel.rda <- rda(spe.hel, grid_env)
hel.R2aj <- RsquareAdj(spe.hel.rda)$adj.r.squared 

# Khi-carré
spe.chisq.rda <- rda(spe.chisq, grid_env)
chisq.R2aj <- RsquareAdj(spe.chisq.rda)$adj.r.squared 

# Total
spe.total.rda <- rda(spe.total, grid_env)
tot.R2aj <- RsquareAdj(spe.total.rda)$adj.r.squared 

# log-chorde
spe.ln.norm.rda <- rda(spe.ln.norm, grid_env)
ln.norm.R2aj <- RsquareAdj(spe.ln.norm.rda)$adj.r.squared

# double square root (Box-cox chord, bc.exp = 0.25)
spe.bc.rda <- rda(spe.bc, grid_env)
boxcox.R2aj <- RsquareAdj(spe.bc.rda)$adj.r.squared 

transfo.R2aj <- data.frame("Transformation" = c("Aucune", "Chorde", "Hellinger",
                                        "Khi-carré", "Profils d'abondances relatives",
                                        "log-Chorde", "Double racine carrée"),
                   "R2ajusté" = c(brute.R2aj, chorde.R2aj, hel.R2aj, chisq.R2aj,
                                   tot.R2aj, ln.norm.R2aj, boxcox.R2aj))


write.csv(transfo.R2aj,"02_outdata/transfo.R2aj.csv", row.names = FALSE)

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
png("03_figs/carte_liens3_logchord.png",
    height = 4, width = 3, # taille
    units = "in", # unités
    res = 300) # résolution
plot(grpWD2cst_constr_hclust, k=3, links=TRUE, main = "Trois groupements", las=1, xlab="x (m)",
     ylab="y (m)", cex=0.6, lwd=1, cex.axis = 0.5, cex.main = 0.75, cex.lab = 0.75)
dev.off()

# Now with k=4 clusters:
png("03_figs/carte_liens4_logchord.png",
    height = 4, width = 3, # taille
    units = "in", # unités
    res = 300) # résolution
plot(grpWD2cst_constr_hclust, k=4, links=TRUE, main = "Quatre groupements", las=1, xlab="x (m)",
     ylab="", yaxt="n", cex=0.6, lwd=1, cex.axis = 0.5, cex.main = 0.75, cex.lab = 0.75)
dev.off()

# Now with k=5 clusters:
png("03_figs/carte_liens5_logchord.png",
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
surfs <- read_csv("02_outdata/surfs.csv")
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
#write_csv(comm3_tbl, "02_outdata/comm3.csv")
#write_csv(comm4_tbl, "02_outdata/comm4.csv")
#write_csv(comm5_tbl, "02_outdata/comm5.csv")

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
ggsave(graph_comm3, file = '03_figs/graph_comm3.png', # nom du fichier avec extension
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
ggsave(graph_comm4, file = '03_figs/graph_comm4.png', # nom du fichier avec extension
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
ggsave(graph_comm5, file = '03_figs/graph_comm5.png', # nom du fichier avec extension
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

# # Importation d'une orthomosaïque
# require(stars)
# ortho <- read_stars('00_rawdata/2021-09-28-sbl-cloutier-z1-P4RTK-MS.tif', 
#                     proxy = TRUE,
#                     NA_value = 0) %>% # car les valeurs manquantes sont 0 dans ce raster
#   st_transform(2950) %>% # on transforme au code 2950
#   st_crop(zone) %>% # on extrait la zone du Lac Croche
#   st_as_stars() # on convertit au format stars 
#    
# # plot(ortho, rgb = 1:3) # pour une image couleur
# 
# st_crs(ortho)$epsg 
# 
# # Carte statique des communautés (3)
# tmap_mode("plot")
# 
# carte_3comm <- tm_shape(ortho, # on charge notre orthomosaïque
#                         bbox = grid_arbres) + # on définit l'étendue à partir des communautés
#   tm_rgba() + # pour image RGB avec transparence (a)
#   tm_shape(grid_arbres) + # on charge les points
#   tm_polygons(col = "comm3_nom",
#               alpha = 0.5,
#               title = "Communautés") +
#   tm_layout(legend.position = c('right', 'bottom'),
#             inner.margins = c(0.05, .05, .05, .05)) + # on augmente les marges
#   tm_compass(type = "8star", position = c("left", "top")) + # rose des vents
#   tm_scale_bar() + # l'échelle
#   tm_graticules(lwd = 0.3,
#                 alpha = 0.6)
# #carte_3comm # on visualise le résultat
# 
# tmap_save(tm = carte_3comm,
#           filename = '03_figs/carte_3comm.png', # on nomme le fichier # l'extension définit le format
#           dpi = 600, # résolution en pixels par pouce
#           width = 8) # largeur de l'image en pouces
# 
# # Carte statique des communautés (4)
# carte_4comm <- tm_shape(ortho, # on charge notre orthomosaïque
#                         bbox = grid_arbres) + # on définit l'étendue à partir des communautés
#   tm_rgba() + # pour image RGB avec transparence (a)
#   tm_shape(grid_arbres) + # on charge les points
#   tm_polygons(col = "comm4_nom",
#               alpha = 0.5,
#               title = "Communautés") +
#   tm_layout(legend.position = c('right', 'bottom'),
#             inner.margins = c(0.05, .05, .05, .05)) + # on augmente les marges
#   tm_compass(type = "8star", position = c("left", "top")) + # rose des vents
#   tm_scale_bar() + # l'échelle
#   tm_graticules(lwd = 0.3,
#                 alpha = 0.6)
# #carte_4comm # on visualise le résultat
# 
# tmap_save(tm = carte_4comm,
#           filename = '03_figs/carte_4comm.png', # on nomme le fichier # l'extension définit le format
#           dpi = 600, # résolution en pixels par pouce
#           width = 8) # largeur de l'image en pouces
# 
# # Carte statique des communautés (5)
# carte_5comm <- tm_shape(ortho, # on charge notre orthomosaïque
#                         bbox = grid_arbres) + # on définit l'étendue à partir des communautés
#   tm_rgba() + # pour image RGB avec transparence (a)
#   tm_shape(grid_arbres) + # on charge les points
#   tm_polygons(col = "comm5_nom",
#               alpha = 0.5,
#               title = "Communautés") +
#   tm_layout(legend.position = c('right', 'bottom'),
#             inner.margins = c(0.05, .05, .05, .05)) + # on augmente les marges
#   tm_compass(type = "8star", position = c("left", "top")) + # rose des vents
#   tm_scale_bar() + # l'échelle
#   tm_graticules(lwd = 0.3,
#                 alpha = 0.6)
# #carte_5comm # on visualise le résultat
# 
# tmap_save(tm = carte_5comm,
#           filename = '03_figs/carte_5comm.png', # on nomme le fichier # l'extension définit le format
#           dpi = 600, # résolution en pixels par pouce
#           width = 8) # largeur de l'image en pouces

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
svglite("03_figs/ACP_3comm.svg")
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
ggsave(acp_graph_sites_sp, file = '03_figs/ACP_3comm.png', # nom du fichier avec extension
       width = 8, # largeur en pouces
       height = 6, # hauteur en pouces
       dpi = 300) # résolution en pixels par pouce

#### Ordination canonique - sélection des variables significatives en RDA ----
source('01_scripts/r_functions/triplot.rda.R')
## Selection of explanatory variables in RDA (Practicals using the R language p.31)
dim(spe.transfo) # Response variables
dim(grid_env) # Explanatory variables

# Sélection des variables avec forward.sel
var.sel = forward.sel(spe.transfo, grid_env)

write.csv(var.sel,"02_outdata/var_sel.csv", row.names = FALSE) # nécessaire?

# Run rda on selected variables
(spe.rda.sel <- rda(spe.transfo ~ ., grid_env[, var.sel$order[1:5]]))

summary(spe.rda.sel)

grid_env.sel <- grid_env[, var.sel$order[1:5]] # On crée un sous-jeu de données pour les variables significatives

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

source('01_scripts/r_functions/triplot.rda.R')
# Triplots with function triplot.rda(), scalings 1 and 2, lc scores

# Figure : "RDA plot with triplot.rda"
png("03_figs/RDA_triplot.rda.lc.png",
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
svglite("03_figs/RDA_triplot.rda.lc.svg",
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
png("03_figs/RDA_triplot.rda.wa.png",
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

