#### EXTRACTION DES DONNÉES LiDAR ####
#### Importation des packages ----
library(tidyverse)
#library(raster)
library(stars)
library(sf)
library(starsExtra)
library(circular)
library(dplyr)

#### Importation des données ----
# Polygones des placettes forestières ----
grid <- read_sf("02_outdata/grid.geojson") %>% # on lit le fichier
  st_transform(2950) # on transforme au code 2950

# Zone d’échantillonnage ----
zone <- read_sf("00_rawdata/Z1.geojson") %>% # on lit le fichier
  st_geometry() %>% # on ne garde que la géométrie
  st_transform(2950) # on transforme au code 2950

# MNT ----
mnt <- read_stars('00_rawdata/mnt_2950.tif',  # on pointe vers le fichier ## package stars
                  proxy = TRUE) %>% # proxy = TRUE pour ne pas charger le fichier au complet dans votre mémoire vive
  st_crop(zone) %>% # on extrait la zone du Lac Croche
  st_as_stars() # on convertit au format stars

st_crs(mnt)$epsg # 2950

png("03_figs/mnt.png",
    height = 8, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
plot(mnt, col = terrain.colors(10), main = 'Modèle numérique de terrain', reset = FALSE) # on inspecte le résultat en couleur)
plot(st_geometry(grid), add = TRUE) # on ajoute les placettes
dev.off()

# Pentes ----
pentes <- read_stars('00_rawdata/pentes_2950.tif',  # on pointe vers le fichier
                     proxy = TRUE) %>% # proxy = TRUE pour ne pas charger le fichier au complet dans votre mémoire vive
  st_crop(zone) %>% # on extrait la zone du Lac Croche
  st_as_stars() # on convertit au format stars

st_crs(pentes)$epsg # 2950

png("03_figs/pentes.png",
    height = 8, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
plot(pentes, col = terrain.colors(10), main = 'Pente LiDAR (%)', reset = FALSE) # on inspecte le résultat en couleur
plot(st_geometry(grid), add = TRUE) # on ajoute les placettes
dev.off()

# Slope (calculé avec R) ----
slope_deg <- slope(mnt) # package starsExtra

plot(slope_deg, col = terrain.colors(10), reset = FALSE) # on inspecte le résultat en couleur
plot(st_geometry(grid), add = TRUE) # on ajoute les placettes

png("03_figs/slope_deg.png",
    height = 8, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
plot(slope_deg, col = terrain.colors(10), main = 'Pente (°)', reset = FALSE) # on inspecte le résultat en couleur
plot(st_geometry(grid), add = TRUE) # on ajoute les placettes
dev.off()

# même si slope peut donner des valeurs de 0-360, ici il reste dans le premier quadrant 0-90,
# donc nous pouvons rapporter ce 0-90 sur 100.

slope_pourc <- (100*slope_deg)/90

png("03_figs/slope_pourc.png",
    height = 8, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
plot(slope_pourc, col = terrain.colors(10), main = 'Pente in R (%)', reset = FALSE) # on inspecte le résultat en couleur
plot(st_geometry(grid), add = TRUE) # on ajoute les placettes
dev.off()

# Orientation de la pente - eastness / northness ----
# transformation en radians en même temps
northness <- cos(aspect(mnt) * pi / 180) 
eastness <- sin(aspect(mnt) * pi / 180)

png("03_figs/northness.png",
    height = 8, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
plot(northness, col = heat.colors(20), main = "northness", reset = FALSE) # on inspecte le résultat en couleur
plot(st_geometry(grid), add = TRUE) # on ajoute les placettes
dev.off()

png("03_figs/eastness.png",
    height = 8, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
plot(eastness, col = heat.colors(20), main = "eastness", reset = FALSE) # on inspecte le résultat en couleur
plot(st_geometry(grid), add = TRUE) # on ajoute les placettes
dev.off()

# Orientation de la pente - eastness exposure / northness exposure ----
slope_rad <- (slope_deg*pi)/180

#plot(slope_rad, col = terrain.colors(10), main = 'Pente in R (rad)', reset = FALSE) # on inspecte le résultat en couleur
#plot(st_geometry(grid), add = TRUE) # on ajoute les placettes

northness_exposure = sin(slope_rad)*northness
eastness_exposure = sin(slope_rad)*eastness

png("03_figs/northness_exposure.png",
    height = 8, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
plot(northness_exposure, col = heat.colors(20), main = "northness exposure", reset = FALSE) # on inspecte le résultat en couleur
plot(st_geometry(grid), add = TRUE) # on ajoute les placettes
dev.off()

png("03_figs/eastness_exposure.png",
    height = 8, width = 8, # taille
    units = "in", # unités
    res = 300) # résolution
plot(eastness_exposure, col = heat.colors(20), main = "eastness exposure", reset = FALSE) # on inspecte le résultat en couleur
plot(st_geometry(grid), add = TRUE) # on ajoute les placettes
dev.off()

# Indice d'humidité topographique ----
twi <- read_stars('00_rawdata/twi_2950.tif',  # on pointe vers le fichier
                     proxy = TRUE) %>% # proxy = TRUE pour ne pas charger le fichier au complet dans votre mémoire vive
  st_crop(zone) %>% # on extrait la zone du Lac Croche
  st_as_stars() # on convertit au format stars

st_crs(twi)$epsg # 2950

png("03_figs/twi.png",
    height = 8, width = 8, # taille
    units = "in", # unités
    res = 500) # résolution
plot(twi, col = terrain.colors(10), main = "Indice d'humidité topographique", reset = FALSE) # on inspecte le résultat en couleur
plot(st_geometry(grid), add = TRUE) # on ajoute les placettes
dev.off()

#### Extraire des données pour plusieurs polygones ####
source("01_scripts/r_functions/extr_rast_stat.R") # inclut le package circular
# MNT ----
mnt_moy <- vector() # On crée un vecteur vide pour stocker les résultats de la boucle
for (i in 1:nrow(grid)) { # se lit: pour chaque élément i allant de 1 jusqu'à n = nombre de placettes)
  mnt_moy[i] <- extr_rast_stat(grid[i, ], # l'élément i de mnt est le résultat de la fonction extra_rast_stat() la placette i
                               rast = mnt, # en utilisant le raster mnt pour l'extraction
                               stat = "mean") # et on calcule la valeur moyenne de hauteur pour chaque placette
}
mnt_moy

# Pentes moyennes (LiDAR) ----
pentes_moy <- vector() # On crée un vecteur vide pour stocker les résultats de la boucle
for (i in 1:nrow(grid)) { # se lit: pour chaque élément i allant de 1 jusqu'à n = nombre de placettes)
  pentes_moy[i] <- extr_rast_stat(grid[i, ], # l'élément i de pente_moy est le résultat de la fonction extra_rast_stat() la placette i
                                 rast = pentes, # en utilisant le raster pentes pour l'extraction
                                 stat = "mean") # et on calcule la valeur moyenne de pente pour chaque placette
}
pentes_moy

# Pentes moyennes (calculées avec R) ----
slope_moy <- vector() # On crée un vecteur vide pour stocker les résultats de la boucle
for (i in 1:nrow(grid)) { # se lit: pour chaque élément i allant de 1 jusqu'à n = nombre de placettes)
  slope_moy[i] <- extr_rast_stat(grid[i, ], # l'élément i de pente_moy est le résultat de la fonction extra_rast_stat() la placette i
                                 rast = slope_pourc, # en utilisant le raster pentes pour l'extraction
                                 stat = "mean") # et on calcule la valeur moyenne de pente pour chaque placette
}
slope_moy

# Indice d'humidité topographique ----
twi_moy <- vector() # On crée un vecteur vide pour stocker les résultats de la boucle
for (i in 1:nrow(grid)) { # se lit: pour chaque élément i allant de 1 jusqu'à n = nombre de placettes)
  twi_moy[i] <- extr_rast_stat(grid[i, ], # l'élément i de mnt est le résultat de la fonction extra_rast_stat() la placette i
                               rast = twi, # en utilisant le raster mnt pour l'extraction
                               stat = "mean") # et on calcule la valeur moyenne de hauteur pour chaque placette
}
twi_moy

# Orientation de la pente / Eastness ----
east_moy <- vector() # On crée un vecteur vide pour stocker les résultats de la boucle
for (i in 1:nrow(grid)) { # se lit: pour chaque élément i allant de 1 jusqu'à n = nombre de stations)
  east_moy[i] <- extr_rast_stat(grid[i, ], # l'élément i de ori_pentes est le résultat de la fonction extra_rast_stat() la station i
                                rast = eastness, # en utilisant le raster ori_pentes pour l'extraction
                                stat = "mean") # et on calcule la valeur moyenne d'orientation pour chaque parcelle
}
east_moy

# Orientation de la pente / Northness ----
north_moy <- vector() # On crée un vecteur vide pour stocker les résultats de la boucle
for (i in 1:nrow(grid)) { # se lit: pour chaque élément i allant de 1 jusqu'à n = nombre de stations)
  north_moy[i] <- extr_rast_stat(grid[i, ], # l'élément i de ori_pentes est le résultat de la fonction extra_rast_stat() la station i
                                 rast = northness, # en utilisant le raster ori_pentes pour l'extraction
                                 stat = "mean") # et on calcule la valeur moyenne d'orientation pour chaque parcelle
}
north_moy

# Orientation de la pente / Eastness exposure ----
east_exp_moy <- vector() # On crée un vecteur vide pour stocker les résultats de la boucle
for (i in 1:nrow(grid)) { # se lit: pour chaque élément i allant de 1 jusqu'à n = nombre de stations)
  east_exp_moy[i] <- extr_rast_stat(grid[i, ], # l'élément i de ori_pentes est le résultat de la fonction extra_rast_stat() la station i
                                    rast = eastness_exposure, # en utilisant le raster ori_pentes pour l'extraction
                                    stat = "mean") # et on calcule la valeur moyenne d'orientation pour chaque parcelle
}
east_exp_moy

# Orientation de la pente / Northness exposure ----
north_exp_moy <- vector() # On crée un vecteur vide pour stocker les résultats de la boucle
for (i in 1:nrow(grid)) { # se lit: pour chaque élément i allant de 1 jusqu'à n = nombre de stations)
  north_exp_moy[i] <- extr_rast_stat(grid[i, ], # l'élément i de ori_pentes est le résultat de la fonction extra_rast_stat() la station i
                                     rast = northness_exposure, # en utilisant le raster ori_pentes pour l'extraction
                                     stat = "mean") # et on calcule la valeur moyenne d'orientation pour chaque parcelle
}
north_exp_moy

#### Compilation des données dans une seule table ----
grid_topo <- bind_cols(grid, # on part de la couche grid
                       pentes = pentes_moy,
                       slope = slope_moy,
                       elev_mnt = mnt_moy,
                       twi = twi_moy,
                       eastness = east_moy,
                       northness = north_moy,
                       eastness_exp = east_exp_moy,
                       northness_exp = north_exp_moy)
grid_topo # on inspecte le tout

grid_topo$grid_id <- as.character(grid_topo$grid_id) # On transforme les grid_id en character

#### Visualisation des variables topographiques ----
plot(grid_topo["pentes"])
plot(grid_topo["slope"])


plot(grid_topo["elev_mnt"])
plot(grid_topo["twi"])
plot(grid_topo["eastness"])
plot(grid_topo["northness"])
plot(grid_topo["eastness_exp"])
plot(grid_topo["northness_exp"])

# OU : 
par(mfrow = c(2,2))
plot(grid_topo[4:8])
plot(grid_topo[9:12])

par(mfrow = c(1,1)) # on remet à 1:1
# Création des cartes choroplètes
library(tmap)

carte_slope <- tm_shape(grid_topo) +
  tm_fill(col = "slope",
          title = "Pentes avec R (%)") +
  tm_borders() 

carte_eastness <- tm_shape(grid_topo) +
  tm_fill(col = "eastness",
          title = "Eastness") +
  tm_borders() 

carte_northness <- tm_shape(grid_topo) +
  tm_fill(col = "northness",
          title = "Northness") +
  tm_borders() 

carte_eastness_exp <- tm_shape(grid_topo) +
  tm_fill(col = "eastness_exp",
          title = "Eastness exposure") +
  tm_borders() 

carte_northness_exp <- tm_shape(grid_topo) +
  tm_fill(col = "northness_exp",
          title = "Northness exposure") +
  tm_borders() 

carte_elev_mnt <- tm_shape(grid_topo) +
  tm_fill(col = "elev_mnt",
          title = "Élévation (m)") +
  tm_borders()

carte_twi <- tm_shape(grid_topo) +
  tm_fill(col = "twi",
          title = "Indice d'humidité
  topographique") +
  tm_borders()

(carte_mnt_twi <- tmap_arrange(carte_elev_mnt, carte_twi, ncol = 2))
(carte_east_north <- tmap_arrange(carte_eastness, carte_northness, ncol = 2))
(carte_east_north_exp <- tmap_arrange(carte_eastness_exp, carte_northness_exp, ncol = 2))

tmap_save(tm = carte_slope, # objet tmap qu'on veut enregistrer
          filename = "03_figs/carte_slope.png", # nom et chemin du fichier qu'on enregistre, incluant l'extension
          height = 4, # hauteur de la carte
          width = 4, # largeur de la carte
          units = "in", # unités utilisées pour la hauteur et la largeur
          dpi = 300) # par défault dpi = 300, on augmente pour une meilleure résolution

tmap_save(tm = carte_mnt_twi, # objet tmap qu'on veut enregistrer
          filename = "03_figs/carte_mnt_twi.png", # nom et chemin du fichier qu'on enregistre, incluant l'extension
          height = 4, # hauteur de la carte
          width = 8, # largeur de la carte
          units = "in", # unités utilisées pour la hauteur et la largeur
          dpi = 300) # par défault dpi = 300, on augmente pour une meilleure résolution

tmap_save(tm = carte_east_north, # objet tmap qu'on veut enregistrer
          filename = "03_figs/carte_east_north.png", # nom et chemin du fichier qu'on enregistre, incluant l'extension
          height = 4, # hauteur de la carte
          width = 8, # largeur de la carte
          units = "in", # unités utilisées pour la hauteur et la largeur
          dpi = 300) # par défault dpi = 300, on augmente pour une meilleure résolution

tmap_save(tm = carte_east_north_exp, # objet tmap qu'on veut enregistrer
          filename = "03_figs/carte_east_north_exp.png", # nom et chemin du fichier qu'on enregistre, incluant l'extension
          height = 4, # hauteur de la carte
          width = 8, # largeur de la carte
          units = "in", # unités utilisées pour la hauteur et la largeur
          dpi = 300) # par défault dpi = 300, on augmente pour une meilleure résolution

#### Exportation de grid_topo en GEOJSON ----
st_write(grid_topo, dsn = "02_outdata/grid_topo.geojson", driver = "GeoJSON") ## package sf