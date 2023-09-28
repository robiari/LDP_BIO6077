#### PRÉTRAITEMENT ARBRES ET PLACETTES FORESTIÈRES ####
# 1) script arbres/grid; 2) script traits; 3) script environnement; 4) script Travaildesession (analyses)
# and export clean data between each one.

#### Importation des packages ----
library(tidyverse)
library(sf)
#library(esri2sf)


#### Importation et nettoyage des données d'arbres (type sf) ----
Z1_polygones <- read_sf("00_rawdata/Z1_polygones.geojson") %>% # on lit le fichier
  st_transform(2950) %>% # on transforme au code 2950
  select(OBJECTID, Label, surf_m2) %>% # on sélectionne les variables pertinentes
  filter(Label != 'Mort' & Label != 'Acer' & Label != 'Populus' & Label != 'Feuillus' & Label != 'Conifere' & Label != 'Betula' & Label != 'PRPE') # on enlève les classes non représentatives

# Regroupement de tous les Picea dans une même classe
Z1_polygones$Label[Z1_polygones$Label %in% c("PIGL","PIRU","PIMA")] <- "Picea" 

#### Exploration de la table de données d'arbres (type sf) ----
head(Z1_polygones)                 # Aperçu des 6 premières lignes
nrow(Z1_polygones)                 # Nombre d'individus
ncol(Z1_polygones)                 # Nombre de variables
dim(Z1_polygones)                  # Dimensions de la table de données (rows, columns)
colnames(Z1_polygones)             # Nom des variables

# Aperçu du nombre d'individus par espèces dans toute la SBL
arbres_SBL <- st_drop_geometry(Z1_polygones) %>%
  group_by(Label) %>%
  count(.)
arbres_SBL

#### Création des placettes forestières de 20 m x 20 m ----
# Création de la grille (grid) représentant les placettes
grid_Z1 = st_make_grid(
  Z1_polygones,
  cellsize = 20,
  crs = if (missing(Z1_polygones)) NA_crs_ else st_crs(Z1_polygones),
  what = "polygons",
  square = TRUE)

#### Sélection des placettes comportant des annotations ----
grid_contain <- st_contains(grid_Z1, Z1_polygones)
sel_logical = lengths(grid_contain) > 0 # returns of logical vector for those cells we keep
grid_sel <- grid_Z1[sel_logical] # we filter only those cells

#### Sélection des placettes comportant SUFFISAMMENT d'annotations ----
# Ajout d'une variable grid_id
grid_sel_id <- grid_sel %>% 
  st_as_sf() %>% 
  mutate(grid_id = 1:nrow(.))

# Intersection des deux couches (grille et annotations)
polygones_grid <- st_intersection(Z1_polygones, grid_sel_id)

# Détermination de la proportion de l'aire « annotée » dans chaque placette
# ? ATTENTION: Répétition des calculs de superficie avec la section : Superficies de la canopée
species_area_grid <- polygones_grid %>% 
  select(grid_id) %>% 
  mutate(area = st_area(.)) %>% 
  group_by(grid_id) %>% 
  summarise_at(vars(area), sum) %>% 
  mutate(area_prop = as.numeric(area) / 400)

# Extraction des placettes avec plus de 40% d'annotations (40% de surface annotée)
thresh <- 0.40
tmp <- species_area_grid %>% 
  dplyr::filter(area_prop >= thresh) %>% 
  select(grid_id) %>% 
  st_drop_geometry()

# Sélectionner les placettes qui en ont assez
grid <- right_join(grid_sel_id, tmp)

grid_id <- grid$grid_id # on extrait les grid_id
rownames(grid) <- grid_id # on utilise les grid_id comme noms de lignes

# Enlever les grid 81 et 114 qui sont exclues du reste, car dans le zone des chalets
grid <- grid[!(grid$grid_id == 81 | grid$grid_id == 114),]

grid_id <- grid$grid_id # on extrait les grid_id sans 81 et 114
rownames(grid) <- grid_id # on utilise les grid_id comme noms de lignes

#### Matrice des coordonnées des placettes ----
grid_coord = st_centroid(grid) %>% # On extrait les centroïdes des placettes
  st_coordinates() # On extrait uniquement la matrice des coordonnées (pas d'objet sf)

grid_coord <-  data.frame(X = grid_coord[,1],
                            Y = grid_coord[,2],
                            grid_id = grid_id)
rownames(grid_coord) <- grid_id # on utilise les grid_id comme noms de lignes

grid <- grid %>%
  left_join(grid_coord)

#### Exportation de grid en GEOJSON et de grid_coord en csv ----
# à décommenter lorsqu'on roule le code pour la première fois
st_write(grid, dsn = "02_outdata/grid.geojson", driver = "GeoJSON") ## package sf
#st_write(grid, dsn = "02_outdata/grid.shp") ## package sf
write.csv(grid_coord,"02_outdata/grid_coord.csv", row.names = FALSE)

#### Sélection des arbres se trouvant dans les placettes sélectionnées ----
# Filtrage les données d'arbres pour conserver uniquement les individus dans les placettes sélectionnées ----
arbres_cent_sel <- Z1_polygones %>% 
  st_centroid() %>% # Transformation de la couche de polygones en couche de points en fonction des coordonnées du centroïde
  st_filter(grid, .pred = st_contains) %>% # Sélection des arbres qui ont leur centroïde dans une placette sélectionnée
  st_join(grid, join = st_within) %>% # Jointure spatiale avec la grille pour conserver la variable grid_id
  select(OBJECTID, grid_id) %>% # On crée une table avec des FID des arbres sélectionnés
  st_drop_geometry()

# Jointure de table avec Z1_polygones pour conserver uniquement les arbres dont le centroïde se trouve à l'intérieur d'une placette
arbres_poly_sel <- right_join(Z1_polygones, arbres_cent_sel)
head(arbres_poly_sel)
#st_write(arbres_poly_sel, dsn = "02_outdata/arbres_poly_sel.geojson", driver = "GeoJSON") ## package sf

# Aperçu du nombre d'individus par espèces dans les placettes sélectionnées
arbres_SBL_sel <- st_drop_geometry(arbres_poly_sel) %>%
  group_by(Label) %>%
  count(.)
arbres_SBL_sel

write.csv(arbres_SBL_sel,"02_outdata/arbres_SBL_sel.csv", row.names = FALSE)

#### Graphique des placettes forestières sélectionnées avec polygones ----
# plot(st_geometry(grid), border = 'red')
# plot(arbres_poly_sel['Label'], add = TRUE) # Les couleurs sont en fonction des Label.
# plot(st_geometry(grid), border = 'red', add= TRUE)

#### Conversion au format large ----
arbres_large <- st_drop_geometry(arbres_poly_sel) %>%
  group_by(grid_id, Label) %>%
  count(.) %>% # décompte du nombre d'individus de chaque espèce dans chacune des placettes
  select(grid_id, Label, n) %>% # on ne conserve que ces variables pertinentes
  pivot_wider(id_cols = grid_id, # convertir en format large ## package tidyverse
              names_from = Label, # utiliser Label pour former des colonnes
              values_from = n, # remplir les cellules (valeurs) avec le nombre d'individus par espèce
              values_fill = 0) # remplir les valeurs manquantes avec des zéros

arbres_large$grid_id <- as.character(arbres_large$grid_id) # On transforme les grid_id en character
head(arbres_large)

spe <- arbres_large %>% # on part des données en format large
  ungroup() %>% # on enlève les groupes formés précédemment
  select(-grid_id) %>% # on enlève la colonne grid_id
  as.matrix()

rownames(spe) <- grid_id # on utilise les grid_id comme noms de lignes
head(spe) # format "large" corrigé et transformé en matrix avec nom

write.csv(spe,"02_outdata/spe.csv", row.names = TRUE)

#### Superficies de la canopée ----
# ? ATTENTION: Répétition des calculs de superficie avec la section : Sélection des placettes comportant SUFFISAMMENT d'annotations

# Superficie de la canopée par espèce par placette
surfs <- arbres_poly_sel %>% # on part des données initiales d'arbres
  st_drop_geometry() %>%
  group_by(grid_id, Label) %>% # on groupe par station et par espèce d'arbre
  summarise(sum_surf_m2 = sum(surf_m2)) %>% # on calcule la somme des superficies des canopées par station
  mutate(sum_surf_m2_ha = (sum_surf_m2 *25)) # on reporte le tout en m2 ha-1, par espèce
surfs # on examine le résultat

# Superficie de canopée relative par espèce par placette
surfs <- surfs %>% # on part du tableau précédent
  ungroup() %>%
  group_by(grid_id) %>% # on groupe par placette
  mutate(surf_rel = sum_surf_m2_ha / sum(sum_surf_m2_ha)) # on calcule la superficie de canopée relative de chaque espèce par placette
surfs$grid_id <- as.character(surfs$grid_id)

write.csv(surfs,"02_outdata/surfs.csv", row.names = FALSE)

#On peut aussi calculer la superficie de la canopée présente, toutes espèces confondues.
#Une forêt mature peut avoir jusqu'à 100% de sa surface recouverte par une canopée.

surf_grid <- surfs %>% 
  ungroup() %>%  # on enlève les groupes précédemment définis
  group_by(grid_id) %>% 
  summarise(surf_tot = sum(sum_surf_m2_ha),
            surf_tot_pourc = (sum(sum_surf_m2_ha)/10000)*100)
surf_grid