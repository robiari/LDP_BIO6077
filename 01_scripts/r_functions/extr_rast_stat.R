## FONCTION QUI PERMET D'EXTRAIRE DES VALUES MATRICIELLES DE PLUSIEURS POLYGONES EN MÊME TEMPS
library(circular)
extr_rast_stat <- function(vect, # 1er argument = la couche vectorielle sf
                           rast, # 2e argument = la couche matricielle
                           stat = c("mean", # la statistique à calculer
                                    "mean.circular", # pour orientation pente
                                    "min",
                                    "max",
                                    "var",
                                    "sd")) { # début fonction
  extr <- st_crop(rast, vect) %>% # extraction matricielle à l'aide d'une couche vectorielle
    st_as_stars() %>% # on lit les données matricielles si stars_proxy
    .[[1]] # on extrait les valeurs de la bande unique
  
  # Calcul différent pour chaque statistique
  
  # Moyenne
  if (stat == 'mean') {
    res = mean(extr, na.rm = TRUE)
    
    # Moyenne circulaire
  } else if (stat == 'mean.circular') {
    extr_circ <- circular(as.numeric(extr),
                          units = "degrees",
                          template = "geographics",
                          modulo = "2pi")
    res = mean(extr_circ, na.rm = TRUE)
    
    # Minimum  
  } else if (stat == 'min') {
    res = min(extr, na.rm = TRUE)
    
    # Maximum
  } else if (stat == 'max') {
    res = max(extr, na.rm = TRUE)
    
    # Variance
  } else if (stat == 'var') {
    res = var(as.numeric(extr), na.rm = TRUE)
    
  } else if (stat == 'sd') {
    res = sd(as.numeric(extr), na.rm = TRUE)
    
    # Écart-type
  } else res = NA # sinon, valeur manquante
  
  # On retourne le résultat
  return(res)
} # fin de la fonction
  