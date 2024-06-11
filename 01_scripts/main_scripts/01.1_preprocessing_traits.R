#### PRÉTRAITEMENT TRAITS FOLIAIRES ####

# project : 2022-Gravel-MSc-UdeM, 2018-BeauchampRioux-MSc-UdeM, 2019-Blanchard-MSc-UdeM, 2017-Dessain-MSc, CABO-General, 2019-CABO-General
# site_id : SBLUdeM, SBLUdeM-a, SBLUdeM-f, sblac

# Logique: 1) tri en fct du projet; 2) tri en fct du site; 3) tri en fonction des dates; 4) tri en fonction des espèces.

#### Chargement des packages ----
library(dplyr)
library(tidyverse)

#### Importation des fichiers de données de FULCRUM ----
# Tous les fichiers ont déjà été filtrés par site et projet à l'aide des filtres de Fulcrum avant le téléchargement des données.

# fichiers plants et bulks
plants_raw <- read_csv("00_rawdata/traits/plants.csv") # use plants for species ID
bulksamples_raw <- read_csv("00_rawdata/traits/bulk_leaf_samples.csv")

# fichiers de traits
Cfractions_raw_data <- read_csv("00_rawdata/traits/carbon_fractions_bags.csv")
CN_raw_data <- read_csv("00_rawdata/traits/cn_leaf_concentrations.csv")
SLA_raw_data <- read_csv("00_rawdata/traits/leaf_area_and_water_samples.csv")
pigments_raw_data <- read_csv("00_rawdata/traits/pigments_extracts.csv")

# fichiers de espèces d'arbres présents à la SBL
especes <- read.table('00_rawdata/especes.txt', header = TRUE, sep = ";")

#### Nettoyage des données de FULCRUM ----
# Plants file ----
plants <- plants_raw %>%
  select(plant_id,
         scientific_name,
         '_project',
         site_id, # déjà filtré par site et projet dans Fulcrum
         plant_remarks) %>% # comments
  arrange(plant_id) %>% #ascending order
  relocate(plant_id, .before = c(scientific_name, '_project')) %>% #to change column order
  rename(project = '_project')

plants #verify result

# Il n'y a pas de commentaires importants, supprimer la colonne de commentaires
plants <- plants %>%
  select(-plant_remarks)

plants #verify result

# Garder uniquement les plant_id qui sont reliés à des espèces d'arbres de la SBL
plants_sbl <- inner_join(plants, especes, by = "scientific_name")
plants_sbl #check result

# Regroupement de tous les Picea dans une même classe
plants_sbl$codesp[plants_sbl$codesp %in% c("PIGL","PIRU","PIMA")] <- "Picea" 

# Bulk samples file ----
bulksamples <- bulksamples_raw %>%
  select(plant_id,
         sample_id,
         date_sampled,
         sample_remarks) %>%
  arrange(plant_id) %>% #ascending order
  relocate(plant_id, .before = sample_id) #to change column order

bulksamples #verify result

# Il n'y a pas de commentaires importants, supprimer la colonne de commentaires
bulksamples <- bulksamples %>%
  select(-sample_remarks)

bulksamples #verify result

# Sélection des échantillons ayant été récoltés entre mi-juin et mi-août seulement ----
# Aucune donnée en 2020 et 2021
bulk_dates <- bulksamples[bulksamples$date_sampled >= "2017-06-15" & bulksamples$date_sampled <= "2017-08-15"
                          | bulksamples$date_sampled >= "2018-06-15" & bulksamples$date_sampled <= "2018-08-15"
                          | bulksamples$date_sampled >= "2019-06-15" & bulksamples$date_sampled <= "2019-08-15"
                          | bulksamples$date_sampled >= "2022-06-15" & bulksamples$date_sampled <= "2022-08-15", ]

# Carbon and Nitrogen leaf concentrations file ----
CN_data <- CN_raw_data %>%
  select(sample_id,
         c_perc, # pourcentage de carbone
         n_perc, # pourcentage d'azote
         analysis_remarks) %>%
  mutate('C_N_ratio' = c_perc / n_perc) # calcul du ratio C:N

CN_data #verify result

# Il y a un duplicat. Il faudra donc en faire la moyenne.
# Pour l'instant, il est possible de supprimer la colonne de commentaires.
CN_data <- CN_data %>%
  select(-analysis_remarks)

CN_data #verify result

# Mean CN_data for 138769323 (same sample, measured twice)
mean_138769323 <- CN_data %>%
  filter(sample_id %in% "138769323") %>% #Select both samples
  select(-1) %>% #remove column that are numerical
  colMeans() %>% #to have the mean of every column (every trait)
  as.data.frame.list() %>% #change to dataframe
  mutate(sample_id = 138769323)

mean_138769323 #check result

CN_data <- CN_data %>%
  filter(sample_id != "138769323") %>% #remove both samples
  rbind(mean_138769323) #add row with mean of both samples

CN_data #check result

# Leaf area and water samples file ----
SLA_data <- SLA_raw_data %>%
  select(sample_id,
         specific_leaf_area_m2_kg, # SLA (inverse de LMA)
         leaf_mass_per_area_g_m2, # LMA
         leaf_dry_matter_content_mg_g, # LDMC (inverse de LWC)
         leaf_water_content_mg_g, # LWC
         sample_remarks) # comments

SLA_data #verify result

# Après avoir vérifié les commentaires, enlever les données erronées
SLA_data <- SLA_data %>%
  subset(is.na(sample_remarks)) %>% # enlever les lignes avec des commentaires (contrôle de qualité)
  select(-sample_remarks) # retirer les colonnes de commentaires

SLA_data #verify result   

# Carbon fractions file ----
Cfractions_data <- Cfractions_raw_data %>%
  select(bottle_id, # sample_id
         sample_type,
         hemicellulose_perc,
         cellulose_perc,
         lignin_perc,
         quality_flag_bag,
         sample_remarks, # comments
         quality_comments_bag) %>% # comments
  filter(sample_type == "sample", quality_flag_bag == "good") %>% #remove blanks and bad samples
  rename(sample_id = bottle_id) %>% #change the name of the variable bottle_id for sample_id
  select(-sample_type,-quality_flag_bag) #remove column sample_type and quality_flag_bag

Cfractions_data #verify result  

# Après avoir vérifié les commentaires, enlever les données erronées
Cfractions_data <- Cfractions_data %>%
  subset(is.na(sample_remarks)) %>% # enlever les lignes avec des commentaires (contrôle de qualité)
  select(-sample_remarks, -quality_comments_bag) # retirer les colonnes de commentaires  

Cfractions_data #verify result  

# Pigments extracts file ----
pigments_data <- pigments_raw_data %>%
  select(vial_id, # sample_id
         'sample_code',
         chla_mg_g_disk_mass,
         chlb_mg_g_disk_mass,
         carot_mg_g_disk_mass,
         chl_a_chl_b_ratio,
         quality_flag_extract,
         quality_comments_extract) %>% # comments
  rename(sample_id = vial_id) %>% #change the name of the variable vial_id for sample_id
  filter(sample_code != c("S1", "S2"), quality_flag_extract == "good") %>% #remove S1 and S2 (not samples) and keep only good extracts
  select(-sample_code, -quality_flag_extract) %>% #remove column sample_code and quality_flag_extract
  arrange(sample_id) #ascending order 

pigments_data #verify result 

# Il n'y a pas de commentaires importants, supprimer la colonne de commentaires
pigments_data <- pigments_data %>%
  select(-quality_comments_extract)

pigments_data #verify result

#### Vérification du choix de clés primaires ----
## When the result is ? 0 lines ? : means that each value comes once, good choice of primary key
#primary key: plants
plants %>% 
  count(plant_id) %>% 
  filter(n > 1)

#primary key: bulksamples
bulksamples %>% 
  count(sample_id) %>% 
  filter(n > 1)

bulksamples %>% 
  count(plant_id) %>% 
  filter(n > 1)

#primary key: CN
CN_data %>% 
  count(sample_id) %>% 
  filter(n > 1)

#primary key: SLA
SLA_data %>% 
  count(sample_id) %>% 
  filter(n > 1)

#primary key: Cfractions
Cfractions_data %>% 
  count(sample_id) %>% 
  filter(n > 1)

#primary key: pigments
pigments_data %>% 
  count(sample_id) %>% 
  filter(n > 1)

# Enlever les duplicats et les colonnes vides du fichier pigments_data
pigments_data <- filter(pigments_data, rowSums(is.na(pigments_data)) != ncol(pigments_data)) %>% # enlever les colonnes remplies de NA
  filter(sample_id != "11914908" & sample_id != "13991813" & sample_id != "28686796") # enlever les duplicats

pigments_data #verify result

# Revérifier les doublons du fichier pigments_data
pigments_data %>% 
  count(sample_id) %>% 
  filter(n > 1)

#### Vérification de la structure des fichiers ----
str(plants_sbl) #Here plant_id is a character
str(bulksamples) #Here plant_id and sample_id are characters
str(CN_data) #Here sample_id is a character
str(SLA_data) #Here sample_id is a number !
str(Cfractions_data) #Here sample_id is a character
str(pigments_data) #Here sample_id is a number !

# Changer les sample_id qui sont en format numérique
SLA_data$sample_id <- as.character(SLA_data$sample_id)
pigments_data$sample_id <- as.character(pigments_data$sample_id)

str(SLA_data) #be sure it worked
str(pigments_data) #be sure it worked

#### Jointures de tables avec Bulksamples ----
CN <- inner_join(bulk_dates, CN_data, by = "sample_id")
CN #check result

SLA <- inner_join(bulk_dates, SLA_data, by = "sample_id")
SLA #check result

Cfractions <- inner_join(bulk_dates, Cfractions_data, by = "sample_id")
Cfractions #check result

pigments <- inner_join(bulk_dates, pigments_data, by = "sample_id")
pigments #check result

# Moyenne des 2 samples du même arbre ----
# Mean SLA for 139284436 & 139293925 (same tree, Populus)
SLA_sel <- SLA %>%
  filter(sample_id == 139284436)
SLA_sel

mean_populus <- SLA %>%
  filter(sample_id %in% c("139284436", "139293925")) %>% #Select both samples
  select(-(1:3)) %>% #remove column that are numerical
  colMeans() %>% #to have the mean of every column (every trait)
  as.data.frame.list() %>% #change to dataframe
  mutate(plant_id = 139284299, sample_id = 139284436, date_sampled = "2022-07-21")

mean_populus #check result

SLA <- SLA %>%
  filter(sample_id != "139284436" & sample_id != "139293925") %>% #remove both samples
  rbind(mean_populus) #add row with mean of both samples

#duplicated(SLA$sample_id) #Check if there is duplicates in column plant_id

# Mean Cfractions for 139284436 & 139293925 (same tree, Populus)
mean_populus2 <- Cfractions %>%
  filter(sample_id %in% c("139284436", "139293925")) %>% #Select both samples
  select(-(1:3)) %>% #remove column that are numerical
  colMeans() %>% #to have the mean of every column (every trait)
  as.data.frame.list() %>% #change to dataframe
  mutate(plant_id = 139284299, sample_id = 139284436, date_sampled = "2022-07-21")

mean_populus2 #check result

Cfractions <- Cfractions %>%
  filter(sample_id != "139284436" & sample_id != "139293925") %>% #remove both samples
  rbind(mean_populus2) #add row with mean of both samples

#duplicated(Cfractions$sample_id) #Check if there is duplicates in column plant_id

# Mean pigments for 139284436 & 139293925 (same tree, Populus)
mean_populus3 <- pigments %>%
  filter(sample_id %in% c("139284436", "139293925")) %>% #Select both samples
  select(-(1:3)) %>% #remove column that are numerical
  colMeans() %>% #to have the mean of every column (every trait)
  as.data.frame.list() %>% #change to dataframe
  mutate(plant_id = 139284299, sample_id = 139284436, date_sampled = "2022-07-21")

mean_populus3 #check result

pigments <- pigments %>%
  filter(sample_id != "139284436" & sample_id != "139293925") %>% #remove both samples
  rbind(mean_populus3) #add row with mean of both samples

#duplicated(pigments$sample_id) #Check if there is duplicates in column plant_id

# Vérification des doublons ----
# L'idée ici est de vérifier s'il y a parfois plus d'un échantillon pour un même arbre.
# Dans un tel cas, on ferait une moyenne des valeurs, comme il l'a été fait plus haut

# plants
plants_sbl %>% 
  count(plant_id) %>% 
  filter(n > 1)

# bulksamples
bulksamples %>% 
  count(plant_id) %>% 
  filter(n > 1)

# CN
CN %>% 
  count(plant_id) %>% 
  filter(n > 1)

# SLA
SLA %>% 
  count(plant_id) %>% 
  filter(n > 1)

# Cfractions 
Cfractions %>% 
  count(plant_id) %>% 
  filter(n > 1)

# pigments
pigments %>% 
  count(plant_id) %>% 
  filter(n > 1)

# Tout est beau, nous pouvons donc procéder à une jointure avec le fichier Plants pour avoir l'ID

#### Jointures de tables avec Plants ----
# Regroupement de tous les Picea dans une même classe
CN_clean <- inner_join(plants_sbl, CN, by = "plant_id")
# Regroupement de tous les Picea dans une même classe
CN_clean$codesp[CN_clean$codesp %in% c("PIGL","PIRU","PIMA")] <- "Picea"
CN_clean #check result

SLA_clean <- inner_join(plants_sbl, SLA, by = "plant_id")
# Regroupement de tous les Picea dans une même classe
SLA_clean$codesp[SLA_clean$codesp %in% c("PIGL","PIRU","PIMA")] <- "Picea"
SLA_clean #check result

Cfractions_clean <- inner_join(plants_sbl, Cfractions, by = "plant_id")
# Regroupement de tous les Picea dans une même classe
Cfractions_clean$codesp[Cfractions_clean$codesp %in% c("PIGL","PIRU","PIMA")] <- "Picea"
Cfractions_clean #check result

pigments_clean <- inner_join(plants_sbl, pigments, by = "plant_id")
# Regroupement de tous les Picea dans une même classe
pigments_clean$codesp[pigments_clean$codesp %in% c("PIGL","PIRU","PIMA")] <- "Picea"
pigments_clean #check result

# Calcul du nombre d'échantillons par espèce et par trait ----
n_CN <- CN_clean %>%
  group_by(codesp) %>%
  count(.)
n_CN

n_SLA <- SLA_clean %>%
  group_by(codesp) %>%
  count(.)
n_SLA

n_Cfractions <- Cfractions_clean %>%
  group_by(codesp) %>%
  count(.)
n_Cfractions

n_pigments <- pigments_clean %>%
  group_by(codesp) %>%
  count(.)
n_pigments

# Jointures de tables pour obtenir une seule grande table de traits ----
traits <- full_join(CN_clean, SLA_clean, by = c("plant_id", "sample_id", "scientific_name", "codesp", "project", "site_id", "date_sampled"))
traits <- full_join(traits, Cfractions_clean, by = c("plant_id", "sample_id", "scientific_name", "codesp", "project", "site_id", "date_sampled"))
traits <- full_join(traits, pigments_clean, by = c("plant_id", "sample_id", "scientific_name", "codesp", "project", "site_id", "date_sampled"))

# On enlève les variables qui ne sont pas pertinentes pour l'analyse
traits <- traits %>%
  select(-c_perc, -n_perc, # C et N ne sont pas pertinents lorsqu'ils sont considérés individuellement.
         -specific_leaf_area_m2_kg, -leaf_dry_matter_content_mg_g) %>% # SLA et LMA, LDMC et LWC sont des inverses.
  mutate('Carot:Chl_tot' = carot_mg_g_disk_mass / (chla_mg_g_disk_mass+chlb_mg_g_disk_mass)) %>% # Carotenoids : tot Chl est préférable que les pigments individuellement
  select(-chla_mg_g_disk_mass, -chlb_mg_g_disk_mass, -carot_mg_g_disk_mass)
head(traits)

n_traits <- traits %>%
  group_by(codesp) %>%
  count(.)
n_traits

str(traits)

# Création d'un fichier de données plus élagué (seulement les espèces et les traits) ----
traits_clean <- traits %>%
  select(-plant_id, -sample_id, -scientific_name, -project, -site_id, -date_sampled)

traits_clean

#### Moyenne des traits par espèce ----
library(doBy)
traits_mean <- summaryBy(traits_clean[2:9] ~ codesp, traits_clean, FUN = mean, na.rm = TRUE)
head(traits_mean)

#### Exportation des données (fichiers) ----
write.csv(CN_clean,"02_outdata/traits/CN.csv", row.names = FALSE)
write.csv(SLA_clean,"02_outdata/traits/SLA.csv", row.names = FALSE)
write.csv(Cfractions_clean,"02_outdata/traits/Cfractions.csv", row.names = FALSE)
write.csv(pigments_clean,"02_outdata/traits/pigments.csv", row.names = FALSE)
write.csv(traits,"02_outdata/traits/traits.csv", row.names = FALSE)
write.csv(traits_clean,"02_outdata/traits/traits_clean.csv", row.names = FALSE)
write.csv(traits_mean,"02_outdata/traits/traits_mean.csv", row.names = FALSE)

traits_short <- traits_mean %>%
  rename('C:N' = 'C_N_ratio.mean',
         LMA = 'leaf_mass_per_area_g_m2.mean',
         LWC = 'leaf_water_content_mg_g.mean',
         Hemicellulose = 'hemicellulose_perc.mean',
         Cellulose = 'cellulose_perc.mean',
         Lignine = 'lignin_perc.mean',
         'Carot:Chltot' = 'Carot:Chl_tot.mean',
         'Chla:Chlb' = 'chl_a_chl_b_ratio.mean')
head(traits_short)

write.csv(traits_short,"02_outdata/traits/traits_short.csv", row.names = FALSE)

#### Statistiques sur les traits ----
x <- traits_short[2:9]

traits_stat <- data.frame('Moyenne' = apply(x, 2, mean, na.rm = TRUE),
                          'Écart-type' = apply(x, 2, sd, na.rm = TRUE),
                          'Étendue' = apply(x, 2, max, na.rm = TRUE)-apply(x, 2, min, na.rm = TRUE),
                          'Cov' = (apply(x, 2, sd, na.rm = TRUE)/apply(x, 2, mean, na.rm = TRUE))*100) %>%
  t() %>% # transposer le tableau
  as.data.frame()

write.csv(traits_stat,"02_outdata/traits/traits_stat.csv", row.names = TRUE)
