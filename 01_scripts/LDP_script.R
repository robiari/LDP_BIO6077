## 2023-09-07
## Ariane Roberge

# Load the data file
install.packages("remotes")
remotes::install_github("lter/lterdatasampler")

if(!"tinytex" %in% rownames(installed.packages()))
  install.packages("tinytex")

tinytex::install_tinytex()

if(!"remotes" %in% rownames(installed.packages()))
  install.packages("remotes")
remotes::install_github("crsh/prereg")
