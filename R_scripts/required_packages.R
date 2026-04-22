# general function
install.packages("LatticeKrig")
install.packages("fields")
install.packages("spam64")
install.packages("here")
install.packages("sf")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(
  "rhdf5",
  ask    = FALSE,    
  update = FALSE 
)

# for timing 
install.packages("tictoc")

# for visuals 
install.packages("maps")
install.packages("cmocean")
install.packages("scico")
install.packages("RColorBrewer")
install.packages("zoo") 
install.packages("mgcv") 
install.packages("dplyr")

# for hpc
install.packages("doParallel")
install.packages("foreach")
install.packages("sp")
install.packages("gstat")

