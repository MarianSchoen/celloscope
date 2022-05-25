# written by Marian Schoen
# this script installs all necessary bioconductor packages for 'celloscope'

# all bioconductor packages:
bioc.packages <- c(
  "rhdf5" 
)

all.installed.packages <- installed.packages()

for(package in bioc.packages){
  print(paste0("Checking ", package))
  if(package %in% all.installed.packages[, "Package"]){
    print("--- already installed, next")
    next  
  } else {
    BiocManager::install(package, ask = FALSE)  
    print(paste0(package, "--- installed"))
  }
}
