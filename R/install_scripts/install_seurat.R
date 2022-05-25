# written by Marian Schoen
# this script installs all necessary CRAN packages for 'celloscope'

r <- getOption("repos")
r["CRAN"] <- "http://cran.us.r-project.org"
options(repos = r, timeout = 3600);


# all CRAN packages: 
cran.packages <- c(
  "Seurat", "SeuratObject"
)

all.installed.packages <- installed.packages()

for(package in cran.packages){
  print(paste0("Checking ", package))
  if(package %in% all.installed.packages[, "Package"]){
    print("--- already installed, next")
    next  
  } else {
    install.packages(package)    
    print(paste0(package, "--- installed"))
  }
}
