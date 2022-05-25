# written by Marian Schoen
# this script installs all necessary packages from 
# 'gitlab.spang-lab.de' for 'celloscope'


# all packages from our gitlab server
# For those: 
# (1) `git clone` the repository in the Dockerfile
# (2) add the path here, and install them from local file
spang.gitlab.packages <- c(
  "dtdpipeline", "diggeR"
)

all.installed.packages <- installed.packages()

for(package in spang.gitlab.packages){
  print(paste0("Checking ", package))
  if(package %in% all.installed.packages[, "Package"]){
    print("--- already installed, next")
    next  
  } else {
    install.packages(
      pkgs = package
      , type = "source"
      , repos = NULL
    )
    print(paste0(package, "--- installed"))
  }
}