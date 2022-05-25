# written by Marian Schoen
# this script installs all necessary CRAN packages for 'celloscope'

r <- getOption("repos")
r["CRAN"] <- "http://cran.us.r-project.org"
options(repos = r, timeout = 3600);

n.cpu <- system(
  command = "echo $(grep -c processor /proc/cpuinfo)"
  , intern = TRUE
  , ignore.stderr = FALSE
  )

if(!is.na(as.numeric(n.cpu))){
  n.cpu <- as.numeric(n.cpu)
} else {
  n.cpu <- 1
}
print(paste0("I assume that there are ", n.cpu, " useable cpus."))

cran.packages <- c(
  "jsonlite", "viridis", "viridisLite", "ggplot2",
  "digest", "tibble", "config", "rjson", "httr", "shinyalert", "shinyWidgets", 
  "shiny", "stats", "graphics", "grDevices", "utils", "datasets", "methods", 
  "base", "devtools", "BiocManager", "shinyjs", "shiny.router"
)

all.installed.packages <- installed.packages()

for(package in cran.packages){
  print(paste0("Checking ", package))
  if(package %in% all.installed.packages[, "Package"]){
    print("--- already installed, next")
    next  
  } else {
    install.packages(
      pkgs = package
      , Ncpus = n.cpu
      , verbose = FALSE
      )
    print(paste0(package, "--- installed"))
    }
}
