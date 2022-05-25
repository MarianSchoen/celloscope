# this script set's up the environment, to start up the celloscope app: 

library(shiny)

# First I check if all required packages are installed: 
starting <- TRUE
missing.packages <- c()

required.packages <- c(
  "shinyWidgets", "diggeR", "shinyalert", "httr", "rjson", "ggplot2", 
  "viridis", "dtdpipeline", "Seurat", "rhdf5", "shinyjs", "shiny.router"
)

all.installed.packages <- installed.packages()
missing.packages <- which(!required.packages %in% all.installed.packages)
if(length(missing.packages) > 0){
  warning(
    paste0(
      "There are missing packages: ", 
      paste0(required.packages[missing.packages], collapse = ", ")
    )
  )
  starting <- FALSE
} else {
  for(package in required.packages){
    library(package, character.only = TRUE)
  }
}

if(!starting){
  fluidPage(
    h3("the following packages are not installed:"),
    h3(paste0(required.packages[missing.packages], collapse = ", ")),
    h3("Please install them and restart the App")
  )
}

source("R/helper.R")

useShinyalert(force =TRUE)