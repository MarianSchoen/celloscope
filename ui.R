library(shiny.router)

celloscope.version <- read.table("VERSION")[1,1]

# source("20_celloscope_page.R")
source("R/10_start_up.R")
source("routing.R")
source("R/helper.R")

addResourcePath("www", "www")

# set up the User Interface: 
shinyUI(
  fluidPage(
    # for session information, I need some cookies:
    tags$head(
      tags$script(
        src = paste0(
          "https://cdn.jsdelivr.net/npm/js-cookie@rc/",
          "dist/js.cookie.min.js"
        )
      ),
      tags$link(
        rel = "stylesheet"
        , type = "text/css"
        , href = "celloscope_style.css" 
        ), 
      tags$script(src = "www/script.js"), 
      tags$script(src = "www/fileInput_text.js")
    ),
    
    tags$ul(
      tags$li(
        a(
          href = route_link("/"), 
          paste0("Celloscope - v", celloscope.version)
          )),
      tags$li(a(href = route_link("contact"), "Contact"))
    ),
    router$ui, 
    
    p(
      "This website is provided by the ",
      actionLink(
        inputId = "institute-link"
        , label = "Institute of Functional Genomics"
        , onclick = "window.open('https://www.spang-lab.de', '_blank')"
      )
      , " - Statistical Bioinformatics Department", tags$br(),
      actionLink(
        inputId = "ur-link"
        , label = "University of Regensburg"
        , onclick = "window.open('https://www.uni-regensburg.de', '_blank')"
      ), "© 2022 Version 1.0", tags$br(),
      a(href = route_link("contact"), "Contact"), 
      " --- ",
      a(href = route_link("disclaimer"), "Legal Notice"),  
      " --- ",
      a(href = route_link("privacy"), "Privacy Policy") 
      , align = "center"
    )
  )
)