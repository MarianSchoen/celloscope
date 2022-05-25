
source("helper.R")
useShinyalert(force =TRUE)
path.to.data <- "./data/shinyUseableData/"

available.datasets <- list_datasets(
  path = path.to.data
)
# list_datasets returns names, here I do not need them. 
names(available.datasets) <- NULL

home_page <- div(
    # I decided to use a sidebar on the left with all the sliders, 
    # and a mainPanel with the plots, and text (see below)
    sidebarLayout( 
      # organized the left side (all the settings) for a better structure in 5 Tabs:
      # Loading; Modify Plot; Genes (consists of Pull Genes and calculate Pathway;
      # probably need a better name for that); Saving; and Glacier
      sidebarPanel(
        # a button for debugging
        actionButton(
          inputId = "brow"
          , label = "browser"
        )
      ),
      mainPanel()
    )
)

settings_page <- div(
  titlePanel("Settings"),
  p("This is a settings page")
)

contact_page <- div(
  titlePanel("Contact"),
  p("This is a contact page")
)

router <- make_router(
  route("/", home_page),
  route("settings", settings_page),
  route("contact", contact_page)
)

ui <- fluidPage(
  tags$ul(
    tags$li(a(href = route_link("/"), "Dashboard")),
    tags$li(a(href = route_link("settings"), "Settings")),
    tags$li(a(href = route_link("contact"), "Contact"))
  ),
  router$ui
)

server <- function(input, output, session) {
  router$server(input, output, session)
}

shinyApp(ui, server)