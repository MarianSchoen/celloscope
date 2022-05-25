source("R/helper.R")
#TODO: I don't know whether the 'path.to.data' may ever change, but maybe then 
# I should have it in some kind of config file instead of here? 
path.to.data <- "./data/shinyUseableData/"

available.datasets <- list_datasets(
  path = path.to.data
)
# list_datasets returns names, here I do not need them. 
names(available.datasets) <- NULL

celloscope.page <- shinyUI(
  fluidPage(
    # I decided to use a sidebar on the left with all the sliders, 
    # and a mainPanel with the plots, and text (see below)
    sidebarLayout( 
      # organized the left side (all the settings) for a better structure in 5 Tabs:
      # Loading; Modify Plot; Genes (consists of Pull Genes and calculate Pathway;
      # probably need a better name for that); Saving; and Glacier
      sidebarPanel(
        # # a button for debugging
        # actionButton(
        #   inputId = "brow"
        #   , label = "browser"
        # ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            title = "Loading",
            h3("Load a dataset"),
            selectInput(
              inputId = "dataset", label = "Which dataset?",
              choices = available.datasets,
              selected = available.datasets[1]
            ), 
            fluidRow(
              column(
                4,
                actionButton(inputId = "reload.dataset", label = "Load the data")
              ),
              column(
                8,
                textOutput(outputId = "loaded.dataset")
              ),
              
            ),
            fluidRow(
              column(
                6,
                actionButton(inputId = "fill.pheno",
                             label = "Fill dataset with QC-metrics"),
              ),
              column(
                6,
                textOutput(outputId = "filled.dataset")
              ),
            ),
            
            h5("(dataset doesn't need to be loaded!)"),
            fileInput(
              inputId = "dragdropupload"
              , label = "Drag and drop here"
              , multiple = TRUE
              , buttonLabel = "Browse..."
              , placeholder = "No file selected"
            ),
            h3("Load parameters for the plot"),
            selectInput(
              inputId = "parameters", label = "Load saved parameters",
              choices = c("None"), selected = "None"
            )
          ),
          tabPanel(
            title = "Modifying Plots",
            # this Tab is again organized in 3 Tabs:
            # one for subsetting data and two for the two scatter plots
            tabsetPanel(
              type = "tabs",
              tabPanel(title = "Subsetting",
                       h3("reduce dataset to 10 % samples"),
                       selectInput(
                         label = "deactivated; load a dataset first",
                         inputId = "reduce.col",
                         choices = c("None"),
                         selected = "None"
                       ),
                       switchInput(inputId = "reduce.dataset"
                                   , labelWidth = 0
                                   , value = FALSE
                                   , size = "small"
                                   , onLabel = "reduced"
                                   , offLabel = "not reduced"),
                       h3("Data Subsetting"),
                       actionButton(inputId ="reload.all",
                                    label = "Reload all plots"),
                       # we have up to 4 sliders with a SelectInput 
                       # for selecting the pheno entry,
                       # if the pheno has enough numerical entries.
                       # but they will load as soon as
                       # the dataset is loaded
                       selectInput(
                         label = "deactivated; load a dataset first",
                         inputId = "first.slider.col",
                         choices = c("None"),
                         selected = "None"
                       ),
                       sliderInput(
                         inputId = "first.slider"
                         , label = "deactivated; load a dataset first"
                         , min = 0
                         , max = 0
                         , value = c(0,0)
                       ),
                       selectInput(
                         label = "deactivated; load a dataset first",
                         inputId = "second.slider.col",
                         choices = c("None"),
                         selected = "None"
                       ),
                       sliderInput(
                         inputId = "second.slider"
                         , label = "deactivated; load a dataset first"
                         , min = 0
                         , max = 0
                         , value = c(0,0)
                       ),
                       selectInput(
                         label = "deactivated; load a dataset first",
                         inputId = "third.slider.col",
                         choices = c("None"),
                         selected = "None"
                       ),
                       sliderInput(
                         inputId = "third.slider"
                         , label = "deactivated; load a dataset first"
                         , min = 0
                         , max = 0
                         , value = c(0,0)
                       ),
                       selectInput(
                         label = "deactivated; load a dataset first",
                         inputId = "fourth.slider.col",
                         choices = c("None"),
                         selected = "None"
                       ),
                       sliderInput(
                         inputId = "fourth.slider"
                         , label = "deactivated; load a dataset first"
                         , min = 0
                         , max = 0
                         , value = c(0,0)
                       ),
                       # next are up to two selectize Input
                       # (depending on the available pheno entries)
                       # for the categorical subsetting
                       # both selectize inputs have a select input
                       # for selecting the pheno entry, a switch,
                       # wether to in- or exclude the selected 
                       # features and a button to clear all selected features
                       # because the App doesn't register excluding the last
                       # selected feature
                       selectInput(
                         inputId = "first.selectize.col",
                         label = "deactivated; load a dataset first",
                         choices = c("None"),
                         selected = "None"
                       ),
                       fluidRow(
                         column(
                           width = 8,
                           switchInput(
                             inputId = "in.ex.first.selectize",
                             onLabel = "include selected features",
                             offLabel = "exclude selected features",
                             labelWidth = 0,
                             value = FALSE,
                             size = "small"
                           ),
                         ),
                         column(
                           width = 4,
                           actionButton(
                             inputId = "clear.first.selectize",
                             label = "clear selected"
                           )
                         ),
                         
                       ),
                       selectizeInput(
                         inputId = "first.selectize",
                         label = "deactivated; load a dataset first",
                         choices = NULL,
                         multiple = TRUE
                       ),
                       selectInput(
                         inputId = "second.selectize.col",
                         label = "deactivated; load a dataset first",
                         choices = c("None"),
                         selected = "None"
                       ),
                       
                       fluidRow(
                         column(
                           width = 8,
                           switchInput(
                             inputId = "in.ex.second.selectize",
                             onLabel = "include selected features",
                             offLabel = "exclude selected features",
                             labelWidth = 0,
                             value = FALSE,
                             size = "small"
                           ),
                         ),
                         column(
                           width = 4,
                           actionButton(
                             inputId = "clear_second_selectize",
                             label = "clear selected"
                           )
                         ),
                         
                       ),
                       selectizeInput(
                         inputId = "second.selectize",
                         label = "deactivated; load a dataset first",
                         choices = NULL,
                         multiple = TRUE
                       ),
              ),
              # this panel is for modifying
              # the first scatter plot
              tabPanel(
                title = "Upper Scatter Plot",
                h3("Modifications Upper Scatter Plot"),
                actionButton(inputId = "first.reload.plot",
                             label = "Reload the upper plot"), 
                actionButton(
                  inputId = "first.reload.brush"
                  , label = "Reload the brush histogram"
                ),
                actionButton(
                  inputId = "first.show.brush.rect"
                  , label = "Show brushed rectangle"
                ),
                actionButton(
                  inputId = "first.remove.brush.rect"
                  , label = "remove brushed rectangle"
                ),
                actionButton(
                  inputId = "first.zoom.in"
                  , label = "Zoom in"
                ),
                actionButton(
                  inputId = "first.zoom.out"
                  , label = "Zoom out"
                ),
                selectInput(
                  inputId = "first.colour.scheme", label = "Which colour map?"
                  , choices = c("viridis", "rainbow", "ggplot-default", "heat", "reversed-heat")
                  , selected = "rainbow"
                ), 
                selectInput(
                  inputId = "first.colouring", label = "Which colouring?",
                  choices = c("None"),
                  selected = "None"
                ),
                selectInput(
                  inputId = "first.x.axis", label = "x axis?",
                  choices = c("None"),
                  selected = "None"
                ),
                selectInput(
                  inputId = "first.y.axis", label = "y axis?",
                  choices = c("None"),
                  selected = "None"
                )
              ),
              # for modifying the second scatter plot
              tabPanel(
                title = "Lower Scatter Plot",
                h3("Modifications Lower Scatter Plot"),
                actionButton(inputId = "second.reload.plot",
                             label = "Reload the lower plot"), 
                actionButton(
                  inputId = "second.reload.brush"
                  , label = "Reload the brush histogram"
                ),
                actionButton(
                  inputId = "second.show.brush.rect"
                  , label = "Show brushed rectangle"
                ),
                actionButton(
                  inputId = "second.remove.brush.rect"
                  , label = "remove brushed rectangle"
                ),
                actionButton(
                  inputId = "second.zoom.in"
                  , label = "Zoom in"
                ),
                actionButton(
                  inputId = "second.zoom.out"
                  , label = "Zoom out"
                ),
                selectInput(
                  inputId = "second.colour.scheme", label = "Which colour map?"
                  , choices = c("viridis", "rainbow", "ggplot-default", "heat", "reversed-heat")
                  , selected = "rainbow"
                ), 
                selectInput(
                  inputId = "second.colouring", label = "Which colouring?",
                  choices = c("None"),
                  selected = "None"
                ),
                selectInput(
                  inputId = "second.x.axis", label = "x axis?",
                  choices = c("None"),
                  selected = "None"
                ),
                selectInput(
                  inputId = "second.y.axis", label = "y axis?",
                  choices = c("None"),
                  selected = "None"
                )
              )
            ),
            verbatimTextOutput(outputId = "currently.subsetting")
          ),
          
          
          tabPanel(
            title = "Genes",
            # there is the possibility to "Pull genes", which means 
            # that the expression and the z-score of this gene will be available 
            # for the axis, the coloring, the sliders, etc.
            h3("Pull Genes"),
            # as there are way to many genes, this subsetting is done in 
            # R (server-side). 
            # see this link: https://shiny.rstudio.com/articles/selectize.html
            selectizeInput(
              inputId = "pullable.genes", label = "Want to pull a gene?",
              choices = NULL, multiple = TRUE
            ),
            textInput(inputId = "pull.genes.cut",
                      label = paste0("cut-value (either a quantil ", 
                                     "(ending with '%'), a expr-value, or 'mean')"),
                      value = "50%"
            ),
            actionButton(inputId = "pull.genes", label = "Pull the genes"), 
            verbatimTextOutput("pulled.a.gene"),
            
            h3("Calculate Pathway"),
            h5(
              "Select genes either by the selectize input, or in the text field (',' - separated!!!):"
            ),
            
            selectizeInput(
              inputId = "pullable.genes.pathway", label = "Genes in your pathway?",
              choices = NULL, multiple = TRUE
            ),
            textInput(inputId = "genes.pathway.text"
                      , label = "Genes in your pathway (comma seperated):"
                      , value = ""
            ),
            textInput(inputId = "pathway.name"
                      , label = "Name of your pathway:"
                      , value = ""
            ),
            actionButton(
              inputId = "calculate.pathways", label = "Calculate pathway enrichment"
            ), 
            verbatimTextOutput("pathway.output")
          ),
          
          
          # this panel is for saving the subsetted dataset or the
          # settings (sliders etc...)
          tabPanel(
            title = "Saving",
            h3("Save the dataset"),
            actionButton(inputId = "save.dataset",
                         label = "Save current data frame"),
            textOutput(outputId = "saved.dataset"),
            br(),
            
            textInput(inputId = "file.prefix.dataset"
                      , label = "file name for saving the dataset"
                      , value = gsub(pattern = " "
                                     , replacement = "__"
                                     , x = gsub(pattern = "-|:"
                                                , x = format(Sys.time()
                                                             , "%Y-%m-%d %H:%M"),
                                                replacement = "_"))
            ),
            actionButton(inputId = "get.Sys.time.dataset",
                         label = "get Sys.time as file name for saving"),
            
            h3("Save the parameters"),
            actionButton(inputId = "save.parameters",
                         label = "save subsetting parameters"),
            textOutput(outputId = "saved.parameters"),
            br(),
            textInput(inputId = "file.prefix.parameters"
                      , label = "file name for saving parameters"
                      , value = gsub(pattern = " "
                                     , replacement = "__"
                                     , x = gsub(pattern = "-|:"
                                                , x = format(Sys.time()
                                                             , "%Y-%m-%d %H:%M"),
                                                replacement = "_"))
            ),
            actionButton(inputId = "get.Sys.time.parameters",
                         label = "get Sys.time as file name for saving"),
            
            
            
            
            h3("Save the plots"),
            actionButton(inputId = "save.plots",
                         label = "save the selected plots"),
            textOutput(outputId = "saved.plots"),
            textInput(inputId = "file.prefix.plots"
                      , label = "file name for saving plots"
                      , value = paste0("plots_",gsub(pattern = " "
                                                     , replacement = "__"
                                                     , x = gsub(pattern = "-|:"
                                                                , x = format(Sys.time()
                                                                             , "%Y-%m-%d %H:%M"),
                                                                replacement = "_")))),
            actionButton(inputId = "get.Sys.time.plots",
                         label = "get Sys.time as file name for saving"),
            checkboxGroupInput(
              inputId = "scatter.save"
              , label = "scatter plots:"
              , choices = c("None")
              , selected = NULL
            ),
            checkboxGroupInput(
              inputId = "brush.save"
              , label = "brush histograms:"
              , choices = c("None")
              , selected = NULL
            ),
            checkboxGroupInput(
              inputId = "hist.save"
              , label = " slider histograms:"
              , choices = c("None")
              , selected = NULL
            )
            
          ), 
          # this panel is for the diggeR workaround
          tabPanel(
            title = "diggeR", 
            h3("diggeR workaround"), 
            textInput(
              inputId = "diggeR.token"
              , label = "diggeR access token"
            ),
            actionButton(
              inputId = "diggeR.save.token"
              , label = "save access token"
            ),
            tabsetPanel(
              type = "tabs",
              tabPanel(
                title = "basic search", 
                actionButton(
                  inputId = "show.tags"
                  , label = "Show common tags"
                ),
                actionButton(
                  inputId = "close.tags"
                  , label = "close tags"
                ),
                verbatimTextOutput(outputId = "diggeR.tags.text"),
                textInput(
                  inputId = "diggeR.search.tags"
                  , label = "Tags (';'-separated)"
                ), 
                textInput(
                  inputId = "diggeR.search.name"
                  , label = "dataset name" 
                )
              ),
              tabPanel(
                title = "advanced search", 
                actionButton(
                  inputId = "show.tags.ad"
                  , label = "Show common tags"
                ),
                actionButton(
                  inputId = "close.tags.ad"
                  , label = "close tags"
                ),
                verbatimTextOutput("diggeR.tags.text.ad"),
                textInput(
                  inputId = "diggeR.search.tags.ad"
                  , label = "Tags (';'-separated)"
                ), 
                textInput(
                  inputId = "diggeR.search.name.ad"
                  , label = "dataset name" 
                ),
                textInput(
                  inputId = "diggeR.search.author"
                  , label = "name of author" 
                ),
                textInput(
                  inputId = "diggeR.search.hash"
                  , label = "dataset hash" 
                ),
                textInput(
                  inputId = "diggeR.search.parent"
                  , label = "parent of given hash" 
                ),
                textInput(
                  inputId = "diggeR.search.commit.tags"
                  , label = "commit tags (';'-separated)" 
                ),
                textInput(
                  inputId = "diggeR.search.project"
                  , label = "name of the project within the dataset was generated" 
                )
              )
            ),
            actionButton(
              inputId = "diggeR.search"
              , label = "diggeR search"
            ),
            selectInput(
              inputId = "diggeR.datasets"
              , label = "Select datasets to download"
              , choices = NULL
              , multiple = TRUE
              , selected = NULL
            ),
            actionButton(
              inputId = "diggeR.download"
              , label = "diggeR download"
            ),
          )  
        )
      ),
      # this is now the main panel with all plots
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            title = "Visualizations", 
            verbatimTextOutput(outputId = "first.brush.text"),
            # first of all the two scatter plots
            plotOutput(
              outputId = "first.scatter.plot"
              , brush = "first.scatter.brush"
            ),
            verbatimTextOutput(outputId = "second.brush.text"),
            plotOutput(
              outputId = "second.scatter.plot"
              , brush = "second.scatter.brush"
            ),
            # then the histograms for the brushes
            fluidRow(
              splitLayout(
                cellWidths = c("50%", "50%")
                , plotOutput(outputId = "first.selected.hist.plot")
                , plotOutput(outputId = "first.not.selected.hist.plot")
              )
            ),
            fluidRow(
              splitLayout(
                cellWidths = c("50%", "50%")
                , plotOutput(outputId = "second.selected.hist.plot")
                , plotOutput(outputId = "second.not.selected.hist.plot")
              )
            ),
            # and then the histograms for the sliders
            fluidRow(
              splitLayout(
                cellWidths = c("50%", "50%")
                , plotOutput(
                  outputId = "hist.first.slider"
                  , brush = "first.hist.brush"
                )
                , plotOutput(
                  outputId = "hist.second.slider"
                  , brush = "second.hist.brush"
                )
              )
            ),
            fluidRow(
              splitLayout(
                cellWidths = c("50%", "50%")
                , plotOutput(
                  outputId = "hist.third.slider"
                  , brush = "third.hist.brush"
                )
                , plotOutput(
                  outputId = "hist.fourth.slider"
                  , brush = "fourth.hist.brush"
                )
              )
            )
          ), 
          tabPanel(
            title = "diggeR", 
            dataTableOutput(
              outputId = "diggeR.results"
            )
          )
        )
      )
    )
  )
)
