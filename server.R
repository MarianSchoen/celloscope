# written by Marian Schoen and Daniel Maar

celloscope.version <- read.table("VERSION")[1,1]

source("R/10_start_up.R")

source("ui.R")

source("R/global.R")

shinyServer(
  function(input, output, session) {
  
  router$server()

  # to be able to load large datasets via the drag and drop input
  options(shiny.maxRequestSize=1024*1024^2)
  # data stores all values which can change during runtime
  data <- reactiveValues()
  # statics holds all values which won't change
  statics <- list()

  # a text which holds the normally used digger tags
  statics$diggeR.tags <- paste0(
    "\n- sortedimmunecells for sorted RNA-seq data."
    , "\n- Optionally / if known, we further attach the sc technology: illumina10x, ...."
    , "\n- accession numbers: GSE11567, SRS123124 (in this format:"
    , "capital letters, no separator between GSE and accession number,...)"
    , "\n- surname of first author (e.g., Tirosh, first character is capital)"
    , "\n- celltypelabelled if there are celltype labels available"
    , "\n- quality_controlled if one has run dtdpipeline::quality_control(...)"
    , "on the data, and stored the resulting QC metrics"
    , "\n- disease type (melanoma, ccRCC...)"
  )

  # the path to the directory where all available datasets are stored
  statics$path.to.data <- "./data/shinyUseableData/"

  datasets <- list_datasets(
    path = statics$path.to.data
  )

  data$available.datasets <- datasets
  data$dataset.paths <- names(datasets)


  # wether the App should print out everything its doing into the console,
  # for debug reasons
  statics$verbose <- TRUE


  #if theres currently a dataset loaded
  data$loaded <- FALSE

  # every hour, I check whether there are datasets to delete. Therefore:
  delete.datasets <- reactiveTimer(intervalMs = 1000*60*60)
  observe(
    {
      delete.datasets()

      all.datasets <- read_datasetlist_csv(path = statics$path.to.data)

      valid.datasets <- all.datasets

      current.time <- format(Sys.time(), "%Y-%m-%d,%H:%M")
      for(file.idx in 1:nrow(all.datasets)){
        if(all.datasets$file.name[file.idx] == "mock_data"){
          next
        }
        file.time <- strptime(
          x = gsub(" ", "", all.datasets$uploaded[file.idx])
          , format = "%Y-%m-%d,%H:%M"
        )
        time.difference <- as.numeric(
          gsub(
            pattern = "Time difference of || days"
            , replacement = ""
            , difftime(
              time1 = current.time
              , time2 = file.time
              , units = "days"
            )
          )
        )
        if(time.difference > 1){
          valid.datasets <- valid.datasets[-file.idx, ]
          unlink(
            paste0(
              statics$path.to.data, all.datasets$file.name[file.idx], ".h5"
            )
          )
        }
        write_datasetlist_csv(
          path = statics$path.to.data
          , the.files = valid.datasets
        )
      }
    }
  )

  observeEvent(
    # to reduce a dataset to 10% of the samples
    eventExpr = input$reduce.dataset
    , handlerExpr = {
      if(statics$verbose){
        cat(file = stderr(),"\nin reduce.dataset")
      }

      if(data$loaded){
        start_loading_screen(session = session)

        if(input$reduce.dataset){
          if(statics$verbose){
            cat(file = stderr(),"\nreducing")
          }

          if(input$reduce.col == "completly random"){
            if(statics$verbose){
              cat(file = stderr(),"\nrandom")
            }
            # the random chosen samples
            chos.samples <- sample(
              x = rownames(data$slider_selectizes$booleans$selected)
              , size = round(nrow(data$slider_selectizes$booleans$selected))*0.1
              , replace = FALSE
            )

          } else {
            if(statics$verbose){
              cat(file = stderr(),paste0("\nwith col '", input$reduce.col, "'"))
            }

            col <- data$general$data.set[, input$reduce.col]
            names(col) <- data$general$data.set$colnames
            chos.samples <- c()
            # samples 10% of every category from the selected column
            for(entry in unique(col)){
              profiles <- names(col)[col == entry]
              chos.samples <- c(sample(
                x = profiles
                , size = ceiling(length(profiles)*0.1)
                , replace = FALSE
              ),
              chos.samples)
            }

          }

          # add.col is an additional column for the selected dataframe
          # add.col holds for every samples whether it is included in the
          # reduced dataset or not
          if(statics$verbose){
            cat(file = stderr(),"\nconstructing add.col")
          }
          add.col <- rep(FALSE, times = nrow(data$slider_selectizes$booleans$selected))
          names(add.col) <- rownames(data$slider_selectizes$booleans$selected)
          add.col[chos.samples] <- TRUE

          if(statics$verbose) cat(file = stderr(),"\ncbinding")
          data$slider_selectizes$booleans$selected$in.reduced <- add.col
          # currently not used
          data$use.old.subset <- FALSE
        } else {
          data$slider_selectizes$booleans$selected$in.reduced <- TRUE
        }
        end_loading_screen(session = session)
      }
      if(statics$verbose) cat(file = stderr(),"\nDONE with reduce.dataset")
    }
  )
  observeEvent(
    #saving the current subsetted dataset
    eventExpr = input$save.dataset
    , handlerExpr = {
      if(statics$verbose){
        cat(file = stderr(),"\nin save.dataset")
      }
      start_loading_screen(session = session)
      #the current loaded filev
      f.name <- data$dataset.paths[data$available.datasets == data$general$dataset.name]

      #getting all selected profiles and then saving the data
      selected.profiles <- select_profiles(
        data = data
        , statics = statics
      )
      if(length(selected.profiles) != 0){
        #reading it, so we can store everything afterwards
        dataset <- dtdpipeline::read_data(f.name)
        sc.pheno <- data$general$data.set
        dataset$additional.info$celltypes <- as.array(
          x = unique(dataset$sc.pheno$celltype)
        )
        dataset$additional.info$C.entries <- NULL
        write_data(
          sc.counts = dataset$sc.counts[, selected.profiles]
          , sc.pheno = sc.pheno[selected.profiles, ]
          , sc.feature.id.type = dataset$sc.feature.id.type
          , sc.normalization.type = dataset$sc.normalization
          , additional.info = dataset$additional.info
          , filename = paste0(
            statics$path.to.data
            , input$file.prefix.dataset
            , ".h5"
          )
          , force = TRUE
        )

        if(statics$verbose) cat(file = stderr(),"\nupdating SelectInput")
        # since there's now a new dataset available we need to include it
        # into the available datasets and update the selectInput
        dataset.choices <- list.files(
          path = statics$path.to.data,
          pattern = ".h5"
        )
        available.datasets <-gsub(
          pattern = ".h5"
          , replacement = ""
          , x = dataset.choices
        )
        updateSelectInput(
          inputId = "dataset",
          choices = available.datasets
        )
        data$output$texts$saved.dataset <- paste0(
          "saved dataset in: ", gsub(
            pattern = statics$path.to.data
            , x = f.name
            , replacement = paste0(statics$path.to.data, input$variables$file.prefix.dataset, "_"))
        )
        end_loading_screen(session = session)
      } else {
        shinyalert(
          title = "No cells selected!"
          , text = "No cells are selected via the sliders and the selectizes. Please select some cells."
          , closeOnEsc = TRUE
          , closeOnClickOutside = FALSE
          , showConfirmButton = TRUE
          , session = session
          , timer = 0
          , animation = FALSE
          , size = "m"
          , immediate = TRUE
        )
      }
      if(statics$verbose) cat(file = stderr(),"\nDONE with save.dataset")
    }
  )
  #The next observe event are just to synchronize both digger searches
  observeEvent(
    input$diggeR.search.name,
    {
      updateTextInput(
        session = session
        , inputId = "diggeR.search.name.ad"
        , value = input$diggeR.search.name
      )
    }
  )
  observeEvent(
    input$diggeR.search.name.ad,
    {
      updateTextInput(
        session = session
        , inputId = "diggeR.search.name"
        , value = input$diggeR.search.name.ad
      )
    }
  )
  observeEvent(
    input$diggeR.search.tags,
    {
      updateTextInput(
        session = session
        , inputId = "diggeR.search.tags.ad"
        , value = input$diggeR.search.tags
      )
    }
  )
  observeEvent(
    input$diggeR.search.tags.ad,
    {
      updateTextInput(
        session = session
        , inputId = "diggeR.search.tags"
        , value = input$diggeR.search.tags.ad
      )
    }
  )
  observeEvent(
    input$show.tags,
    {
      show_tags(
        statics = statics
        , data = data
      )
    }
  )
  observeEvent(
    input$show.tags.ad,
    {
      show_tags(
        statics = statics
        , data = data
      )
    }
  )
  observeEvent(
    input$close.tags,
    {
      close_tags(data = data)
    }
  )
  observeEvent(
    input$close.tags.ad,
    {
      close_tags(data = data)
    }
  )
  observeEvent(
    # to search datatomb
    input$diggeR.search,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin diggeR.search\npreparing search")
      }
      if(is.null(data$general$diggeR.token.av)){
        data$general$diggeR.token.av <- FALSE
      }
      if(!data$general$diggeR.token.av){
        shinyalert(
          title = "missing access token"
          , text = paste0(
            "You didn't submit an access token. You can't"
            , " search for and download datasets without an access"
            , "token!")
          , closeOnEsc = TRUE
          , closeOnClickOutside = TRUE
          , showConfirmButton = TRUE
          , session = session
          , animation = FALSE
          , size = "m"
        )
      } else{
        start_loading_screen(session = session)
        # getting the search parameters
        if(input$diggeR.search.author == ""){
          author <- NA
        }
        if(input$diggeR.search.tags == ""){
          tags <- NULL
        } else {
          # getting the tags in a digger readable format
          tags <- strsplit(
            x = input$diggeR.search.tags
            , split = ";")
          tags <- gsub(
            x = unlist(tags)
            , pattern = " "
            , replacement = ""
          )
          if("" %in% tags) tags <- tags[-which(tags == "")]
        }
        if(input$diggeR.search.name == ""){
          dname <- NA
        }
        if(input$diggeR.search.hash == ""){
          hash <- NA
        }
        if(input$diggeR.search.parent == ""){
          parent <- NA
        }
        if(input$diggeR.search.commit.tags == ""){
          commit.tags <- NA
        } else {
          commit.tags <- strsplit(
            x = input$diggeR.search.commit.tags
            , split = ";")
          commit.tags <- gsub(
            x = unlist(commit.tags)
            , pattern = " "
            , replacement = ""
          )
          if("" %in% commit.tags) commit.tags <- commit.tags[-which(commit.tags == "")]
        }
        if(input$diggeR.search.project == ""){
          project <- NA
        }
        if(statics$verbose) cat(file = stderr(),"\nsearching...")
        results <- diggeR::search(
          tags = unique(c(tags, "sc"))
          , name = dname
          , author = author
          , parentOf = parent
          , hash = hash
          , commit = unique(commit.tags)
          , projectname = project
          ,
        )
        if(nrow(results) == 0){
          shinyalert(
            title = "no results"
            , text = paste0(
              "Your search request didn't match any datatsets"
              , ". Be aware that the uppercases and lowercases"
              , " lead to different results!"
            )
            , size = "m"
            , session = session
            , immediate = TRUE
            , closeOnEsc = TRUE
            , closeOnClickOutside = TRUE
            , showConfirmButton = TRUE
            , animation = FALSE
          )
        } else {
          # preparing the results
          all.datasets <- strsplit(x = results$name, split = "/")
          dataset.names <- c()
          if(statics$verbose) cat(file = stderr(),"\npreparing results")

          for(dataset in all.datasets){
            dataset.names <- c(dataset.names, dataset[length(dataset)])

          }
          dataset.names <- gsub(pattern = ".h5", replacement = "", x = dataset.names)
          results$name <- dataset.names
          hashes <- results$hash
          names(hashes) <- dataset.names
          data$diggeR$hashes <- hashes
          all.metas <- list()
          all.tags <- c()
          all.C.entries <- c()
          all.celltypes <- c()
          all.normalizations <- c()
          # preparing every found datasets separatly
          for(i in 1:length(hashes)){
            meta <- metadata(hashes[i])
            all.metas <- c(all.metas, list(meta))
            names(all.metas)[length(all.metas)] <- names(hashes[i])
            # we don't need this part of the metadata, since it is to complex
            # to show reasonable
            meta$data$Rsessioninfo <- NULL
            tags <- Reduce(
              f = function(x,y){return(paste(x,y, sep = ", "))}
              , x = meta$tags)
            if(is.null(tags)){
              all.tags[i] <- "None"
            } else {
              all.tags[i] <-tags
            }
            C.entries <- Reduce(
              f = function(x,y){return(paste(x,y, sep = ", "))}
              , x = meta$data$C.entries)
            if(is.null(C.entries)){
              all.C.entries[i] <- "None"
            } else {
              all.C.entries[i] <- C.entries
            }
            celltypes <- Reduce(
              f = function(x,y){return(paste(x,y, sep = ", "))}
              , x = meta$data$celltypes)
            if(is.null(celltypes)){
              all.celltypes[i] <- "None"
            } else {
              all.celltypes[i] <-celltypes
            }
            data.normalization <- unlist(meta$data$data.normalization)
            dnorm <- ""
            for(j in 1:length(data.normalization)){
              dnorm <- paste0(dnorm
                              , names(data.normalization)[j], ": "
                              , data.normalization[j])
              if(j != length(data.normalization)) dnorm <- paste0(dnorm, "; ")
            }
            all.normalizations[i] <- dnorm
          }

          results <- data.frame(
            name = results$name
            , description = results$description
            , C.entries = all.C.entries
            , celltypes = all.celltypes
            , tags = all.tags
            , data.normalization = all.normalizations
          )
          data$output$frames$diggeR$info <- results
          updateSelectInput(
            inputId = "diggeR.datasets"
            , label = "Select datasets to download"
            , choices = dataset.names
            , selected = NULL
          )
          output$diggeR.results <- renderDataTable(
            expr = {data$output$frames$diggeR$info}
          )
          end_loading_screen(session = session)
        }
      }

      if(statics$verbose) cat(file = stderr(),"\nDONE with diggeR.search")
    }
  )
  observeEvent(
    input$diggeR.save.token,
    # to save a token for diggeR
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin diggeR.save.token")
      }
      reset_brushes(session = session)
      # try to set the token
      try_err <- try(set_token(input$diggeR.token)
                     , silent = TRUE)
      if(class(try_err) == "try-error" | try_err == "anonymous"){
        shinyalert(
          title = "invalid token"
          , text = paste0("Your submitted token is invalid. Please submit a valid one! "
                          , "Access tokens are handed out by the auth server, in our case: "
                          , "https://auth.spang-lab.de")
          , closeOnEsc = TRUE
          , closeOnClickOutside = TRUE
          , showCancelButton = FALSE
          , session = session
          , showConfirmButton = TRUE
          , animation = FALSE
          , size = "m"
        )
      } else {
        shinyalert(
          title = "Token successfully submitted"
          , text = paste0("Your submitted token is valid. Now you can search for"
                          , " and download datasets!")
          , closeOnEsc = TRUE
          , closeOnClickOutside = TRUE
          , showCancelButton = FALSE
          , session = session
          , showConfirmButton = TRUE
          , animation = FALSE
          , size = "m"
        )
        data$general$diggeR.token.av = TRUE
      }
      updateTextInput(
        session = session
        , inputId = "diggeR.token"
        , value = ""
      )
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with diggeR.save.token")
      }
    }
  )
  observeEvent(
    input$diggeR.download,
    # to download a dataset
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin diggeR.download")
      }

      reset_brushes(session = session)
      if(!is.null(input$diggeR.datasets)){
        start_loading_screen(session = session)
        for(dataset in input$diggeR.datasets){
          file_name <- download(
            data$diggeR$hashes[dataset]
            , file = paste0(
              statics$path.to.data
              , names(data$diggeR$hashes[dataset])
              , ".h5"))

          if(statics$verbose) cat(file = stderr(),"\nupdating available datasets")
          data$available.datasets <- c(
            data$available.datasets
            , names(data$diggeR$hashes[dataset])
          )
          data$dataset.paths <- c(
            data$dataset.paths
            , paste0(
              statics$path.to.data
              , names(data$diggeR$hashes[dataset])
              , ".h5"
            )
          )
        }
        updateSelectInput(
          inputId = "dataset"
          , choices = data$available.datasets
        )
        end_loading_screen(session = session)
        if(statics$verbose) cat(file = stderr(),"\nDONE with diggeR.download")
      }
    })

  # for debug reasons uncommen this
  # observeEvent(input$brow,{browser()})

  observeEvent(
    input$first.show.brush.rect,
    {
      if(data$loaded){
        show_brush_rect(
          statics = statics,
          data = data,
          id = "first",
          session = session
        )
        reset_brushes(session = session)

        output$first.scatter.plot <- renderPlot(({data$output$plots$scatter$first}))
      }
    }
  )
  observeEvent(
    input$first.remove.brush.rect,{
      if(data$loaded){
        remove_brush_rect(
          id = "first",
          data = data,
          statics = statics,
          session = session
        )
        reset_brushes(session = session)

        output$first.scatter.plot <- renderPlot(({data$output$plots$scatter$first}))
      }
    }
  )
  observeEvent(
    input$second.show.brush.rect,
    {
      if(data$loaded){
        show_brush_rect(statics = statics,
                        data = data,
                        id = "second",
                        session = session)
        reset_brushes(session = session)
        output$second.scatter.plot <- renderPlot(({data$output$plots$scatter$second}))
      }
    }
  )
  observeEvent(
    input$second.remove.brush.rect,{
      if(data$loaded){
        remove_brush_rect(id = "second",
                          data = data,
                          statics = statics,
                          session = session)
        reset_brushes(session = session)
        output$second.scatter.plot <- renderPlot(({data$output$plots$scatter$second}))
      }
    }
  )
  observeEvent(
    # reloads both scatterplots
    input$reload.all,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin reaload.all")
      }
      reset_brushes(session = session)
      if(data$loaded){
        start_loading_screen(session = session)
        if(statics$verbose){
          cat(file = stderr(),"\nupdateing first plot")
        }
        plot_scatter_plot(
          data = data, statics = statics
          , input = input
          , id = "first"
          , session = session
        )
        if(statics$verbose){
          cat(file = stderr(),"\nupdating second plot")
        }
        plot_scatter_plot(
          data = data, statics = statics
          , input = input
          , session = session
          , id = "second"
        )
        if(statics$verbose){
          cat(file = stderr(),"\nupdating histograms")
        }
        plot_histograms(
          data = data, statics = statics
          , input = input
          , session = session
        )
        output$hist.first.slider <- renderPlot(({data$output$plots$slider.histograms$first}))
        output$hist.second.slider <- renderPlot(({data$output$plots$slider.histograms$second}))
        output$hist.third.slider <- renderPlot(({data$output$plots$slider.histograms$third}))
        output$hist.fourth.slider <- renderPlot(({data$output$plots$slider.histograms$fourth}))
        output$first.scatter.plot <- renderPlot(({data$output$plots$scatter$first}))
        output$second.scatter.plot <- renderPlot(({data$output$plots$scatter$second}))

        end_loading_screen(session = session)
      }

      if(statics$verbose){
        cat(file = stderr(),"\n\nDONE with reaload.all")
      }

    }
  )
  observeEvent(
    # fills the pheno of the selected dataset with the qc metrics through executing
    # the quality control function
    input$fill.pheno,
    {
      start_loading_screen(session = session)
      if(statics$verbose){
        reset_brushes(session = session)
        cat(file = stderr(),"\n\nin fill pheno")
        cat(file = stderr(),"\nexecuting quality_control")
      }
      # getting the correct dataset name
      if(!grepl(pattern = ".h5$", x = data$general$dataset.name)){
        f.name <- paste0(data$general$dataset.name, ".h5")
      }
      if(!grepl(pattern = statics$path.to.data, x = f.name)){
        f.name <- paste0(statics$path.to.data, f.name)
      }
      # the quality_control-script can calculate the qc-metrics and
      # than store the result again
      dtdpipeline::quality_control(
        file.name = f.name,
        verbose = FALSE,
        print.output = FALSE,
        problem.fun = warning,
        restore.data = TRUE,
        calculate.qc.metrics = TRUE
      )
      data$output$texts$filled.dataset <- paste0(
        "filled dataset: '", data$general$dataset.name,
        " with QC-metrics. Please load it again to update the pheno"
      )
      output$filled.dataset <- renderText({paste0(data$output$texts$filled.dataset)})
      end_loading_screen(session = session)
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with fill_pheno")
      }
    })

  observeEvent(
    #loading the selected dataset
    eventExpr = input$reload.dataset,
    handlerExpr = {
      if(statics$verbose){
        cat(file = stderr(),"\nin reload_dataset\nresetting data")
      }

      all.data.sets <- list_datasets(
        path = statics$path.to.data
        , owner.cookie = input$cookies$`session-number`
      )

      data$loaded <- TRUE
      data$output <- NULL
      data$plots <- NULL
      data$genes <- NULL
      data$slider_selectizes <- NULL
      data$diggeR <- NULL
      start_loading_screen(session = session)
      data$general$dataset.name <- input$dataset
      reset_brushes(session = session)
      # getting the pheno data
      pheno.data <- load_data_set(
        load.expression = FALSE
        , session = session
        , verbose = statics$verbose
        , path.to.data = data$dataset.paths[data$available.datasets == data$general$dataset.name]
      )
      if(!is.null(pheno.data)){
        data$general$data.set <- pheno.data
        #since the pullable genes are just the rownames of the expr. matris
        #we're just loading it via load_data_set
        pullable.genes <- load_data_set(
          load.expression = FALSE
          , verbose = statics$verbose
          , session = session
          , load.row.names.exp = TRUE
          , path.to.data = data$dataset.paths[data$available.datasets == data$general$dataset.name]
        )
        if(!is.null(pullable.genes)){
          HGNC.genes <- as.vector(unique(
            HGNC_ENSGID$hgnc_symbol[HGNC_ENSGID$ensembl_gene_id %in% pullable.genes]
          ))
          data$genes$pullable.genes <- HGNC.genes
          pheno.data <- data$general$data.set


          if(statics$verbose){
            cat(file = stderr(),"\ngetting the slider data\ngetting is.slider")
          }
          # is.slider indicates whether a column of the pheno data
          # can be used as a numeric slider (is numeric)
          is.slider <- rep(x = F, times = ncol(pheno.data))
          for(i in 1:ncol(pheno.data)) {
            is.slider[i] <- is.numeric(pheno.data[ ,i])
          }
          names(is.slider) <- colnames(pheno.data)
          # since we need this multiple times we store it in data
          data$slider_selectizes$booleans$is.numeric <- is.slider
          # whether the switchs for the categorical subsetting (selectizeInput)
          # are TRUE or FALSE
          # The switches themself store this information but as you will see later
          # we need to store this separatly (mostly to save some lines of code =))...
          data$slider_selectizes$booleans$selectize.switch <- c(FALSE,FALSE)

          if(statics$verbose) cat(file = stderr(),"\ngetting ranges and slider.values")
          # the ranges list holds the min, the max and the stepsize for every entry
          # of the pheno data which can be a slider
          ranges <- list()
          # slider.values stores the current values of every entry
          # of the pheno data which can be a slider
          slider.values <- list()
          #filling slider.values and ranges
          for(col in 1:ncol(pheno.data)){
            if(is.slider[col]){
              if(statics$verbose) cat(file = stderr(),paste0("\nexecuting min_max_step for ", names(is.slider)[col]))
              ranges$tmp <- min_max_step(pheno.data[ ,col])
              slider.values$tmp <- ranges$tmp[1:2]
              names(ranges)[length(ranges)] <- names(is.slider[col])
              names(slider.values)[length(slider.values)] <- names(is.slider[col])
            }
          }
          data$slider_selectizes$numerics$ranges <- ranges
          data$slider_selectizes$numerics$slider.values <- slider.values
          ids <- c(
            "first.colouring", "second.colouring", "first.x.axis", "second.x.axis",
            "first.y.axis", "second.y.axis"
          )
          if(statics$verbose){
            cat(file = stderr(),"\nupdating plot options")
          }
          i <- 1
          for(input.id in ids){

            if(statics$verbose){
              cat(file = stderr(),paste0("\nupdating ", input.id))
            }

            # these are our "standards" for the axis and the colouring
            if(grepl(pattern = "colouring", x = input.id) & "celltype" %in% colnames(data$general$data.set)){
              input.val <- "celltype"
            } else if(grepl(pattern = "x.axis", x = input.id) & any(grepl(
              pattern = "UMAP"
              , x = colnames(data$general$data.set)))){
              UMAPS <- colnames(data$general$data.set)[grepl(
                pattern = "UMAP"
                , x = colnames(data$general$data.set))]
              input.val <- UMAPS[grepl(
                pattern = "1",
                x = UMAPS
              )]
            } else if(grepl(pattern = "y.axis", x = input.id) & any(grepl(
              pattern = "UMAP"
              , x = colnames(data$general$data.set)))){
              UMAPS <- colnames(data$general$data.set)[grepl(
                pattern = "UMAP"
                , x = colnames(data$general$data.set))]

              input.val <- UMAPS[grepl(
                pattern = "2",
                x = UMAPS
              )]
            } else {
              input.val <- colnames(data$general$data.set)[i]
              i <- i + 1
            }
            if(length(input.val) == 0){
              input.val <- colnames(data$general$data.set)[i]
              i <- i + 1
            }
            if(grepl(pattern = "axis", x = input.id)){
              if(all(is.numeric(data$general$data.set[, input.val]))){

                if(grepl("x", x = input.id)){
                  data$plots[[gsub(x = input.id, pattern = "...axis", replacement = ".scatter.plot")]]$x.lims <- data$slider_selectizes$numerics$ranges[[input.val]][1:2]

                } else {
                  data$plots[[gsub(x = input.id, pattern = "...axis", replacement = ".scatter.plot")]]$y.lims <- data$slider_selectizes$numerics$ranges[[input.val]][1:2]

                }
              } else {
                if(grepl("x", x = input.id)){
                  data$plots[[gsub(x = input.id, pattern = "...axis", replacement = ".scatter.plot")]]$x.lims <- unique(data$general$data.set[, input.val])
                } else {
                  data$plots[[gsub(x = input.id, pattern = "...axis", replacement = ".scatter.plot")]]$y.lims <- unique(data$general$data.set[, input.val])
                }
              }
            }
            updateSelectInput(
              session = session
              , inputId =  input.id
              , choices =  sort(colnames(data$general$data.set))
              , selected = input.val
            )
          }
          if(statics$verbose){
            cat(file = stderr(),"\nupdating gene stuff")
          }
          # updating the SelectizeInputs
          ids <- c(
            "pullable.genes", "pullable.genes.pathway"
          )
          for(input.id in ids){
            updateSelectizeInput(
              session = session
              , inputId = input.id
              , choices = data$genes$pullable.genes
              , server = TRUE
            )
          }


          # data$slider_selectizes$numerics$n.categorical stores the number of categorical columns of the
          # pheno.data
          data$slider_selectizes$numerics$n.categorical <- length(which(!is.slider))
          # data$slider_selectizes$booleans$selected stores for every column of the pheno.data
          # all selected elements, which got selected via selectized input and
          # the switches. Therefore e.g. if a column contains the 3 categories
          # c(a,b,d) and selectizeInput has selected a and d, but the switch is
          # on FALSE (=exclude selected values) then every entry in selected
          # for this column which has the entry "b" in the pheno.data will be TRUE,
          # the remaining will be FALSE
          data$slider_selectizes$booleans$selected <- as.data.frame(matrix(
            data = TRUE
            , nrow = nrow(pheno.data)
            , ncol = length(which(!is.slider))
          ))
          rownames(data$slider_selectizes$booleans$selected) <- data$general$data.set$colnames
          colnames(data$slider_selectizes$booleans$selected) <- names(which(!is.slider))
          # update_sliders and update_selectizes initializes the sliders
          # and selectizes with the values from pheno.data
          update_sliders(
            data = data, statics = statics,
            session = session,
            input = input
          )
          update_selectizes(
            data = data, statics = statics,
            session = session,
            input = input
          )

          update_saves(
            data = data, statics = statics
            , id = "clear"
            , session = session
            , input = input
            , type = "scatter"
          )
          update_saves(
            data = data, statics = statics
            , id = "clear"
            , session = session
            , input = input
            , type = "brush"
          )
          data$output$texts$loaded.dataset <- paste0(
            "Currently loaded dataset: '", data$general$dataset.name, "'"
          )
          update_subsetting_output(
            data = data, statics = statics,
            input = input
          )
          if(!dir.exists(paste0(statics$path.to.data, "parameters"))) {
            dir.create(paste0(statics$path.to.data, "parameters"))
          }
          if(!dir.exists(paste0(statics$path.to.data, "/parameters/"))) dir.create(paste0(statics$path.to.data, "/parameters/"))
          available.parameters <- c("None", gsub(pattern = paste0("_selectize__",
                                                                  data$general$dataset.name,
                                                                  ".csv"),
                                                 replacement = "",
                                                 x = list.files(path = paste0(statics$path.to.data, "/parameters/"),
                                                                pattern = paste0("_selectize__", data$general$dataset.name))))
          updateSelectInput(session = session
                            , inputId = "parameters"
                            , choices =  available.parameters
                            , selected = available.parameters[1]
                            , label = "Load saved parameters"
          )
          reset_brushes(session = session)
          end_loading_screen(session = session)
        }

      }

      output$first.scatter.plot <- renderPlot(({data$output$plots$scatter$first}))
      output$second.scatter.plot <- renderPlot(({data$output$plots$scatter$second}))
      output$loaded.dataset <- renderText({paste0(data$output$texts$loaded.dataset)})
      output$filled.dataset <- renderText({paste0(data$output$texts$filled.dataset)})
      output$hist.first.slider <- renderPlot(({data$output$plots$slider.histograms$first}))
      output$hist.second.slider <- renderPlot(({data$output$plots$slider.histograms$second}))
      output$hist.third.slider <- renderPlot(({data$output$plots$slider.histograms$third}))
      output$hist.fourth.slider <- renderPlot(({data$output$plots$slider.histograms$fourth}))
      output$second.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$second.selected}))
      output$second.not.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$second.not.selected}))
      output$first.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$first.selected}))
      output$first.not.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$first.not.selected}))
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with reload dataset\n")
      }
    }
  )

  observeEvent(
    #reload the first scatter plot
    input$first.reload.plot,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin first.reload.plot")
      }
      reset_brushes(session = session)
      plot_scatter_plot(
        data = data, statics = statics
        , input = input
        , session = session
        , id = "first"
      )
      plot_histograms(
        data = data, statics = statics
        , input = input
        , session = session
      )
      output$hist.first.slider <- renderPlot(({data$output$plots$slider.histograms$first}))
      output$hist.second.slider <- renderPlot(({data$output$plots$slider.histograms$second}))
      output$hist.third.slider <- renderPlot(({data$output$plots$slider.histograms$third}))
      output$hist.fourth.slider <- renderPlot(({data$output$plots$slider.histograms$fourth}))
      output$first.scatter.plot <- renderPlot(({data$output$plots$scatter$first}))
      output$first.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$first.selected}))
      output$first.not.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$first.not.selected}))
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with first.reload.plot\n")
      }
    }
  )
  observeEvent(
    #reload the second scatter plot
    input$second.reload.plot,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin second.reload.plot")
      }
      reset_brushes(session = session)
      plot_scatter_plot(
        data = data, statics = statics
        , session = session
        , input = input
        , id = "second"
      )
      plot_histograms(
        data = data, statics = statics
        , input = input
        , session = session
      )
      output$hist.first.slider <- renderPlot(({data$output$plots$slider.histograms$first}))
      output$hist.second.slider <- renderPlot(({data$output$plots$slider.histograms$second}))
      output$hist.third.slider <- renderPlot(({data$output$plots$slider.histograms$third}))
      output$hist.fourth.slider <- renderPlot(({data$output$plots$slider.histograms$fourth}))
      output$second.scatter.plot <- renderPlot(({data$output$plots$scatter$second}))
      output$second.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$second.selected}))
      output$second.not.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$second.not.selected}))
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with second.reload.plot\n")
      }
    }
  )

  observeEvent(
    # pull genes
    input$pull.genes,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin pull.genes")
      }
      start_loading_screen(session = session)
      reset_brushes(session = session)
      cut <- input$pull.genes.cut
      # if we should cut by the mean, cut will be 'mean'
      # otherwise it has to be a numeric with possible a '%' at the end
      b.perc <- grepl(pattern = "%", x = cut)
      if(cut != "mean"){
        cut <- as.numeric(gsub(pattern = "%", x = cut, replacement = ""))
      }
      # if it isn't a numeric (or 'mean') cut will now be NA
      if(is.na(cut)){
        if(statics$verbose){
          cat(file = stderr(),"\nno cut value => doing nothing")
        }
        data$output$texts$pulled.genes <- paste0(
          "please provide either a cut expression value, a quantil, or 'mean'"
        )
      } else {
        # if the cut is a quantil it should be between 0 and 100
        if((b.perc & (cut>0 | cut<100)) | (!b.perc)){

          if(statics$verbose){
            cat(file = stderr(),"\nloading the data")
          }
          the.genes <- input$pullable.genes
          # those genes for which a zscore value should be calculated
          zscore.genes <- c()
          # those genes for which a cut value should be calculated
          cut.genes <- c()
          if(b.perc){
            cut.name <- paste0(cut, "%")
          }
          else{
            cut.name <- cut
          }
          for(gene in the.genes){
            # checking if we pulled this gene previously
            if(gene %in% names(data$genes$pulled.genes)){
              # checking if the pulled genes got pulled with this cut value
              if(!cut.name %in% data$genes$pulled.genes[[gene]]){
                data$genes$pulled.genes[[gene]] <- c(data$genes$pulled.genes[[gene]], cut.name)
                cut.genes <- c(cut.genes, gene)
              }
            } else {
              data$genes$pulled.genes[[gene]] <- cut.name
              zscore.genes <- c(zscore.genes, gene)
              cut.genes <- c(cut.genes, gene)
            }
          }
          # then checking if there are genes selected actually
          if(!(length(zscore.genes) == 0 & length(cut.genes) == 0)){
            the.expression <- load_data_set(
              load.expression = TRUE
              , session = session
              , path.to.data = data$dataset.paths[data$available.datasets == data$general$dataset.name]
              , data = data
              , verbose = statics$verbose
              , subset.data = FALSE
            )
            # if there is a problem with the dataset, load_data_set will return NULL
            if(!is.null(the.expression)){
              if(statics$verbose){
                cat(file = stderr(),"\ngetting the genes")
              }
              # in case a quantil isn't well defined, we will need the following parameters
              # but since we checked currently no quantil, every quantil is well definded
              quantil.error <- FALSE
              data$genes$quantil.error <- list("lower.quantils" = c(),
                                               "upper.quantils" = c(),
                                               "genes" = c(),
                                               "expr.vecs" = matrix(nrow=0, ncol=0))

              # converting HGNC symbols to ENGS to convert the count matrix from ENGS to HGNC
              pulled.ens <- as.vector(HGNC_ENSGID$ensembl_gene_id[HGNC_ENSGID$hgnc_symbol %in% the.genes])
              pulled.ens <- pulled.ens[pulled.ens %in% rownames(the.expression)]
              pulled.counts <- rownames_mapper_wrapper(
                count.matrix = the.expression[pulled.ens, , drop = FALSE]
                , input.feature.type = "ensembl_gene_id"
                , verbose = statics$verbose
                , output.feature.type = "hgnc_symbol")
              pulled.counts <- pulled.counts$count.matrix
              rm(the.expression)
              for(gene in the.genes){
                expr.vec <- pulled.counts[gene, ]
                gene <- gsub(pattern = "-", replacement = ".", gene)

                # we transfer the gene expression to zscore.
                if(gene %in% zscore.genes){
                  # if a gene is constant, zscore can not be computed
                  if(all(expr.vec == mean(expr.vec))){
                    data$general$data.set[[paste0("zscore.", gene)]] <- expr.vec
                  } else {
                    # scaling data
                    data$general$data.set[[paste0("zscore.", gene)]] <- (expr.vec- mean(expr.vec))/(sd(expr.vec)+1)
                  }
                  # adding the raw expression as well
                  data$general$data.set[[paste0("raw.", gene)]] <- expr.vec
                }
                if(gene %in% cut.genes){
                  # if we cut by the mean
                  if(cut == "mean"){
                    # inserting new column in the pheno, holding the info
                    # wether this gene of the cell is lower or higher
                    # than the mean
                    tmp.mean <- mean(x = expr.vec)
                    data$general$data.set[[paste0("above.mean.", gene)]] <- ifelse(
                      test = expr.vec > tmp.mean
                      , yes = "> mean"
                      , no = "<= mean"
                    )
                    pull.genes.helper(data = data, statics = statics,
                                      gene = gene,
                                      text = "above.mean.")
                  } else if(b.perc){
                    if(statics$verbose){
                      cat(file = stderr(),"\ncut by quantil")
                    }
                    # if we cut by a certain quantil
                    # sorting the expression
                    sorted.expr.vec <- expr.vec[order(expr.vec)]
                    quantil <- round(length(sorted.expr.vec)*cut/100)
                    # if the quantil isn't well defined
                    if(length(which(sorted.expr.vec == sorted.expr.vec[quantil]))>1){
                      quantil.error <- TRUE
                      if(statics$verbose){
                        cat(file = stderr(),"\nquantil not well defined")
                      }
                      # getting the next upper/lower quantil
                      # the next lower quantil is the maximum (highest in the order in case
                      # this quantil is also not well defined) of all values which are lower
                      # than the value of the original quantil
                      # Since a quantil is defined trough all values which are <= than the quantil value
                      # the original quantil, the "real" original quantil (aka the upper quantil) is just
                      # the maxmimum of all values which are lower or equal as the original value
                      # similar for the higher quantil
                      lower.quantil <- max(which(sorted.expr.vec < sorted.expr.vec[quantil]))
                      upper.quantil <- max(which(sorted.expr.vec <= sorted.expr.vec[quantil]))
                      data$genes$quantil.error$lower.quantils <- c(data$genes$quantil.error$lower.quantils,
                                                                   lower.quantil)
                      data$genes$quantil.error$upper.quantils <- c(data$genes$quantil.error$upper.quantils,
                                                                   upper.quantil)
                      data$genes$quantil.error$genes <- c(data$genes$quantil.error$genes,
                                                          gene)
                      data$genes$quantil.error$expr.vecs <- combine_counts(
                        data$genes$quantil.error$expr.vecs,
                        matrix(data = expr.vec, ncol = 1, dimnames = list(names(expr.vec), gene)))
                    } else {
                      # getting the value of the quantil
                      quantil <- sorted.expr.vec[quantil]
                      data$general$data.set[[paste0(cut, "%_quantil.", gene)]] <- ifelse(
                        test = expr.vec > quantil
                        , yes = paste0("> ", cut, "%_quantil.")
                        , no = paste0("<= ", cut, "%_quantil.")
                      )
                      pull.genes.helper(data = data, statics = statics,
                                        gene = gene,
                                        text = paste0(cut, "%_quantil."))
                    }
                  } else {
                    # if the cut is a expr. value
                    data$general$data.set[[paste0(cut, "_expr.", gene)]] <- ifelse(
                      test = expr.vec > cut
                      , yes = paste0("> expr. ", cut)
                      , no = paste0("<= expr. ", cut)
                    )
                    pull.genes.helper(data = data, statics = statics,
                                      gene = gene,
                                      text = paste0(cut, "_expr."))
                  }
                }
              }
              if(quantil.error){
                # if a quantil was not well defined I ask the user to select a well defined
                # quantil
                # getting the informations for this quantil
                next_quantil_error(
                  data = data
                  , cut = cut
                  , session
                )
              } else {
                # updating the text for the pulled gene
                data$output$texts$pulled.genes <- paste0(
                  "Pulled the gene(s) '", paste0(input$pullable.genes, collapse = ", ")
                  , "'"
                )
                if(statics$verbose){
                  cat(file = stderr(),"\nupdating the Select input")
                }
                # updating the selectInputs for the plots
                update_plot_options(session = session,
                                    data = data,
                                    statics = statics,
                                    input = input)

                update_sliders(
                  data = data, statics = statics
                  , input = input
                  , session = session
                )
                update_selectizes(
                  data = data, statics = statics
                  , input = input
                  , session = session
                )
                end_loading_screen(session = session)
              }
              rm(pulled.counts)
              rm(expr.vec)

            } else {
              if(statics$verbose){
                cat(file = stderr(),"\nno genes selected => doing nothing")
              }
              data$output$texts$pulled.genes <- "please select a viable gene to pull (which didn't got pulled in this setting)"
            }
          }

        } else {
          if(statics$verbose){
            cat(file = stderr(),"\nwrong quantil => doing nothing")
          }
          data$output$texts$pulled.genes <- "the quantil needs to be between 0 and 100"
        }
      }

      if(statics$verbose) cat(file = stderr(),"\nDONE with pull.genes\n")
    }
  )

  observeEvent(
    # if quantils by pull genes aren't well definded Celloscope will
    # use shinyalerts to ask which next well defined quantil the user wants to use
    input$quantil.error,{
      if(statics$verbose){
        cat(file = stderr(),"\n\nin quantil.error")
      }
      if(input$quantil.error){
        quantil <- data$genes$quantil.error$upper.quantils[1]
      } else {
        quantil <- data$genes$quantil.error$lower.quantils[1]
      }
      if(is.finite(quantil)){
        gene <- data$genes$quantil.error$genes[1]
        expr.vec <- data$genes$quantil.error$expr.vecs[, 1]
        sorted.expr.vec <- expr.vec[order(expr.vec)]
        cut <- round(quantil/length(sorted.expr.vec)*100, digits =3)
        quantil <- sorted.expr.vec[quantil]
        data$general$data.set[[paste0(cut, "%_quantil.", gene)]] <- ifelse(
          test = expr.vec > quantil
          , yes = paste0("> ", cut, "%_quantil.")
          , no = paste0("<= ", cut, "%_quantil.")
        )
        pull.genes.helper(data = data, statics = statics,
                          gene = gene,
                          text = paste0(cut, "%_quantil."))
        data$genes$quantil.error$lower.quantils <- data$genes$quantil.error$lower.quantils[-1]
        data$genes$quantil.error$upper.quantils <- data$genes$quantil.error$upper.quantils[-1]
        data$genes$quantil.error$genes <- data$genes$quantil.error$genes[-1]
        data$genes$quantil.error$expr.vecs <- data$genes$quantil.error$expr.vecs[, -1, drop = FALSE]
        if(statics$verbose){
          cat(file = stderr(),"\nreducing data$genes$quantil.error")
        }
      }
      if(length(data$genes$quantil.error$genes) == 0){
        # if this was the last not defined quantil
        if(statics$verbose){
          cat(file = stderr(),"\nupdating text")
        }
        # updating the text for the pulled gene
        data$output$texts$pulled.genes <- paste0(
          "Pulled the gene(s) '", paste0(input$pullable.genes, collapse = ", ")
          , "'"
        )
        if(statics$verbose){
          cat(file = stderr(),"\nupdating the Select input")
        }
        # updating the selectInputs for the plots
        update_plot_options(
          session = session
          , data = data
          , statics = statics
          , input = input)

        update_sliders(
          data = data, statics = statics
          , input = input
          , session = session
        )
        update_selectizes(
          data = data, statics = statics
          , input = input
          , session = session
        )
      } else {
        # if there are again some not defined qantils
        if(statics$verbose){
          cat(file = stderr(),"\npreparing next message")
        }
        next_quantil_error(
          data = data
          , cut = cut
          , session
        )
      }
    })


  observeEvent(
    #some values got selected trough the first slider
    input$first.slider,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin first.slider")
      }
      variable.name <- input$first.slider.col
      #if the slider is activated (if we have enough pheno entries)
      #we update the slider.values
      if(variable.name != "None"){
        data$slider_selectizes$numerics$slider.values[[variable.name]] <- input$first.slider
      }
      #update the Text output
      update_subsetting_output(
        data = data
        , statics = statics
        , input = input
      )
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with first slider\n")
      }
    }
  )

  observeEvent(
    #some values got selected trough the second slider
    input$second.slider,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin second.slider")
      }
      variable.name <- input$second.slider.col
      if(variable.name != "None"){
        data$slider_selectizes$numerics$slider.values[[variable.name]] <- input$second.slider
      }
      update_subsetting_output(
        data = data
        , statics = statics
        , input = input)
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with second slider\n")
      }
    }
  )

  observeEvent(
    #some values got selected trough the third slider
    input$third.slider,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin third slider")
      }
      variable.name <- input$third.slider.col
      if(variable.name != "None"){
        data$slider_selectizes$numerics$slider.values[[variable.name]] <- input$third.slider
      }
      update_subsetting_output(
        data = data
        , statics = statics
        , input = input)
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with third slider\n")
      }
    }
  )

  observeEvent(
    #some values got selected trough the fourth slider
    input$fourth.slider,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin fourth slider")
      }
      variable.name <- input$fourth.slider.col
      if(variable.name != "None"){
        data$slider_selectizes$numerics$slider.values[[variable.name]] <- input$fourth.slider
      }
      update_subsetting_output(
        data = data
        , statics = statics
        , input = input
        )
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with fourth slider\n")
      }
    }
  )

  observeEvent(
    #unselected everything from the first selectize
    input$clear.first.selectize,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\ninclear.first.selectize")
      }
      if(data$loaded){
        pheno.data <- data$general$data.set
        updateSelectizeInput(
          inputId = "first.selectize",
          label = input$first.selectize.col,
          choices = unique(pheno.data[, input$first.selectize.col]),
          session = session,
          server = TRUE,
          selected = NULL
        )
        #since we have now no entries selected via the selectize; we have to
        #put all entries on TRUE if the in.ex.first.selectize is on FALSE
        # (no entries excluded) or on FALSE if the in.ex.first.selectize is on
        # TRUE (no entries included)
        data$slider_selectizes$booleans$selected[, input$first.selectize.col] <- !input$in.ex.first.selectize
        update_subsetting_output(
          data = data, statics = statics,
          input = input)
      }
      if(statics$verbose){
        cat(file = stderr(),"\n\nDONE with clear.first.selectize")
      }
    }
  )
  observeEvent(
    # draged and dropped a private dataset
    input$dragdropupload,{
      if(statics$verbose){
        cat(file = stderr(),"\n\nin dragdropupload")
      }
      start_loading_screen(session = session)
      
      file.end.pos <- regexpr("\\.([[:alnum:]]+)$", input$dragdropupload$name)
      
      if(file.end.pos > -1){
        file.type <- substring(input$dragdropupload$name, file.end.pos + 1L)
      } else {
        file.type <- ""
      }
      
      if(file.type != "h5"){
        end_loading_screen(session = session)
        warning.message <- paste0(
          "celloscope detected a '", file.type, "'. Currently, only 'h5' ", 
          "files are supported"
        )
        shinyalert(
          title = "Unsupported file type"
          , text = warning.message
          , type = "warning"
          , closeOnEsc = TRUE
          , closeOnClickOutside = TRUE
          , showCancelButton = TRUE
          , animation = FALSE
          , size = "m"
          , immediate = TRUE
          , session = session
          , showConfirmButton = FALSE
        )
      } else {
        new.dataset <- gsub(
          pattern = ".h5"
          , replacement = ""
          , x = input$dragdropupload$name
        )

        file.copy(
          from = input$dragdropupload$datapath
          , to = paste0(statics$path.to.data, "/", input$dragdropupload$name)
        )

        data$available.datasets <- c(
          data$available.datasets
          , new.dataset
        )
        data$dataset.paths <- c(
          data$dataset.paths
          , paste0(statics$path.to.data, "/", input$dragdropupload$name)
        )

        all.files <- read_datasetlist_csv(path = statics$path.to.data)

        all.files <- rbind(
          all.files,
          c(digest::digest(
            object = input$cookies$`session-number`
            , algo = "sha256"),
            input$dragdropupload$name,
            format(Sys.time(), "%Y-%m-%d, %H:%M")
          )
        )

        write_datasetlist_csv(
          path = statics$path.to.data
          , the.files = all.files
        )


        all.datasets <- list_datasets(
          path = statics$path.to.data
          , owner.cookie = input$cookies$`session-number`
        )

        data$dataset.paths <- names(all.datasets)

        names(all.datasets) <- NULL

        data$available.datasets <- all.datasets

        updateSelectInput(
          session = session
          , inputId = "dataset"
          , choices = data$available.datasets
          , selected = new.dataset
        )

        end_loading_screen(session = session)
        if(statics$verbose){
          cat(file = stderr(),"\n\nDONE with dragdropupload")
        }
      }
    })
  observeEvent(
    #unselected everythiFng from the second selectize
    input$clear_second_selectize,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\ninclear.second.selectize")
      }
      pheno.data <- data$general$data.set
      updateSelectizeInput(
        inputId = "second.selectize",
        label = input$second.selectize.col,
        choices = unique(pheno.data[, input$second.selectize.col]),
        session = session,
        server = TRUE,
        selected = NULL
      )
      data$slider_selectizes$booleans$selected[, input$second.selectize.col] <- !input$in.ex.second.selectize
      update_subsetting_output(
        data = data
        , statics = statics
        , input = input)
      if(statics$verbose){
        cat(file = stderr(),"\n\nDONE with clear.second.selectize")
      }
    }
  )

  observeEvent(
    #some entries got selected trough the first selectize
    input$first.selectize,
    {
      pheno.data <- data$general$data.set
      if(statics$verbose){
        cat(file = stderr(),"\n\nin first.selectize")
      }
      #all rows which are selected via this selectize
      selected.ids <- c()
      for(selected in input$first.selectize){
        selected.ids <- c(
          selected.ids
          , which(pheno.data[, input$first.selectize.col] == selected)
        )
      }
      data$slider_selectizes$booleans$selected[selected.ids, input$first.selectize.col] <- input$in.ex.first.selectize
      #if selected ids is null the next line would through an error
      if(!is.null(selected.ids)){
        data$slider_selectizes$booleans$selected[-selected.ids, input$first.selectize.col] <- !input$in.ex.first.selectize
      } else {
        data$slider_selectizes$booleans$selected[, input$first.selectize.col] <- !input$in.ex.first.selectize
      }
      update_subsetting_output(
        data = data
        , statics = statics
        , input = input)
      if(statics$verbose){
        cat(file = stderr(),"\n\nDONE with first.selectize")
      }
    }

  )

  observeEvent(
    # some entries got selected trough the second selectize
    input$second.selectize,
    {
      pheno.data <- data$general$data.set
      if(statics$verbose){
        cat(file = stderr(),"\n\nin second.selectize")
      }
      selected.ids <- c()
      for(selected in input$second.selectize){
        selected.ids <- c(
          selected.ids
          , which(pheno.data[, input$second.selectize.col] == selected)
        )
      }
      data$slider_selectizes$booleans$selected[selected.ids, input$second.selectize.col] <- input$in.ex.second.selectize
      if(!is.null(selected.ids)){
        data$slider_selectizes$booleans$selected[-selected.ids, input$second.selectize.col] <- !input$in.ex.second.selectize
      }
      else {
        data$slider_selectizes$booleans$selected[, input$second.selectize.col] <- !input$in.ex.second.selectize
      }
      update_subsetting_output(
        data = data
        , statics = statics
        , input = input)
      if(statics$verbose){
        cat(file = stderr(),"\n\nDONE with second selectize")
      }
    }

  )

  observeEvent(
    #if the switch from the first selectize got switched
    input$in.ex.first.selectize,
    {

      if(statics$verbose){
        cat(file = stderr(),"\n\nin in.ex.first.selectize")
      }
      if(data$loaded){
        # if this function gets activated but there is nothing to
        # do (input$in.ex.first.selectize == data$slider_selectizes$booleans$selectize.switch[1])
        # we won't do something
        # if we don't have categorical entries in our pheno, input$first.selectize.col == NULL
        # and therefore the line
        # data$slider_selectizes$booleans$selected[, input$first.selectize.col] <- !data$slider_selectizes$booleans$selected[, input$first.selectize.col]
        # would result in an error
        if(data$slider_selectizes$numerics$n.categorical > 0 & input$in.ex.first.selectize != data$slider_selectizes$booleans$selectize.switch[1]){
          if(statics$verbose){
            cat(file = stderr(),"\nchanging selection")
          }
          #just reverse the selection
          data$slider_selectizes$booleans$selected[, input$first.selectize.col] <- !data$slider_selectizes$booleans$selected[, input$first.selectize.col]
          data$slider_selectizes$booleans$selectize.switch[1] <- input$in.ex.first.selectize
        }
        update_subsetting_output(
          data = data
          , statics = statics
          , input = input)
      }
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with in.ex.first.selectize")
      }
    })

  observeEvent(
    #if the switch from the second selectize got switched
    input$in.ex.second.selectize,
    {

      if(statics$verbose){
        cat(file = stderr(),"\n\nnin in.ex.second.selectize")
      }
      if(data$loaded){
        if(data$slider_selectizes$numerics$n.categorical > 0 & input$in.ex.second.selectize != data$slider_selectizes$booleans$selectize.switch[2]){
          data$slider_selectizes$booleans$selected[, input$second.selectize.col] <- !data$slider_selectizes$booleans$selected[, input$second.selectize.col]
          data$slider_selectizes$booleans$selectize.switch[2] <- input$in.ex.second.selectize
        }
      }
      update_subsetting_output(
        data = data
        , statics = statics
        , input = input)
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with in.ex.second.selectize")
      }
    })

  observeEvent(
    {
      #if the category of one slider gets changed
      input$first.slider.col
      input$second.slider.col
      input$third.slider.col
      input$fourth.slider.col
    },
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin slider.cols")
      }
      reset_brushes(session = session)
      if(data$loaded){
        update_sliders(
          data = data
          , statics = statics
          , session = session
          , input = input
        )
        plot_histograms(
          data = data
          , statics = statics
          , input = input
          , session = session
        )
      }
      output$hist.first.slider <- renderPlot(({data$output$plots$slider.histograms$first}))
      output$hist.second.slider <- renderPlot(({data$output$plots$slider.histograms$second}))
      output$hist.third.slider <- renderPlot(({data$output$plots$slider.histograms$third}))
      output$hist.fourth.slider <- renderPlot(({data$output$plots$slider.histograms$fourth}))
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with input.sliders")
      }
    }
  )

  observeEvent(
    #if the category of the first selectize gets changes
    input$first.selectize.col,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin first.selectize.col")
      }
      if(data$loaded){
        #this part gets triggered if one changes the column (of the pheno matrix) for the first selectize
        #We want to have all entries selected at the start (if we first choose this column), but if
        #a category is very large we don't want them written down in the selectize, because it would
        #look kind of overloaded and crowded therefore we start with excluding all
        #selected entries (first selectize switch is FALSE) and don't select any
        data$slider_selectizes$booleans$selectize.switch[1] <- FALSE
        update_selectizes(
          data = data, statics = statics,
          session = session,
          input = input
        )
      }
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with input.sliders")
      }
    }

  )

  observeEvent(
    #if the category of the second selectize gets changes
    input$second.selectize.col,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin second.selectize.col")
      }
      if(data$loaded){
        data$slider_selectizes$booleans$selectize.switch[2] <- FALSE
        update_selectizes(
          data = data
          , statics = statics
          , session = session
          , input = input
        )
      }
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with input.sliders")
      }
    }

  )

  observeEvent(
    #to reload the brush on the first scatterplot
    input$first.reload.brush,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin first.reload.brush")
      }
      plot_brush_histograms(
        data = data
        , statics = statics
        , input = input
        , session = session
        , id = "first"
      )
      reset_brushes(session = session)
      output$first.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$first.selected}))
      output$first.not.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$first.not.selected}))
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with first.reload.brush\n")
      }
    }
  )
  observeEvent(
    # to reload the brush on the second scatterplot
    input$second.reload.brush,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin first.reload.brush")
      }
      #this is the brush in the first big plot
      plot_brush_histograms(
        data = data
        , statics = statics
        , input = input
        , session = session
        , id = "second"
      )
      reset_brushes(session = session)
      output$second.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$second.selected}))
      output$second.not.selected.hist.plot <- renderPlot(({data$output$plots$brush.histograms$second.not.selected}))
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with first.reload.brush\n")
      }
    }
  )

  observeEvent(
    # to save the parameters of the current selection
    input$save.parameters, {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin save.paramters")
      }
      start_loading_screen(session = session)
      #converting the slider.values-list to a matrix so we can write it as a csv
      slider.values <- data$slider_selectizes$numerics$slider.values
      slider.values.m <- matrix(ncol=length(slider.values), nrow = 2)
      colnames(slider.values.m) <- names(slider.values)
      for(i in 1:length(slider.values)){
        slider.values.m[, i] <- slider.values[[i]]
      }
      write.csv(
        x = slider.values.m
        , file = paste0(
          statics$path.to.data
          , "parameters/"
          , input$file.prefix.parameters
          , "_slider_values__"
          , data$general$dataset.name
          , ".csv"
        )
      )
      write.csv(
        x = data$slider_selectizes$booleans$selected
        , file = paste0(
          statics$path.to.data
          , "parameters/"
          , input$file.prefix.parameters
          , "_selectize__"
          , data$general$dataset.name
          , ".csv"
        )
      )
      data$output$texts$saved.parameters <- paste0(
        "saved parameters in: "
        , statics$path.to.data
        , input$file.prefix.parameters
      )
      if(!dir.exists(paste0(statics$path.to.data, "/parameters/"))){
        dir.create(
          paste0(
            statics$path.to.data
            , "/parameters/"
          )
        )
      }
      available.parameters <- c(
        "None"
        , gsub(pattern = paste0(
          "_selectize__"
          , data$general$dataset.name
          , ".csv"
        )
        , replacement = ""
        , x = list.files(
          path = paste0(
            statics$path.to.data
            , "/parameters/"
          )
          , pattern = paste0(
            "_selectize__"
            , data$general$dataset.name
          )
        )
        )
      )
      updateSelectInput(
        session = session
        , inputId = "parameters"
        , choices =  available.parameters
        , selected = available.parameters[1]
        , label = "Load saved parameters"
      )
      reset_brushes(session = session)
      end_loading_screen(session = session)
      if(statics$verbose){
        cat(file = stderr(),"\n\nDONE with save.paramters")
      }
    }
  )

  observeEvent(
    # write the current sys.time as a prefix for saving the parameters
    input$get.Sys.time.parameters, {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin get.Sys.time.parameters")
      }
      # setting the file prefix to the Sys.time
      # this just converts '2021-10-20 15:48:51 CEST'
      # to '2021_10_20__15_48'
      # since we don't want to have any blanks, '-', ':',
      # or the timezone in our filename
      updateTextInput(
        session = session
        , inputId = "file.prefix.parameters"
        , value = gsub(
          pattern = " "
          , replacement = "__"
          , x = gsub(pattern = "-|:"
                     , x = format(
                       Sys.time()
                       , "%Y-%m-%d %H:%M"
                     )
                     , replacement = "_"
          )
        )
      )
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with get.Sys.time.parameters\n")
      }
    }
  )
  observeEvent(
    # write the current sys.time as a prefix for saving the plots
    input$get.Sys.time.plots, {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin get.Sys.time.plots")
      }
      # setting the file prefix to the Sys.time
      # this just converts '2021-10-20 15:48:51 CEST'
      # to '2021_10_20__15_48'
      # since we don't want to have any blanks, '-', ':',
      # or the timezone in our filename
      updateTextInput(
        session = session
        , inputId = "file.prefix.plots"
        , value = paste0(
          "plots_",gsub(
            pattern = " "
            , replacement = "__"
            , x = gsub(
              pattern = "-|:"
              , x = format(
                Sys.time()
                , "%Y-%m-%d %H:%M"
              )
              , replacement = "_"
            )
          )
        )
      )
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with get.Sys.time.plots\n")
      }
    }
  )
  observeEvent(
    input$save.plots, {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin save.plots")
      }
      start_loading_screen(session = session)
      file.name <- paste0(statics$path.to.data, "plots/", input$file.prefix.plots)
      if(!grepl(".pdf$", file.name)){
        file.name <- paste0(file.name, ".pdf")
      }
      pdf(file.name)
      for(id in c("first", "second", "third", "fourth")){
        if(input[[paste0(id, ".slider.col")]] %in% input$hist.save) print(data$output$plots$slider.histograms[[id]])
      }
      if("first scatter plot" %in% input$scatter.save) print(data$output$plots$scatter$first)
      if("second scatter plot" %in% input$scatter.save) print(data$output$plots$scatter$second)
      if("first brush histogram" %in% input$brush.save){
        print(data$output$plots$brush.histograms$first.selected)
        print(data$output$plots$brush.histograms$first.not.selected)
      }
      if("second brush histogram" %in% input$brush.save){
        print(data$output$plots$brush.histograms$second.selected)
        print(data$output$plots$brush.histograms$second.not.selected)
      }
      dev.off()
      end_loading_screen(session = session)
      if(statics$verbose){
        cat(file = stderr(),"\n\nDONE with save.plots")
      }
    }
  )
  observeEvent(
    # write the curren sys.time as a prefix for saving the current dataset
    input$get.Sys.time.dataset, {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin get.Sys.time.dataset")
      }
      #setting the file prefix to the Sys.time
      updateTextInput(
        session = session
        , inputId = "file.prefix.dataset"
        , value = gsub(
          pattern = " "
          , replacement = "__"
          , x = gsub(
            pattern = "-|:"
            , x = format(
              Sys.time()
              , "%Y-%m-%d %H:%M"
            )
            , replacement = "_"
          )
        )
      )
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with get.Sys.time.dataset\n")
      }
    }
  )

  observeEvent(
    # to load saved parameters
    input$parameters, {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin paramters")
      }
      start_loading_screen(session = session)
      if(input$parameters != "None"){
        reset_brushes(session = session)
        slider.file <- paste0(
          statics$path.to.data
          , "parameters/"
          , input$parameters
          , "_slider_values__"
          , data$general$dataset.name
          , ".csv"
        )
        # since read.csv inserts and additional column "X" at the beginning we remove it
        slider.values.m <- read.csv(slider.file)[,-1]
        # converting the matrix back to our list
        for(i in 1:ncol(slider.values.m)){
          data$slider_selectizes$numerics$slider.values[[colnames(slider.values.m)[i]]] <- slider.values.m[, i]
        }
        # loading the selectize
        selectize.file <- paste0(
          statics$path.to.data, "parameters/"
          , input$parameters
          , "_selectize__"
          , data$general$dataset.name
          , ".csv"
        )
        selectizes <- data.frame(read.csv(selectize.file))
        rownames(selectizes) <- selectizes$X
        selectizes$X <- NULL
        selectizes <- data.frame(selectizes == "TRUE")
        data$slider_selectizes$booleans$selected <- selectizes
        # updating everthing
        update_sliders(
          data = data
          , statics = statics
          , session = session
          , input = input
        )
        update_selectizes(
          data = data
          , statics = statics
          , session = session
          , input = input
        )
        update_subsetting_output(
          data = data
          , statics = statics
          , input = input
        )
        plot_histograms(
          data = data
          , statics = statics
          , input = input
          , session = session
        )
      }
      output$hist.first.slider <- renderPlot(({data$output$plots$slider.histograms$first}))
      output$hist.second.slider <- renderPlot(({data$output$plots$slider.histograms$second}))
      output$hist.third.slider <- renderPlot(({data$output$plots$slider.histograms$third}))
      output$hist.fourth.slider <- renderPlot(({data$output$plots$slider.histograms$fourth}))
      end_loading_screen(session = session)
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with paramets\n")
      }
    }
  )


  observeEvent(
    #calculating pathways
    eventExpr = input$calculate.pathways
    , handlerExpr = {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin calculate.pathways\ngetting the seleceted genes")
      }
      start_loading_screen(session = session)
      the.genes.selectize <- input$pullable.genes.pathway
      the.genes.text <- gsub(
        pattern = " "
        , replacement = ""
        , x = unlist(
          strsplit(
            x = input$genes.pathway.text
            , split = ","
          )
        )
      )

      the.genes <- unique(c(the.genes.selectize, the.genes.text))
      the.expression <- load_data_set(
        load.expression = TRUE
        , path.to.data = data$dataset.paths[data$available.datasets == data$general$dataset.name]
        , verbose = statics$verbose
        , data = data
        , session = session
        , subset.data = FALSE
      )
      if(!is.null(the.expression)){
        if(statics$verbose){
          cat(file = stderr(),"\nconverting to hgnc")
        }
        mapper <- HGNC_ENSGID
        mapper <- mapper[, c(
          which(colnames(mapper)=="hgnc_symbol")
          ,which(colnames(mapper)=="ensembl_gene_id"))]
        pulled.ens <- mapper[mapper$hgnc_symbol %in% the.genes, ]
        rownames(pulled.ens) <- 1:nrow(pulled.ens)
        # try.res <- try(pulled.ens <- map_feature_names(the.genes))
        # while(class(try.res)=="try-error"){
        #   try.res <- try(pulled.ens <- map_feature_names(the.genes))
        # }
        useable.genes <- intersect(rownames(the.expression), pulled.ens$ensembl_gene_id)
        if(statics$verbose){
          cat(file = stderr(),"\nmaking pathway name usable")
        }
        useable.name <- gsub(" ", "", gsub("-", ".", input$pathway.name))
        if(statics$verbose){
          cat(file = stderr(),"\ncalculating the scores")
        }
        the.expression <- the.expression[useable.genes,]
        add.frame <- calc.pathway(
          the.counts = the.expression
          , a.set.of.genes = useable.genes
          , verbose = statics$verbose
          , data = data
        )
        if(class(add.frame) == "logical"){
          if(add.frame){
            shinyalert(
              title = "To less genes"
              , text = "Your pathway contains to less genes; Please select more genes"
              , type = "warning"
              , closeOnEsc = TRUE
              , closeOnClickOutside = TRUE
              , showCancelButton = TRUE
              , animation = FALSE
              , size = "m"
              , immediate = TRUE
              , session = session
              , showConfirmButton = FALSE
            )
          }
        } else {
          add.frame <- add.frame[data$general$data.set$colnames, ]
          data$general$data.set <- cbind(data$general$data.set, add.frame)
          if(any(useable.name == gsub(
            x = colnames(data$general$data.set)
            , pattern = "cluster_|mean_z_score_"
            , replacement = ""
          ))){
            n.names <- length(which(grepl(
              pattern = paste0("cluster_", useable.name)
              , x = colnames(data$general$data.set))))
            useable.name <- paste0(useable.name, n.names)
          }
          colnames(data$general$data.set)[
            (ncol(data$general$data.set)-1):ncol(data$general$data.set)
          ] <- c(paste0("cluster_", useable.name), paste0("mean_z_score_", useable.name))



          if(length(useable.genes) > 0){
            if(statics$verbose){
              cat(file = stderr(),"\nupdating is.slider")
            }
            pheno.data <- data$general$data.set
            is.slider <- data$slider_selectizes$booleans$is.numeric
            is.slider <- c(is.slider,FALSE,TRUE)
            names(is.slider)[(length(is.slider)-1):length(is.slider)] <- c(
              paste0("cluster_", useable.name)
              , paste0("mean_z_score_", useable.name)
            )
            data$slider_selectizes$booleans$is.numeric <- is.slider
            if(statics$verbose){
              cat(file = stderr(),"\nupdating ranges and slider.values")
            }

            # update ranges and slider values
            ranges <- data$slider_selectizes$numerics$ranges
            slider.values <- data$slider_selectizes$numerics$slider.values
            if(statics$verbose){
              cat(file = stderr(),"\nexecuting min_max_step")
            }
            ranges$tmp <- min_max_step(pheno.data[, paste0("mean_z_score_", useable.name)])
            slider.values$tmp <- ranges$tmp[1:2]
            names(ranges)[length(ranges)] <- paste0("mean_z_score_", useable.name)
            names(slider.values)[length(slider.values)] <- paste0("mean_z_score_", useable.name)
            data$slider_selectizes$numerics$ranges <- ranges
            data$slider_selectizes$numerics$slider.values <- slider.values

            selected <- data$slider_selectizes$booleans$selected
            selected <- cbind(selected, rep(x = TRUE, times = nrow(selected)))
            colnames(selected)[ncol(selected)] <- paste0("cluster_", useable.name)
            data$slider_selectizes$booleans$selected <- selected
            update_plot_options(
              session = session
              , data = data
              , statics = statics
              , input = input
            )

            update_sliders(
              data = data
              , statics = statics
              , input = input
              , session = session
            )
            update_selectizes(
              data = data
              , statics = statics
              , input = input
              , session = session
            )

            if(statics$verbose){
              cat(file = stderr(),"\nupdating SelectInput")
            }
            updateSelectInput(
              session = session
              , inputId = "first.colouring"
              , choices = colnames(data$general$data.set)
              , selected = input$first.colouring
            )

            updateSelectInput(
              session = session
              , inputId = "second.colouring"
              , choices = colnames(data$general$data.set)
              , selected = input$second.colouring
            )
            updateSelectInput(
              session = session
              , inputId = "second.x.axis"
              , choices = colnames(data$general$data.set)
              , selected = input$second.x.axis
            )
            updateSelectInput(
              session = session
              , inputId = "first.x.axis"
              , choices = colnames(data$general$data.set)
              , selected = input$first.x.axis
            )
            updateSelectInput(
              session = session
              , inputId = "first.y.axis"
              , choices = colnames(data$general$data.set)
              , selected = input$first.y.axis
            )
            updateSelectInput(
              session = session
              , inputId = "second.y.axis"
              , choices = colnames(data$general$data.set)
              , selected = input$second.y.axis
            )
          }


          data$output$texts$calculated.pathway <- paste0(
            "For pathway '"
            , useable.name, "' ",
            round(digits = 2, length(useable.genes)/length(pulled.ens$ensembl_gene_id))*100
            , "% out of "
            , length(the.genes)
            , " genes are found in the rownames of the dataset"
          )
          end_loading_screen(session = session)

        }

        rm(the.expression)
      }

      if(statics$verbose){
        cat(file = stderr(),"\nDONE with calculate.pathways\n")
      }
    }
  )
  #the slider are getting updated trough the brushes
  observeEvent(
    #brush on the first histogram; then the correspoding slider needs to be updated
    input$first.hist.brush,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin first.hist.brush")
      }
      updateSliderInput(
        session = session
        , inputId = "first.slider"
        , value = c(input$first.hist.brush$xmin, input$first.hist.brush$xmax)
      )
      session$resetBrush("first.hist.brush")
      if(statics$verbose){
        cat(file = stderr(),"\nDone with first.hist.brush\n")
      }
    }
  )

  observeEvent(
    #brush on the second histogram; then the correspoding slider needs to be updated
    input$second.hist.brush,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin second.hist.brush")
      }
      updateSliderInput(
        session = session
        , inputId = "second.slider"
        , value = c(input$second.hist.brush$xmin, input$second.hist.brush$xmax)
      )
      session$resetBrush("second.hist.brush")
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with second.hist.brush\n")
      }
    }
  )

  observeEvent(
    #brush on the third histogram; then the correspoding slider needs to be updated
    input$third.hist.brush,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin third.hist.brush")
      }
      updateSliderInput(
        session = session
        , inputId = "third.slider"
        , value = c(input$third.hist.brush$xmin, input$third.hist.brush$xmax)
      )
      session$resetBrush("third.hist.brush")
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with third.brush\n")
      }
    }
  )

  observeEvent(
    #brush on the fourth histogram; then the correspoding slider needs to be updated
    input$fourth.hist.brush,
    {
      if(statics$verbose){
        cat(file = stderr(),"\n\nin fourth.hist.brush")
      }
      updateSliderInput(
        session = session
        , inputId = "fourth.slider"
        , value = c(input$fourth.hist.brush$xmin, input$fourth.hist.brush$xmax)
      )
      session$resetBrush("fourth.hist.brush")
      if(statics$verbose){
        cat(file = stderr(),"\nDONE with fourth.hist.brush\n")
      }
    }
  )


  observeEvent(
    input$first.x.axis,
    {
      if(data$loaded){
        # if the data is numerical
        if(all(is.numeric(data$general$data.set[,input$first.x.axis]))){
          data$plots$first.scatter.plot$x.lims <- data$slider_selectizes$numerics$ranges[[input$first.x.axis]][1:2]
        } else {
          data$plots$first.scatter.plot$x.lims <- unique(data$general$data.set[, input$first.x.axis])
        }
      }
    })
  observeEvent(
    input$first.y.axis,
    {
      if(data$loaded){
        if(all(is.numeric(data$general$data.set[,input$first.y.axis]))){
          data$plots$first.scatter.plot$y.lims <- data$slider_selectizes$numerics$ranges[[input$first.y.axis]][1:2]
        } else {
          data$plots$first.scatter.plot$y.lims <- unique(data$general$data.set[, input$first.y.axis])
        }
      }
    })
  observeEvent(
    input$second.x.axis,
    {
      if(data$loaded){
        if(all(is.numeric(data$general$data.set[,input$second.x.axis]))){
          data$plots$second.scatter.plot$x.lims <- data$slider_selectizes$numerics$ranges[[input$second.x.axis]][1:2]
        } else {
          data$plots$second.scatter.plot$x.lims <- unique(data$general$data.set[, input$second.x.axis])
        }
      }
    })
  observeEvent(
    input$second.y.axis,
    {
      if(data$loaded){
        if(all(is.numeric(data$general$data.set[,input$second.y.axis]))){
          data$plots$second.scatter.plot$y.lims <- data$slider_selectizes$numerics$ranges[[input$second.y.axis]][1:2]
        } else {
          data$plots$second.scatter.plot$y.lims <- unique(data$general$data.set[, input$second.y.axis])
        }
      }
    })
  observeEvent(
    input$first.zoom.in,
    {
      if(data$loaded){
        zoom_in(
          data = data, statics = statics
          , input = input
          , session = session
          , id = "first"
        )
      }
      plot_scatter_plot(
        data = data
        , statics = statics
        , input = input
        , id = "first"
        , session = session
      )

      output$first.scatter.plot <- renderPlot(({data$output$plots$scatter$first}))
    }
  )
  observeEvent(
    input$second.zoom.in,
    {
      if(data$loaded){
        zoom_in(
          data = data
          , statics = statics
          , input = input
          , session = session
          , id = "second")
        plot_scatter_plot(
          data = data
          , statics = statics
          , input = input
          , session = session
          , id = "second"
        )
        output$second.scatter.plot <- renderPlot(({data$output$plots$scatter$second}))
        reset_brushes(session = session)
      }
    }
  )
  observeEvent(
    input$first.zoom.out,
    {
      if(data$loaded){
        zoom_out(
          id = "first"
          , data = data
          , session = session
          , input = input
          , output = output
          , statics = statics
        )
      }
    }
  )
  observeEvent(
    input$second.zoom.out,
    {
      if(data$loaded){
        zoom_out(
          id = "second"
          , data = data
          , session = session
          , input = input
          , statics = statics
          , output = output
        )
      }
    }
  )
  #rendering the output
  output$first.brush.text <- renderText({paste0(data$output$texts$first.brush)})
  output$second.brush.text <- renderText({paste0(data$output$texts$second.brush)})
  output$pathway.output <- renderText(paste0(data$output$texts$calculated.pathway))
  output$pulled.a.gene <- renderText(paste0(data$output$texts$pulled.genes))
  output$saved.dataset <- renderText(paste0(data$output$texts$saved.dataset))
  output$diggeR.tags.text <- renderText(paste0(data$output$texts$diggeR$tags))
  output$diggeR.tags.text.ad <- renderText(paste0(data$output$texts$diggeR$tags.ad))
  output$saved.parameters <- renderText(paste0(data$output$texts$saved.parameters))
  output$currently.subsetting <- renderText(paste0(data$output$texts$currently.subsetting))
  }
)