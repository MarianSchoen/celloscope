#' read the dataset list csv file 
#'
#' @param path string, path to the file
#'
#' @return
#' @export
read_datasetlist_csv <- function(
  path = statics$path.to.data
) {
  owner.information <- read.csv(
    file = paste0(path, "datasetlist.csv")
    , header = TRUE
  )
  return(owner.information)
}

#' write the dataset list csv file 
#'
#' @param path string, path to the file
#'
#' @return
#' @export
write_datasetlist_csv <- function(
  path
  , the.files 
){
  write.csv(
    x = the.files
    , file = paste0(path, "datasetlist.csv")
    , row.names = FALSE
  )
}

#' list all available datasets 
#' 
#' @param path string to the data folder (usually "./data/shinyUseableData/")
#'
#' @return named vector with all available datasets
#' @export
#'
#' @examples
list_datasets <- function(
  path = statics$path.to.data
  , owner.cookie = NULL
){
  
  dataset.choices <- list.files(
    path = path
    , pattern = ".h5"
    , full.names = FALSE
  )
  
  available.datasets <- gsub(
    pattern = ".h5"
    , replacement = ""
    , x = dataset.choices
  )
  
  owner.information <- read_datasetlist_csv(path = path)
  
  useable.idx <- which(owner.information$owner.hash == "everybody")
  
  if(!is.null(owner.cookie)){
    owner.hash <- digest::digest(
      object = owner.cookie
      , algo = "sha256"
      )
    useable.idx <- c(
      useable.idx, 
      which(owner.information$owner.hash == owner.hash)
    )
  } 
    
  the.files <- gsub(
    " ", "", owner.information$file.name[useable.idx]
  )
  the.files <- gsub(
    ".h5", "", the.files 
  )
  useable.files <- intersect(
    available.datasets, the.files
  )
  names(useable.files) <- paste0(path, useable.files, ".h5")
  
  return(useable.files)
}

#' function, thats loads either the expression (load.expression = TRUE), the pheno-data 
#' (load.expression = FALSE) or the row names of the expression matrix 
#' (load.row.names.exp = TRUE & load.expression = FALSE)
#' if f.name is NULL, the function expects path.to.data to hold the
#' complete path of the data (inclusive the file name)
loading <- function(
  f.name = NULL
  , load.expression = FALSE
  , load.row.names.exp = FALSE
  , verbose = TRUE
  , path.to.data
  , data = NULL
  , subset.data = FALSE
){
  if(verbose){
    cat(file = stderr(), "\nin load_data_set")
  }
  if(is.null(f.name)){
    f.name <- path.to.data
  }
  if(!grepl(pattern = ".h5$", x = f.name)){
    f.name <- paste0(f.name, ".h5")
  }
  if(!(grepl(pattern = path.to.data, x = f.name) | f.name == path.to.data)){
    f.name <- paste0(path.to.data, f.name)
  }
  if(!file.exists(f.name)){
    stop( 
      paste0(
        "In 'load_data_set':", f.name, "' doesn't exist"
        )
    )
  } 
  #getting the content of the h5 file in order to look if the required part is there
  content <- rhdf5::h5ls(
    f.name
    , recursive = TRUE
  )
  if(load.expression){
    if(subset.data){
      # if the reducing of the dataset should also effect the counts for 
      # e.g. calculating a pathway; currently an deactivated feature (subset.data is always false)
      if(data$use.old.subset){
        sc.counts <- data$reduced.counts
      } else {
        
        #checking if sc.counts is there
        all.sc <- content[which(content$group == "/singlecell"),]
        if(! "data" %in% all.sc$name) stop(paste0(
          "In 'load_data_set': (dataset: '", f.name, "') the provided dataset has no sc.counts "
        ))
        sc.counts <- h5read(f.name, "singlecell/data")
        rownames(sc.counts) <- as.vector(h5read(f.name, "singlecell/geneids"))
        colnames(sc.counts) <- as.vector(h5read(f.name, "singlecell/cellids"))
        reduced.samples <- rownames(data$slider_selectizes$booleans$selected)[which(data$slider_selectizes$booleans$selected$in.reduced)]
        data$reduced.counts <-sc.counts[, reduced.samples]
        sc.counts <- sc.counts[, reduced.samples]
        data$use.old.subset <- TRUE
      }
    } else {
      
      #checking if sc.counts is there
      all.sc <- content[which(content$group == "/singlecell"),]
      if(! "data" %in% all.sc$name){
        stop(
          paste0(
            "In 'load_data_set': (dataset: '"
            , f.name
            , "') the provided dataset has no sc.counts "
            )
          )
      }
      sc.counts <- h5read(f.name, "singlecell/data")
      rownames(sc.counts) <- as.vector(h5read(f.name, "singlecell/geneids"))
      colnames(sc.counts) <- as.vector(h5read(f.name, "singlecell/cellids"))
    }
    if(verbose){
      cat(file = stderr(), "\nDONE with load_data_set")
    }
    if(is.numeric(sc.counts) & is.matrix(sc.counts)){
      return(sc.counts)
    } else {
      stop(
        "In 'load_data_set': sc.counts is no numeric matrix"
      )
    }
  }
  if(!load.expression){
    #checking if sc.pheno is there
    all.sc <- content[which(content$group == "/singlecell"),]
    if(! "pheno" %in% all.sc$name) stop(paste0(
      "In 'load_data_set': (dataset: '", f.name, "') the provided dataset has no sc.counts "
    ))
    sc.pheno <- data.frame()
    sc.pheno.idx <- which(content$group == "/singlecell/pheno")
    for(id in sc.pheno.idx){
      add.pheno <- data.frame(h5read(f.name, paste0("singlecell/pheno/", content$name[id])))
      colnames(add.pheno) <- content$name[id]
      if(ncol(sc.pheno>0)){
        sc.pheno <- cbind(sc.pheno, add.pheno)
      } else{
        sc.pheno <- add.pheno
      }
    }
    if(load.row.names.exp){
      #checking if the rownames of sc.counts are there
      all.sc <- content[which(content$group == "/singlecell"),]
      if(! "geneids" %in% all.sc$name){
        stop(
          paste0(
            "In 'load_data_set': (dataset: '"
            , f.name
            , "') the provided dataset has no sc.counts "
            )
          )
      }
      if(verbose){
        cat(file = stderr(), "\nDONE with load_data_set")
      }
      ret <- as.vector(h5read(f.name, "singlecell/geneids"))
      if(!is.null(ret)){
        return(ret)
      } else {
        stop(
          "In 'load_data_set': rownames of counts are NULL"
        )
      }
    }
    if(verbose){
      cat(file = stderr(), "\nDONE with load_data_set")
    }
    if(is.data.frame(sc.pheno)){
      return(sc.pheno)
    } else {
      stop(
        "In 'load_data_set': sc.pheno is no dataframe"
      )
    }
  }
}


# wrapper function for loading()
load_data_set <- function(
  f.name = NULL
  , load.expression = FALSE
  , load.row.names.exp = FALSE
  , path.to.data
  , data = NULL
  , subset.data = FALSE
  , session
  , verbose
){
  result <- try(loading(
    f.name = f.name
    , load.expression = load.expression
    , load.row.names.exp = load.row.names.exp
    , verbose = verbose
    , path.to.data = path.to.data
    , data = data
    , subset.data = subset.data
  ))
  if(class(result)=="try-error"){
    shinyalert(
      title = "Dataset Error"
      , text = "There is something wrong with your dataset. Please check it and try again"
      , closeOnEsc = TRUE
      , closeOnClickOutside = FALSE
      , session = session
      , showCancelButton = FALSE
      , showConfirmButton = TRUE
      , animation = FALSE
      , size = "m"
      , immediate = TRUE
    )
  } else {
    return(result)
  }
  return(NULL)
}

# this function calculates the average z-score per profile of 
# a given count marix only using the genes in gene.vec
avg_z_score <- function(
  gene.vec, 
  count.matrix, 
  log2.transform = TRUE
){
  
  if(max(as.vector(count.matrix)) < 20 && log2.transform){
    count.matrix <- 2**count.matrix - 1
  }
  
  useable.genes <- intersect(gene.vec, rownames(count.matrix))
  
  if(length(useable.genes) == 0){
    stop(
      "In 'avg_z_score': 'length(useable.genes)' is zero."
      )
  }
  sd.vec <- apply(
    X = count.matrix[useable.genes, , drop = FALSE]
    , MARGIN = 1
    , FUN = sd
  )
  mean.vec <- apply(
    X = count.matrix[useable.genes, , drop = FALSE]
    , MARGIN = 1
    , FUN = mean
  )
  
  # Attention: 'counts.z.score' is profiles x genes (other way around)
  counts.z.score <- sapply(
    X = useable.genes
    , FUN = function(gene){
      tmp <- count.matrix[gene, ] - mean.vec[gene]
      tmp <- tmp/(sd.vec[gene]+1)
      
      return(tmp)
    }
  )
  
  activation.per.profile <- rowSums(
    x = counts.z.score**2
  )
  normed.activation.per.profile <- sqrt(activation.per.profile/length(gene.vec))
  
  return(normed.activation.per.profile)
}

#function to zoom out
zoom_out <- function(id, data, session, input, statics, output){
  if(all(is.numeric(data$general$data.set[,input[[paste0(id , ".x.axis")]]]))){
    data$plots[[paste0(id , ".scatter.plot")]]$x.lims <- data$slider_selectizes$numerics$ranges[[input[[paste0(id , ".x.axis")]]]][1:2]
  } else {
    data$plots[[paste0(id , ".scatter.plot")]]$x.lims <- unique(data$general$data.set[, input[[paste0(id , ".x.axis")]]])
  }
  if(all(is.numeric(data$general$data.set[,input[[paste0(id , ".y.axis")]]]))){
    data$plots[[paste0(id , ".scatter.plot")]]$y.lims <- data$slider_selectizes$numerics$ranges[[input[[paste0(id , ".y.axis")]]]][1:2]
  } else {
    data$plots[[paste0(id , ".scatter.plot")]]$y.lims <- unique(data$general$data.set[, input[[paste0(id , ".y.axis")]]])
  }
  data$plots[[paste0(id , ".scatter.plot")]]$n.samples <- nrow(data$general$data.set)
  plot_scatter_plot(
    data = data
    , statics = statics
    , input = input
    , id = id
    , session = session
  )
  
  output[[paste0(id , ".scatter.plot")]] <- renderPlot(({data$output$plots$scatter[[id]]}))
  reset_brushes(session = session)
  
}


# setting the next quantil error message
next_quantil_error <- function(data, cut, session){
  gene <- data$genes$quantil.error$genes[1]
  upper.quantil <- data$genes$quantil.error$upper.quantils[1]
  lower.quantil <- data$genes$quantil.error$lower.quantils[1]
  expr.vec <- data$genes$quantil.error$expr.vecs[, 1]
  sorted.expr.vec <- expr.vec[order(expr.vec)]
  quantil <- round(length(sorted.expr.vec)*cut/100)
  lower.text <- paste0("Lower quantil: ",
                       round(lower.quantil/length(sorted.expr.vec), digits = 3)*100,
                       "% (expr.val: ", sorted.expr.vec[lower.quantil], ")")
  upper.text <- paste0("Upper quantil: ",
                       round(upper.quantil/length(sorted.expr.vec), digits = 3)*100,
                       "% (expr.val: ", sorted.expr.vec[upper.quantil], ")")
  # if one of those quantils isn't defined, because the
  # original value was the highest or lowest value
  if(is.infinite(lower.quantil) & is.infinite(upper.quantil)){
    text <- paste0("There are several cells which have the same value as your choosen",
                   " quantil for the gene ", gene,". But this gene is constant over all quantils, so no",
                   " quantil is defineable.")
    show.cancel <- FALSE
    show.confirm <- TRUE
    upper.text <- "OK"
  }
  if(is.infinite(lower.quantil)){
    text <- paste0("There are several cells which have the same value as your choosen",
                   " quantil for the gene ", gene,". The next lower quantil would be 0%.",
                   "Please select the upper quantil."
    )
    show.cancel <- FALSE
    show.confirm <- TRUE
  } else if(is.infinite(upper.quantil)){
    text <- paste0("There are several cells which have the same value as your choosen",
                   " quantil for the gene ", gene,". The next upper quantil would be 100%.",
                   "Please select the lower quantil."
    )
    show.cancel <- TRUE
    show.confirm <- FALSE
  } else {
    text <- paste0("There are several cells which have the same value as your choosen",
           " quantil. Please select one of the suggested quantils, which ",
           "are the next upper or lower quantil")
    show.cancel <- TRUE
    show.confirm <- TRUE
  }
  # and asking the user; if quantils for more than one genes aren't well defined
  # the observe event for input$quantil.error will ask for all other quantils
  shinyalert(
    title = "quantil not defined",
    text = text,
    closeOnEsc = FALSE,
    showCancelButton = show.cancel,
    showConfirmButton = show.confirm,
    confirmButtonText = upper.text,
    confirmButtonCol = "#d0d0d0",
    cancelButtonText = lower.text,
    animation = FALSE,
    size = "m",
    inputId = "quantil.error",
    immediate = T,
    session = session
  )
}

# function to show all widely used tags
show_tags <- function(statics, data){
  data$output$texts$diggeR$tags <- paste0(
    "These tags are widely used and are our standards:"
    , statics$diggeR.tags
  )
  data$output$texts$diggeR$tags.ad <- paste0(
    "These tags are widely used and are our standards:"
    , statics$diggeR.tags
  )
} 

# function to close the shown tags
close_tags <- function(data){
  data$output$texts$diggeR$tags <- ""
  data$output$texts$diggeR$tags.ad <- ""
}


#function to calculates the pathway scores
calc.pathway <- function(
  a.set.of.genes
  , the.counts
  , verbose
  , data
){
  if(verbose){
    print(
      paste0(
        "calculating avg.z.score at ", Sys.time()
      )
    )
  }
  
  avg.z.score <- avg_z_score(
    gene.vec = a.set.of.genes
    , count.matrix = the.counts
  )
  celltype.column <- which(
    grepl(
      pattern = "celltype|celltypes"
      , x = colnames(data$general$data.set)
      )
    )
  celltypes <- data$general$data.set[, celltype.column]
  names(celltypes) <- data$general$data.set$colnames
  
  #will be the meta data for the seurat object
  meta.data<- data.frame(
    "celltypes" = celltypes
    , "z.score" = avg.z.score[names(celltypes)]
    )
  rownames(meta.data) <- data$general$data.set$colnames
  the.data <- CreateSeuratObject(
    counts = the.counts
    , meta.data = meta.data
    )
  complete.return.frame <- data.frame(
    "pathway_cluster" = "cluster",
    "mean_cluster_z_score" = 0
    )
  for(celltype in unique(meta.data$celltypes)){
    if(verbose){
      print(
        paste0(
          "starting with celltype "
          , celltype
          , " ("
          , which(unique(meta.data$celltypes) == celltype)
          , " of "
          , length(unique(meta.data$celltypes))
          , ")"
          )
        )
    }
    # calculating the pathway scores
    sub.data<- subset(x = the.data, subset = celltypes == celltype)
    the.counts <- as.matrix(GetAssayData(sub.data, assay = "RNA", slot = "data"))
    try.res <- try(sub.data <- FindVariableFeatures(sub.data, verbose = verbose))
    if(class(try.res)=="try-error"){
      return(TRUE)
    }
    sub.data <- ScaleData(sub.data, verbose = verbose)
    sub.data <- RunPCA(sub.data, verbose = verbose)
    sub.data <- FindNeighbors(
      sub.data
      , dims = 1:min(ncol(sub.data@reductions$pca), 10)
      , verbose = verbose
      )
    # 2 = using louvain
    sub.data <- FindClusters(sub.data, algorithm = 2, verbose = verbose)
    z.score.added <- rep(0, times = length(sub.data$orig.ident))
    names(z.score.added) <- rownames(sub.data@meta.data)
    for(cluster in unique(sub.data$seurat_clusters)){
      a <- mean(sub.data$z.score[sub.data$seurat_clusters==cluster])
      z.score.added[which(sub.data$seurat_clusters == cluster)] <- a
    }
    
    return.frame <- data.frame(
      "pathway_cluster" = paste0(
        "cluster"
        , sub.data$seurat_clusters
        , "_"
        , celltype
        )
      , "mean_cluster_z_score" = z.score.added
      , stringsAsFactors = TRUE
    )
    if(verbose){
      print(
        paste0(
          "found "
          , length(unique(return.frame$pathway_cluster))
          , " cluster"
          )
        )
    }
    complete.return.frame <- rbind(complete.return.frame, return.frame)
  }
  complete.return.frame <- complete.return.frame[-1, ]
  return(complete.return.frame)
}

# like floor; only with floating number
floor_f <- function(x, digits = 0){
  x <- x*10^digits
  return(floor(x)/10^digits)
}

# like ceiling; only with floating number
ceiling_f <- function(x, digits = 0){
  x <- x*10^digits
  return(ceiling(x)/10^digits)
}

#' function that extracts a minimum, a maximum and a suggested step size, which
#' is (max-min)/length(vector)/10 out of an data vector. It also rounds every returned
#' element to a logically usefull amount of digits (min gets rounded down, max gets
#' rounded up)
min_max_step <- function(data.vec){
  #if every element has the same value one just returns an unrounded an unrounded min and max
  # and a stepsize of 0
  if(all((data.vec - data.vec[1])==0)){
    return(c(data.vec[1], data.vec[1], 0))
  }
  #sorting it and computing the stepsize from one element in the sorted vector to the next element
  data.vec <- sort(data.vec)
  steps <- data.vec[2:length(data.vec)] - data.vec[1:(length(data.vec)-1)]
  #excluding the stepsizes of zero
  steps <- steps[which(steps != 0 )]
  # now rounding it to the mean number of digits of the stepsizes

  
  step <- round(x = (max(data.vec) - min(data.vec))/length(data.vec),
        digits = ceiling(abs(log10(1/mean(steps)))))/10
  # if the mean(steps) is 1 the following code would be an endless loop 
  if(mean(steps) == 1){
    step <- 1
  }
  # if the step is 0 the digits argument in round is to low;
  # in the following I increase it as long as step is zero
  if(step == 0){
    mul <- 2
    while(step == 0){
      step <- round(x = (max(data.vec) - min(data.vec))/length(data.vec),
                    digits = ceiling(abs(log10(1/mean(steps)))*mul))/10
      mul <- mul +1
    }
    
  }
  
  return(
    c(
      floor_f(
        x = min(data.vec)
        , digits = ceiling(abs(log10(1/mean(steps)))))
      , ceiling_f(
        x = max(data.vec)
        , digits = ceiling(abs(log10(1/mean(steps)))))
      , step
      )
    )
}

#' function which updates the selectizes after e.g. a different category 
#' got selected via a corresponding selectInput also works as an 
#' initialiation for all selectizes
update_selectizes <- function(data,
                              statics,
                              input,
                              session){
  if(statics$verbose){
    cat(file = stderr(), "\nin update_selectizes")
  }
  pheno.data <- data$general$data.set
  is.slider <- data$slider_selectizes$booleans$is.numeric
  # there are only those columns available for a selectize which aren't currently selected
  available.selectize <- names(which(!is.slider))
  
  if(statics$verbose){
    cat(file = stderr(), "\nreducing available.selectize")
  }
  available.selectize <- available.selectize[
    which(available.selectize != input$first.selectize.col &
            available.selectize != input$second.selectize.col)
    ]
  # a mapper to map 1 -> first, 2 -> second
  word.mapper <- c("first", "second")
  # getting all selected entries to afterwards update the selectizes correctly
  selected <- list()
  for(id in c(1,2)){
    # first of all we check if there is a pheno entry available for this selectize
    if(statics$verbose){
      cat(file = stderr(), paste0("id: ", id))
    }
    if(length(which(!is.slider)) > (id-1)){
      # if there is one available we check if the currently selected entry is a
      # part of the pheno
      # if one switches datasets or loads a dataset at the first time the currently selected
      # entry is one of the previous loaded dataset and not of the dataset which gets 
      # loaded here; therefore we need to check this and if its no part of the currently
      # selected dataset we set it to the id. (first or second) available pheno entry
      cur_sel <- input[[paste0(word.mapper[id], '.selectize.col')]]
      if(!(cur_sel %in% names(which(!is.slider)))){
        cur_sel <- names(which(!is.slider))[id]
      }
      # then we update the SelectInput
      updateSelectInput(
        inputId = paste0(word.mapper[id], ".selectize.col")
        , label = paste0("category for ", word.mapper[id], " selectize Input")
        , choices = sort(c(available.selectize, cur_sel))
        , session = session
        , selected = cur_sel
        )
      # gettiung all selected entries to afterwards update the selectizes correctly
      if(statics$verbose){
        cat(file = stderr(), "\ngetting selected entries")
      }
      if(data$slider_selectizes$booleans$selectize.switch[id]){
        selected[[id]] <- unique(pheno.data[data$slider_selectizes$booleans$selected[, cur_sel], cur_sel])
      }
      else{
        selected[[id]] <- unique(pheno.data[!data$slider_selectizes$booleans$selected[, cur_sel], cur_sel])
      }
      updateSelectizeInput(
        inputId = paste0(word.mapper[id], ".selectize")
        , label = cur_sel
        , choices = sort(unique(pheno.data[, cur_sel]))
        , session = session
        , server = TRUE
        , selected = selected[[id]]
      )
      updateSwitchInput(
        session = session
        , inputId = paste0("in.ex.", word.mapper[id], ".selectize")
        , value = data$slider_selectizes$booleans$selectizw.switch[id]
      )
    } else {
      # if there aren't enough pheno entries
      updateSelectInput(
        label = "deactivated; to less pheno entries"
        , inputId = paste0(word.mapper[id], ".selectize.col")
        , choices = c("None")
        , selected = "None"
        , session = session
      )
      updateSelectizeInput(
        inputId = paste0(word.mapper[id], ".selectize")
        , label = "deactivated; to less pheno entries"
        , choices = NULL
        , session = session
      )
      updateSwitchInput(
        inputId = paste0("in.ex.", word.mapper[id], ".selectize")
        , onLabel = "deactivated"
        , offLabel = "deactivated"
        , session = session
      )
    }
  }
  if(statics$verbose){
    cat(file = stderr(), "\nudpate data subset options")
  }
  updateSelectInput(
    inputId = "reduce.col"
    , label = "select column from which uniform data sampling"
    , choices = c("completly random", names(which(!is.slider)))
    , session = session
    , selected = "completly random"
  )
  
  if(statics$verbose){
    cat(file = stderr(), "\nDONE with update_selectizes")
  }
}




#' function which updates the sliders after e.g. a different category 
#' got selected via a corresponding selectInput also works as an 
#' initialiation for all sliders
update_sliders <- function(data,
                           statics,
                           input,
                           session){
  if(statics$verbose){
    cat(file = stderr(), "\nin update_sliders")
  }
  pheno.data <- data$general$data.set
  is.slider <- data$slider_selectizes$booleans$is.numeric
  if(statics$verbose){
    cat(file = stderr(), "\ngetting available.sliders")
  }
  available.sliders <- names(which(is.slider))
  # maps 1->first, 2->second etc.
  word.mapper <- c("first", "second", "third", "fourth")
  # get all columns of the pheno which aren't currently selected
  # since we don't want to open the posibility to have multiple
  # sliders from the same entry
  if(statics$verbose){
    cat(file = stderr(), "\nreducing available.sliders")
  }

  available.sliders <- available.sliders[which(available.sliders != input$first.slider.col &
                                                 available.sliders != input$second.slider.col &
                                                 available.sliders != input$third.slider.col &
                                                 available.sliders != input$fourth.slider.col)]
  slider.values <- data$slider_selectizes$numerics$slider.values
  ranges <- data$slider_selectizes$numerics$ranges
  if(statics$verbose){
    cat(file = stderr(), "\nupdating sliders and select input")
  }
  # updating sliders and selectInput
  for(id in 1:4){
    if(length(which(is.slider)) > id){
      # if some of the input hasn't changed during a switch of the dataset,
      # it can be that the selected input of the SelectInput isn't a pheno
      # entry and therefore the code will crashed afterwards
      if(input[[paste0(word.mapper[id], '.slider.col')]] %in% names(which(is.slider))){
        selected <- input[[paste0(word.mapper[id], '.slider.col')]]
      } else{
        selected <-  available.sliders[1]
        available.sliders <- available.sliders[-1]
      }
      updateSelectInput(
        session = session
        , label = paste0("value of the ", word.mapper[id], " slider")
        , inputId = paste0(word.mapper[id], ".slider.col")
        , choices = sort(c(available.sliders, selected))
        , selected = selected
      )
      updateSliderInput(
        session = session
        , inputId = paste0(word.mapper[id], ".slider")
        , label = selected
        , min = ranges[[selected]][1]
        , max = ranges[[selected]][2]
        , value = slider.values[[selected]]
        , step = ranges[[selected]][3]
      )
    
  }
    else{
      updateSliderInput(
        session = session
        , inputId = paste0(word.mapper[id], ".slider")
        , label = "deactivated; to less pheno entries"
        , min = 0
        , max = 0
        , value = c(0,0)
      )
      updateSelectInput(
        session = session
        , label = "deactivated; to less pheno entries"
        , inputId = paste0(word.mapper[id], ".slider.col")
        , choices = c("None")
        , selected = "None"
      )
    }
  }
  if(statics$verbose){
    cat(file = stderr(), "\nDONE with update_sliders")
  }
}
  
#' function to plot the 4 histograms of the sliders at the bottom
#' of the app
plot_histograms <- function(data,
                            statics,
                            input,
                            session){ 
  if(statics$verbose){
    cat(file = stderr(), "\nin plot_histograms")
  }
  timer.new <- Sys.time()
  if(is.null(data$timer$plot_histograms)){
    data$timer$plot_histograms <- Sys.time()
    timer.new <- data$timer$plot_histograms + 2
  }
  time.dif <- as.numeric(timer.new - data$timer$plot_histograms)
  data$timer$plot_histograms <- Sys.time()
  if(time.dif > 1){
    pheno.data <- data$general$data.set
    # get all selected profiles
    selected.profiles <- select_profiles(data = data, statics = statics)
    
    if(statics$verbose){
      cat(file = stderr(), "\nupdating the histograms")
    }
    all.names <- c()
    for(id in c("first", "second", "third", "fourth")){
      if(statics$verbose) cat(file = stderr(), paste0("\nid: ", id))
      # the selected column through the first select input
      variable.name <- input[[paste0(id, '.slider.col')]]
      if(statics$verbose){
        cat(file = stderr(), paste0("\nupdating histogram with ", variable.name))
      }
      if(variable.name != "None"){
        # getting now the subsetted data
        subsetted <- pheno.data[
          selected.profiles, 
          variable.name
          , drop = FALSE]
        # plotting it
        colnames(subsetted) <- "tmp"
        data$output$plots$slider.histograms[[id]] <- ggplot(
          data = subsetted
          , mapping = aes(x = tmp)
        ) + geom_histogram(bins = 100) + 
          xlim(range(pheno.data[, variable.name])) +
          xlab(variable.name)
        all.names <- c(all.names, variable.name)
      }
    }
    updateCheckboxGroupInput(
      session = session
      , inputId = "hist.save"
      , choices = all.names
      , selected = all.names
    )
  }
  if(statics$verbose){
    cat(file = stderr(), "\nDONE with plot_histograms")
  }
}


# function that returns all selected profiles
select_profiles <- function(data,
                            statics){
  if(statics$verbose){
    cat(file = stderr(), "\nin select_profiles")
  }
  pheno.data <- data$general$data.set
  # at the beginning all profiles are selected
  selected.profiles <- rep(TRUE, times = nrow(pheno.data))
  # iterating through all columns of the pheno...
  for(cn in colnames(pheno.data)){
    # ... if they are sliders
    if(data$slider_selectizes$booleans$is.numeric[cn]){
      
      # ... selected profiles only is true for the data points which
      # lies between the selected slider values
      selected.profiles <- pheno.data[, cn] >= data$slider_selectizes$numerics$slider.values[[cn]][1] &
        pheno.data[, cn] <= data$slider_selectizes$numerics$slider.values[[cn]][2] &
        selected.profiles
    }
  }
  
  
  # since data$slider_selectizes$booleans$selected is a boolean matrix with profiles in the rows
  # and pheno columns in the columns we can just use an apply
  selected.profiles <- selected.profiles & apply(X = data$slider_selectizes$booleans$selected, MARGIN = 1, FUN = all)
  if(statics$verbose){
    cat(file = stderr(), "\nDONE with select_profiles")
  }
  return(which(selected.profiles))
}

# plots one of the two big scatter plots, which corresponds 
# to id (either 'first' or 'second')
plot_scatter_plot <- function(data,
                              statics,
                              input,
                              id,
                              session){
  
  if(statics$verbose){
    cat(file = stderr(), "\n\nin plot_scatter_plot\nchecking timer")
  }
  timer.new <- Sys.time()
  if(is.null(data$timer$plot_scatter_plot[[id]])){
    data$timer$plot_scatter_plot[[id]] <- Sys.time()
    timer.new <- data$timer$plot_scatter_plot[[id]] + 2.5
  }
  time.dif <- as.numeric(timer.new - data$timer$plot_scatter_plot[[id]])
  data$timer$plot_scatter_plot[[id]] <- Sys.time()
  if(time.dif > 1){
    if(statics$verbose){
      cat(file = stderr(), "\nclear brushes")
    }
    session$resetBrush(paste0(id, ".scatter.brush"))
    
    pheno.data <- data$general$data.set
    if(statics$verbose){
      cat(file = stderr(), "\ngetting the selected profiles")
    }
    # getting the selected profiles
    selected.profiles <- select_profiles(data = data, statics = statics)
    data$plots[[paste0(id, ".scatter.plot")]]$selected.profiles <- selected.profiles
    if(statics$verbose){
      cat(file = stderr(), "\nupdating the big plot")
    }
    coloring <- input[[paste0(id, ".colouring")]]
    if(!coloring %in% colnames(pheno.data)) coloring <- colnames(pheno.data)[1]
    # we need the choosen x.axis, y.axis and the colouring
    data$plots[[paste0(id, ".scatter.plot")]]$coloring <- coloring
    if(statics$verbose) cat(file = stderr(), "\nselecting the data")
    plot.data <- pheno.data[
      selected.profiles, 
      c(input[[paste0(id, ".x.axis")]], input[[paste0(id, ".y.axis")]], coloring)
    ]
    colnames(plot.data) <- c("x", "y", "coloring")
    # if the data is categorical, we now need to calculate the zooming in 
    # if its numerical it later gets calculated in
    if(!all(is.numeric(data$general$data.set[, input[[paste0(id, ".y.axis")]]]))){
      choosen.cats <- data$plots[[paste0(id, ".scatter.plot")]]$y.lims
      zoomed <- which(plot.data$y %in% choosen.cats)
      plot.data <- plot.data[zoomed, ]
    }
    if(!all(is.numeric(data$general$data.set[, input[[paste0(id, ".x.axis")]]]))){
      choosen.cats <- data$plots[[paste0(id, ".scatter.plot")]]$x.lims
      zoomed <- which(plot.data$x %in% choosen.cats)
      plot.data <- plot.data[zoomed, ]
    }
    if(statics$verbose) cat(file = stderr(), "\nplotting it")
    # plotting this
    if(!is.null(data$plots[[paste0(id, ".scatter.plot")]]$n.samples)){
      if(length(selected.profiles) > data$plots[[paste0(id, ".scatter.plot")]]$n.samples){
        edit <- paste0("; zoomed in on: ", data$plots[[paste0(id, ".scatter.plot")]]$n.samples)
      } else {
        edit <- ""
      }
    } else {
      edit <- ""
    }
  
    data$output$plots$scatter[[id]] <- ggplot(
      data = plot.data, 
      aes_string(x = "x", y = "y", col = "coloring")
    ) + 
      geom_point() + 
      ggtitle(
        paste0("With current subsetting: ", length(selected.profiles), " selected"
               , edit)
      ) + xlab(
        input[[paste0(id, ".x.axis")]]
        ) + ylab(
          input[[paste0(id, ".y.axis")]]
          ) + labs(
            color = coloring
            ) 
    # if both axis aren't numeric we plot a jitter plot
    if(
      !is.numeric(plot.data[, "x"]) && 
      !is.numeric(plot.data[, "y"]) 
    ){
      data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + 
        geom_jitter(width = 0.2, height = 0.2)  
    }
    
    #if one axis isn't numeric we additionally plot a violin plot
    if(
      !is.numeric(plot.data[, "x"]) ||
      !is.numeric(plot.data[, "y"]) 
    ){
      data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + 
        geom_violin(scale = "width")  + geom_point()
    }
    
    if(statics$verbose){
      cat(file = stderr(), "\ncoloring it")
    }
    # for the colouring the plot we have to distinguish between categorical, 
    # and continous variables: 
    if(is.numeric(plot.data[, "coloring"])){ 
      # This means the the colouring is continous...
      
      # i don't know what number to pick here. 
      # this is the number of colours in the colour palete: 
      n.colours <- 100
      if(input[[paste0(id, ".colour.scheme")]] == "viridis"){
        data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + scale_color_viridis_c()
      } 
      if(input[[paste0(id, ".colour.scheme")]] == "rainbow"){
        data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + 
          scale_color_gradientn(colours = rainbow(n.colours))
      }
      
      if(input[[paste0(id, ".colour.scheme")]] == "ggplot-default"){
        # default stays as it is: 
        data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] 
      }
      if(input[[paste0(id, ".colour.scheme")]] == "heat"){
        data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + 
          scale_color_gradientn(colours = heat.colors(n.colours))
      }
      if(input[[paste0(id, ".colour.scheme")]] == "reversed-heat"){
        
        data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + scale_colour_gradient2(
          low = "#fcffc2",
          high = "#000000",
          midpoint = mean(plot.data[, "coloring"])+1.8*sd(plot.data[, "coloring"]),
          mid = "#e60000",
          limits = c(min(plot.data[, "coloring"]), max(plot.data[, "coloring"]))
        ) 
      }
    } else {
      # this means the colouring is categorical/discrete: 
      n.colours <- length(unique(plot.data[, "coloring"]))
      
      if(input[[paste0(id, ".colour.scheme")]] == "viridis"){
        data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + scale_color_viridis_d()
      } 
      if(input[[paste0(id, ".colour.scheme")]] == "rainbow"){
        data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + 
          scale_color_manual(values = rainbow(n.colours))
      }
      
      if(input[[paste0(id, ".colour.scheme")]] == "ggplot-default"){
        # default stays as it is: 
        data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] 
      }
      if(input[[paste0(id, ".colour.scheme")]] == "heat"){
        data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + 
          scale_color_manual(values = heat.colors(n.colours))
      }
      if(input[[paste0(id, ".colour.scheme")]] == "reversed-heat"){
        data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + 
          scale_color_manual(values = rev(heat.colors(n.colours)))
      }
    }
    if(statics$verbose){
      cat(file = stderr(), "\nscaling the plot")
    }
    
    if(all(is.numeric(data$general$data.set[, input[[paste0(id, ".x.axis")]]]))){
      data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + 
        xlim(data$plots[[paste0(id, ".scatter.plot")]]$x.lims)
    }
    if(all(is.numeric(data$general$data.set[, input[[paste0(id, ".y.axis")]]]))){
      data$output$plots$scatter[[id]] <- data$output$plots$scatter[[id]] + 
        ylim(data$plots[[paste0(id, ".scatter.plot")]]$y.lims)
    }
    update_saves(
      data = data
      , statics = statics
      , id = id
      , session = session
      , input = input
      , type = "scatter"
    )
    data$output$texts[[paste0(id, ".brush")]] <- "Via brush selected profiles: None"
    data$plots[[paste0(id, ".scatter.plot")]]$active <- TRUE
    data$plots[[paste0(id, ".scatter.plot")]]$brush.values <- NULL
  }
  if(statics$verbose){
    cat(file = stderr(), "\nDONE with plot_scatter_plot")
  }
}
# remove a rectangle which highlights a brush
remove_brush_rect <- function(id,
                              data,
                              statics,
                              session){
  if(statics$verbose){
    cat(file = stderr(), "\n\nin remove_brush_rect")
  }
  
  session$resetBrush(paste0(id, ".scatter.brush"))
  sc.plot <- data$output$plots$scatter[[id]]
  if("brush.rect" %in% names(sc.plot$layers)){
    sc.plot$layers[[which(names(sc.plot$layers) == "brush.rect")]] <- NULL
    data$output$plots$scatter[[id]] <- sc.plot
  }
  if(statics$verbose){
    cat(file = stderr(), "\nDONE with remove_brush_rect")
  }
}
# updating which plot can be saved
update_saves <- function(data,
                         statics,
                          type,
                          id,
                          session,
                          input){

  if(statics$verbose){
    cat(file = stderr(), "\n\ninupdate_saves")
  }
  # if a new dataset is loaded, no scatterplot will be printed and therefore 
  # one shouldn't be able to select a scatterplot to save
  if(id == "clear"){
    updateCheckboxGroupInput(
      session = session
      , inputId = paste0(type, ".save")
      , choices = "None"
      , selected = NULL
    )
  } else {
    data$variables[[paste0(type, ".save")]] <- unique(c(data$variables[[paste0(type, ".save")]], id))
    if(type == "brush") edit <- " histogram"
    if(type == "scatter") edit <- " plot"
    updateCheckboxGroupInput(
      session = session
      , inputId = paste0(type, ".save")
      , choices = paste0(data$variables[[paste0(type, ".save")]], paste0(" ", type, edit))
      , selected = unique(c(input[[paste0(type, ".save")]], paste0(id, paste0(" ", type, edit))))
    )
  }
}
# function to zoom into a scatter plot
zoom_in <- function(data,
                    statics,
                    input,
                    session,
                    id){
  if(statics$verbose){
    cat(file = stderr(), "\n\nin zoom_in")
  }
  timer.new <- Sys.time()
  if(is.null(data$timer$zoom_in[[id]])){
    data$timer$zoom_in[[id]] <- Sys.time()
    timer.new <- data$timer$zoom_in[[id]] + 4
  }
  time.dif <- as.numeric(timer.new - data$timer$zoom_in[[id]])
  data$timer$zoom_in[[id]] <- Sys.time()
  if(time.dif > 1){
    # if the data is numerical
    pheno.data <- data$general$data.set
    # getting all profiles which are currently selected via the sliders and
    # selectizes
    selected.profiles <- select_profiles(data = data, statics = statics)
    #reducing the data to these profiles
    if(statics$verbose) cat(file = stderr(), "\nreducing it to brushed profiles")
    # we need the choosen x.axis, y.axis
    selected.data <- pheno.data[
      selected.profiles, 
      c(input[[paste0(id, ".x.axis")]], input[[paste0(id, ".y.axis")]])
    ]
    colnames(selected.data) <- c("x", "y")
    brushed.data <- brushedPoints(
      df = selected.data
      , brush = input[[paste0(id, ".scatter.brush")]]
    )
    if(all(is.numeric(data$general$data.set[,input$first.x.axis]))){
      xlims <- c(
        input[[paste0(id, ".scatter.brush")]]$xmin, input[[paste0(id, ".scatter.brush")]]$xmax
      )
    } else {
      xlims <- unique(brushed.data$x)
    }
    if(!is.null(xlims)) data$plots[[paste0(id, ".scatter.plot")]]$x.lims <- xlims
    if(all(is.numeric(data$general$data.set[,input$first.y.axis]))){
      ylims <- c(
        input[[paste0(id, ".scatter.brush")]]$ymin, input[[paste0(id, ".scatter.brush")]]$ymax
      )
    } else {
      ylims <- unique(brushed.data$y)
    }
    if(!is.null(ylims)) data$plots[[paste0(id, ".scatter.plot")]]$y.lims <- ylims
    data$plots[[paste0(id, ".scatter.plot")]]$n.samples <- nrow(brushed.data)
    if(statics$verbose){
      cat(file = stderr(), "\nDONE with zoom_in")
    }
  }
  
  

}

# reset the brushes
reset_brushes <- function(session){
  
  for(id in c("first", "second", "third", "fourth")){
    session$resetBrush(paste0(id, ".hist.brush"))
  }
  for(id in c("first", "second")){
    session$resetBrush(paste0(id, ".scatter.brush"))
  }
}


show_brush_rect <- function(statics,
                            data,
                            id,
                            session){
  if(statics$verbose){
    cat(file = stderr(), "\n\nin show.brush.rect")
  }
  
  timer.new <- Sys.time()
  if(is.null(data$timer$show_brush_rect[[id]])){
    data$timer$show_brush_rect[[id]] <- Sys.time()
    timer.new <- data$timer$show_brush_rect[[id]] + 4
  }
  time.dif <- as.numeric(timer.new - data$timer$show_brush_rect[[id]])
  data$timer$show_brush_rect[[id]] <- Sys.time()
  if(time.dif > 1){
    session$resetBrush(paste0(id, ".scatter.brush"))
    brush.values <- data$plots[[paste0(id, ".scatter.plot")]]$brush.values
    scatter.plot <- data$output$plots$scatter[[id]]
    if(!is.null(brush.values) & !("brush.rect" %in% names(scatter.plot$layers))){
      scatter.plot <- scatter.plot +
        geom_rect(
          data = brush.values
          , inherit.aes = FALSE
          , mapping = aes(
            xmin = xmin
            , xmax = xmax
            , ymin = ymin
            , ymax = ymax
          )
          , color = "#66b3ff"
          , alpha = 0.07
          , fill = "darkblue"
          , size = 0.1
        )
      names(scatter.plot$layers)[length(scatter.plot$layers)] <- "brush.rect"
      data$output$plots$scatter[[id]] <- scatter.plot 
    }
  }
  
  if(statics$verbose){
    cat(file = stderr(), "\n\nDONE with show.brush.rect")
  }
}

# plots the 2 histograms of the brushed data (first: all
# profiles which got selected via brush, second: the
# unselected profiles)
plot_brush_histograms <- function(data,
                                  statics,
                                  input,
                                  session,
                                  id){
  if(statics$verbose){
    cat(file = stderr(), "\n\nin plot_brush_histograms")
  }
  # first checking wether there is a plot
  if(is.null(data$plots[[paste0(id, ".scatter.plot")]]$active)) data$plots[[paste0(id, ".scatter.plot")]]$active <- FALSE
  if(data$plots[[paste0(id, ".scatter.plot")]]$active){
    timer.new <- Sys.time()
    if(is.null(data$timer$plot_brush_histograms[[id]])){
      data$timer$plot_brush_histograms[[id]] <- Sys.time()
      timer.new <- data$timer$plot_brush_histograms[[id]] + 4
    }
    time.dif <- as.numeric(timer.new - data$timer$plot_brush_histograms[[id]])
    data$timer$plot_brush_histograms[[id]] <- Sys.time()
    if(time.dif > 1){
      scatter.plot <- data$output$plots$scatter[[id]]
      scatter.plot <- ggplot_build(scatter.plot)
      plot.data <- data$plots[[paste0(id, ".scatter.plot")]]
      dataset <- data.frame(
        "plot_entry" = scatter.plot$plot$data$coloring
        , "x" = scatter.plot$plot$data$x
        , "y" = scatter.plot$plot$data$y
      )
      
      # if no samples are plotted; there will be no colnames; therefore I check if the needed
      # colnames are there
      if("x" %in% colnames(dataset) & !is.null(input[[paste0(id, ".scatter.brush")]]) & "y" %in% colnames(dataset)){
        is.numerical <- is.numeric(dataset$plot_entry)
        if(!is.numerical){
          coloring <- unique(scatter.plot$data[[1]]$colour)
          names(coloring) <- unique(scatter.plot$plot$data$coloring)
        }
        brushed.dataset <- brushedPoints(
          df = dataset
          , brush = input[[paste0(id, ".scatter.brush")]]
          , allRows = TRUE
        )
        selected.data <- brushed.dataset[brushed.dataset$selected_, ]
        not.selected.data <- brushed.dataset[!brushed.dataset$selected_, ]
        if(!is.numerical){
          not.selected.data$plot_entry <- factor(not.selected.data$plot_entry
                                                 , levels = unique(not.selected.data$plot_entry))
          selected.data$plot_entry <- factor(selected.data$plot_entry
                                             , levels = unique(selected.data$plot_entry))
        }
        #plotting the two histograms
        if(statics$verbose) cat(file = stderr(), "\nplotting the histograms\nplotting the not.selected.hist.plot")
        data$output$plots$brush.histograms[[paste0(id, ".not.selected")]] <- ggplot(
          data = not.selected.data
          , mapping = aes(x = plot_entry)
        ) + theme(axis.text.x = element_text(angle = 75, hjust = 1))
        
        if(is.numerical){
          
          data$output$plots$brush.histograms[[paste0(id, ".not.selected")]] <- data$output$plots$brush.histograms[[paste0(id, ".not.selected")]] + 
            geom_histogram() + 
            xlim(range(dataset$plot_entry))
        } else {
          count.max <- max(
            c(
              table(selected.data[, "plot_entry"]), 
              table(not.selected.data[, "plot_entry"])
            )
          )
          # this rounds count max to the next full 100 multiple
          count.max <- count.max + 100 - count.max %% 10
          data$output$plots$brush.histograms[[paste0(id, ".not.selected")]] <- data$output$plots$brush.histograms[[paste0(id, ".not.selected")]] + 
            geom_histogram(stat = "count", fill = coloring[as.vector(unique(not.selected.data$plot_entry))])+ 
            ylim(c(0, count.max))
        }
        
        data$output$plots$brush.histograms[[paste0(id, ".not.selected")]] <- data$output$plots$brush.histograms[[paste0(id, ".not.selected")]] + 
          xlab(data$plots[[paste0(id, ".scatter.plot")]]$coloring) + 
          ggtitle(paste0("Not selected profiles (",id, " plot)")) 
        if(statics$verbose) cat(file = stderr(), "\nplotting the selected.hist.plot")
        
        
        data$output$plots$brush.histograms[[paste0(id, ".selected")]] <- ggplot(
          data = selected.data
          , mapping = aes(x = plot_entry)
        ) + theme(axis.text.x = element_text(angle = 75, hjust = 1))
        
        if(is.numerical){
          data$output$plots$brush.histograms[[paste0(id, ".selected")]] <- data$output$plots$brush.histograms[[paste0(id, ".selected")]] + 
            geom_histogram() + 
            xlim(range(dataset$plot_entry))
        } else {
          data$output$plots$brush.histograms[[paste0(id, ".selected")]] <- data$output$plots$brush.histograms[[paste0(id, ".selected")]] + 
            geom_histogram(stat = "count", fill = coloring[as.vector(unique(selected.data$plot_entry))]) + 
            ylim(c(0, count.max))
        }
        
        data$output$plots$brush.histograms[[paste0(id, ".selected")]] <- data$output$plots$brush.histograms[[paste0(id, ".selected")]] + 
          xlab(data$plots[[paste0(id, ".scatter.plot")]]$coloring) + 
          ggtitle(paste0("Via brush selected profiles (",id, " plot)"))
        data$output$texts[[paste0(id, ".brush")]] <- paste0("Via brush selected profiles: ", nrow(selected.data))
        if(is.numerical){
          selected.build <- ggplot_build(data$output$plots$brush.histograms[[paste0(id, ".selected")]])
          not.selected.build <- ggplot_build(data$output$plots$brush.histograms[[paste0(id, ".not.selected")]])
          
          maxy <- max(selected.build$data[[1]]$y, not.selected.build$data[[1]]$y)+1
          data$output$plots$brush.histograms[[paste0(id, ".selected")]] <- data$output$plots$brush.histograms[[paste0(id, ".selected")]]+
            ylim(c(0, maxy))
          data$output$plots$brush.histograms[[paste0(id, ".not.selected")]] <- data$output$plots$brush.histograms[[paste0(id, ".not.selected")]] +
            ylim(c(0, maxy))
        }
        update_saves(
          data = data
          , statics = statics
          , id = id
          , session = session
          , input = input
          , type = "brush"
        )
        brush.values <- input[[paste0(id, ".scatter.brush")]]
        brush.values <- data.frame(
          xmin = brush.values$xmin
          , xmax = brush.values$xmax
          , ymin = brush.values$ymin
          , ymax = brush.values$ymax)
        data$plots[[paste0(id, ".scatter.plot")]]$brush.values <- brush.values
        session$resetBrush(paste0(id, ".scatter.brush"))
        remove_brush_rect(
          id = id
          , data = data
          , statics = statics
          , session = session
          )
      } else {
        shinyalert(
          title = "No cells select"
          ,text = paste0("There are no cells selected via the sliders and select input or via brush",
          ". Please select some cells")
          , closeOnEsc = TRUE
          , closeOnClickOutside = TRUE
          , showCancelButton = FALSE
          , showConfirmButton = TRUE
          , session = session
          , animation = FALSE
          , immediate = TRUE
          , size = "m"
        )
      }
    }
      }

  if(statics$verbose){
    cat(file = stderr(), "\nDONE with plot_brush_histograms")
  }
}

# this function updates the text output, which shows,
# which profiles are currently selected
update_subsetting_output <- function(data,
                                     statics,
                                     input){
  if(statics$verbose){
    cat(file = stderr(), "\nin update_subsetting_output")
  }
  sliders.text <- "all sliders:\n"
  selectize.text <- "all selectize:\n"
  for(col.name in names(data$slider_selectizes$booleans$is.numeric)){
    if(data$slider_selectizes$booleans$is.numeric[col.name]){
      if(!all(data$slider_selectizes$numerics$slider.values[[col.name]] == data$slider_selectizes$numerics$ranges[[col.name]][1:2])){
        
        add.text <- gsub(x = Reduce(f = paste,
                                    x = data$slider_selectizes$numerics$slider.values[[col.name]]),
                         replacement = " - ",
                         pattern = " ")
        sliders.text <- paste0(sliders.text,
                               col.name,
                               ": ",
                               add.text,
                               "\n")
      }
    } else{
      if(!all(data$slider_selectizes$booleans$selected[, col.name])){
        selected <- unique(data$general$data.set[data$slider_selectizes$booleans$selected[, col.name], col.name])
        if(length(selected) == 0) selected <- "NONE"
        add.text <- gsub(x = Reduce(f = paste,
                                    x = selected),
                         replacement = ", ",
                         pattern = " ")
        selectize.text <- paste0(selectize.text,
                                 col.name,
                                 ": ",
                                 add.text,
                                 "\n")
      }
    }
  }
  sliders.text <- paste0(sliders.text,
                         "+ everything from all remaining sliders")
  selectize.text <- paste0(selectize.text,
                           "+ everything from all remaining selectizes")
  text.out <- paste0("currently subbsetting:\n\n",
                     sliders.text,
                     "\n\n",
                     selectize.text)
  data$output$texts$currently.subsetting <- text.out
  if(statics$verbose){
    cat(file = stderr(), "\nDONE with update_subsetting_output")
  }
}

# this function updates the plot options
update_plot_options <- function(session,
                                data,
                                statics,
                                input){
  if(statics$verbose){
    cat(file = stderr(), "\n\nin update_plot_options")
  }
  ids <- c(
    "first.colouring"
    , "second.colouring"
    , "first.x.axis"
    , "second.x.axis"
    , "first.y.axis"
    , "second.y.axis"
  )
  i <- 1
  for(input.id in ids){
    if(input[[input.id]] == "None"){
      if(grepl(pattern = "colouring", x = input.id) & "celltype" %in% colnames(data$general$data.set)){
        input.val <- "celltype"
      } else if(grepl(pattern = "x.axis", x = input.id) & "UMAP_1" %in% colnames(data$general$data.set)){
        input.val <- "UMAP_1"
      } else if(grepl(pattern = "y.axis", x = input.id) & "UMAP_2" %in% colnames(data$general$data.set)){
        input.val <- "UMAP_2"
      } else {
        input.val <- colnames(data$general$data.set)[i]
        i <- i + 1
      }
    } else{
      input.val <- input[[input.id]]
    }
    
    updateSelectInput(
      session = session
      , inputId =  input.id
      , choices =  sort(colnames(data$general$data.set))
      , selected = input.val
    )
  }
  if(statics$verbose){
    cat(file = stderr(), "\n\nDONE with update_plot_options")
  }
}

# updating everything (slider selectizes plotable feature) after a gene got pulled
pull.genes.helper <- function(data,
                              statics,
                              gene,
                              text){
  if(statics$verbose){
    cat(file = stderr(), "\n\nin pull.genes.helper\nupdating is.slider")
  }
  pheno.data <- data$general$data.set
  is.slider <- data$slider_selectizes$booleans$is.numeric
  # add two new logicals, one for the zscore, one for the categorical:
  is.slider <- c(is.slider, TRUE, TRUE, FALSE)
  names(is.slider)[(length(is.slider)-2):length(is.slider)] <- c(
    paste0("zscore.", gene), paste0("raw.", gene), paste0(text, gene)
  )
  data$slider_selectizes$booleans$is.numeric <- is.slider
  if(statics$verbose){
    cat(file = stderr(), "\nupdating ranges and slider.values")
  }
  
  # update ranges and slider values
  ranges <- data$slider_selectizes$numerics$ranges
  slider.values <- data$slider_selectizes$numerics$slider.values
  if(statics$verbose){
    cat(file = stderr(), "\nexecuting min_max_step")
  }
  ranges$tmp <- min_max_step(pheno.data[ ,paste0("zscore.", gene)])
  slider.values$tmp <- ranges$tmp[1:2]
  names(ranges)[length(ranges)] <- paste0("zscore.", gene)
  names(slider.values)[length(slider.values)] <- paste0("zscore.", gene)
  
  ranges$tmp <- min_max_step(pheno.data[ ,paste0("raw.", gene)])
  slider.values$tmp <- ranges$tmp[1:2]
  names(ranges)[length(ranges)] <- paste0("raw.", gene)
  names(slider.values)[length(slider.values)] <- paste0("raw.", gene)
  data$slider_selectizes$numerics$ranges <- ranges
  data$slider_selectizes$numerics$slider.values <- slider.values
  selected <- data$slider_selectizes$booleans$selected
  selected <- cbind(selected, rep(x = TRUE, times = nrow(selected)))
  colnames(selected)[ncol(selected)] <- paste0(text, gene)
  data$slider_selectizes$booleans$selected <- selected
  if(statics$verbose){
    cat(file = stderr(), "\nDONE with pull.genes.helper")
  }
}
start_loading_screen <- function(session){
  shinyalert(
    title = ""
    , text = "loading please wait"
    , closeOnEsc = FALSE
    , closeOnClickOutside = FALSE
    , showCancelButton = FALSE
    , session = session
    , showConfirmButton = FALSE
    , timer = 0
    , animation = FALSE
    , size = "m"
  )
}
end_loading_screen <- function(session){
  shinyalert(
    title = ""
    , text = "loading please wait"
    , closeOnEsc = FALSE
    , closeOnClickOutside = FALSE
    , session = session
    , showCancelButton = FALSE
    , showConfirmButton = FALSE
    , timer = 1
    , animation = FALSE
    , size = "m"
    , immediate = TRUE
  )
}