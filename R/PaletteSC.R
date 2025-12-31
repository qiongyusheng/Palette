###################################################################
##### Palette core functions
###################################################################
# Written by Qiongyu Sheng on Dec 26, 2025.
#'
#' @import Matrix
#' @importClassesFrom Matrix dgCMatrix
NULL
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))

#' The Palette.Assay Class
#'
#' The Palette.Assay class is an intermediate data storage class that
#' stores a data matrix and other related information needed for
#' performing downstream analyses.
#'
#' @slot Modal Data modality information, such as RNA, ATAC.
#' @slot Sample Data sample/batch information.
#' @slot Raw Raw data matrix (features by cells).
#' @slot data Normalized data matrix (features by cells).
#' @slot RepresentativeCell List of information for representative cells, including cell names and cluster infromation.
#' @slot embedding Intra-modal joint dimensionality reduction results.
#'
#' @exportClass Palette.Assay
#' @name Palette.Assay-class
#' @rdname Palette.Assay-class
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods as new slot slot<-
#' @useDynLib PaletteSC
#'
Palette.Assay <- methods::setClass(
  Class = "Palette.Assay",
  slots = c(
    Modal = "character",
    Sample = "character",
    Raw = "AnyMatrix",
    data = "AnyMatrix",
    RepresentativeCell = "list",
    embedding = "AnyMatrix"
  )
)

#' The Palette Object
#'
#' The Palette Object is a S4 classes for single-cell multimodal data
#' and associated information, such as dimensionality reduction embeddings,
#' mosaic integration results.
#'
#' @slot Assays List of Palette.Assay.
#' @slot Meta A data.frame containing meta informatoin for each Palette.Assay.
#' @slot Data.index List of information for each Palette.Assay, including modality and sample/batch information.
#' @slot Co.Features Identified high variable features shared across batches.
#' @slot pre.clusters List of clustering information for each data matrix.
#' @slot Int.result List of integrated low-dimensional embedding.
#'
#' @exportClass Palette
#' @name Palette-class
#' @rdname Palette-class
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass
#'
Palette <- methods::setClass(
  Class = "Palette",
  slots = c(
    Assays = "list",
    Meta = "data.frame",
    Data.index = "list",
    Co.Features = "list",
    pre.clusters = "list",
    Int.result = "list"
  )
)

#' Create a new Palette object
#'
#' @param data.list List of multimodal data matrix. If the names of the list are not \code{NULL}, then the names information will be used to name each Palette.Assay. If the names of the list are \code{NULL}, then the combined sample and modal information will be used to name each Palette.Assay.
#' @param samples A character vector of up to the same length as data.list, with each element represent the sample/batch information for each data matrix.
#' @param modals A character vector of up to the same length as data.list, with each element represent the modality information for each data matrix.
#' @param filter Whether to filter the data.
#' @param min.cells Include features detected in at least this many cells. Use only if \code{filter = TRUE}.
#' @param min.features Include cells where at least this many features are detected. Use only if \code{filter = TRUE}.
#' @param modal.order Orders of modalities for filter condition
#' @param sparse Whether use sparse format ('dgCMatirx' class). (default TRUE)
#'
#' @return \code{Palette} object.
#'
#' @import Matrix
#' @importFrom rlang %||%
#' @importFrom Rcpp evalCpp
#'
#' @export
Create.Palette.Object <- function(data.list,
                                  samples,
                                  modals,
                                  filter = FALSE,
                                  min.cells = 0,
                                  min.features = 0,
                                  modal.order,
                                  sparse = TRUE){

  if(missing(data.list) || missing(samples) || missing(modals)){
    stop("Missing required data. data.list, samples and modals are the data that must be provided.")
  }

  if(length(data.list) != length(samples) || length(data.list) != length(modals) || length(samples) != length(modals)){
    stop("Wrong match information. The lengths of data.list, samples and modals must be equal.")
  }

  if(filter){
    modal.order <- modal.order %||% unique(modals)
    a <- setdiff(modals,modal.order)
    if(length(a) != 0){
      stop("Incomplete modal information in modal.order")
    }
    b <- setdiff(modal.order,modals)
    if(length(b) != 0){
      stop("Redundant modal information in modal.order")
    }

    min.cells <- min.cells %||% as.numeric(rep(0,length(modal.order)))
    min.features <- min.features %||% as.numeric(rep(0,length(modal.order)))
    if(length(min.cells) != length(modal.order) || length(min.features) != length(modal.order)){
      stop("The lengths of min.cells, min.features and modal.order must be equal.")
    }
  }

  if(Reduce("+",grepl('&', modals)) > 0){
    warning("character '&' in modal name will replaced by '-")
    modals <- gsub('&', '-', modals)
  }
  if(Reduce("+",grepl('_', modals)) > 0){
    warning("character '_' in modal name will replaced by '-")
    modals <- gsub('_', '-', modals)
  }
  if(Reduce("+",grepl('&', samples)) > 0){
    warning("character '&' in sample name will replaced by '-")
    samples <- gsub('&', '-', samples)
  }
  if(Reduce("+",grepl('_', samples)) > 0){
    warning("character '_' in sample name will replaced by '-")
    samples <- gsub('_', '-', samples)
  }

  if(is.null(names(data.list)) || any(names(x = data.list) == '')){
    names(data.list) <- paste0(modals,"_",samples)
  }else{
    if(Reduce("+",grepl('&', names(data.list))) > 0){
      warning("character '&' in data name will replaced by '-")
      names(data.list) <- gsub('&', '-', names(data.list))
    }
    if(Reduce("+",grepl('_', names(data.list))) > 0){
      warning("character '_' in data name will replaced by '-")
      names(data.list) <- gsub('_', '-', names(data.list))
    }
  }
  data.name <- names(data.list)

  data.list <- lapply(data.list, check.data, sparse = sparse)

  if(filter){
    min.cells <- min.cells[match(modals, modal.order)]
    min.features <- min.features[match(modals, modal.order)]

    filter.cell.list <- lapply(1:length(data.list), function(x){
      if(min.features[x] > 0){
        nfeatures <- Matrix::colSums(x = data.list[[x]] > 0)
        x <- colnames(data.list[[x]])[which(nfeatures >= min.features[x])]
      }else{
        x <- colnames(data.list[[x]])
      }
      x
    })

    filter.feature.list <- lapply(1:length(data.list), function(x){
      if(min.cells[x] >0){
        num.cells <- Matrix::rowSums(x = data.list[[x]] > 0)
        x <- rownames(data.list[[x]])[which(num.cells >= min.cells[x])]
      }else{
        x <- rownames(data.list[[x]])
      }
      x
    })
  }else{
    filter.cell.list <- lapply(data.list, colnames)
    filter.feature.list <- lapply(data.list, rownames)
  }

  for(j in unique(modals)){
    idx.modal <- which(modals == j)
    feature.list <- filter.feature.list[idx.modal]
    if(length(idx.modal) > 1){
      modal.feature <- Reduce(intersect, feature.list)
      data.list[idx.modal] <- lapply(1:length(idx.modal), function(x){
        x <- data.list[idx.modal][[x]][modal.feature, ]
        x
      })
    }
  }

  dim.vector <- unlist(lapply(data.list,dim))
  if(is.numeric(Reduce("*",as.numeric(dim.vector))) && Reduce("*",as.numeric(dim.vector)) != 0){
    names(data.list) <- data.name
    assay <- lapply(1:length(data.list),function(x){
      x <- Create.Palette.Assay(count = data.list[[x]],
                                modal = modals[x],
                                sample = samples[x],
                                filter = FALSE,
                                sparse = sparse)
      x
    })
    names(assay) <- data.name
    data.idx <- list(assay.idx = data.frame(sample = samples,
                                            modal = modals,class = rep("original",length(samples))))
    cell.name <- unlist(lapply(unique(samples), function(x){
      idx <- which(samples == x)[1]
      x <- colnames(data.list[[idx]])
      x
    }))
    Palette <- methods::new(
      Class = "Palette",
      Assays = assay,
      Meta = data.frame(cell.name),
      Data.index = data.idx
    )
  }else{
    stop("Check whether the filter conditions or the input information match.")
  }
  return(Palette)
}

check.data <- function(counts,
                       sparse = TRUE){

  if(missing(x = counts)){
    stop("Must provide 'counts'")
  }else if(!missing(x = counts)){
    if(anyDuplicated(x = rownames(x = counts))){
      warning(
        "Non-unique features (rownames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = counts) <- make.unique(names = rownames(x = counts))
    }
    if(anyDuplicated(x = colnames(x = counts))){
      warning(
        "Non-unique cell names (colnames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      colnames(x = counts) <- make.unique(names = colnames(x = counts))
    }
    if (is.null(x = colnames(x = counts))) {
      stop("No cell names (colnames) names present in the input matrix")
    }
    if (any(rownames(x = counts) == '')) {
      stop("Feature names of counts matrix cannot be empty", call. = FALSE)
    }
    if (nrow(x = counts) > 0 && is.null(x = rownames(x = counts))) {
      stop("No feature names (rownames) names present in the input matrix")
    }
    if(sparse){
      if (!inherits(x = counts, what = 'dgCMatrix')){
        counts <- as(as.matrix(counts), 'dgCMatrix')
      }
    }
  }
  return(counts)
}

#' Create a new Palette.Assay class
#'
#' @param count A data matrix with feature by cell.
#' @param filter Whether to filter the data.
#' @param sample A character represents the sample/batch information for the data matrix.
#' @param modal A character represents the modality information for the data matrix.
#' @param filter Whether to filter the data.
#' @param min.cells Include features detected in at least this many cells. Use only if \code{filter = TRUE}.
#' @param min.features Include cells where at least this many features are detected. Use only if \code{filter = TRUE}.
#' @param slot Which slot ('Raw' or 'data') to put the \code{count}.
#' @param sparse Whether to use sparse format ('dgCMatirx' class). (default TRUE)
#'
#' @return \code{Palette.Assay} object.
#'
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @keywords internal
#'
Create.Palette.Assay <- function(count,
                                 sample,
                                 modal,
                                 filter = FALSE,
                                 min.cells = 0,
                                 min.features = 0,
                                 slot = "Raw",
                                 sparse = TRUE){
  count <- check.data(counts = count, sparse = sparse)
  if(filter){
    if (min.cells > 0) {
      num.cells <- Matrix::rowSums(x = count > 0)
      count <- count[which(x = num.cells >= min.cells), ]
    }
    if (min.features > 0) {
      nfeatures <- Matrix::colSums(x = count > 0)
      count <- count[, which(x = nfeatures >= min.features)]
    }
  }
  if (!is.vector(x = rownames(x = count))) {
    rownames(x = count) <- as.vector(x = rownames(x = count))
  }
  if (!is.vector(x = colnames(x = count))) {
    colnames(x = count) <- as.vector(x = colnames(x = count))
  }

  if(slot == "Raw"){
    assay <- methods::new(
      Class = "Palette.Assay",
      Modal = modal,
      Sample = sample,
      Raw = count
    )
  }else if(slot == "data"){
    assay <- methods::new(
      Class = "Palette.Assay",
      Modal = modal,
      Sample = sample,
      Raw = count,
      data = count
    )
  }else{
    stop("Unknown slot.")
  }

  return(assay)
}

#' Add a Palette.Assay to a Palette object
#'
#' @param assay A Palette.Assay.
#' @param Palette.object A Palette object.
#' @param check Whether to check the feature or cell information between \code{assay} and \code{Palette.object}.
#' @param class The data category carried by the assay, where ‘original’ represents real data and ‘infer’ represents inferred data。
#' @param name The name of \code{assay}.
#'
#' @return \code{Palette} object.
#'
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @keywords internal
Add.Assay <- function(assay, Palette.object, check = TRUE, class = "original",name = NULL){

  modal <- assay@Modal
  sample <- assay@Sample
  modals <- Palette.object@Data.index[["assay.idx"]][,"modal"]
  samples <- Palette.object@Data.index[["assay.idx"]][,"sample"]
  classes <- Palette.object@Data.index[["assay.idx"]][,"class"]

  modal_logic <- modals == modal
  sample_logic <- samples == sample

  if(Reduce("+",modal_logic) == 0 && Reduce("+",sample_logic) == 0){
    warning("No modal or sample match between Palette object and new assay")
    new.name <- colnames(assay@Raw)
    if(dim(Palette.object@Meta)[2] > 1){
      Palette.object@Meta <- Palette.object@Meta[,"cell.name"]
      warning("Length mismatched meta data removed")
    }
    cell.name <- c(Palette.object@Meta[,"cell.name"],new.name)
    Palette.object@Meta <- data.frame(cell.name)
  }else{
    if(check){
      if(Reduce("+",modal_logic) != 0){
        idx.m <- which(modal_logic == 1)
        feature.list <- list(old.feature = rownames(Palette.object@Assays[[idx.m[1]]]@Raw),
                             new.feature = rownames(assay@Raw))
        co.feature <- Reduce(intersect, feature.list)
        for(i in idx.m){
          Palette.object@Assays[[i]]@Raw <- Palette.object@Assays[[i]]@Raw[co.feature,]
        }
        assay@Raw <- assay@Raw[co.feature,]
      }
      if(Reduce("+",sample_logic) != 0){
        idx.s <- which(sample_logic == 1)
        cell.list <- list(old.cell = colnames(Palette.object@Assays[[idx.s[1]]]@Raw),
                          new.cell = colnames(assay@Raw))
        co.cell <- Reduce(intersect, cell.list)
        for(i in idx.s){
          Palette.object@Assays[[i]]@Raw <- Palette.object@Assays[[i]]@Raw[ , co.cell]
        }
        assay@Raw <- assay@Raw[,co.cell]
      }
    }
  }
  name <- name %||%paste0(modal,"_",sample)
  Palette.object@Assays[[name]] <- assay
  modals <- c(modals,modal)
  samples <- c(samples,sample)
  classes <- c(classes,class)
  data.idx <- data.frame(sample = samples,
                         modal = modals,
                         class = classes)
  Palette.object@Data.index[["assay.idx"]] <- data.idx
  return(Palette.object)
}

#' Normalize the count data present in a given assay.
#'
#' @param object A \code{Palette} object.
#' @param modal A character vector of modality names to be normalized.
#' @param normal.method Method for normalization.('LogNormalize' for transcriptome, and 'CLR' for proteome)
#' @param scale.factor The scale factor used for each cell. (default 10000)
#' @param margin If performing CLR normalization, normalize across features (1) or cells (2)
#'
#' @return \code{Palette} object.
#'
#' @import Seurat
#'
#' @export
#'
Normalize.Data <- function(object,
                           modal = NULL,
                           normal.method = "LogNormalize",
                           scale.factor = 10000,
                           margin = 1){
  if(!(inherits(modal,what = "character"))){
    stop("Parameter modal should input as character vector.")
  }
  if(!(inherits(normal.method,what = "character"))){
    stop("Parameter normal.method should input as character vector.")
  }

  if(!(inherits(margin,what = "numeric"))){
    stop("Parameter margin should input as numeric vector.")
  }

  if(length(modal)!=length(normal.method) || length(normal.method) != length(margin)){
    stop("Unmatched information.")
  }
  if(!all(modal %in% object@Data.index[["assay.idx"]]$modal)){
    stop("Unknown modality information.")
  }
  if(!all(normal.method %in% c("LogNormalize","CLR"))){
    stop("Unknown normalization method.")
  }
  if(!all(margin %in% c(1,2))){
    stop("The margin parameter must be selected between 1 and 2.")
  }
  for(i in 1:length(modal)){
    modal_logic <- ifelse(object@Data.index[["assay.idx"]]$modal == modal[i], TRUE, FALSE)
    idx <- which(modal_logic)
    if(normal.method[i]=="CLR"){
      for(j in idx){
        object@Assays[[j]]@data <- as(Seurat::NormalizeData(object@Assays[[j]]@Raw,
                                                     normalization.method = normal.method[i],
                                                     scale.factor = scale.factor,
                                                     margin = margin[i],
                                                     verbose = FALSE),"dgCMatrix")
      }
    }else{
      for(j in idx){
        object@Assays[[j]]@data <- Seurat::NormalizeData(object@Assays[[j]]@Raw,
                                                  normalization.method = normal.method[i],
                                                  scale.factor = scale.factor,
                                                  margin = margin[i],
                                                  verbose = FALSE)
      }
    }
  }
  return(object)
}

#' This function performs inverse document frequency (IDF)–based scaling followed by log normalization on the raw assay data of a specified modality (and sample) within a Palette object.
#'
#' @param object A \code{Palette} object.
#' @param modal A character of modality.
#' @param sample A character of sample/batch.
#' @param scale.factor The scale factor used for each cell. (default 10000)
#' @param idf A precomputed IDF vector to use. If NULL, compute based on the input data matrix.
#' @param verbose Whether to show progress bar for calculations
#'
#' @return \code{Palette} object.
#'
#' @import Seurat
#' @importFrom rlang %||%
#'
#' @export
#'
IDFLog.Data <- function(object,
                        modal = NULL,
                        sample = NULL,
                        scale.factor = 10000,
                        idf = NULL,
                        verbose = TRUE){

  if(missing(object) || missing(modal) || length(modal) != 1){
    stop("Must provide a Palette object and a modal information")
  }
  if(!(modal %in% object@Data.index[["assay.idx"]]$modal)){
    stop("Unknown modal information.")
  }
  modal_logic <- object@Data.index[["assay.idx"]]$modal %in% modal

  if(!is.null(sample) && !(sample %in% object@Data.index[["assay.idx"]]$sample)){
    stop("Unknown samople information.")
  }
  sample <- sample %||% unique(object@Data.index[["assay.idx"]]$sample)
  sample_logic <- object@Data.index[["assay.idx"]]$sample %in% sample
  logic <- modal_logic * sample_logic

  if (verbose) {
    message("Performing IDF and log normalization")
  }

  if(is.null(idf)){
    rsums <- list()
    col_sum <- 0
    for(i in 1:length(logic)){
      if(logic[i]){
        rsums <- c(rsums,list(Matrix::rowSums(x = object@Assays[[i]]@Raw)))
        col_sum <- col_sum + ncol(object@Assays[[i]]@Raw)
      }
    }
    rsums.vec <- Reduce("+",rsums)
    idf <- col_sum / rsums.vec
  }else{
    idx <- which(logic)[1]
    if(length(idf) != nrow(x = object@Assays[[idx]]@Raw)){
      stop("Length of supplied IDF vector does not match",
           " number of rows in input matrix")
    }
  }

  for(j in 1:length(logic)){
    if(logic[j]){
      data <- object@Assays[[j]]@Raw
      if (inherits(x = data, what = "data.frame")) {
        data <- as.matrix(x = data)
      }
      if (!inherits(x = data, what = "CsparseMatrix")) {
        data <- as(object = data, Class = "CsparseMatrix")
      }
      norm.data <- Matrix::Diagonal(n = length(x = idf), x = idf)%*%data
      norm.data <- Seurat::NormalizeData(norm.data,
                                         normalization.method = "LogNormalize",
                                         scale.factor = scale.factor,
                                         margin = 1,
                                         block.size = NULL,
                                         verbose = FALSE)
      colnames(x = norm.data) <- colnames(x = data)
      rownames(x = norm.data) <- rownames(x = data)
      # set NA values to 0
      vals <- slot(object = norm.data, name = "x")
      vals[is.na(x = vals)] <- 0
      slot(object = norm.data, name = "x") <- vals
      object@Assays[[j]]@data <- norm.data
    }
  }
  return(object)
}

#' Identifies highly variable features (HVFs).
#'
#' @param object A \code{Palette} object.
#' @param modal A character of modality.
#' @param nfeatures The number of selected features in each sample/batch.
#' @param verbose Whether to show progress bar for calculations
#'
#' @return \code{Palette} object.
#'
#' @import Seurat
#'
#' @export
#'
Find.HVFs <- function(object,
                      modal,
                      nfeatures = 2000,
                      verbose = FALSE){
  idx <- which(object@Data.index[["assay.idx"]]$modal == modal)
  HVF.set <- c()
  for(i in idx){
    HVFs <- rownames(Seurat::FindVariableFeatures(object@Assays[[i]]@Raw,
                                  nfeatures = nfeatures,
                                  verbose = verbose))
    HVF.set <- c(HVF.set, HVFs)
    rm(HVFs)
  }
  co_feature <- unique(HVF.set)
  object@Co.Features[[modal]] <- co_feature

  return(object)
}

#' Get the feature name information of a modality
#'
#' @param object A \code{Palette} object.
#' @param modal A character of modality.
#' @param sample A character of sample/batch.
#' @param check Whether to traverse all batches contained in the modality and take the feature union.
#'
#' @return \code{Palette} object.
#' @keywords internal
Get.Features <- function(object = NULL,
                         modal = NULL,
                         sample = NULL,
                         check = FALSE){
  if(is.null(object) || is.null(modal)){
    stop("Must provide Palette object and modal information")
  }else{
    if(is.null(sample)){
      modal_logic <- ifelse(object@Data.index[["assay.idx"]]$modal == modal, TRUE, FALSE)
      if(check){
        feature.set <- c()
        for(i in 1:length(modal_logic)){
          if(modal_logic[i]){
            feature <- rownames(object@Assays[[i]]@Raw)
            feature.set <- c(feature.set, feature)
            rm(feature)
          }
        }
        out <- unique(feature.set)
      }else{
        idx <- which(modal_logic == 1)[1]
        out <- rownames(object@Assays[[idx]]@Raw)
      }
    }else{
      modal_logic <- ifelse(object@Data.index[["assay.idx"]]$modal == modal, TRUE, FALSE)
      sample_logic <- ifelse(object@Data.index[["assay.idx"]]$sample == sample, TRUE, FALSE)
      idx <- which(modal_logic*sample_logic == 1)
      if(length(idx) != 1){
        stop("Please check the input modal and sample information")
      }else{
        out <- rownames(object@Assays[[idx]]@Raw)
      }
    }
  }
  return(out)
}

#' Specify a set of features as HVFs for a modality in Palette object.
#'
#' @param object A \code{Palette} object.
#' @param modal A character of modality.
#' @param features A character vector specifying the feature names to use as the HVF.If \code{features = NULL}, all shared features of the modality are used as HVFs.
#'
#' @return A character vector of feature names.
#'
#' @export
#'
Add.HVFs <- function(object,
                     modal = NULL,
                     features = NULL){
  if(missing(object) || is.null(modal)){
    stop("Must specify an object and provide modality information")
  }
  available.feautres <- Get.Features(object = object,
                                     modal = modal,
                                     check = TRUE)
  if(is.null(features)){
    message("No features provide, use all shared features by default.")
    input.features <- available.feautres
  }else{
    input.features <- intersect(features, available.feautres)
    if(length(input.features)==0){
      stop("Provided features are not included in this modality.")
    }
  }
  object@Co.Features[[modal]] <- input.features
  return(object)
}

#' This function performs multi-dataset singular value decomposition (SVD) to obtain a shared low-dimensional representation, using either covariance aggregation or concatenation strategies depending on data size.
#'
#' @param x A list of matrices to be jointly decomposed.
#' @param n.dims Number of dimensions (singular vectors) to retain in the reduced representation.
#' @param thresh Row-size threshold determining whether to use standard covariance-based SVD or large-scale matrix methods.
#'
#' @import bigstatsr
#' @importFrom irlba irlba
#' @keywords internal
Run_Multi_SVD <- function(x, n.dims = 50, thresh = 10000L){

  if(nrow(x[[1]]) < thresh){

    d <- lapply(x, function(a){
      a <- tcrossprod(a)/(ncol(a)-1)
      a
    })

    d <- Reduce("+", d)
    V <- irlba::irlba(d, nv = n.dims, work = min(1000,n.dims*3))[["u"]] %>% as.matrix()
    r <- lapply(x, function(b){
      crossprod(V, b)
    })
    return(r)

  }else{
        d <- lapply(x, function(a){
          a <- a/sqrt(ncol(a)-1)
          a
        })
        d <- Reduce(cbind, d)
        V <- RSpectra::svds(d, k = n.dims)[["u"]] %>% as.matrix()
        r <- lapply(x, function(b){
          crossprod(V, b)
        })
        return(r)
    }
}

#' This function performs singular value decomposition (SVD) across multiple matrices to generate a shared low-dimensional representation and projections.
#'
#' @param x A list of matrices to be jointly decomposed.
#' @param n.dims Number of singular vectors (dimensions) to retain.
#' @param thresh Row-size threshold that determines whether covariance-based SVD or concatenated-matrix SVD is applied.
#'
#' @importFrom RSpectra svds
#' @keywords internal
Run_Multi_SVD_s <- function(x, n.dims = 50, thresh = 10000L){

  if(nrow(x[[1]]) < thresh){

    d <- lapply(x, function(a){
      a <- tcrossprod(a)/(ncol(a)-1)
      a
    })

    d <- Reduce("+", d)
    V <- RSpectra::svds(d, k = n.dims)[["u"]] %>% as.matrix()
    r <- lapply(x, function(b){
      crossprod(V, b)
    })
    return(list(r,V))

  }else{
    d <- lapply(x, function(a){
      a <- a/sqrt(ncol(a)-1)
      a
    })
    d <- Reduce(cbind, d)
    V <- RSpectra::svds(d, k = n.dims)[["u"]] %>% as.matrix()
    r <- lapply(x, function(b){
      crossprod(V, b)
    })
    return(list(r,V))
  }
}


# FastSparseRowScaleWithKnownStats <- function(mat, mu, sigma, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
#   .Call('_Seurat_FastSparseRowScaleWithKnownStats', PACKAGE = 'Seurat', mat, mu, sigma, scale, center, scale_max, display_progress)
# }

#' This function dentifies cell clusters within a Palette object using dimension reduction (PCA, LSI or SVD) and neighbor graph–based clustering, and optionally selects representative cells from each cluster based on silhouette scores and quantile thresholds.
#'
#' @param object A Palette object containing assay data.
#' @param modal Modalities to be clustered (e.g., RNA, ATAC, ADT).
#' @param sample Specific sample(s) to include in clustering.
#' @param dim.reduce Whether to perform clustering on reduced-dimensional data.
#' @param dims Dimensions to use for clustering.
#' @param method Dimensionality reduction method ('PCA', 'LSI' or 'SVD').
#' @param joint Whether to cluster in a shared low-dimensional space across modalities.
#' @param resolution Resolution parameter controlling clustering granularity.
#' @param nn.k Number of neighbors for kNN/SNN graph construction.
#' @param nn.method Neighbor search method (default "annoy").
#' @param annoy.metric Distance metric for annoy-based neighbor search.
#' @param prune.SNN Pruning parameter for shared nearest-neighbor graph.
#' @param downsample Whether to select representative cells within each cluster.
#' @param quantile.cutoff Quantile threshold for selecting representative cells.
#' @param min.cell Minimum number of cells to retain per cluster.
#' @param verbose Whether to print progress messages.
#'
#' @import Matrix
#' @import Seurat
#' @import bigstatsr
#' @importFrom rlang %||%
#' @import irlba
#' @importFrom Rcpp evalCpp
#' @importFrom stats quantile
#'
#' @export
#'

Find.Cluster <- function(object,
                         modal = NULL,
                         sample = NULL,
                         dim.reduce = TRUE,
                         dims = list(1:50),
                         method = c("PCA","LSI",'SVD'),
                         joint = FALSE,
                         resolution = 1,
                         nn.k = 20,
                         nn.method = "annoy",
                         annoy.metric = "euclidean",
                         prune.SNN = 1/15,
                         downsample = TRUE,
                         quantile.cutoff = 0.75,
                         min.cell = 20,
                         verbose = TRUE){
  if(missing(object)){
    stop("Palette object missing.")
  }

  modal <- modal %||% unique(object@Data.index[["assay.idx"]]$modal)
  sample <- sample %||% unique(object@Data.index[["assay.idx"]]$sample)

  if(length(resolution) != length(modal)){
    resolution <- rep(resolution[1],length(modal))
  }

  if(dim.reduce){

    if(length(setdiff(method,c("PCA","LSI",'SVD')))){
      stop("Unknown dimensionality reduction method.")
    }

    if(!(inherits(x = dims, what = 'list'))){
      dims <- list(dims)
    }

    if(length(dims) != length(modal)){
      dims <- lapply(1:length(modal),function(x){
        x <- dims[[1]]
      })
    }
    names(dims) <- modal

    if(length(method) == 1){
      method <- rep(method,length(modal))
    }else if(length(method) != length(modal)){
      stop("Unmatched length between parameters method and modal.")
    }
    names(method) <- modal
  }

  modal_logic <- object@Data.index[["assay.idx"]]$modal %in% modal
  sample_logic <- object@Data.index[["assay.idx"]]$sample %in% sample
  logic <- modal_logic * sample_logic

  if(Reduce("+",logic) == 0){
    stop("Palette.Assay Not found")
  }

  cluster.name <- names(object@Assays)
  for(i in 1:length(modal)){

    idx_i <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
    cluster.name_i <- cluster.name[idx_i]
    features <- object@Co.Features[[modal[i]]]
    if(is.null(features)){
      stop("HVFs Not found!")
    }
    data.list_i <- lapply(idx_i, function(x) object@Assays[[x]]@data[features,,drop = FALSE])
    names(data.list_i) <- cluster.name_i
    if(dim.reduce){

      method_i <- method[i]
      dim_i <- dims[[i]]

      if(method_i == "PCA"){
        scale.list <- lapply(data.list_i,function(x){
          x <- Seurat::ScaleData(object = x,
                          features = features,
                          do.scale = TRUE,
                          do.center = TRUE,
                          scale.max = 10,
                          block.size = 1000,
                          min.cells.to.block = 3000,
                          verbose = FALSE)
          x
        })
        names(scale.list) <- cluster.name_i
        if(joint){
          input.data <- Run_Multi_SVD(scale.list,
                                      max(dim_i))
        }else{
          input.data <- lapply(1:length(scale.list),function(x){
            V <- irlba::irlba(scale.list[[x]], nv = max(dim_i), work = max(dim_i)*3)[["u"]] %>% as.matrix()
            x <- crossprod(V, scale.list[[x]])
            x
          })
        }
        input.data <- lapply(input.data, function(x) x[dim_i,,drop = FALSE])

        names(input.data) <- cluster.name_i
        rm(data.list_i,scale.list)
        gc(verbose = FALSE)

      }else if(method_i %in% c("LSI",'SVD')){

        if(joint){
          input.data <- Run_Multi_SVD(data.list_i,
                                      max(dim_i))
        }else{

          if(nrow(data.list_i[[1]]) > 10000L){
            input.data <- lapply(1:length(data.list_i),function(x){
              d <- bigstatsr::FBM(nrow(data.list_i[[x]]),
                                  ncol(data.list_i[[x]]),
                                  init = 0)
              ind_nozero <- which(data.list_i[[x]] != 0, arr.ind = TRUE)
              d[ind_nozero] <- data.list_i[[x]][ind_nozero]
              V <- bigstatsr::big_SVD(d,k = max(dim_i))[["u"]] %>% as.matrix()
              x <- crossprod(V, data.list_i[[x]])
              x
            })

          }else{
            input.data <- lapply(1:length(data.list_i),function(x){
              V <- irlba::irlba(data.list_i[[x]], nv = max(dim_i), work = max(dim_i)*3)[["u"]] %>% as.matrix()
              x <- crossprod(V, data.list_i[[x]])
              x
            })
          }

        }
        input.data <- lapply(input.data, function(x) x[dim_i,,drop = FALSE])
        names(input.data) <- cluster.name_i
        rm(data.list_i)
        gc(verbose = FALSE)
      }

    }else{
      input.data <- data.list_i
      rm(data.list_i)
      gc(verbose = FALSE)
    }

    out_i <- list()
    for(j in 1:length(input.data)){
      snn <- Seurat::FindNeighbors(t(input.data[[j]]),
                            k.param = nn.k,
                            nn.method = nn.method,
                            prune.SNN = prune.SNN,
                            annoy.metric = annoy.metric,
                            verbose = verbose)$snn
      out_i[[j]] <- safe_find_clusters(snn,
                                  resolution = resolution[i],
                                  verbose = verbose)

      Cluster <- as.character(out_i[[j]][,1])

      if(downsample == TRUE){
        ASW_score <- silhouette_cpp(as.numeric(out_i[[j]][,1]),as.matrix(input.data[[j]]))
        Cluster.uni <- as.character(unique(Cluster))
        flag <- c()
        if(quantile.cutoff > 0){
          for(jj in 1:length(Cluster.uni)){
            Index <- which(Cluster == Cluster.uni[jj])
            asw_jj <- ASW_score[Index]
            if(length(asw_jj) < min.cell){
             flag <- c(flag,Index)
            }else{
              thresh <- stats::quantile(asw_jj,quantile.cutoff)
              index <- which(asw_jj > thresh)
              if (length(index) > min.cell) {
                flag <- c(flag, Index[index])
              }else{
                index <- rank(asw_jj,ties.method = "max")[seq(min.cell)]
                flag <- c(flag, Index[index])
              }
            }
          }
        }else{
          flag <- seq(ncol(input.data[[j]]))
        }

        object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["CellName"]] <- colnames(input.data[[j]])[flag]
        object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["Cluster"]] <- Cluster[flag]

      }else{
        flag <- seq(ncol(input.data[[j]]))

        object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["CellName"]] <- colnames(input.data[[j]])[flag]
        object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["Cluster"]] <- Cluster[flag]
      }
    }
    names(out_i) <- cluster.name_i
    current.cluster <- names(object@pre.clusters)
    if(!is.null(current.cluster) && length(current.cluster) > 0){
      replace.idx <- which(current.cluster %in% cluster.name_i)
      if(length(replace.idx) > 0){
        object@pre.clusters <- object@pre.clusters[-replace.idx]
      }
    }
    object@pre.clusters <- c(object@pre.clusters,out_i)
  }

  return(object)
}

clust_labels_to_int = function(clust_labels) {
  uniq_labels = unique(clust_labels)
  clust_labels_int = rep(NA, length(clust_labels))

  for (ii in 1:length(clust_labels)) {
    clust_labels_int[ii] = match(clust_labels[ii], uniq_labels)
  }

  return(clust_labels_int)
}

#' This function DownSample selects representative cells from clusters within a Palette object, either in a supervised manner (using known group labels from metadata) or unsupervised (using pre-computed clusters), optionally leveraging dimensionality reduction (PCA, LSI or SVD) and silhouette-based filtering.
#'
#' @param object A Palette object containing assay data, clustering results, and metadata.
#' @param modal Modalities to include (e.g., RNA, ATAC, ADT).
#' @param sample Specific sample(s) to include from the object.
#' @param supervised Logical; whether to perform supervised downsampling using predefined groups.
#' @param group Metadata column name specifying labels for supervised sampling.
#' @param dim.reduce Logical; whether to perform sampling in reduced-dimensional space.
#' @param dims Dimensions to use for representation (can be a list per modality).
#' @param method Dimensionality reduction method; "PCA" for RNA/ADT, "LSI" or "SVD" for ATAC.
#' @param joint Logical; whether to compute clustering in a shared reduced space across modalities.
#' @param quantile.cutoff Quantile threshold for selecting top-scoring representative cells within clusters.
#' @param min.cell Minimum number of representative cells to retain per cluster.
#'
#' @import Matrix
#' @import Seurat
#' @importFrom rlang %||%
#' @importFrom irlba irlba
#' @importFrom RSpectra svds
#' @importFrom Rcpp evalCpp
#' @importFrom stats quantile
#'
#' @export
#'
DownSample <- function(object,
                       modal = NULL,
                       sample = NULL,
                       supervised = TRUE,
                       group = NULL,
                       dim.reduce = TRUE,
                       dims = list(1:50),
                       method = c("PCA","LSI",'SVD'),
                       joint = FALSE,
                       quantile.cutoff = 0.75,
                       min.cell = 20){

  if(missing(object)){
    stop("Palette object missing.")
  }

  if(supervised && is.null(group)){
    stop('Must selcet meta data under supervised case.')
  }
  if(supervised){
    if(length(intersect(group,colnames(object@Meta))) == 0){
      stop('Group information Not found.')
    }
  }

  modal <- modal %||% unique(object@Data.index[["assay.idx"]]$modal)
  sample <- sample %||% unique(object@Data.index[["assay.idx"]]$sample)

  if(dim.reduce){

    if(length(setdiff(method,c("PCA","LSI",'SVD')))){
      stop("Unknown dimensionality reduction method.")
    }

    if(!(inherits(x = dims, what = 'list'))){
      dims <- list(dims)
    }

    if(length(dims) != length(modal)){
      dims <- lapply(1:length(modal),function(x){
        x <- dims[[1]]
      })
    }
    names(dims) <- modal

    if(length(method) == 1){
      method <- rep(method,length(modal))
    }else if(length(method) != length(modal)){
      stop("Unmatched length between 'method' and 'modal'.")
    }
    names(method) <- modal
  }

  modal_logic <- object@Data.index[["assay.idx"]]$modal %in% modal
  sample_logic <- object@Data.index[["assay.idx"]]$sample %in% sample
  logic <- modal_logic * sample_logic

  if(Reduce("+",logic) == 0){
    stop("Palette.Assay Not found!")
  }

  if(supervised){
    meta_cluster <- as.character(clust_labels_to_int(as.character(object@Meta[,group])))
    idx <- which(!is.na(meta_cluster))
    meta_cluster <- meta_cluster[idx]
    meta_cluster <- factor(meta_cluster,
                           levels = as.character(1:length(unique(meta_cluster))))
    names(meta_cluster) <- object@Meta[,'cell.name'][idx]
    for(i in 1:length(modal)){
      idx_i <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
      for(j in idx_i){
        clust_name <- names(object@Assays)[j]
        tmp_name <- colnames(object@Assays[[j]]@Raw)
        tmp_name <- intersect(tmp_name,names(meta_cluster))

        tmp_clust <- as.data.frame(meta_cluster[tmp_name])
        rownames(tmp_clust) <- tmp_name
        object@pre.clusters[[clust_name]] <- tmp_clust
      }
    }
  }else{
    for(i in 1:length(modal)){
      idx_i <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
      for(j in idx_i){
        clust_name <- paste0(modal[i],'_',object@Data.index[["assay.idx"]]$sample[j])
        if(is.null(object@pre.clusters[[clust_name]])){
          stop('Clustering information must be provided.')
        }
      }
    }
  }

  cluster.name <- names(object@Assays)
  for(i in 1:length(modal)){

    idx_i <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
    cluster.name_i <- cluster.name[idx_i]
    features <- object@Co.Features[[modal[i]]]
    if(is.null(features)){
      stop("HVFs not find!")
    }
    data.list_i <- lapply(idx_i, function(x) object@Assays[[x]]@data[features,])

    idx <- which(names(object@pre.clusters) %in% cluster.name_i)
    idx <- which(names(object@pre.clusters) %in% cluster.name_i)
    Cluster_list <- lapply(idx,function(x){
      x <- object@pre.clusters[[x]]
      x
    })
    names(Cluster_list) <- names(object@pre.clusters)[idx]
    Cluster_list <- Cluster_list[cluster.name_i]

    names(data.list_i) <- cluster.name_i
    if(dim.reduce){

      method_i <- method[i]
      dim_i <- dims[[i]]

      if(method_i == "PCA"){
        scale.list <- lapply(data.list_i,function(x){
          x <- Seurat::ScaleData(object = x,
                          features = features,
                          do.scale = TRUE,
                          do.center = TRUE,
                          scale.max = 10,
                          block.size = 1000,
                          min.cells.to.block = 3000,
                          verbose = FALSE)
          x
        })
        names(scale.list) <- cluster.name_i
        if(joint){
          input.data <- Run_Multi_SVD(scale.list,
                                      max(dim_i))
        }else{
          input.data <- lapply(1:length(scale.list),function(x){
            V <- irlba::irlba(scale.list[[x]], nv = max(dim_i), work = max(dim_i)*3)[["u"]] %>% as.matrix()
            x <- crossprod(V, scale.list[[x]])
            x
          })
        }
        input.data <- lapply(input.data, function(x) x[dim_i,])

        names(input.data) <- cluster.name_i
        rm(data.list_i,scale.list)
        gc(verbose = FALSE)

      }else if(method_i %in% c("LSI",'SVD')){

        if(joint){
          input.data <- Run_Multi_SVD(data.list_i,
                                      max(dim_i))
        }else{

          input.data <- lapply(1:length(data.list_i),function(x){
            V <- RSpectra::svds(data.list_i[[x]], k = max(dim_i))[["u"]] %>% as.matrix()
            x <- crossprod(V, data.list_i[[x]])
            x
          })

        }
        input.data <- lapply(input.data, function(x) x[dim_i,])
        names(input.data) <- cluster.name_i
        rm(data.list_i)
        gc(verbose = FALSE)
      }

    }else{
      input.data <- data.list_i
      rm(data.list_i)
      gc(verbose = FALSE)
    }

    for(j in 1:length(input.data)){

      Cluster <- as.character(Cluster_list[[j]][,1])
      ASW_score <- silhouette_cpp(as.numeric(Cluster_list[[j]][,1]),as.matrix(input.data[[j]]))
      Cluster.uni <- as.character(unique(Cluster))
      flag <- c()
      if(quantile.cutoff > 0){
        for(jj in 1:length(Cluster.uni)){
          Index <- which(Cluster == Cluster.uni[jj])
          asw_jj <- ASW_score[Index]
          if(length(asw_jj) < min.cell){
            flag <- c(flag,Index)
          }else{
            thresh <- stats::quantile(asw_jj,quantile.cutoff)
            index <- which(asw_jj > thresh)
            if (length(index) > min.cell) {
              flag <- c(flag, Index[index])
            }else{
              index <- rank(asw_jj,ties.method = "max")[seq(min.cell)]
              flag <- c(flag, Index[index])
            }
          }
        }
      }else{
        flag <- seq(ncol(input.data[[j]]))
      }

      object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["CellName"]] <- colnames(input.data[[j]])[flag]
      object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["Cluster"]] <- Cluster[flag]
    }
    current.cluster <- names(object@pre.clusters)
    if(!is.null(current.cluster) && length(current.cluster) > 0){
      replace.idx <- which(current.cluster %in% cluster.name_i)
      if(length(replace.idx) > 0){
        object@pre.clusters <- object@pre.clusters[-replace.idx]
      }
    }
    object@pre.clusters <- c(object@pre.clusters,Cluster_list)
  }

  return(object)

}

#' This function constructs a global kernel (similarity) matrix across batches or modalities by integrating cluster structures, using dimensionality reduction, L2 normalization, and cosine similarity–based filtering to identify reliable inter-cluster links.
#'
#' @param data.list A list of matrices containing input data from different batches/modalities.
#' @param cluster A list of clustering assignments corresponding to \code{data.list}.
#' @param dim.reduce Logical; whether to perform dimensionality reduction before kernel construction.
#' @param cos.dims Dimensions to use for cosine similarity calculation.
#' @param Angle.var Angular variance (in degrees) used as a threshold for filtering weak links.
#' @param max.Angle Maximum angular threshold (in degrees) for acceptable inter-cluster similarity.
#' @param joint Logical; whether to perform joint pairwise dimensionality reduction across datasets.
#' @param seed Random seed to ensure reproducibility.
#'
#' @import Matrix
#' @import Seurat
#' @importFrom rlang %||%
#' @importFrom irlba irlba
#' @importFrom RSpectra svds
#' @importFrom Rcpp evalCpp
#' @importFrom stats quantile
#' @keywords internal
Build.Kernel.cluster <- function(data.list,
                                  cluster,
                                  dim.reduce = TRUE,
                                  cos.dims = NULL,
                                  Angle.var = 10,
                                  max.Angle = 35,
                                  joint = TRUE,
                                  seed = 123L){

  set.seed(seed)
  cluster.uni <- list()
  idx <- list()
  for(i in 1:length(cluster)){
    cluster_i <- cluster[[i]][,1]
    clusteri.uni <- as.character(unique(cluster_i))
    ni <- length(clusteri.uni)
    idxi <- lapply(1:ni,function(x) which(cluster_i == clusteri.uni[x]))
    cluster.uni[[i]] <- clusteri.uni
    idx[[i]] <- idxi
  }

  Angle.var <- Angle.var*pi/180
  max.Angle <- max.Angle*pi/180
  gc(verbose = F)

  if(joint && dim.reduce){
    message('Performing pair-wise dimensionality reduction and L2 normalization, and identifying batch-specific cluster median cosine similarity. \n')
    subKernel <- list()
    scaler <- c()
    for(i in 1:(length(cluster)-1)){
      for(j in (i+1):length(cluster)){
        tmp.list <- data.list[c(i,j)]
        svd.list <- Run_Multi_SVD(tmp.list, max(cos.dims))
        svd.list <- lapply(svd.list, function(x) as.matrix(x[cos.dims,]))
        # L2.list <- lapply(svd.list, function(x) L2Norm_seurat(t(x)))
        L2.list <- lapply(svd.list, function(x) t(cosineNorm(x)))
        rm(svd.list)

        tmpi <- FindMedian_self(tcrossprod(L2.list[[1]]),idx[[i]])
        rownames(tmpi) <- colnames(tmpi) <- cluster.uni[[i]]

        tmpj <- FindMedian_self(tcrossprod(L2.list[[2]]),idx[[j]])
        rownames(tmpj) <- colnames(tmpj) <- cluster.uni[[j]]
        gc(verbose = F)

        tmp <- FindMedian(L2.list[[1]]%*%t(L2.list[[2]]),idx[[i]],idx[[j]])
        filterLink <- FindGroupLink(tmp,
                                    tmpi,
                                    tmpj,
                                    Angle.var,
                                    max.Angle)

        link <- PairwiseKernel_norm(L2.list[[1]]%*%t(L2.list[[2]]),
                                    filterLink,
                                    idx[[i]],
                                    idx[[j]],
                                    seed = seed)

        subKernel <- c(subKernel,link[['subK']])
        scaler <- c(scaler,link[["scaler"]])

      }
    }
    message('Done !')
  }else{
    if(dim.reduce){
      message("Performing dimensionality reduction and L2 normalization. \n")
      svd.list <- Run_Multi_SVD(data.list, max(cos.dims))
      svd.list <- lapply(svd.list, function(x) as.matrix(x[cos.dims,]))
      # L2.list <- lapply(svd.list, function(x) L2Norm_seurat(t(x)))
      L2.list <- lapply(svd.list, function(x) t(cosineNorm(x)))
      message("Dimensionality reduction and L2 normalization done! \n")
      rm(svd.list)
    }else{
      message("Performing L2 normalization. \n")
      # L2.list <- lapply(data.list, function(x) L2Norm_seurat(t(x)))
      L2.list <- lapply(svd.list, function(x) t(cosineNorm(x)))
      message("L2 normalization done! \n")
    }
    gc(verbose = F)
    self.sim <- list()

    message("Identifying batch-specific cluster median cosine similarity. \n")
    for(i in 1:length(cluster)){
      tmp <- FindMedian_self(tcrossprod(L2.list[[i]]),idx[[i]])
      rownames(tmp) <- colnames(tmp) <- cluster.uni[[i]]
      self.sim[[names(data.list)[i]]] <- tmp
    }
    gc(verbose = F)

    subKernel <- list()
    scaler <- c()
    for(i in 1:(length(cluster)-1)){
      for(j in (i+1):length(cluster)){
        tmp <- FindMedian(L2.list[[i]]%*%t(L2.list[[j]]),idx[[i]],idx[[j]])
        filterLink <- FindGroupLink(tmp,
                                    self.sim[[i]],
                                    self.sim[[j]],
                                    Angle.var,
                                    max.Angle)
        link <- PairwiseKernel_norm(L2.list[[i]]%*%t(L2.list[[j]]),
                                    filterLink,
                                    idx[[i]],
                                    idx[[j]],
                                    seed = seed)

        subKernel <- c(subKernel,link[['subK']])
        scaler <- c(scaler,link[["scaler"]])
      }
    }
    message('Done !')
  }

  gc(verbose = F)
  max_s <- max(scaler)
  N <- Reduce("+",lapply(data.list,ncol))
  ZK <- Matrix::Matrix(0, nrow = N, ncol = N)

  K <- Insert_submat(ZK,subKernel,sapply(data.list,ncol),scale = T,scaler = max_s)
  K <- Matrix::drop0(K)

  Name.list <- lapply(data.list,colnames)
  TName <- Reduce("c",Name.list)
  colnames(K) <- rownames(K) <- TName
  gc(verbose = F)

  return(K)
}

#' Build simplified kernel matrix across clusters
#'
#' @description
#' `Build.Kernel.cluster_S` constructs a simplified global kernel matrix (`K`) and
#' corresponding block matrix (`B`) across clusters, based on input data matrices
#' and cluster assignments, using normalized kernel construction.
#'
#' @param data.list A list of matrices containing input data from different batches or modalities.
#' @param cluster A list of clustering assignments corresponding to each dataset in `data.list`.
#' @param seed An integer random seed to ensure reproducibility (default: 123L).
#'
#' @return A list containing two sparse matrices:
#' \item{K}{The constructed global kernel (similarity) matrix.}
#' \item{B}{The corresponding block (batch-specific) kernel matrix.}
#' @keywords internal
Build.Kernel.cluster_S <- function(data.list,
                                    cluster,
                                    seed = 123L){

  set.seed(seed)
  cluster.uni <- unique(unlist(cluster))
  nn <- length(cluster.uni)
  idx <- list()
  for(i in 1:length(cluster)){
    cluster_i <- cluster[[i]][,1]
    idxi <- lapply(1:nn,function(x){
      tmp_idx <- which(cluster_i == cluster.uni[x])
      if(length(tmp_idx) == 0){
        tmp_idx <- -1
      }
      x <- tmp_idx
      x
    })
    idx[[i]] <- idxi
  }

  N_list <- sapply(data.list,ncol)
  N <- Reduce("+",N_list)
  K_list  <- sKernel_norm(N_list,
                     N,nn,
                     idx)
  K <- K_list[['K']]
  B <- K_list[['B']]
  rm(K_list)
  K <- Matrix::drop0(K)
  B <- Matrix::drop0(B)

  Name.list <- lapply(data.list,colnames)
  TName <- Reduce("c",Name.list)
  colnames(K) <- rownames(K) <- TName
  colnames(B) <- rownames(B) <- TName

  return(list(K,B))
}

#' Find subspace embedding for a given modality
#'
#' @description
#' `Find_SubSpace` computes a low-dimensional subspace representation for cells
#' in a specified modality by constructing a kernel matrix, performing dimensionality
#' reduction or random projection, and extracting subspace embeddings.
#'
#' @param object An analysis object containing assay data, metadata, and representative cells.
#' @param modal A character string specifying which modality to process.
#' @param sample A character string specifying which sample(s) to include (default: NULL, use all).
#' @param lambda_i A numeric regularization parameter applied when adjusting the kernel diagonal.
#' @param dim.reduce Logical; whether to perform dimensionality reduction (default: TRUE).
#' @param sub_dims Integer vector specifying which subspace dimensions to retain.
#' @param cos_dims Integer vector specifying the dimensions used for cosine similarity computation.
#' @param Angle_var Numeric; variance angle (in degrees) used for filtering links (default: 10).
#' @param max_Angle Numeric; maximum allowable angle (in degrees) for filtering links (default: 35).
#' @param joint Logical; whether to compute joint SVD across batches (default: TRUE).
#' @param proj_scale Logical; whether to scale features during projection (default: TRUE).
#' @param proj_center Logical; whether to center features during projection (default: TRUE).
#' @param L2.norm Logical; whether to apply L2 normalization to the final embedding (default: TRUE).
#' @param seed Integer random seed for reproducibility (default: 123L).
#' @param pre.reduce Logical; whether to perform pre-dimensionality reduction before random projection (default: TRUE).
#' @param pre_dims List of integer vectors specifying dimensions to keep in pre-reduction (default: list(1:300)).
#' @param RP_thresh Integer; number of rows above which random projection is applied (default: 10000).
#'
#' @return A list containing:
#' \item{emb}{A matrix of subspace embeddings for cells.}
#' \item{Index}{Indices of the assays used.}
#' \item{Subspace}{The learned subspace basis matrix.}
#' @keywords internal
Find_SubSpace <- function(object,
                          modal = NULL,
                          sample = NULL,
                          lambda_i = NULL,
                          dim.reduce = TRUE,
                          sub_dims = NULL,
                          cos_dims = NULL,
                          Angle_var = 10,
                          max_Angle = 35,
                          joint = TRUE,
                          proj_scale = T,
                          proj_center = T,
                          L2.norm = TRUE,
                          seed = 123L,
                          pre.reduce = TRUE,
                          pre_dims = list(1:300),
                          RP_thresh = 10000){

  set.seed(seed)

  idx <- which(object@Data.index[["assay.idx"]]$modal == modal)
  features <- object@Co.Features[[modal]]
  data.input <- lapply(idx,function(x) object@Assays[[x]]@data[features,object@Assays[[x]]@RepresentativeCell[[1]]])
  cluster.input <- lapply(idx,function(x) as.data.frame(object@Assays[[x]]@RepresentativeCell[[2]]))
  names(data.input) <- names(object@Assays)[idx]

  K <- Build.Kernel.cluster(data.input,
                                  cluster.input,
                                  joint = joint,
                                  dim.reduce = TRUE,
                                  cos.dims = cos_dims,
                                  Angle.var = Angle_var,
                                  max.Angle = max_Angle,
                                  seed = seed)

  data.mat <- Reduce(cbind,lapply(idx,function(x) object@Assays[[x]]@data[features,]))
  if(nrow(data.mat) > RP_thresh){
    if(pre.reduce){
      message(paste0("Performing pre-dimensionality reduction on modal ",modal,". \n"))
      ncells <- sum(unlist(lapply(idx,function(x) ncol(object@Assays[[x]]@data))))
      pre_dims[[1]] <- as.integer(intersect(1:ncells,pre_dims[[1]]))
      data_res <- Run_Multi_SVD(lapply(idx,function(x) object@Assays[[x]]@data[features,]),
                                n.dims = max(pre_dims[[1]]))

      data_res <- as.matrix(Reduce(cbind,data_res))[as.integer(pre_dims[[1]]),]
      colnames(data_res) <- colnames(data.mat)
      rownames(data_res) <- paste0("dim_",seq(nrow(data_res)))
      rm(data.mat)
      gc(verbose = F)
    }else{
      # RP
      message(paste0("Performing random projection on modal ",modal,". \n"))
      RP_mat <- JL_Trans(nrow(data.mat),seed = seed)
      data_FBM <- bigstatsr::FBM(nrow(data.mat), ncol(data.mat), init = 0)
      ind_nozero <- which(data.mat != 0, arr.ind = TRUE)
      data_FBM[ind_nozero] <- data.mat[ind_nozero]
      data_res <- bigstatsr::big_apply(data_FBM,function(X,ind){
        big_cprodMat(RP_mat, X[, ind])
      },a.combine = "cbind",block.size = 1000)[]

      colnames(data_res) <- colnames(data.mat)
      rownames(data_res) <- paste0("dim_",seq(nrow(data_res)))
      rm(data.mat)
      gc(verbose = F)
    }
  }else{
    data_res <- data.mat
    rm(data.mat)
    gc(verbose = F)
  }
  scaled.data <- Seurat::ScaleData(object = data_res,
                            features = rownames(data_res),
                            do.scale = proj_scale,
                            do.center = proj_center,
                            scale.max = 10,
                            block.size = 1000,
                            min.cells.to.block = 3000,
                            verbose = FALSE)

  data.use <- scaled.data[,colnames(K)]
  DVec <- rowSums(K)
  fill <- mean(DVec[which(DVec != 0)])
  Zidx <- which(diag(K) == 0)
  diag(K)[Zidx] <- fill*(1-lambda_i)

  r <- GetS(data.use, K, DVec, lambda_i)
  Subspace <- suppressWarnings(RSpectra::eigs_sym(calAv_cpp,k=max(sub_dims),which = "LA",n = nrow(r),args = r))[["vectors"]] %>% as.matrix() # 使用函数，保证数值稳定
  L_dim <- as.integer(intersect(seq(ncol(Subspace)),sub_dims))
  embedding <- crossprod(as.matrix(scaled.data), Subspace)
  if(L2.norm){
    embedding <- t(cosineNorm(t(embedding)))
  }
  rownames(embedding) <- colnames(scaled.data)
  colnames(embedding) <- paste0("dim_",seq(ncol(embedding)))
  return(list(emb = embedding, Index = idx))
}

#' Find simplified subspace embedding for a given modality
#'
#' @description
#' `Find_SubSpace_s` computes a simplified low-dimensional subspace representation
#' for cells in a given modality by constructing kernel matrices, applying dimensionality
#' reduction or random projection, and extracting subspace embeddings.
#' Compared to `Find_SubSpace`, this version uses `Build.Kernel.cluster_S`
#' and `Run_Multi_SVD_s` for efficiency.
#'
#' @param object An analysis object containing assay data, metadata, and representative cells.
#' @param modal A character string specifying the modality to process.
#' @param sample A character string specifying which sample(s) to include (default: NULL, use all).
#' @param lambda_i A numeric regularization parameter applied during kernel adjustment.
#' @param dim.reduce Logical; whether to perform dimensionality reduction (default: TRUE).
#' @param sub_dims Integer vector specifying which subspace dimensions to retain.
#' @param proj_scale Logical; whether to scale features during projection (default: TRUE).
#' @param proj_center Logical; whether to center features during projection (default: TRUE).
#' @param L2.norm Logical; whether to apply L2 normalization to the final embedding (default: TRUE).
#' @param seed Integer random seed for reproducibility (default: 123L).
#' @param pre.reduce Logical; whether to perform pre-dimensionality reduction before random projection (default: TRUE).
#' @param pre_dims List of integer vectors specifying dimensions to keep in pre-reduction (default: list(1:300)).
#' @param RP_thresh Integer; number of rows above which random projection is applied (default: 10000).
#'
#' @return A list containing:
#' \item{emb}{A matrix of subspace embeddings for cells.}
#' \item{Index}{Indices of the assays used.}
#' \item{Subspace}{The learned subspace basis matrix.}
#' @keywords internal
Find_SubSpace_s <- function(object,
                             modal = NULL,
                             sample = NULL,
                             lambda_i = NULL,
                             dim.reduce = TRUE,
                             sub_dims = NULL,
                             proj_scale = T,
                             proj_center = T,
                             L2.norm = TRUE,
                             seed = 123L,
                             pre.reduce = TRUE,
                             pre_dims = list(1:300),
                             RP_thresh = 10000){

  set.seed(seed)

  idx <- which(object@Data.index[["assay.idx"]]$modal == modal)
  features <- object@Co.Features[[modal]]
  data.input <- lapply(idx,function(x) object@Assays[[x]]@data[features,object@Assays[[x]]@RepresentativeCell[[1]]])
  cluster.input <- lapply(idx,function(x) as.data.frame(object@Assays[[x]]@RepresentativeCell[[2]]))
  names(data.input) <- names(object@Assays)[idx]


  K_list <- Build.Kernel.cluster_S(data.input,
                             cluster.input,
                             seed = seed)

  data.mat <- Reduce(cbind,lapply(idx,function(x) object@Assays[[x]]@data[features,]))
  if(nrow(data.mat) > RP_thresh){
    if(pre.reduce){
      message(paste0("Performing pre-dimensionality reduction on modal ",modal,". \n"))
      ncells <- sum(unlist(lapply(idx,function(x) ncol(object@Assays[[x]]@data))))
      pre_dims[[1]] <- as.integer(intersect(1:ncells,pre_dims[[1]]))
      res <- Run_Multi_SVD_s(lapply(idx,function(x) object@Assays[[x]]@data[features,]),
                                n.dims = max(pre_dims[[1]]))
      data_res <- res[[1]]
      V_pre <- res[[2]]
      rm(res)
      rownames(V_pre) <- features
      data_res <- as.matrix(Reduce(cbind,data_res))[as.integer(pre_dims[[1]]),]
      V_pre <- V_pre[,as.integer(pre_dims[[1]])]
      colnames(data_res) <- colnames(data.mat)
      rownames(data_res) <- paste0("dim_",seq(nrow(data_res)))
      colnames(V_pre) <- paste0("dim_",seq(nrow(data_res)))
      rm(data.mat)
      gc(verbose = F)
    }else{
      # RP
      message(paste0("Performing random projection on modal ",modal,". \n"))
      RP_mat <- JL_Trans(nrow(data.mat),seed = seed)
      data_FBM <- bigstatsr::FBM(nrow(data.mat), ncol(data.mat), init = 0)
      ind_nozero <- which(data.mat != 0, arr.ind = TRUE)
      data_FBM[ind_nozero] <- data.mat[ind_nozero]
      # data_res <- bigstatsr::big_cprodMat(RP_mat,data_FBM)[]
      data_res <- bigstatsr::big_apply(data_FBM,function(X,ind){
        big_cprodMat(RP_mat, X[, ind])
      },a.combine = "cbind",block.size = 1000)[]

      colnames(data_res) <- colnames(data.mat)
      rownames(data_res) <- paste0("dim_",seq(nrow(data_res)))
      rm(data.mat)
      gc(verbose = F)
    }
  }else{
    data_res <- data.mat
    rm(data.mat)
    gc(verbose = F)
  }

  scaled.data <- Seurat::ScaleData(object = data_res,
                            features = rownames(data_res),
                            do.scale = proj_scale,
                            do.center = proj_center,
                            scale.max = 10,
                            block.size = 1000,
                            min.cells.to.block = 3000,
                            verbose = FALSE)

  K <- K_list[[1]]
  B <- K_list[[2]]
  rm(K_list)
  data.use <- scaled.data[,colnames(K)]

  r <- GetS_new(data.use, K, B, lambda_i)
  Subspace <- suppressWarnings(RSpectra::eigs_sym(calAv_cpp,k=max(sub_dims),which = "LA",n = nrow(r),args = r))[["vectors"]] %>% as.matrix() # 使用函数，保证数值稳定
  L_dim <- as.integer(intersect(seq(ncol(Subspace)),sub_dims))
  embedding <- crossprod(as.matrix(scaled.data), Subspace)

  if(L2.norm){
    embedding <- t(cosineNorm(t(embedding)))
  }
  rownames(embedding) <- colnames(scaled.data)
  colnames(embedding) <- paste0("dim_",seq(ncol(embedding)))
  return(list(emb = embedding, Index = idx))
}

#' Find subspace embeddings across modalities
#'
#' @description
#' `Find.Subspace` computes subspace embeddings for one or more modalities in the input object.
#' Depending on the `supervised` flag, it will call either the supervised version
#' or the unsupervised version function. The resulting embeddings are stored in the corresponding assays.
#'
#' @param object An analysis object containing assay data, metadata, and representative cells.
#' @param modal A character vector specifying modalities to process (default: all in object).
#' @param sample A character vector specifying which sample(s) to include (default: all in object).
#' @param lambda A numeric vector of regularization parameters, one per modality (default: 0.9).
#' @param dim.reduce Logical; whether to perform dimensionality reduction before subspace learning (default: TRUE).
#' @param sub.dims A list of integer vectors specifying which subspace dimensions to retain per modality.
#' @param cos.dims A list of integer vectors specifying dimensions used in cosine kernel computation.
#' @param Angle.var A numeric vector controlling angular variance in kernel construction (default: 10).
#' @param max.Angle A numeric vector specifying the maximum angle threshold for kernel construction (default: 35).
#' @param joint Logical; whether to compute joint subspaces across batches (default: TRUE).
#' @param supervised Logical; whether to run supervised version using `Find_SubSpace_s` (default: FALSE).
#' @param proj.scale Logical; whether to scale features during projection (default: TRUE).
#' @param proj.center Logical; whether to center features during projection (default: TRUE).
#' @param L2.norm Logical; whether to apply L2 normalization to the final embedding (default: TRUE).
#' @param seed Integer random seed for reproducibility (default: 123L).
#' @param pre.reduce Logical; whether to perform pre-dimensionality reduction before Bi-sPCA (default: FALSE).
#' @param pre_dims List of integer vectors specifying dimensions to keep in pre-reduction (default: list(1:3000)).
#' @param RP_thresh Integer; number of rows above which random projection is applied (default: 10000).
#'
#' @return The input object with updated subspace embeddings stored in each assay.
#'
#' @import Matrix
#' @import Seurat
#' @importFrom rlang %||%
#' @importFrom irlba irlba
#' @importFrom RSpectra svds
#' @importFrom Rcpp evalCpp
#'
#' @export
#'
Find.Subspace <- function(object,
                          modal = NULL,
                          sample = NULL,
                          lambda = c(0.8),
                          dim.reduce = TRUE,
                          sub.dims = list(1:40),
                          cos.dims = list(1:50),
                          Angle.var = 15,
                          max.Angle = 50,
                          joint = FALSE,
                          supervised = FALSE,
                          proj.scale = TRUE,
                          proj.center = TRUE,
                          L2.norm = FALSE,
                          seed = 123L,
                          pre.reduce = FALSE,
                          pre_dims = list(1:100),
                          RP_thresh = 10000){
  set.seed(seed)

  if(missing(object)){
    stop("Must provide object")
  }
  sample <- sample %||% unique(object@Data.index[["assay.idx"]]$sample)
  modal <- modal %||% unique(object@Data.index[["assay.idx"]]$modal)

  drop.idx <- c()
  for(i in 1:length(modal)){
    tmp_idx <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
    if(length(tmp_idx) == 1){
      drop.idx <- c(drop.idx,i)
    }
  }

  if(length(drop.idx) > 0){
    modal <- modal[-drop.idx]
  }

  if(!(inherits(sub.dims,what = "list"))){
    sub.dims <- list(sub.dims)
  }
  if(!(inherits(cos.dims,what = "list"))){
    cos.dims <- list(cos.dims)
  }

  if(length(dim.reduce) != length(modal)){
    dim.reduce <- rep(dim.reduce[1],length(modal))
  }
  if(length(sub.dims) != length(modal)){
    sub.dims <- replicate(length(modal),sub.dims[[1]], simplify = FALSE)
  }
  if(length(cos.dims) != length(modal)){
    cos.dims <- replicate(length(modal),cos.dims[[1]], simplify = FALSE)
  }
  if(!(inherits(lambda,what = "numeric"))){
    lambda <- as.numeric(lambda)
  }
  if(length(lambda) != length(modal)){
    lambda <- rep(lambda[1],length(modal))
  }

  if(supervised){

    for(i in 1:length(modal)){

      out <- Find_SubSpace_s(object = object,
                            modal = modal[i],
                            sample = sample,
                            lambda_i = lambda[i],
                            dim.reduce = dim.reduce[i],
                            sub_dims = sub.dims[[i]],
                            proj_scale = proj.scale,
                            proj_center = proj.center,
                            L2.norm = L2.norm,
                            seed = seed,
                            pre.reduce = pre.reduce,
                            pre_dims = pre_dims,
                            RP_thresh = RP_thresh)

      embedding <- out[["emb"]]
      idx <- out[["Index"]]

      for(j in 1:length(idx)){
        object@Assays[[idx[j]]]@embedding <- t(embedding[colnames(object@Assays[[idx[j]]]@data),])
      }
    }

  }else{
    if(!(inherits(Angle.var,what = "numeric"))){
      Angle.var <- as.numeric(Angle.var)
    }
    if(!(inherits(max.Angle,what = "numeric"))){
      max.Angle <- as.numeric(max.Angle)
    }
    if(length(Angle.var) != length(modal)){
      Angle.var <- rep(Angle.var[1],length(modal))
    }
    if(length(max.Angle) != length(modal)){
      max.Angle <- rep(max.Angle[1],length(modal))
    }

    for(i in 1:length(modal)){

      out <- Find_SubSpace(object = object,
                            modal = modal[i],
                            sample = sample,
                            lambda_i = lambda[i],
                            dim.reduce = dim.reduce[i],
                            sub_dims = sub.dims[[i]],
                            cos_dims = cos.dims[[i]],
                            Angle_var = Angle.var[i],
                            max_Angle = max.Angle[i],
                            joint = joint,
                            proj_scale = proj.scale,
                            proj_center = proj.center,
                            L2.norm = L2.norm,
                            seed = seed,
                            pre.reduce = pre.reduce,
                            pre_dims = pre_dims,
                            RP_thresh = RP_thresh)

      embedding <- out[["emb"]]
      idx <- out[["Index"]]

      for(j in 1:length(idx)){
        object@Assays[[idx[j]]]@embedding <- t(embedding[colnames(object@Assays[[idx[j]]]@data),])
      }
    }
  }

  return(object)
}

#' Generate a Johnson-Lindenstrauss (JL) random projection matrix
#'
#' This function constructs a random projection matrix using either the
#' "li" method or the "gaussian" method, to reduce dimensionality
#' while approximately preserving pairwise distances.
#'
#' @param Dim Integer. Original data dimensionality.
#' @param eps Numeric. Distortion parameter (0 < eps <= 1) controlling the approximation error.
#' @param max_dim Integer. Maximum allowed projected dimension (default: 15000).
#' @param seed Integer. Random seed for reproducibility (default: NA, no fixed seed).
#' @param method Character. Projection method, one of `"li"` or `"gaussian"`.
#'
#' @return A random projection matrix stored as a `bigstatsr::FBM` object.
#' @importFrom stats rnorm
#' @keywords internal
JL_Trans <- function(Dim,
                     eps = 0.1,
                     max_dim = 15000L,
                     seed = NA_integer_,
                     method = "gaussian"){

  # method should be one of 'li' or 'gaussian'

  if (!is.na(x = seed)) {
    set.seed(seed = seed)
  }
  method <- method[1L]

  if (!is.null(x = eps)) {
    if (eps > 1 || eps <= 0) {
      stop("'eps' must be 0 < eps <= 1")
    }
    ncol <- floor(x = 4 * log(x = Dim) / ((eps ^ 2) / 2 - (eps ^ 3 / 3)))
  }

  ncols <- min(ncol, max_dim)

  if(ncols > Dim){
    message("Using raw dimensions.")
    # random_matrix <- diag(1,Dim)
    random_matrix <- bigstatsr::FBM(Dim,Dim,init = 0)
    diag(random_matrix) <- 1
    return(random_matrix)
  }else{
    if(!(method %in% c("li","gaussian"))){
      stop("Only methods 'li' and 'gaussian' are supported.")
    }

    if(method == "li"){
      s <- ceiling(sqrt(ncols))
      prob <- c((1/(2*s)),(1-(1/s)),(1/(2*s)))
      # random_matrix <- matrix(sample(seq.int(-1L, 1L),size = Dim * ncols,replace = TRUE,prob = prob),nrow = nrow)
      random_matrix <- bigstatsr::FBM(Dim,ncols,init = sample(seq.int(-1L, 1L),size = Dim * ncols,replace = TRUE,prob = prob))
      return(random_matrix)
    }

    if(method == "gaussian"){
      #crandom_matrix <- matrix(rnorm(Dim*ncols,mean=0,sd=(1/sqrt(ncols))),Dim,ncols)
      random_matrix <- bigstatsr::FBM(Dim,ncols,init = rnorm(Dim*ncols,mean=0,sd=(1/sqrt(ncols))))
      return(random_matrix)
    }
  }

}

#' Build a modality-sample linkage graph
#'
#' This function constructs an undirected graph representing connections
#' between modalities and samples based on their co-occurrence in the dataset.
#' It checks whether the resulting network is fully connected.
#'
#' @param object A data object containing `Data.index` with assay modality
#' and sample information.
#'
#' @return An `igraph` object representing the modality-sample network.
#'
#' @import igraph
#' @import Matrix
#' @keywords internal
dataLink <- function(object){

  modal <- object@Data.index[["assay.idx"]]$modal
  sample <- object@Data.index[["assay.idx"]]$sample
  modal.uni <- unique(modal)
  sample.uni <- unique(sample)
  adj_mat <- as(matrix(0,length(modal.uni)+length(sample.uni),length(modal.uni)+length(sample.uni)),"dgCMatrix")
  rownames(adj_mat) <- colnames(adj_mat) <- c(modal.uni,sample.uni)

  for(i in 1:length(modal)){
    adj_mat[modal[i],sample[i]] <- adj_mat[sample[i],modal[i]] <- 1
  }

  g <- graph.adjacency(adj_mat, mode="undirected")
  if (components(g)$no != 1) {
    stop("The MBG is not connected")
  }
  return(g)
}

#' Compute embeddings for isolated modalities
#'
#' This function generates low-dimensional embeddings for isolated modalities
#' using dimensionality reduction methods such as PCA, SVD, or LSI.
#'
#' @param object A data object containing assay and feature information.
#' @param modal Character vector specifying modalities to process;
#' defaults to isolated modalities detected in the object.
#' @param method Character vector indicating the dimensionality reduction
#' method to use (e.g., `"PCA"`, `"SVD"`, `"LSI"`).
#' @param dims List of integer vectors specifying which dimensions to retain
#' for each modality.
#' @param do.scale Logical value (or vector) indicating whether to apply scaling
#' before dimensionality reduction.
#'
#' @return The updated object with embeddings stored in each selected assay.
#'
#' @import Matrix
#' @import Seurat
#' @importFrom rlang %||%
#' @importFrom irlba irlba
#'
#' @export
#'
IsolatedModalEmbedding <- function(object,
                                   modal = NULL,
                                   method = "PCA",
                                   dims = NULL,
                                   do.scale = TRUE){

 modal <- modal %||% check_one(object@Data.index[["assay.idx"]]$modal)
 modal <- check_aga(modal,object@Data.index[["assay.idx"]]$modal)
 if(length(modal) == 0){
   message("There are no isolated modal data.")
   return(object)
 }else{
   if(!(inherits(method,what = "character"))){
     method <- as.character(method)
   }
   if(!(inherits(dims,what = "list"))){
     dims <- list(dims)
   }
   if(length(method) != length(modal)){
     method <- rep(method[1],length(modal))
   }
   if(length(dims) != length(modal)){
     dims <- replicate(length(modal),dims[[1]], simplify = FALSE)
   }
   if(length(do.scale) != length(modal)){
     do.scale <- rep(do.scale[1],length(modal))
   }

   for(i in 1:length(modal)){
     idx <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
     features <- object@Co.Features[[modal[i]]]
     data <- object@Assays[[idx]]@data[features,]
     if(method[i] == "PCA"){
       input <- Seurat::ScaleData(object = data,
                           features = features,
                           do.scale = do.scale[i],
                           do.center = TRUE,
                           scale.max = 10,
                           block.size = 1000,
                           min.cells.to.block = 3000,
                           verbose = FALSE)

     }else if(method[i] %in% c("SVD",'LSI')){

       if(do.scale[i]){
         input <-  Seurat::ScaleData(object = data,
                             features = features,
                             do.scale = do.scale[i],
                             do.center = FALSE,
                             scale.max = 10,
                             block.size = 1000,
                             min.cells.to.block = 3000,
                             verbose = FALSE)
       }else{
         input <- data
       }

     }else{
       stop("Unknown dim reduce method.")
     }
     rm(data)
     V <- irlba::irlba(input, nv = max(dims[[i]]), work = max(dims[[i]])*3)[["u"]] %>% as.matrix()
     V <- V[,dims[[i]]]
     emb <- t(V) %*% as.matrix(input)
     colnames(emb) <- colnames(input)
     rownames(emb) <- colnames(V) <- paste0("dim_",seq(nrow(emb)))
     rownames(V) <- rownames(input)

     object@Assays[[idx]]@embedding <- emb
     gc(verbose = FALSE)
   }
   return(object)
 }
}

#' Run Palette integration framework
#'
#' This function performs cross-modal integration and embedding inference
#' across samples and modalities using the Palette framework.
#'
#' @param object A data object containing assays and embeddings for integration.
#' @param nn Integer. Number of nearest neighbors used in modal inference.
#' @param modal.norm Logical. Whether to normalize embeddings across modalities.
#' @param weight Named list of numeric weights for scaling each modality embedding.
#' @param do.scaler Logical. Whether to apply feature scaling during inference.
#' @param scaler.k Integer. Number of top features used for scaling.
#' @param do.var Logical. Whether to select variable features during inference.
#' @param var.k Integer. Number of variable features to retain.
#' @param lambda Numeric. Regularization parameter for integration.
#' @param nDims Integer. Number of latent dimensions to retain.
#' @param supervised Logical. Whether to use supervised information during integration.
#' @param group Character vector. Optional grouping labels for supervised integration.
#' @param origin.infer Logical. Whether to perform inference in the original feature space.
#' @param center Logical. Whether to center embeddings before scaling.
#' @param scale Logical. Whether to scale embeddings before integration.
#' @param nn.k Integer. Number of nearest neighbors used in the SNN graph.
#' @param k Integer. Number of nearest neighbors for constructing graphs in Bi-sPCA.
# Parameter controlling local connectivity during inference.
#' @param nn.method Character. Method for nearest neighbor search (e.g., `"annoy"`).
#' @param annoy.metric Character. Distance metric for annoy neighbor search.
#' @param L2Norm Logical. Whether to apply L2 normalization to embeddings during MBG-guied indering.
#' @param emb_L2 Logical. Whether to apply L2 normalization on final embeddings.
#' @param sample_ratio Numeric in (0, 1]. Subsampling ratio used for SNN construction.
#' @param snn_full_max_cells Integer. Threshold number of cells for subsampling.
#' @param seed Integer. Random seed for reproducibility.
#'
#'
#' @return An updated object with integrated Palette embeddings.
#' @export
#'
Run.Palette <- function(object,
                        nn = 2,
                        modal.norm = TRUE,
                        weight = NULL,
                        do.scaler = TRUE,
                        scaler.k = 50,
                        do.var = TRUE,
                        var.k = 20,
                        lambda = 0.5,
                        nDims = 20,
                        supervised = FALSE,
                        group = NULL,
                        origin.infer = FALSE,
                        center = TRUE,
                        scale = TRUE,
                        nn.k = 30,
                        k = 5,
                        nn.method = "annoy",
                        annoy.metric = "euclidean",
                        L2Norm = TRUE,
                        emb_L2 = TRUE,
                        sample_ratio = 1,
                        snn_full_max_cells = 200000,
                        seed = 123L){

  update_emb_local <- function(emb_obj, new_cols, new_data) {
    stopifnot(nrow(new_data) == nrow(emb_obj$emb))
    pos <- match(new_cols, emb_obj$cols)

    exist <- which(!is.na(pos))
    if (length(exist) > 0) {
      emb_obj$emb[, pos[exist]] <- new_data[, exist, drop = FALSE]
    }

    add <- which(is.na(pos))
    if (length(add) > 0) {
      emb_obj$emb  <- cbind(emb_obj$emb, new_data[, add, drop = FALSE])
      emb_obj$cols <- c(emb_obj$cols, new_cols[add])
    }

    emb_obj
  }

  align_emb <- function(x, cols_target) {
    out <- matrix(0, nrow = nrow(x$emb), ncol = length(cols_target),
                  dimnames = list(rownames(x$emb), cols_target))
    hit <- match(x$cols, cols_target)
    ok <- which(!is.na(hit))
    if (length(ok) > 0) {
      out[, hit[ok]] <- x$emb[, ok, drop = FALSE]
    }
    out
  }

  message('Constructing MBG.')
  g <- dataLink(object)
  message('MBG Done!')

  modal  <- object@Data.index[["assay.idx"]]$modal
  sample <- object@Data.index[["assay.idx"]]$sample
  modal.uni  <- unique(modal)
  sample.uni <- unique(sample)

  # 预先建立 (sample|modal)->assay index 查表
  assay_idx <- split(seq_along(modal), paste(sample, modal, sep = "|"))

  # 按 sample 收集 cellnames，用于构建 Tname
  CN.list <- list()
  for (s in sample.uni) {
    idx <- which(sample == s)[1]
    CN.list[[s]] <- colnames(object@Assays[[idx]]@data)
  }
  Tname <- Reduce(c, CN.list)

  # ---------------------------
  # 1)  weight_M_list
  # ---------------------------
  weight_M_list  <- list()
  real_modal_mat <- list()
  emb.list       <- list()

  for (m in modal.uni) {
    idx <- which(modal == m)

    if (length(idx) > 1) {
      if (length(idx) <= length(sample.uni)) {
        d <- Reduce(c, lapply(idx, function(x) colnames(object@Assays[[x]]@data)))
      } else {
        next
      }
    } else {
      d <- colnames(object@Assays[[idx]]@data)
    }

    pos <- match(d, Tname)
    M <- Matrix::sparseMatrix(
      i = seq_along(d),
      j = pos,
      x = 1,
      dims = c(length(d), length(Tname)),
      dimnames = list(d, Tname)
    )
    weight_M_list[[m]] <- M

    if (length(idx) > 1) {
      real_modal_mat[[m]] <- do.call(cbind, lapply(idx, function(x) object@Assays[[x]]@embedding))
    } else {
      real_modal_mat[[m]] <- object@Assays[[idx]]@embedding
    }

    emb.list[[m]] <- list(
      emb  = real_modal_mat[[m]],
      cols = colnames(real_modal_mat[[m]])
    )
  }

  # ---------------------------
  # 2) MBG-guided modality inference
  # ---------------------------
  message('Inferring missing modality matrices guided by MBG.\n')
  ss.idx <- c()
  mm.idx <- c()

  for (s in sample.uni) {
    idxs <- which(sample == s)
    loss.modal <- setdiff(modal.uni, modal[idxs])
    if (length(loss.modal) == 0) next

    for (mm in loss.modal) {
      message(paste0("Performing inference for modal \"", mm, "\" of  sample \"", s, "\" in embedding space.\n"))

      paths <- lapply(all_shortest_paths(g, from = s, to = mm)[["res"]], names)

      # -----------------------
      # Case 1: path length == 4
      # -----------------------
      if (length(paths[[1]]) == 4) {

        start.m <- Reduce(c, lapply(paths, function(x) x[2]))
        start.m.uni <- unique(start.m)

        # ----- 1.1 单起点 -----
        if (length(start.m.uni) == 1) {

          idx_start <- assay_idx[[paste(s, start.m.uni, sep = "|")]]
          start.data <- object@Assays[[idx_start]]@embedding
          if (L2Norm) start.data <- cosineNorm(start.data)

          if (length(paths) > 1) {
            through.data <- do.call(cbind, lapply(seq_along(paths), function(p) {
              tmp.s <- paths[[p]][3]
              idx_t <- assay_idx[[paste(tmp.s, start.m.uni, sep = "|")]]
              object@Assays[[idx_t]]@embedding
            }))

            arrive.data <- do.call(cbind, lapply(seq_along(paths), function(p) {
              tmp.s <- paths[[p]][3]
              idx_a <- assay_idx[[paste(tmp.s, mm, sep="|")]]
              object@Assays[[idx_a]]@embedding
            }))
            arrive_name <- colnames(arrive.data)
          } else {
            tmp.s <- paths[[1]][3]
            through.data <- object@Assays[[assay_idx[[paste(tmp.s, start.m.uni, sep="|")]]]]@embedding
            arrive.data  <- object@Assays[[assay_idx[[paste(tmp.s, mm, sep="|")]]]]@embedding
            arrive_name  <- colnames(arrive.data)
          }

          target <- real_modal_mat[[mm]]
          if (L2Norm) through.data <- cosineNorm(through.data)

          infer_result <- modal_infer_one_one(
            start.data, through.data, nn, L2Norm, arrive.data, target,
            do.scaler, scaler.k, do.var, var.k
          )

          weight_M <- infer_result[["M"]]
          colnames(weight_M) <- arrive_name
          rownames(weight_M) <- colnames(start.data)
          weight_M_list[[mm]][arrive_name, colnames(start.data)] <- t(weight_M)

          infer.data <- infer_result[[2]]
          rownames(infer.data) <- rownames(target)
          colnames(infer.data) <- colnames(start.data)

          emb.list[[mm]] <- update_emb_local(
            emb.list[[mm]],
            new_cols = colnames(start.data),
            new_data = infer.data
          )

          ss.idx <- c(ss.idx, s)
          mm.idx <- c(mm.idx, mm)

          if (origin.infer) {
            arrive_o <- do.call(cbind, lapply(seq_along(paths), function(p) {
              tmp.s <- paths[[p]][3]
              idx_o <- assay_idx[[paste(tmp.s, mm, sep="|")]]
              object@Assays[[idx_o]]@data
            }))

            infer_o <- arrive_o %*% t(weight_M)
            rownames(infer_o) <- rownames(arrive_o)
            colnames(infer_o) <- colnames(start.data)

            tmp1 <- Create.Palette.Assay(
              infer_o, filter = FALSE, min.cells = 0, min.features = 0,
              modal = mm, sample = s, slot = "data", sparse = TRUE
            )
            object <- Add.Assay(tmp1, object, check = FALSE, class = "infer")
          }

          rm(start.data, through.data, target, arrive_name, infer_result,
             weight_M, infer.data)
          gc(FALSE)

        } else {
          # ----- 1.2 多起点 -----
          start.list <- lapply(start.m.uni, function(x) {
            idx_start <- assay_idx[[paste(s, x, sep="|")]]
            object@Assays[[idx_start]]@embedding
          })

          mm_path <- Reduce(c, lapply(paths, function(x) x[2]))

          through.list <- lapply(start.m.uni, function(x) {
            idxp <- which(mm_path == x)
            if (length(idxp) == 1) {
              tmp.s <- paths[[idxp]][3]
              idx_t <- assay_idx[[paste(tmp.s, x, sep="|")]]
              object@Assays[[idx_t]]@embedding
            } else {
              do.call(cbind, lapply(idxp, function(p) {
                tmp.s <- paths[[p]][3]
                idx_t <- assay_idx[[paste(tmp.s, x, sep="|")]]
                object@Assays[[idx_t]]@embedding
              }))
            }
          })

          if (L2Norm) {
            start.list   <- lapply(start.list, cosineNorm)
            through.list <- lapply(through.list, cosineNorm)
          }

          arrive.data <- do.call(cbind, lapply(start.m.uni, function(x) {
            idxp <- which(mm_path == x)
            if (length(idxp) == 1) {
              tmp.s <- paths[[idxp]][3]
              idx_a <- assay_idx[[paste(tmp.s, mm, sep="|")]]
              object@Assays[[idx_a]]@embedding
            } else {
              do.call(cbind, lapply(idxp, function(p) {
                tmp.s <- paths[[p]][3]
                idx_a <- assay_idx[[paste(tmp.s, mm, sep="|")]]
                object@Assays[[idx_a]]@embedding
              }))
            }
          }))

          arrive_name <- colnames(arrive.data)
          target <- real_modal_mat[[mm]]

          infer_result <- modal_infer_one_multi(
            start.list, through.list, nn, L2Norm, arrive.data, target,
            do.scaler, scaler.k, do.var, var.k
          )

          weight_M <- infer_result[["M"]]
          colnames(weight_M) <- arrive_name
          rownames(weight_M) <- colnames(start.list[[1]])
          weight_M_list[[mm]][arrive_name, colnames(start.list[[1]])] <- t(weight_M)

          infer.data <- infer_result[[2]]
          rownames(infer.data) <- rownames(target)
          colnames(infer.data) <- colnames(start.list[[1]])

          emb.list[[mm]] <- update_emb_local(
            emb.list[[mm]],
            new_cols = colnames(start.list[[1]]),
            new_data = infer.data
          )

          ss.idx <- c(ss.idx, s)
          mm.idx <- c(mm.idx, mm)

          if (origin.infer) {
            arrive_o <- do.call(cbind, lapply(start.m.uni, function(x) {
              idxp <- which(mm_path == x)
              if (length(idxp) == 1) {
                tmp.s <- paths[[idxp]][3]
                object@Assays[[assay_idx[[paste(tmp.s, mm, sep="|")]]]]@data
              } else {
                do.call(cbind, lapply(idxp, function(p) {
                  tmp.s <- paths[[p]][3]
                  object@Assays[[assay_idx[[paste(tmp.s, mm, sep="|")]]]]@data
                }))
              }
            }))

            infer_o <- arrive_o %*% t(weight_M)
            rownames(infer_o) <- rownames(arrive_o)
            colnames(infer_o) <- colnames(start.list[[1]])

            tmp1 <- Create.Palette.Assay(
              infer_o, filter = FALSE, min.cells = 0, min.features = 0,
              modal = mm, sample = s, slot = "data", sparse = TRUE
            )
            object <- Add.Assay(tmp1, object, check = FALSE, class = "infer")
          }

          rm(start.list, through.list, target, arrive.data, arrive_name,
             infer_result, weight_M, infer.data)
          gc(FALSE)
        }

      } else {
        # -----------------------
        # Case 2: path length != 4
        # -----------------------

        # ----- 2.1 paths 只有一条 -----
        if (length(paths) == 1) {

          idx_start <- assay_idx[[paste(s, paths[[1]][2], sep="|")]]
          start.data <- object@Assays[[idx_start]]@embedding

          l  <- length(paths[[1]])
          l.m <- seq(l / 2) * 2

          road.list <- lapply(1:(length(l.m)-1), function(x) {
            ll <- list()
            tmp_s <- paths[[1]][l.m[x] + 1]
            ll[[1]] <- object@Assays[[assay_idx[[paste(tmp_s, paths[[1]][l.m[x]], sep="|")]]]]@embedding
            ll[[2]] <- object@Assays[[assay_idx[[paste(tmp_s, paths[[1]][l.m[x + 1]], sep="|")]]]]@embedding
            ll
          })

          tmp_s2 <- paths[[1]][l.m[length(l.m)] - 1]
          arrive.data <- object@Assays[[assay_idx[[paste(tmp_s2, paths[[1]][l.m[length(l.m)]], sep="|")]]]]@embedding
          arrive_name <- colnames(arrive.data)

          if (L2Norm) {
            start.data <- cosineNorm(start.data)
            road.list <- lapply(road.list, function(x) {
              lapply(x, cosineNorm)
            })
          }

          target <- real_modal_mat[[mm]]

          infer_result <- modal_infer_multi_one(
            start.data, road.list, nn, L2Norm, arrive.data, target,
            do.scaler, scaler.k, do.var, var.k
          )

          weight_M <- infer_result[["M"]]
          colnames(weight_M) <- arrive_name
          rownames(weight_M) <- colnames(start.data)
          weight_M_list[[mm]][arrive_name, colnames(start.data)] <- t(weight_M)

          infer.data <- infer_result[[2]]
          rownames(infer.data) <- rownames(target)
          colnames(infer.data) <- colnames(start.data)

          emb.list[[mm]] <- update_emb_local(
            emb.list[[mm]],
            new_cols = colnames(start.data),
            new_data = infer.data
          )

          ss.idx <- c(ss.idx, s)
          mm.idx <- c(mm.idx, mm)

          if (origin.infer) {
            arrive_o <- object@Assays[[assay_idx[[paste(tmp_s2, paths[[1]][l.m[length(l.m)]], sep="|")]]]]@data
            infer_o <- arrive_o %*% t(weight_M)
            rownames(infer_o) <- rownames(arrive_o)
            colnames(infer_o) <- colnames(start.data)

            tmp1 <- Create.Palette.Assay(
              infer_o, filter = FALSE, min.cells = 0, min.features = 0,
              modal = mm, sample = s, slot = "data", sparse = TRUE
            )
            object <- Add.Assay(tmp1, object, check = FALSE, class = "infer")
          }

          rm(start.data, road.list, target, arrive.data, arrive_name,
             infer_result, weight_M, infer.data)
          gc(FALSE)

        } else {

          # ----- 2.2 paths 多条 -----
          l  <- length(paths[[1]])
          l.m <- seq(l / 2) * 2

          through.m <- lapply(paths, function(x) x[l.m])
          through.m.uni <- unique(through.m)

          # ----- 2.2.1 through.m.uni 只有一个 -----
          if (length(through.m.uni) == 1) {

            idx_start <- assay_idx[[paste(s, paths[[1]][2], sep="|")]]
            start.data <- object@Assays[[idx_start]]@embedding
            if (L2Norm) start.data <- cosineNorm(start.data)

            # 构造 road.list
            road.list <- lapply(1:(length(l.m)-1), function(x) {
              ll <- list()
              ll[[1]] <- do.call(cbind, lapply(seq_along(paths), function(p) {
                tmp_s <- paths[[p]][l.m[x] + 1]
                idx_t <- assay_idx[[paste(tmp_s, paths[[p]][l.m[x]], sep="|")]]
                object@Assays[[idx_t]]@embedding
              }))
              ll[[1]] <- ll[[1]][, !duplicated(colnames(ll[[1]])), drop=FALSE]

              ll[[2]] <- do.call(cbind, lapply(seq_along(paths), function(p) {
                tmp_s <- paths[[p]][l.m[x] + 1]
                idx_t <- assay_idx[[paste(tmp_s, paths[[p]][l.m[x + 1]], sep="|")]]
                object@Assays[[idx_t]]@embedding
              }))
              ll[[2]] <- ll[[2]][, !duplicated(colnames(ll[[2]])), drop=FALSE]

              ll[[1]] <- ll[[1]][, colnames(ll[[2]]), drop=FALSE]

              if (L2Norm) {
                ll <- lapply(ll, cosineNorm)
              }
              ll
            })

            target <- real_modal_mat[[mm]]

            arrive.data <- do.call(cbind, lapply(seq_along(paths), function(p) {
              tmp_s <- paths[[p]][l.m[length(l.m)] - 1]
              idx_a <- assay_idx[[paste(tmp_s, paths[[p]][l.m[length(l.m)]], sep="|")]]
              object@Assays[[idx_a]]@embedding
            }))
            arrive.data <- arrive.data[, !duplicated(colnames(arrive.data)), drop=FALSE]
            arrive_name <- colnames(arrive.data)

            target <- target[, arrive_name, drop=FALSE]

            infer_result <- modal_infer_multi_one(
              start.data, road.list, nn, L2Norm, arrive.data, target,
              do.scaler, scaler.k, do.var, var.k
            )

            weight_M <- infer_result[["M"]]
            colnames(weight_M) <- arrive_name
            rownames(weight_M) <- colnames(start.data)
            weight_M_list[[mm]][arrive_name, colnames(start.data)] <- t(weight_M)

            infer.data <- infer_result[[2]]
            rownames(infer.data) <- rownames(target)
            colnames(infer.data) <- colnames(start.data)

            emb.list[[mm]] <- update_emb_local(
              emb.list[[mm]],
              new_cols = colnames(start.data),
              new_data = infer.data
            )

            ss.idx <- c(ss.idx, s)
            mm.idx <- c(mm.idx, mm)

            if (origin.infer) {
              arrive_o <- do.call(cbind, lapply(seq_along(paths), function(p) {
                tmp_s <- paths[[p]][l.m[length(l.m)] - 1]
                idx_o <- assay_idx[[paste(tmp_s, paths[[p]][l.m[length(l.m)]], sep="|")]]
                object@Assays[[idx_o]]@data
              }))
              arrive_o <- arrive_o[, !duplicated(colnames(arrive_o)), drop=FALSE]
              arrive_o <- arrive_o[, arrive_name, drop=FALSE]

              infer_o <- arrive_o %*% t(weight_M)
              rownames(infer_o) <- rownames(arrive_o)
              colnames(infer_o) <- colnames(start.data)

              tmp1 <- Create.Palette.Assay(
                infer_o, filter = FALSE, min.cells = 0, min.features = 0,
                modal = mm, sample = s, slot = "data", sparse = TRUE
              )
              object <- Add.Assay(tmp1, object, check = FALSE, class = "infer")
            }

            rm(start.data, road.list, target, arrive.data, arrive_name,
               infer_result, weight_M, infer.data)
            gc(FALSE)

          } else {
            # ----- 2.2.2 through.m.uni 多个 -----

            idxmap <- lapply(through.m.uni, function(x) {
              which(vapply(through.m, function(z) identical(z, x), logical(1)))
            })
            names(idxmap) <- vapply(through.m.uni, function(x) paste(x, collapse="|"), character(1))

            # start.list
            start.list <- lapply(through.m.uni, function(x) {
              idx_start <- assay_idx[[paste(s, x[1], sep="|")]]
              object@Assays[[idx_start]]@embedding
            })
            if (L2Norm) start.list <- lapply(start.list, cosineNorm)

            # road.multi.list
            road.multi.list <- lapply(idxmap, function(z) {
              lapply(1:(length(l.m)-1), function(x) {
                ll <- list()
                ll[[1]] <- do.call(cbind, lapply(z, function(p) {
                  tmp_s <- paths[[p]][l.m[x] + 1]
                  idx_t <- assay_idx[[paste(tmp_s, paths[[p]][l.m[x]], sep="|")]]
                  object@Assays[[idx_t]]@embedding
                }))
                ll[[1]] <- ll[[1]][, !duplicated(colnames(ll[[1]])), drop=FALSE]

                ll[[2]] <- do.call(cbind, lapply(z, function(p) {
                  tmp_s <- paths[[p]][l.m[x] + 1]
                  idx_t <- assay_idx[[paste(tmp_s, paths[[p]][l.m[x + 1]], sep="|")]]
                  object@Assays[[idx_t]]@embedding
                }))
                ll[[2]] <- ll[[2]][, !duplicated(colnames(ll[[2]])), drop=FALSE]

                ll[[1]] <- ll[[1]][, colnames(ll[[2]]), drop=FALSE]

                if (L2Norm) ll <- lapply(ll, cosineNorm)
                ll
              })
            })

            target <- real_modal_mat[[mm]]

            arrive.data <- do.call(
              cbind,
              lapply(road.multi.list, function(x) x[[2]][[2]])
            )

            arrive_name <- colnames(arrive.data)
            target <- target[,colnames(arrive.data[, !duplicated(colnames(arrive.data)), drop=FALSE]), drop=FALSE]

            infer_result <- modal_infer_multi_multi(
              start.list, road.multi.list, nn, L2Norm, arrive.data, target,
              do.scaler, scaler.k, do.var, var.k
            )

            weight_M <- infer_result[["M"]]
            colnames(weight_M) <- arrive_name

            if (length(unique(arrive_name)) < length(arrive_name)) {
              unique_arrive_name <- unique(arrive_name)
              map_matrix <- sapply(unique_arrive_name, function(x) arrive_name == x)
              map_matrix <- matrix(as.numeric(map_matrix), nrow = length(arrive_name))
              weight_M <- weight_M %*% map_matrix
              colnames(weight_M) <- unique_arrive_name
            }

            rownames(weight_M) <- colnames(start.list[[1]])
            weight_M_list[[mm]][arrive_name, colnames(start.list[[1]])] <- t(weight_M)

            infer.data <- infer_result[[2]]
            colnames(infer.data) <- colnames(start.list[[1]])

            emb.list[[mm]] <- update_emb_local(
              emb.list[[mm]],
              new_cols = colnames(start.list[[1]]),
              new_data = infer.data
            )

            ss.idx <- c(ss.idx, s)
            mm.idx <- c(mm.idx, mm)

            if (origin.infer) {
              arrive_o <- do.call(cbind, lapply(idxmap, function(z) {
                do.call(cbind, lapply(z, function(p) {
                  tmp_s <- paths[[p]][max(l.m) - 1]
                  idx_o <- assay_idx[[paste(tmp_s, mm, sep="|")]]
                  object@Assays[[idx_o]]@data
                }))
              }))

              infer_o <- arrive_o %*% t(weight_M)
              rownames(infer_o) <- rownames(arrive_o)
              colnames(infer_o) <- colnames(start.list[[1]])

              tmp1 <- Create.Palette.Assay(
                infer_o, filter = FALSE, min.cells = 0, min.features = 0,
                modal = mm, sample = s, slot = "data", sparse = TRUE
              )
              object <- Add.Assay(tmp1, object, check = FALSE, class = "infer")
            }

            rm(start.list, road.multi.list, target, arrive.data, arrive_name,
               infer_result, weight_M, infer.data)
            gc(FALSE)
          }
        }
      }
    }
  }

  message('MBG-guided inferring Done.\n')

  # ---------------------------
  # 3)  ScaleData
  # ---------------------------
  emb.list <- lapply(emb.list, function(x) {
    x$emb <- Seurat::ScaleData(
      object = x$emb,
      features = rownames(x$emb),
      do.scale = scale,
      do.center = center,
      scale.max = 1e6,
      block.size = 1000,
      min.cells.to.block = 3000,
      verbose = FALSE
    )
    x
  })

  # ---------------------------
  # 4) modal normalization
  # ---------------------------
  if (modal.norm) {
    C.list <- lapply(emb.list, function(x){
      sum(sqrt(diag(tcrossprod(x$emb))))
    })
    maxc <- max(unlist(C.list))
    emb.list <- mapply(function(x, cval){
      x$emb <- x$emb * maxc / cval
      x
    }, emb.list, C.list, SIMPLIFY = FALSE)
  }

  for (m in names(emb.list)) {
    rownames(emb.list[[m]]$emb) <- paste0(m, rownames(emb.list[[m]]$emb))
  }

  # ---------------------------
  # 5) weight
  # ---------------------------
  if (!is.null(weight)) {
    embname <- names(emb.list)
    weight <- weight[embname]

    if (any(unlist(weight) == 0)) {
      zidx <- which(unlist(weight) == 0)
      emb.list <- emb.list[-zidx]
      weight   <- weight[-zidx]
      embname  <- names(emb.list)
    }

    emb.list <- mapply(function(x, w){
      x$emb <- x$emb * w
      x
    }, emb.list, weight, SIMPLIFY = FALSE)

    names(emb.list) <- embname
  }

  # ---------------------------
  # 6) Cross-batch alignment + SNN
  # ---------------------------
  message('Cross-batch alignment.\n')

  if (sample_ratio <= 0 || sample_ratio > 1) {
    stop("sample_ratio must be in (0, 1].")
  }

  use_sampling_snn <- (length(Tname) > snn_full_max_cells) && (sample_ratio < 1)

  # 6.1 确定用于构建 SNN 的 cell 集
  if (use_sampling_snn) {
    set.seed(seed)
    Tsample <- unlist(lapply(sample.uni, function(s){
      cn <- CN.list[[s]]
      ns <- max(1, floor(length(cn) * sample_ratio))
      s <- intersect(cn,sample(cn, ns, replace = FALSE))
      s
    }))
    Tsample <- unique(Tsample)
  } else {
    Tsample <- Tname
  }

  if (supervised){
    if (is.null(group)) stop('Must selcet meta data under supervised case.')
    if (length(intersect(group, colnames(object@Meta))) == 0) stop('Group information Not found.')

    meta_cluster <- as.character(object@Meta[, group])
    names(meta_cluster) <- object@Meta[, 'cell.name']
    meta_cluster <- meta_cluster[Tsample]

    idx <- which(!is.na(meta_cluster))
    Tsample <- Tsample[idx]
    meta_cluster <- meta_cluster[Tsample]
  }

  # 6.2 对齐采样子集并拼 emb_sample
  emb.list.sample <- lapply(emb.list, function(x){
    list(emb = align_emb(x, Tsample), cols = Tsample)
  })
  emb_sample <- Reduce(rbind, lapply(emb.list.sample, `[[`, "emb"))
  emb_sample <- emb_sample[, Tsample, drop = FALSE]

  # 6.3 在 Tsample 上构建 SNN
  SNN_sample <- Seurat::FindNeighbors(
    t(emb_sample),
    k.param = nn.k,
    nn.method = nn.method,
    prune.SNN = 0,
    annoy.metric = annoy.metric,
    verbose = FALSE
  )$snn
  SNN_sample <- as(SNN_sample, "dgCMatrix")

  # 6.4 index.list（采样子集上）
  index.list.sample <- lapply(CN.list, function(cn){
    match(intersect(cn, Tsample), Tsample)
  })

  # 6.5 supervised 过滤（采样子集上）
  if (supervised) {
    meta_list.sample <- lapply(CN.list, function(cn){
      meta_cluster[intersect(cn, Tsample)]
    })

    SNN_sample <- filter_SNN(SNN_sample, index.list.sample, meta_list.sample)
    colnames(SNN_sample) <- rownames(SNN_sample) <- Tsample
  }

  # 6.6 在采样子集上求 Subspace
  K_list <- K_SNN(SNN_sample, index.list.sample, k, lambda)
  K <- K_list[["K"]]
  rownames(K) <- colnames(K) <- Tsample

  B <- Matrix::Matrix(0, nrow = nrow(SNN_sample), ncol = ncol(SNN_sample),
                      dimnames = list(Tsample, Tsample))
  diag(B) <- K_list[["B"]]

  Zname <- which(rowSums(K) == 0)
  diag(K)[Zname] <- K_list[["value"]]

  mat <- Getmat(emb_sample, K, B)
  rownames(mat) <- colnames(mat) <- rownames(emb_sample)

  Subspace <- suppressWarnings(
    RSpectra::eigs_sym(mat, k = nDims, which = "LA")
  )[["vectors"]] %>% as.matrix()

  # 6.7 对全量细胞一次性对齐并投影到 Subspace
  emb.list.full <- lapply(emb.list, function(x){
    list(emb = align_emb(x, Tname), cols = Tname)
  })
  emb_full <- Reduce(rbind, lapply(emb.list.full, `[[`, "emb"))
  emb_full <- emb_full[, Tname, drop = FALSE]

  embedding <- crossprod(as.matrix(emb_full), Subspace)

  if (emb_L2) {
    embedding <- t(cosineNorm(t(embedding)))
  }

  message('Palette integration Done! \n')
  rownames(embedding) <- Tname
  colnames(embedding) <- colnames(Subspace) <- paste0("dim_", seq_len(ncol(embedding)))
  object@Int.result[["bind"]] <- embedding

  return(object)
}

#' Transfer labels from reference to query using k-NN
#'
#' This function predicts labels for a query dataset based on a reference dataset
#' using k-nearest neighbors classification.
#'
#' @param ref A matrix or data frame of reference features (rows = samples, columns = features).
#' @param query A matrix or data frame of query features to be labeled.
#' @param ref.label A data frame or vector of labels corresponding to the reference samples.
#' @param k Integer. Number of nearest neighbors to use in the k-NN classifier.
#'
#' @return A data frame with predicted labels for the query samples.
#'
#' @export
#'
Trans.Knowledge <- function(ref,
                            query,
                            ref.label = NULL,
                            k = 5){
  pred_label <- class::knn(ref,
                           query,
                           as.character(ref.label[,1]),
                           k = k, prob = TRUE) %>% as.character()
  out <- data.frame(predict_label = pred_label)
  return(out)
}

#' Compute a latent representation of query data based on a reference and latent space
#'
#' This function projects query data into a latent space derived from reference data,
#' optionally using bi-directional sparse PCA (Bi_sPCA) for dimensionality reduction and alignment.
#'
#' @param query A numeric matrix representing the query dataset (features by cells).
#' @param ref A numeric matrix representing the reference dataset (features by cells).
#' @param latent A numeric matrix representing the latent space of the reference.
#' @param Bi_sPCA Logical. Whether to apply Bi-supervised Principal Component Analysis (Bi-sPCA) for alignment. Default is TRUE.
#' @param dims Integer. Number of dimensions to use for PCA subspace. Default is 30.
#' @param lambda Numeric. Regularization parameter for Bi-sPCA, between 0 and 1. Default is 0.5.
#' @param do.scale Logical. Whether to scale features to unit variance. Default is TRUE.
#' @param do.center Logical. Whether to center features to zero mean. Default is TRUE.
#' @param nn.k Integer. Number of nearest neighbors for SNN graph construction. Default is 20.
#' @param nn.method Character. Nearest neighbor method to use ("annoy", "rann", etc.). Default is "annoy".
#' @param annoy.metric Character. Distance metric for Annoy ("euclidean", "manhattan", etc.). Default is "euclidean".
#' @param prune.SNN Numeric. Threshold to prune edges in the SNN graph (fraction). Default is 1/15.
#' @param verbose Logical. Whether to display messages during computation. Default is TRUE.
#'
#' @return A numeric matrix representing the query projected into the latent space (latent features x query samples).
#'
#' @importFrom MASS ginv
#' @keywords internal
get_representation <- function(query,
                               ref,
                               latent,
                               Bi_sPCA = TRUE,
                               dims = 30,
                               lambda = 0.5,
                               do.scale = TRUE,
                               do.center = TRUE,
                               nn.k = 20,# knn snn参数
                               nn.method = "annoy",# knn snn参数
                               annoy.metric = "euclidean",
                               prune.SNN = 1/15,
                               verbose = TRUE){

  if(all(colnames(ref) %in% colnames(latent))){
    latent <- latent[,colnames(ref)]
  }else{
    stop('Unmatched reference and latent')
  }

  if(nrow(query) != nrow(ref)){
    stop('Unmatched dim between reference and query')
  }

  if(Bi_sPCA){
    if(dims >= nrow(ref)){
      Bi_sPCA = FALSE
      message('Input dimension is lower than the preset reduced dimension, set Bi_sPCA to FALSE.')
    }

    if(lambda < 0){
      lambda = 0
      message('lambda should be between 0 and 1, set lambda to 0.')
    }

    if(lambda > 1){
      lambda = 1
      message('lambda should be between 0 and 1, set lambda to 1.')
    }
  }

  latent <- latent[,colnames(ref)]

  if(Bi_sPCA){
    K1 <- Seurat::FindNeighbors(t(as.matrix(latent)),
                         k.param = nn.k,
                         nn.method = nn.method,
                         prune.SNN = prune.SNN,
                         annoy.metric = annoy.metric,
                         verbose = verbose)$snn

    degrees = colSums(K1)
    degrees[which(degrees == 0)] <- 1e-8
    D <- Matrix::Diagonal(x = 1 / sqrt(degrees))
    K1 <- D %*% K1 %*% D

    K2 <- Seurat::FindNeighbors(t(as.matrix(ref)),
                         k.param = nn.k,
                         nn.method = nn.method,
                         prune.SNN = prune.SNN,
                         annoy.metric = annoy.metric,
                         verbose = verbose)$snn

    degrees = colSums(K2)
    degrees[which(degrees == 0)] <- 1e-8
    D <- Matrix::Diagonal(x = 1 / sqrt(degrees))
    K2 <- D %*% K2 %*% D

    rm(degrees,D)

    obj <- cbind(query,ref)
    rownames(obj) <- paste0('dim',1:nrow(obj))
    rownames(ref) <- rownames(query) <- rownames(obj)

    data <- Seurat::ScaleData(object = obj,
                       features = rownames(ref),
                       do.scale = do.scale,
                       do.center = do.center,
                       scale.max = 10,
                       block.size = 1000,
                       min.cells.to.block = 3000,
                       verbose = FALSE)
    rm(obj)

    colnames(data) <- c(colnames(query),colnames(ref))
    s_ref <- data[,colnames(ref)]
    s_query <- data[,colnames(query)]

    mat <- s_ref%*%(K1-lambda*K2)%*%t(s_ref)

    rownames(mat) <- colnames(mat) <- rownames(ref)
    Subspace <- suppressWarnings(RSpectra::eigs_sym(mat,k = dims,which = "LA"))[["vectors"]] %>% as.matrix()

    embedding <- t(Subspace) %*% data
    colnames(embedding) <- c(colnames(query),colnames(ref))
    rownames(embedding) <- paste0('dim',1:nrow(embedding))
    ref <- embedding[,colnames(ref)]
    query <- embedding[,colnames(query)]

  }else{
    ref <- as.matrix(ref)
    query <- as.matrix(query)
  }

  latent <- as.matrix(latent)

  out <- latent %*% MASS::ginv(ref) %*% query
  colnames(out) <- colnames(query)
  rownames(out) <- rownames(latent)
  return(out)
}

#' Run fastMNN integration for multiple batches
#'
#' This function integrates multiple single-cell batches using the fastMNN algorithm
#' and returns a low-dimensional representation for all cells.
#'
#' @param batches A list of expression matrices, each representing a batch.
#' @param cell_name A character vector of cell names corresponding to all cells across batches.
#' @param is.normalize Logical, whether to normalize data before integration.
#' @param method Normalization method. Options include 'LogNormalize', 'CLR', or 'TF-IDF'.
#' @param nfeatures Number of features to use for variable feature selection (currently unused, can be NULL).
#' @param vargs Optional variable genes to use (currently unused, can be NULL).
#' @param k Number of nearest neighbors for the fastMNN algorithm.
#' @param out.npcs Number of principal components to compute in the output embedding.
#'
#' @return A matrix of integrated low-dimensional embeddings with cells as columns and dimensions as rows.
#'
# #' @importFrom Signac RunTFIDF
# #' @importFrom Signac RunSVD
#'
#' @export
#'
run_fastMNN = function(batches,
                       cell_name,
                       is.normalize = TRUE,
                       method = 'LogNormalize',
                       nfeatures = NULL,
                       vargs = NULL,
                       k = 20,
                       out.npcs = NA){
  .check_batchelor()
  seurat.list <- sapply(batches, function(x){
    Seurat::CreateSeuratObject(x)
  }, simplify = FALSE)
  # modified rownames
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    # VariableFeatures(x) = rownames(batches[[1]])
    x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = Inf)
  })
  # VariableFeatures(batch_seurat) = vargs
  vargenes <- seurat.list[[1]]@assays[["RNA"]]@var.features

  if (is.normalize == TRUE) {
    if(method == 'LogNormalize'){
      seurat.list <- lapply(X = seurat.list, FUN = function(x) {
        x <- Seurat::NormalizeData(x)
        x
      })
    }
    if(method == 'CLR'){
      seurat.list <- lapply(X = seurat.list, FUN = function(x) {
        x <- Seurat::NormalizeData(x, normalization.method = 'CLR', margin = 2)
        x
      })
    }

  }

  if(method == 'TF-IDF'){
    .check_signac()
    obj <- Seurat::CreateSeuratObject(Reduce(cbind,batches))
    obj <-  Signac::RunTFIDF(obj)
    obj[['lsi']] <- Signac::RunSVD(obj@assays[[1]]@data)

    data_list = lapply(1:length(batches), function(i) t(obj[['lsi']]@cell.embeddings[colnames(batches[[i]]),2:(out.npcs+1)]))
    rm(obj)

  }else{

    data_list = lapply(1:length(seurat.list), function(i) seurat.list[[i]]@assays[["RNA"]]@data[vargenes,])

  }

  rm(seurat.list)
  ##########################################################
  # run fastMNN
  t1 = Sys.time()
  out_mnn_total = do.call(batchelor::fastMNN, c(data_list, k = k, d = unname(out.npcs)))
  t2 = Sys.time()
  print(t2-t1)

  fastmnn_res = t(out_mnn_total@assays@data@listData[["reconstructed"]]@seed@components)
  colnames(fastmnn_res) = cell_name
  rownames(fastmnn_res) <- paste0('dim',1:nrow(fastmnn_res))
  return(fastmnn_res)
}

run_Seurat = function(batches,
                      cell_name,
                      is.normalize = TRUE,
                      method = 'LogNormalize',
                      nfeatures = NULL,
                      vargs = NULL,
                      k = 20,
                      dims = NULL,
                      out.npcs = 30){
  # preprocessing

  seurat.list <- sapply(batches, function(x){
    Seurat::CreateSeuratObject(x)
  }, simplify = FALSE)
  # modified rownames
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    # VariableFeatures(x) = rownames(batches[[1]])
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = Inf)
  })
  # VariableFeatures(batch_seurat) = vargs
  vargenes <- seurat.list[[1]]@assays[["RNA"]]@var.features

  if (is.normalize == TRUE) {
    if(method == 'LogNormalize'){
      seurat.list <- lapply(X = seurat.list, FUN = function(x) {
        x <- Seurat::NormalizeData(x)
      })
    }
    if(method == 'CLR'){
      seurat.list <- lapply(X = seurat.list, FUN = function(x) {
        x <- Seurat::NormalizeData(x, normalization.method = 'CLR', margin = 2)
      })
    }

    if(method == 'TF-IDF'){
      .check_signac()
      seurat.list <- lapply(X = seurat.list, FUN = function(x) {
        x <- Signac::RunTFIDF(x)
      })
    }
  }

  # features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = Inf)

  if(method == 'TF-IDF'){
    seurat.list <- lapply(X = seurat.list, FUN = function(x) {
      # x <- RunSVD(x,reduction.name = 'lsi')
      x[['lsi']] <-  Signac::RunSVD(x@assays[[1]]@data)
      x
    })
    cell_anchors = Seurat::FindIntegrationAnchors(object.list = seurat.list,
                                                  k.filter = 200,
                                                  reduction = 'rlsi',
                                                  anchor.features = vargenes,
                                                  dims = 2:(out.npcs+1))
    rm(seurat.list)
    comb <- Seurat::CreateSeuratObject(do.call(cbind, batches))
    Seurat::VariableFeatures(comb) = rownames(batches[[1]])
    comb <- Signac::RunTFIDF(comb)
    comb <- Signac::RunSVD(comb,reduction.name = 'lsi',n = out.npcs)
    gc()
    batch_correct = Seurat::IntegrateEmbeddings(anchorset = cell_anchors,
                                                reductions = comb[["lsi"]],
                                                new.reduction.name = "integrated_lsi",
                                                dims.to.integrate = 1:out.npcs)
    rm(comb,cell_anchors)
    seurat_res <- t(as.matrix(batch_correct[["integrated_lsi"]]@cell.embeddings))

    rm(batch_correct)
  }else{
    cell_anchors = Seurat::FindIntegrationAnchors(object.list = seurat.list,
                                                  k.filter = 200,
                                                  anchor.features = vargenes,
                                                  reduction = 'cca')
    rm(seurat.list)
    if(is.na(out.npcs)){

      if(is.null(dims)){
        dims = 1:length(vargenes)
      }

      batch_correct = Seurat::IntegrateData(anchorset = cell_anchors,dims = dims)
      seurat_res <- as.matrix(batch_correct@assays[["integrated"]]@data)
    }else{
      batch_correct = Seurat::IntegrateData(anchorset = cell_anchors)
      Seurat::DefaultAssay(batch_correct) = "integrated"
      batch_correct = Seurat::ScaleData(object = batch_correct)
      batch_correct = Seurat::RunPCA(object = batch_correct, npcs = out.npcs, verbose = FALSE)
      seurat_res = t(as.matrix(batch_correct@reductions$pca@cell.embeddings))
    }
  }

  colnames(seurat_res) = cell_name
  return(seurat_res)
}

# #' @importFrom harmony RunHarmony
#' @keywords internal
run_Harmony = function(batches, cell_name,
                       is.normalize = TRUE,
                       method = 'LogNormalize',
                       theta = 2,
                       out.npcs = 30){
  # preprocessing
  .check_harmony()
  batch_seurat = Seurat::CreateSeuratObject(do.call(cbind, batches))
  batch_info <- unlist(lapply(1:length(batches),FUN = function(x){
    x <- rep(paste0("bath",x),ncol(batches[[x]]))
    x
  }))
  batch_seurat[["batch"]] <- batch_info
  # modified rownames
  batch_seurat <- Seurat::FindVariableFeatures(batch_seurat, selection.method = "vst", nfeatures = Inf)
  # VariableFeatures(batch_seurat) = rownames(batches[[1]])
  # VariableFeatures(batch_seurat) = vargs
  vargenes <- batch_seurat@assays[["RNA"]]@var.features

  if (is.normalize == TRUE){
    if(method == 'LogNormalize'){
      batch_seurat <- Seurat::NormalizeData(batch_seurat)
      if(is.na(out.npcs)){
        # batch_seurat@assays[['RNA']]@data <- batch_seurat@assays[['RNA']]@counts
        batch_seurat@assays[['RNA']]@scale.data <- as.matrix(batch_seurat@assays[['RNA']]@counts)
        batch_seurat[['pca']] <- Seurat::CreateDimReducObject(embeddings = as.matrix(t(batch_seurat@assays[['RNA']]@data)),key = 'PC_')
      }else{
        batch_seurat <- Seurat::ScaleData(batch_seurat)
        batch_seurat <- Seurat::RunPCA(batch_seurat,npcs = out.npcs)
      }
    }
    if(method == 'CLR'){
      batch_seurat <- Seurat::NormalizeData(batch_seurat, normalization.method = 'CLR', margin = 2)
      if(is.na(out.npcs)){
        # batch_seurat@assays[['RNA']]@data <- batch_seurat@assays[['RNA']]@counts
        batch_seurat@assays[['RNA']]@scale.data <- as.matrix(batch_seurat@assays[['RNA']]@counts)
        batch_seurat[['pca']] <- Seurat::CreateDimReducObject(embeddings = as.matrix(t(batch_seurat@assays[['RNA']]@data)),key = 'PC_')
      }else{
        batch_seurat <- Seurat::ScaleData(batch_seurat)
        batch_seurat <- Seurat::RunPCA(batch_seurat,npcs = out.npcs)
      }
    }

    if(method == 'TF-IDF'){
      .check_signac()
      batch_seurat <- Signac::RunTFIDF(batch_seurat)
      tmp <- Signac::RunSVD(batch_seurat@assays[[1]]@data, n = (out.npcs+1))
      batch_seurat[['pca']] <- Seurat::CreateDimReducObject(embeddings = as.matrix(tmp@cell.embeddings[,2:(out.npcs+1)]),key = 'PC_')
      rm(tmp)
      # batch_seurat <- RunSVD(batch_seurat, n = (out.npcs+1),
      #                        reduction.name = 'pca')
    }
  }else{
    if(is.na(out.npcs)){
      # batch_seurat@assays[['RNA']]@data <- batch_seurat@assays[['RNA']]@counts
      batch_seurat@assays[['RNA']]@scale.data <- as.matrix(batch_seurat@assays[['RNA']]@counts)
      batch_seurat[['pca']] <- Seurat::CreateDimReducObject(embeddings = as.matrix(t(batch_seurat@assays[['RNA']]@data)),key = 'PC_')
    }else{
      batch_seurat <- Seurat::ScaleData(batch_seurat)
      batch_seurat <- Seurat::RunPCA(batch_seurat,npcs = out.npcs)
    }
  }

  #run Harmony
  t1 = Sys.time()
  if(method == 'TF-IDF'){
    batch_seurat = harmony::RunHarmony(object = batch_seurat,
                                      group.by.vars = 'batch',
                                      theta = theta,
                                      dims.use = 1:out.npcs,
                                      # dims.use = 2:(out.npcs+1),
                                      project.dim = F,
                                      plot_convergence = TRUE,
                                      nclust = 50,
                                      max.iter.cluster = 100)
  }else{
    batch_seurat = harmony::RunHarmony(object = batch_seurat,
                                      group.by.vars = 'batch',
                                      theta = theta,
                                      plot_convergence = TRUE,
                                      nclust = 50,
                                      max.iter.cluster = 100)
  }
  t2 = Sys.time()
  print(t2-t1)

  harmony_res = t(as.matrix(batch_seurat@reductions[["harmony"]]@cell.embeddings))
  colnames(harmony_res) = cell_name
  return(harmony_res)
}

#' Get reference indices matching query modalities
#'
#' This function identifies the indices of reference batches that correspond
#' to the modalities present in the query.
#'
#' @param q_modal A character vector of query modalities to match.
#' @param ref_batch A character or factor vector indicating the batch of each reference sample.
#' @param ref_modal A character vector indicating the modality of each reference sample.
#'
#' @return An integer vector of indices in `ref_batch` that match all query modalities. Returns NULL if no match is found.
#' @keywords internal
get_ref_idx <- function(q_modal, ref_batch, ref_modal) {
  q_modal <- unique(q_modal)
  matched_batches <- vector("list", length(q_modal))

  for (i in seq_along(q_modal)) {
    tmp_idx <- which(ref_modal == q_modal[i])
    if (length(tmp_idx) == 0) {
      return(NULL)
    }
    matched_batches[[i]] <- ref_batch[tmp_idx]
  }
  final_batches <- Reduce(intersect, matched_batches)
  ref_idx <- which(ref_batch %in% final_batches)

  return(ref_idx)
}

#' Split modalities into blocks by batch
#'
#' This function groups modalities by batch and identifies unique modality combinations across batches.
#'
#' @param modal A character vector of modalities corresponding to samples.
#' @param batch A character or factor vector indicating the batch of each sample.
#'
#' @return A list where each element contains the indices of batches sharing the same set of modalities.
#' @keywords internal
split_block <- function(modal,batch){
  uni_batch <- unique(batch)
  modal_list <- list()
  for(i in uni_batch){
    tmp_idx <- which(batch == i)
    modal_list[[i]] <- modal[tmp_idx]
  }

  vector_signatures <- sapply(modal_list, function(vec) paste(sort(vec), collapse = ","))

  unique_signatures <- unique(vector_signatures)
  out <- lapply(unique_signatures, function(sig) which(vector_signatures == sig))

  return(out)
}

#' Get the inferred embedding of the query data (Reference-based integration).
#'
#' This function integrates and aligns datasets from different modalities and batches,
#' using specified integration methods and returns the representations for alignment.
#'
#' @param reference A list of expression matrices or data frames for each reference sample.
#' @param query A list of expression matrices or data frames for each query sample.
#' @param modal_ref A character vector indicating the modality of each reference sample.
#' @param modal_query A character vector indicating the modality of each query sample.
#' @param batch_ref A character or factor vector indicating the batch of each reference sample.
#' @param batch_query A character or factor vector indicating the batch of each query sample.
#' @param modal_order Optional vector specifying the desired order of modalities for integration.
#' @param method Integration method for each modality. Options include 'fastMNN', 'Seurat', 'Harmony'.
#' @param is.normalize Logical or vector specifying whether to normalize each modality.
#' @param normalize_method Normalization method(s) for each modality. Options: 'LogNormalize','CLR','TF-IDF'.
#' @param method_out.npcs Number of principal components to output for each modality's integration.
#' @param latent A matrix representing latent features to guide the integration.
#' @param Bi_sPCA Logical indicating whether to use Bi-supervised Principal Component Analysis (Bi-sPCA) for subspace alignment.
#' @param dims Number of dimensions for the Bi-sPCA or subspace projection.
#' @param lambda Weighting parameter for Bi-sPCA subspace combination, between 0 and 1.
#' @param do.scale Logical indicating whether to scale data before subspace calculation.
#' @param do.center Logical indicating whether to center data before subspace calculation.
#' @param nn.k Number of nearest neighbors for constructing graphs in Bi-sPCA.
#' @param nn.method Method for nearest neighbor search. Default is 'annoy'.
#' @param annoy.metric Distance metric for Annoy. Default is 'euclidean'.
#' @param prune.SNN Threshold to prune shared nearest neighbor graph edges.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @return A list of representations for each modality composition block.
#'
#' @export
#'
RQ.Rep <- function(reference = NULL,
                   query = NULL,
                          modal_ref = NULL,
                          modal_query = NULL,
                          batch_ref = NULL,
                          batch_query = NULL,
                          modal_order,
                          method = 'fastMNN',
                          is.normalize = TRUE,
                          normalize_method = c('LogNormalize','CLR','TF-IDF'),
                          method_out.npcs = 40,
                          latent,
                          Bi_sPCA = TRUE,
                          dims = 30,
                          lambda = 0.5,
                          do.scale = TRUE,
                          do.center = TRUE,
                          nn.k = 20,
                          nn.method = "annoy",
                          annoy.metric = "euclidean",
                          prune.SNN = 1/15,
                          verbose = TRUE){

  # idx_q <- which(q_or_ref == 'query')
  # q_modal <- modal[idx_q]
  # q_batch <- batch[idx_q]
  # query <- datalist[idx_q]
  #
  # idx_r <- which(q_or_ref == 'ref')
  # r_modal <- modal[idx_r]
  # r_batch <- batch[idx_r]
  # ref <- datalist[idx_r]
  #
  # rm(datalist,modal,batch,q_or_ref)

  ref <- reference
  r_modal <- modal_ref
  r_batch <- batch_ref

  q_modal <- modal_query
  q_batch <- batch_query

  rm(reference,modal_ref,batch_ref,modal_query,batch_query)

  if(!is.null(modal_order)){
    diff <- setdiff(c(unique(q_modal),modal_order),modal_order)
    if(length(diff) > 0 ){
      modal_order <- unique(q_modal)
      message('The modal_order parameter is incorrect and has been reset')
    }
  }else{
    modal_order <- unique(q_modal)
  }

  if(length(method) != length(modal_order)){
    method <- rep(method[1],length(modal_order))
  }
  names(method) <- modal_order
  if(length(is.normalize) != length(modal_order)){
    is.normalize <- rep(is.normalize[1],length(modal_order))
  }
  names(is.normalize) <- modal_order
  if(length(method_out.npcs) != length(modal_order)){
    method_out.npcs <- rep(method_out.npcs[1],length(modal_order))
  }
  names(method_out.npcs) <- modal_order
  if(length(normalize_method) != length(modal_order)){
    normalize_method <- rep(normalize_method[1],length(modal_order))
  }
  names(normalize_method) <- modal_order

  batch_uni <- unique(q_batch)

  batch_block <- split_block(q_modal,q_batch)
  rep_list <- list()

  for(i in 1:length(batch_block)){

    tmp_block <- batch_uni[batch_block[[i]]]


    tmp_idx <- which(q_batch %in% tmp_block)
    tmp_modal <- q_modal[tmp_idx]
    tmp_data <- query[tmp_idx]

    ref_idx <- get_ref_idx(tmp_modal,r_batch,r_modal)
    if(is.null(ref_idx)){
      stop(paste0('No reference could be found that matched batch ',i,'.'))
    }

    tmp_ref <- ref[ref_idx]
    ref_modal <- r_modal[ref_idx]
    ref_batch <- r_batch[ref_idx]

    int_ref <- list()
    int_query <- list()

    for(j in unique(tmp_modal)){
      method_j <- method[j]
      n_method_j <- normalize_method[j]
      is.normalize_j <- is.normalize[j]
      method_out.npcs_j <- method_out.npcs[j]

      idx_q_j <- which(tmp_modal == j)
      idx_r_j <- which(ref_modal == j)
      data_j <- c(tmp_data[idx_q_j],tmp_ref[idx_r_j])
      cell_name_j <- unlist(sapply(data_j,colnames))
      if(method_j == "Seurat"){
        res <- run_Seurat(data_j,
                          cell_name_j,
                          is.normalize = is.normalize_j,
                          method = n_method_j,
                          out.npcs = method_out.npcs_j)
      }

      if(method_j == "fastMNN"){
        res <- run_fastMNN(data_j,
                           cell_name_j,
                           is.normalize = is.normalize_j,
                           method = n_method_j,
                           out.npcs = unname(method_out.npcs_j))
      }

      if(method_j == "Harmony"){
        res <- run_Harmony(data_j,
                           cell_name_j,
                           is.normalize = is.normalize_j,
                           method = n_method_j,
                           out.npcs = unname(method_out.npcs_j))
      }

      int_ref[[j]] <- res[,unlist(lapply(tmp_ref[idx_r_j],colnames))]
      int_query[[j]] <- res[,unlist(lapply(tmp_data[idx_q_j],colnames))]
    }

    rm(res,data_j)
    if(length(int_ref)>1){
      int_ref_mat <- Reduce(rbind,int_ref)
      int_query_mat <- Reduce(rbind,int_query)
    }else{
      int_ref_mat <- int_ref[[1]]
      int_query_mat <- int_query[[1]]
    }

    rm(int_ref,int_query)
    rep_list[[i]] <- get_representation(int_query_mat,
                                        int_ref_mat,
                                        latent = latent,
                                        Bi_sPCA = Bi_sPCA,
                                        dims = dims,
                                        lambda = lambda,
                                        do.scale = do.scale,
                                        do.center = do.center,
                                        nn.k = nn.k,
                                        nn.method = nn.method,
                                        annoy.metric = annoy.metric,
                                        prune.SNN = prune.SNN,
                                        verbose = verbose)

  }
  # rep_list[['latent']] <- latent
  return(rep_list)
}

#' Cosine normalize data matrix
#' @param x Input data matrix.
#'
#' @import Matrix
#' @keywords internal
cosineNorm <- function(x) {
  l2 <- sqrt(Matrix::colSums(x^2))
  l2 <- pmax(1e-8, l2)
  mat <- scale(x, center = F, scale = l2)
  mat
}

check_one <- function(x){
  x.uni <- unique(x)
  out <- c()
  for(i in x.uni){
    idx <- which(x == i)
    if(length(idx) == 1){
      out <- c(out,i)
    }
  }
  return(out)
}

check_aga <- function(query,ref){
  query <- unique(query)
  out <- c()
  for(i in query){
    idx <- which(ref == i)
    if(length(idx) == 1){
      out <- c(out,i)
    }
  }
  return(out)
}

safe_find_clusters <- function(snn, resolution, verbose = TRUE) {

  result <- tryCatch({
    do.call("FindClusters",
            list(object = snn,
                 resolution = resolution,
                 verbose = verbose),
            envir = asNamespace("Seurat"))
  }, error = function(e) {
    warning(paste("FindClusters failed:", e$message))
    return(NULL)
  })
  return(result)
}

#' @keywords internal
.check_suggested_pkg <- function(pkg, install_hint = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- sprintf("Package '%s' is required for this function but is not installed.", pkg)
    if (!is.null(install_hint)) {
      msg <- paste0(msg, "\n\nInstallation:\n  ", install_hint)
    }
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
.check_batchelor <- function() {
  .check_suggested_pkg(
    "batchelor",
    "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); BiocManager::install('batchelor')"
  )
}

#' @keywords internal
.check_harmony <- function() {
  .check_suggested_pkg(
    "harmony",
    "install.packages('harmony')"
  )
}

#' @keywords internal
.check_signac <- function() {
  .check_suggested_pkg(
    "Signac",
    "install.packages('Signac')"
  )
}

#' @import RcppArmadillo
NULL
