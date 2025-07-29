#' Modified EM test for normalized data
#' 
#' @description Modified version of EM test that works with normalized data instead of raw counts
#' @param x Normalized expression values for a gene
#' @param alpha.ini Initial alpha values
#' @param k0 Maximum iterations
#' @param C Convergence threshold
#' @param labels Cluster labels
#' @param group.num Number of groups
#' @param prior.weight Prior weight
#' @param earlystop Early stopping threshold
#' @return Vector containing p-value and EM statistic
em.stat.norm <- function(x, alpha.ini, k0 = 100, C = 1e-3, labels, group.num, prior.weight = 0.05, earlystop = 1e-4) {
  # Handle edge cases
  if(any(is.na(x)) || any(is.infinite(x))) {
    return(c(1, 0))  # Return no significance if data contains NA or Inf
  }
  
  # Convert labels to numeric
  labels <- as.numeric(labels)
  
  # Initialize parameters
  n <- length(x)
  alpha <- alpha.ini[1,]
  alpha.0 <- alpha.ini[2,]
  
  # Handle the case where all values are the same
  if(sd(x) < 1e-10) {
    return(c(1, 0))
  }
  
  # Calculate initial parameters for normal distribution
  mu <- tapply(x, labels, mean)
  sigma <- tapply(x, labels, sd)
  if(any(is.na(sigma)) || any(sigma < 1e-10)) {
    sigma[is.na(sigma) | sigma < 1e-10] <- 1e-6
  }
  
  tryCatch({
    # EM iterations
    ll.old <- -Inf
    for(k in 1:k0) {
      # E-step: Calculate responsibilities
      r <- matrix(0, n, group.num)
      for(j in 1:group.num) {
        r[,j] <- alpha[j] * dnorm(x, mu[j], sigma[j])
      }
      r.sum <- rowSums(r)
      if(any(r.sum == 0)) {
        r[r.sum == 0,] <- 1/group.num
        r.sum[r.sum == 0] <- 1
      }
      r <- r / r.sum
      
      # M-step: Update parameters
      nk <- colSums(r)
      alpha.new <- (nk + prior.weight * alpha.0) / (n + prior.weight)
      mu.new <- colSums(r * x) / nk
      sigma.new <- sqrt(colSums(r * (x - mu.new)^2) / nk)
      
      # Ensure minimum variance
      sigma.new[sigma.new < 1e-6] <- 1e-6
      
      # Update parameters
      alpha <- alpha.new
      mu <- mu.new
      sigma <- sigma.new
      
      # Check convergence
      ll.new <- sum(log(rowSums(matrix(alpha, n, group.num, byrow = TRUE) * 
                               matrix(dnorm(x, mu, sigma), n, group.num))))
      if(abs(ll.new - ll.old) < earlystop) break
      ll.old <- ll.new
    }
    
    # Calculate test statistic and p-value
    ll.1 <- ll.new
    ll.0 <- sum(dnorm(x, mean(x), sd(x), log = TRUE))
    stat <- 2 * (ll.1 - ll.0)
    
    # Handle numerical issues
    if(is.na(stat) || is.infinite(stat)) {
      return(c(1, 0))
    }
    
    # Ensure non-negative test statistic
    stat <- max(0, stat)
    p.value <- pchisq(stat, df = 2*(group.num-1), lower.tail = FALSE)
    
    return(c(p.value, stat))
  }, error = function(e) {
    return(c(1, 0))  # Return no significance if any error occurs
  })
}

#' Core function for Festem with normalized data
#' 
#' @description Core implementation of Festem algorithm adapted for normalized data
#' @inheritParams FestemCore
#' @return List containing test results and gene rankings
FestemCoreNorm <- function(data, cluster.labels, batch.id,
                          prior.weight = 0.05, prior.weight.filter = 0.9,
                          earlystop = 1e-4, min.percent = 0.01,
                          min.cell.num = 30, seed = 321,
                          num.threads = 2) {
  message("Running Festem with normalized data...")
  set.seed(seed)
  
  # Convert factors
  cluster.labels <- factor(cluster.labels)
  levels(cluster.labels) <- 1:nlevels(cluster.labels)
  batch.id <- factor(batch.id)
  
  # Setup parallel processing
  RNGversion("4.4.0")
  cl <- parallel::makeCluster(getOption("cl.cores", num.threads))
  parallel::clusterSetRNGStream(cl, iseed = seed)
  
  # Process each batch
  em.result <- vector("list", nlevels(batch.id))
  em.result.f <- vector("list", nlevels(batch.id))
  
  for(B in 1:nlevels(batch.id)) {
    # Get batch data
    data.tmp <- data[, batch.id == levels(batch.id)[B]]
    cluster.labels.tmp <- cluster.labels[batch.id == levels(batch.id)[B]]
    
    # Filter small clusters
    cluster.labels.tmp <- factor(cluster.labels.tmp)
    data.tmp <- data.tmp[, summary(cluster.labels.tmp)[cluster.labels.tmp] > 10]
    cluster.labels.tmp <- cluster.labels.tmp[summary(cluster.labels.tmp)[cluster.labels.tmp] > 10]
    cluster.labels.tmp <- factor(cluster.labels.tmp)
    levels(cluster.labels.tmp) <- 1:nlevels(cluster.labels.tmp)
    
    # Filter low-expressed genes
    nonzeros.num <- function(x) sum(abs(x) > 1e-10)
    tmp <- apply(data.tmp, 1, nonzeros.num)
    min.cell <- min(min.percent * ncol(data.tmp), min.cell.num)
    data.tmp <- data.tmp[tmp >= min.cell, ]
    
    # Calculate alpha labels
    alpha.label <- numeric(nlevels(cluster.labels.tmp) - 1)
    for(g in 1:length(alpha.label)) {
      alpha.label[g] <- sum(cluster.labels.tmp == g) / ncol(data.tmp)
    }
    
    # Run EM test
    message(paste0("Batch ", levels(batch.id)[B], ": EM-test ..."))
    if(requireNamespace("pbapply", quietly = TRUE)) {
      em.result[[B]] <- pbapply::pbapply(data.tmp, 1, em.stat.norm,
                                        alpha.ini = rbind(alpha.label, rep(1/nlevels(cluster.labels.tmp), length(alpha.label))),
                                        k0 = 100, C = 1e-3, labels = cluster.labels.tmp,
                                        group.num = nlevels(cluster.labels.tmp),
                                        prior.weight = prior.weight,
                                        earlystop = earlystop, cl = cl)
    } else {
      em.result[[B]] <- parallel::parApply(cl, data.tmp, 1, em.stat.norm,
                                          alpha.ini = rbind(alpha.label, rep(1/nlevels(cluster.labels.tmp), length(alpha.label))),
                                          k0 = 100, C = 1e-3, labels = cluster.labels.tmp,
                                          group.num = nlevels(cluster.labels.tmp),
                                          prior.weight = prior.weight,
                                          earlystop = earlystop)
    }
    
    # Run filtering
    message(paste0("Batch ", levels(batch.id)[B], ": filtering ..."))
    if(requireNamespace("pbapply", quietly = TRUE)) {
      em.result.f[[B]] <- pbapply::pbapply(data.tmp, 1, em.stat.norm,
                                          alpha.ini = rbind(alpha.label, rep(1/nlevels(cluster.labels.tmp), length(alpha.label))),
                                          k0 = 100, C = 1e-3, labels = cluster.labels.tmp,
                                          group.num = nlevels(cluster.labels.tmp),
                                          prior.weight = prior.weight.filter,
                                          earlystop = earlystop, cl = cl)
    } else {
      em.result.f[[B]] <- parallel::parApply(cl, data.tmp, 1, em.stat.norm,
                                            alpha.ini = rbind(alpha.label, rep(1/nlevels(cluster.labels.tmp), length(alpha.label))),
                                            k0 = 100, C = 1e-3, labels = cluster.labels.tmp,
                                            group.num = nlevels(cluster.labels.tmp),
                                            prior.weight = prior.weight.filter,
                                            earlystop = earlystop)
    }
  }
  
  # Combine results
  p.tmp <- matrix(NA, nrow = nrow(data), ncol = length(em.result))
  rownames(p.tmp) <- rownames(data)
  for(i in 1:length(em.result)) {
    p.tmp[colnames(em.result[[i]]), i] <- em.result[[i]][1, ]
  }
  p.tmp <- parallel::parApply(cl, p.tmp, 1, my.min)
  p.tmp <- p.tmp * length(em.result)
  
  em.tmp <- matrix(NA, nrow = nrow(data), ncol = length(em.result))
  rownames(em.tmp) <- rownames(data)
  for(i in 1:length(em.result)) {
    em.tmp[colnames(em.result.f[[i]]), i] <- em.result.f[[i]][2, ]
  }
  em.tmp <- apply(em.tmp, 1, my.max)
  
  parallel::stopCluster(cl)
  
  # Format results
  result.EM <- data.frame(names = rownames(data),
                         p = p.adjust(p.tmp, "BH"),
                         EM = em.tmp)
  
  tmp.na <- result.EM[is.na(result.EM$p), 1]
  result.EM <- result.EM[!is.na(result.EM$p), ]
  gene.names <- result.EM[result.EM$p < 0.05 & result.EM$EM > 0, ]
  gene.names <- gene.names[order(-gene.names$EM), ]
  
  tmp <- result.EM[result.EM$p >= 0.05 & result.EM$EM > 0, ]
  tmp <- tmp[order(-tmp$EM), ]
  gene.names <- rbind(gene.names, tmp)
  tmp <- result.EM[result.EM$EM <= 0, ]
  tmp <- tmp[order(tmp$p, -tmp$EM), ]
  gene.names <- rbind(gene.names, tmp)
  gene.names <- gene.names[, 1]
  gene.names <- c(gene.names, tmp.na)
  EM <- gene.names
  
  list(stat = result.EM, genelist = EM)
}

#' Run Festem with normalized data
#' 
#' @description Run Festem algorithm with Seurat pipeline using normalized data instead of raw counts.
#' This version is specifically designed to work with Seurat objects that only have the "data" layer.
#' 
#' @param object A Seurat object with normalized data in the "data" layer
#' @inheritParams RunFestem
#' @details This function adapts the Festem algorithm to work with normalized data by:
#' 1. Using a modified EM test suitable for continuous normalized data
#' 2. Skipping count-specific preprocessing steps
#' 3. Adjusting filtering criteria for normalized values
#' @return A Seurat object with Festem results added
#' @export
RunFestemData <- function(object, ...) {
  UseMethod("RunFestemData")
}

#' @rdname RunFestemData
#' @export
RunFestemData.Seurat <- function(object, G = NULL, prior = "HVG", batch = NULL,
                                prior.weight = 0.05, prior.weight.filter = 0.9,
                                earlystop = 1e-4, min.percent = 0.01,
                                min.cell.num = 30, seed = 321,
                                num.threads = 1, FDR_level = 0.05,
                                assay = "RNA",
                                prior_parameters = list(HVG_num = 2000,
                                                      PC_dims = 50,
                                                      resolution = 0.7), ...) {
  if(!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Running Festem on a Seurat object requires Seurat")
  }
  
  # Check for data layer
  if(!"data" %in% SeuratObject::Layers(object, assay = assay)) {
    stop("No 'data' layer found in the specified assay")
  }
  
  object@active.assay <- assay
  
  # Generate or get prior labels
  if(prior == "HVG") {
    message("Generating priors ... ")
    object <- Seurat::FindVariableFeatures(object = object,
                                         selection.method = "vst",
                                         nfeatures = prior_parameters[["HVG_num"]] %or% 8000,
                                         verbose = FALSE)
    object <- Seurat::ScaleData(object)
    object <- Seurat::RunPCA(object, verbose = FALSE)
    object <- Seurat::FindNeighbors(object = object,
                                   dims = 1:(prior_parameters[["PC_dims"]] %or% 50))
    
    if(is.null(G)) {
      object <- Seurat::FindClusters(object,
                                   resolution = prior_parameters[["resolution"]],
                                   verbose = FALSE)
      G <- nlevels(object@active.ident)
      prior_label <- object@active.ident
    } else {
      for(k in 1:100) {
        object <- Seurat::FindClusters(object, resolution = 0.01 * k,
                                     verbose = FALSE)
        if(nlevels(object@active.ident) >= G) {
          if(nlevels(object@active.ident) == G) {
            prior_label <- object@active.ident
            break
          } else {
            stop(paste0("Cannot generate a prior with ", G,
                       " clusters. Try setting G as ",
                       nlevels(object@active.ident),
                       " instead and re-run Festem."))
          }
        }
      }
    }
  } else if(prior == "active.ident") {
    prior_label <- object@active.ident
  } else if(length(prior) == 1) {
    if(prior %in% colnames(object@meta.data)) {
      prior_label <- object@meta.data[, prior]
      if(length(unique(prior_label)) != G) {
        warning("Number of clusters in the prior does not equal G. Using prior cluster number.")
        G <- length(unique(prior_label))
      }
    } else {
      stop(paste0(prior, " was not found in metadata"))
    }
  } else if(length(prior) == ncol(object)) {
    prior_label <- factor(prior)
    if(nlevels(prior_label) != G) {
      warning("Number of clusters in the prior does not equal G. Using prior cluster number.")
      G <- nlevels(prior_label)
    }
  } else {
    stop("Invalid prior specification")
  }
  
  # Handle batch information
  if(is.null(batch)) {
    batch_id <- rep(1, ncol(object))
  } else if(length(batch) == 1) {
    if(batch %in% colnames(object@meta.data)) {
      batch_id <- object@meta.data[, batch]
    } else {
      stop(paste0(batch, " was not found in metadata"))
    }
  } else if(length(batch) == ncol(object)) {
    batch_id <- batch
  } else {
    stop("Invalid batch specification")
  }
  
  if(G == 1) {
    stop("Number of components must be larger than 1")
  }
  
  # Get normalized data
  data_use <- SeuratObject::LayerData(object, assay = assay, layer = "data")
  
  # Filter batches
  filter_table <- table(prior_label, batch_id)
  filter_criteria <- apply(filter_table, 2, function(x) {
    y <- sort(x, decreasing = TRUE)
    y[2]
  })
  batch_keep <- names(filter_criteria)[filter_criteria > 10]
  filter_flag <- batch_id %in% batch_keep
  data_use <- data_use[, filter_flag]
  prior_label <- prior_label[filter_flag]
  batch_id <- batch_id[filter_flag]
  
  # Run modified Festem
  result <- FestemCoreNorm(data_use,
                          cluster.labels = prior_label,
                          batch.id = batch_id,
                          prior.weight = prior.weight,
                          prior.weight.filter = prior.weight.filter,
                          earlystop = earlystop,
                          min.percent = min.percent,
                          min.cell.num = min.cell.num,
                          seed = seed,
                          num.threads = num.threads)
  
  # Update Seurat object
  num_features <- sum(result[["stat"]]$p < FDR_level & result[["stat"]]$EM > 0)
  Seurat::VariableFeatures(object) <- result[["genelist"]][1:min(num_features, 6000)]
  
  # Create meta.features data frame
  meta_to_add <- data.frame(
    "Festem_p_adj" = rep(NA, nrow(object[[assay]])),
    "Festem_EM" = rep(NA, nrow(object[[assay]])),
    "Festem_rank" = rep(NA, nrow(object[[assay]]))
  )
  rownames(meta_to_add) <- rownames(object[[assay]])
  meta_to_add[result[["stat"]]$names, "Festem_p_adj"] <- result[["stat"]]$p
  meta_to_add[result[["stat"]]$names, "Festem_EM"] <- result[["stat"]]$EM
  meta_to_add[, "Festem_rank"] <- sapply(rownames(object[[assay]]), function(x) {
    tmp <- which(x == result[["genelist"]])
    if(length(tmp) > 0) tmp else NA
  })
  
  # Add meta.features to the assay
  if(is.null(object[[assay]]@meta.features)) {
    object[[assay]]@meta.features <- meta_to_add
  } else {
    # Merge with existing meta.features
    existing_meta <- object[[assay]]@meta.features
    for(col in colnames(meta_to_add)) {
      existing_meta[[col]] <- meta_to_add[[col]]
    }
    object[[assay]]@meta.features <- existing_meta
  }
  
  return(object)
}
