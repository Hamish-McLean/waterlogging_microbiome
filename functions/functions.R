#' @title Calculate Bray-Curtis distances and run ADONIS
#' @param data Data object
#' @return `data.table` of ANOVA results
#'
adonisModel <- function(data, kingdom, source, design, strata = NULL, nperm = 1000) {
  set.seed(SEED)
  vg <- vegdist(t(counts(data$dds, normalize = T)), method = "bray")
  model <- adonis2(
    update(design, vg ~ .),
    data$colData,
    strata = if (!is.null(strata)) data$colData[[strata]] else NULL,
    # strata = colData(dds)$experiment:colData(dds)$block, # For repeat measures
    permutations = nperm,
    by = "terms"
  )
  results <- data.frame(model)
  total_ss <- results["Total", "SumOfSqs"]
  data.table(
    kingdom = kingdom,
    source  = source,
    factor  = rownames(results),
    df      = results$Df,
    SS      = results$SumOfSqs,
    R2      = results$R2,
    F       = results$F,
    P       = results$Pr..F.,
    stars   = p_stars(results$Pr..F.),
    var     = results$SumOfSqs / total_ss * 100,
    n_perm  = attr(model, "control")$nperm
  )
}


#' @title Analyze alpha diversity
#' @description Perform alpha diversity analysis on each subset of data
#' @param data A list containing countData and colData
#' @param kingdom A string representing the kingdom of the data (e.g., "Fungi" or "Bacteria")
#' @param source A string representing the source of the data (e.g., "root" or "soil")
#' @param metrics A vector of alpha diversity metrics to analyze
#' @param formula A formula for the model to use in the analysis
#' @return A list with the alpha diversity results
#'
alphaModel <- function(data, kingdom, source, metrics, formula) {
  cat(paste0("\nAnalyzing alpha diversity for ", kingdom, ": ", source, "\n"))

  # Run alpha diversity analysis for each metric
  anova_results <- list()
  for (met in metrics) {

    # Ensure the metric column exists
    if (!met %in% names(data$alphaData)) {
      warning("Metric '", met, "' not found in alpha data for ", kingdom, ": ", source)
      next
    }
    
    temp_dt <- copy(data$alphaData)
    temp_dt[, model_metric := rank(get(met))] # Rank of metric
    # Run the permutation model
    res <- tryCatch({
      perm.lmer(update(formula, model_metric ~ .), temp_dt, type = "anova", nperm = 1000)
    }, error = function(e) {
      warning("Error running perm.lmer for ", kingdom, ":", source, ", metric ", met, ": ", e$message)
      NULL # Return NULL on error
    })

    if (!is.null(res)) {
      res$metric <- met
      anova_results[[met]] <- res
    } else {
      anova_results[[met]] <- "Error" # Indicate error
      warning("ANOVA for ", kingdom, ":", source, ", metric ", met, " failed.")
    }
  }
  results <- rbindlist(anova_results)
  # results$kingdom <- kingdom
  # results$source  <- source
  # results$stars   <- p_stars(results$p)
  # results
  data.table(
    kingdom = kingdom,
    source  = source,
    metric  = results$metric,
    factor  = results$Factor,
    df      = results$df,
    LRT     = results$LRT,
    F       = results$F,
    p       = results$p,
    stars   = p_stars(results$p)
  )
}


#' @title Analyze alpha diversity
#' @description Perform alpha diversity analysis on each subset of data
#' @param data A list containing countData and colData
#' @param kingdom A string representing the kingdom of the data (e.g., "Fungi" or "Bacteria")
#' @param source A string representing the source of the data (e.g., "root" or "soil")
#' @param metrics A vector of alpha diversity metrics to analyze
#' @param formula A formula for the model to use in the analysis
#' @return A list with the alpha diversity results
#'
alphaLM <- function(data, kingdom, source, metrics, formula) {
  cat(paste0("\nAnalyzing alpha diversity for ", kingdom, ": ", source, "\n"))

  # Run alpha diversity analysis for each metric
  anova_results <- list()
  for (met in metrics) {

    # Ensure the metric column exists
    if (!met %in% names(data$alphaData)) {
      warning("Metric '", met, "' not found in alpha data for ", kingdom, ": ", source)
      next
    }
    
    temp_dt <- data.table::copy(data$alphaData)
    temp_dt[, model_metric := rank(get(met))] # Rank of metric
    # formula <- update(formula, model_metric ~ .)
    # formula <- update(formula, as.formula(paste0("rank(", met, ") ~ .")))

    rhs_terms <- attr(terms(formula), "term.labels")
    new_formula <- reformulate(
      termlabels = rhs_terms,
      response = "model_metric" #paste0("rank(", met, ")")
    )

    # Create environment with data for robust evaluation
    # f_env <- new.env(parent = environment(new_formula))
    # f_env$temp_dt <- temp_dt
    # environment(new_formula) <- f_env

    # Run the permutation model
    # return(aovp(new_formula, data = temp_dt, seqs = TRUE))

    res <- tryCatch({
      aovp(update(new_formula, model_metric ~ .), data = temp_dt, seqs = TRUE)
    }, error = function(e) {
      warning("Error running aovp for ", kingdom, ":", source, ", metric ", met, ": ", e$message)
      NULL # Return NULL on error
    })
    # return(summary(res))
    if (!is.null(res)) {
      anova_df <- if (length(res) > 3) { # length = 3 for repeat measures anova
        summary(res)[[1]] |> as.data.frame()
      } else {
        bind_rows(
          summary(res)[["Error: plot"]][[1]] |> as.data.frame(), 
          summary(res)[["Error: Within"]][[1]] |> as.data.frame()
        )
      }
      anova_df$metric <- met
      anova_df$factor <- rownames(anova_df) |> trimws()
      anova_df$total_ss <- sum((temp_dt$model_metric - mean(temp_dt$model_metric))^2)
      anova_results[[met]] <- anova_df
    } else {
      anova_results[[met]] <- "Error"
      warning("ANOVA for ", kingdom, ":", source, ", metric ", met, " failed.")
    }
  }
  results <- rbindlist(anova_results)
  # return(results)
  data.table(
    kingdom = kingdom,
    source  = source,
    metric  = results$metric,
    factor  = results$factor,
    df      = results[["Df"]],
    SS      = results[["R Sum Sq"]],
    MS      = results[["R Mean Sq"]],
    Iter    = results[["Iter"]],
    p       = results[["Pr(Prob)"]],
    stars   = p_stars(results[["Pr(Prob)"]]),
    var     = results[["R Sum Sq"]] / results$total_ss * 100
  )
}


#' @title Post-hoc analysis for perm.lmer models
#' @description Splits data by a grouping variable, runs a simplified
#'   perm.lmer model on each subset, and calculates aggregate means.
#' @param data A data.table or data.frame containing the experiment data.
#' @param metric A string representing the name of the column to be analyzed (e.g., "alpha_diversity").
#' @param design A formula object for the fixed effects of the simplified model (e.g., `~ treatment`).
#' @param group_by_var A string representing the column name to split the data by (e.g., "experiment" or "timepoint").
#' @param nperm The number of permutations to use for perm.lmer.
#' @return A list containing two elements:
#'   - `models`: A list of the permutes objects from perm.lmer for each subset.
#'   - `aggregates`: A list of data.tables with the mean of the metric for each treatment group in each subset.
alphaPostHoc <- function(data, metric, design, group_by_var, nperm = 1000) {

  # Split the data into a list of data.tables based on the grouping variable
  split_data <- split(data, by = group_by_var)

  formula <- update(design, metric_rank ~ .)

  # Initialize lists to store the results
  fitted_models <- list()
  aggregate_tables <- list()

  # Iterate through each subset of the data
  for (group_name in names(split_data)) {
    subset_dt <- split_data[[group_name]][, metric_rank := rank(get(metric))]

    message(paste("Running model for", group_by_var, "=", group_name))

    # Run the perm.lmer model on the subset
    model_result <- tryCatch({
      perm.lmer(formula, subset_dt, type = "anova", nperm = nperm)
    }, error = function(e) {
      warning(paste("Error running perm.lmer for", group_by_var, "=", group_name, ":", e$message))
      NULL
    })

    # Store the model result
    fitted_models[[group_name]] <- model_result

    # Calculate and store the aggregate means for the subset
    aggregate_tables[[group_name]] <- aggregate(
      get(metric) ~ treatment,
      data = subset_dt,
      FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    )
  }

  # Return a list of the models and aggregate tables
  return(list(models = fitted_models, aggregates = aggregate_tables))
}


#' @title Generic iterator function
#' @description Applies a function to each subset of data
#' @param data_list A list of data subsets (e.g., FUN or BAC)
#' @param name A string representing the name of the data (e.g, "Fungi" or "Bacteria")
#' @param func A function to apply to each subset
#' @param ... Additional arguments to pass to the function
#' @return A list of results from applying the function to each subset
#'
apply_to_subsets <- function(data_list, func, ...) {
  # Iterate over kingdoms (FUN, BAC) and sources (root, soil)
  kingdom_data <- lapply(names(data_list), function(kingdom) {
    source_data <- lapply(names(data_list[[kingdom]]), function(source) {
      # Select data subset by kingdom and source
      data <- data_list[[kingdom]][[source]]
      # Apply the provided function to each source within each kingdom
      res <- func(data, kingdom, source, ...)
    })
    names(source_data) <- names(data_list[[kingdom]])
    return(source_data)
  })
  names(kingdom_data) <- names(data_list) # Set names for the result list
  return(kingdom_data)
}


#' @title Calculate alpha diversity
#' @description Calculate alpha diversity metrics for each subset of data
#' @param data A list containing countData and colData
#' @return A list with the alpha diversity data added
#' 
calculateAlpha <- function(data, ...) {
  alpha_dt <- plot_alpha(data$countData, data$colData, returnData = TRUE)
  data$alphaData <- alpha_dt[as.data.table(data$colData, keep.rownames = "Samples"), on = "Samples"]
  return(data)
}


#' @title Calculate PCA
#'
calculatePCA <- function(data, name, source, n_pcs = NULL) {
  # Calculate PCA scores
  pca <- des_to_pca(data$dds)

  # Adjust PCA scores with total variance
  data$pca <- pca #t(data.frame(t(pca$x) * pca$percentVar))
  return(data)
}


#' @title Flatten and rowbind a nested list of tables
#' @param data_list A nested list of tables
#' @return A single table of data
#' 
combineTables <- function(data_list) {
  unlist(data_list, recursive = FALSE) |> rbindlist()
}


#' @title Fit LM to copy number data run ANOVA
#' @param data Data object with copy number data
#' 
copyNumberModelLM <- function(data, kingdom, source, design) {
  anova    <- aov(update(design, copy_number ~ .), data$colData)
  anova_df <- summary(anova)[[1]] # anova |> data.frame()
  total_ss <- sum(anova_df[["Sum Sq"]])
  data.table(
    kingdom = kingdom,
    source  = source,
    factor  = rownames(anova_df) |> trimws(),
    df      = anova_df[["Df"]],
    SS      = anova_df[["Sum Sq"]],
    MS      = anova_df[["Mean Sq"]],
    F       = anova_df[["F value"]],
    P       = anova_df[["Pr(>F)"]],
    stars   = p_stars(anova_df[["Pr(>F)"]]),
    var     = anova_df[["Sum Sq"]] / total_ss * 100
  )
}


#' @title Fit LMM to copy number data run ANOVA
#' @param data Data object with copy number data
#' 
copyNumberModelLMM <- function(data, kingdom, source, design) {
  message <- ""
  model <- withCallingHandlers(
    lmer(update(design, copy_number ~ .), data$colData),
    message = function(m) {
      message <<- sub(":.*", "", conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  r2 <- r2_nakagawa(model, verbose = FALSE)
  # part_r2 <- partR2(model, data = data$colData, partvars = c("experiment", "treatment"))
  # print(part_r2)
  results <- Anova(model, type = 3) |> data.frame()
  data.table(
    kingdom = kingdom,
    source  = source,
    factor  = rownames(results),
    df      = results$Df,
    Chisq   = results$Chisq,
    P       = results$Pr..Chisq.,
    R2_marg = r2$R2_marginal,
    R2_cond = r2$R2_conditional,
    stars   = p_stars(results$Pr..Chisq.),
    message = message
  )
}


#' @title DESeq differential abundance analysis
#'
diffAbundance <- function(data, kingdom, source, term) {
  int_term <- paste0(term, "2.treatment1")
  # Extract results for E1 (reference), E2, and interaction
  res_x1  <- results(data$dds, name = "treatment_1_vs_0")
  res_x2  <- results(data$dds, contrast = list(c("treatment_1_vs_0", int_term)))
  res_int <- results(data$dds, name = int_term)

  # Combine into a data.table
  dt <- data.table(
    kingdom = kingdom,
    source = source,
    ASV = rownames(res_x1),
    log2FC_X1 = res_x1$log2FoldChange,
    padj_X1 = res_x1$padj,
    sig_X1 = p_stars(res_x1$padj),
    log2FC_X2 = res_x2$log2FoldChange,
    padj_X2 = res_x2$padj,
    sig_X2 = p_stars(res_x2$padj),
    log2FC_inter = res_int$log2FoldChange,
    padj_inter = res_int$padj,
    sig_inter = p_stars(res_int$padj),
    taxa = data$taxData[rownames(res_x1), ]$rank
  )
  
  # Filter for any significant result
  dt_filt <- dt[
    complete.cases(padj_X1, padj_X2, padj_inter) &
    pmin(padj_X1, padj_X2, padj_inter, na.rm = TRUE) < 0.05
  ]

  # Optionally, sort by log2FC
  dt_filt <- dt_filt[order(pmax(abs(log2FC_X1), abs(log2FC_X2), na.rm = TRUE), decreasing = TRUE)]
  
  return(dt_filt)
}


#' @title DESeq differential abundance analysis
#'
diffAbundanceOverall <- function(data, kingdom, source, term) {
  # Extract results for E1 (reference), E2, and interaction
  res <- results(data$dds, name = "treatment_1_vs_0")

  # Combine into a data.table
  dt <- data.table(
    kingdom     = kingdom,
    source      = source,
    ASV         = rownames(res),
    baseMean    = res$baseMean,
    log2FC      = res$log2FoldChange,
    correlation = ifelse(res$log2FoldChange > 0, "positive", "negative"),
    padj        = res$padj,
    sig         = p_stars(res$padj),
    taxa        = data$taxData[rownames(res), ]$rank
  )
  
  # Filter for any significant result
  dt_filt <- dt[complete.cases(padj) & padj < 0.05]

  # Sort
  dt_filt <- dt_filt[order(abs(baseMean), decreasing = TRUE)]
  
  return(dt_filt)
}


#' @title NMDS plot
#' @description Create a NMDS plot for a given DESeq2 object
#' @param dds DESeq2 object
#' @param shape A string representing factor to use for the shape aesthetic
#' @param filename A string representing the filename to save the plot
#' @return A ggplot object of the NMDS plot
#'
nmdsPlot <- function(dds, shape, filename) {
  vg <- vegdist(t(counts(dds, normalize = T)), method = "bray")
  ord <- metaMDS(vg, trace=0)

  plot <- plotOrd(
    scores(ord), 
    colData(dds), 
    design = "treatment", 
    shape = shape, 
    alpha = 0.75, cbPalette = T
  )
  ggsave(filename, path = FIGURES_DIR)
  return(plot)
}


#' @title Significance stars
#' @description Returns significance stars based on a p-value
#' @param p A p-value or list of p-values
#' @return A string of stars representing significance
#'
p_stars <- function(p) {
  # Define cutoffs and corresponding symbols
  cuts  <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
  syms  <- c("***", "**", "*", ".", "")
  idx   <- findInterval(p, cuts, rightmost.closed = TRUE, left.open = TRUE)
  stars <- syms[idx]
  stars[is.na(p)] <- ""
  stars
}


#' @title PCA linear model anova
#' 
pcaLM <- function(data, kingdom, source, design, type, n_pcs = NULL) {
  pcs      <- 1:ifelse(!is.null(n_pcs), n_pcs, ncol(data$pca$x) - 1)
  pca_data <- as.data.frame(data$pca$x[, pcs])
  variance <- data$pca$percentVar[pcs]
  pc_names <- colnames(pca_data)

  results  <- lapply(seq_along(pcs), function(i) {
    pc_name  <- pc_names[i]
    pc_dat   <- pca_data[[pc_name]]
    df       <- cbind(PC = pc_dat, data$colData)
    # model    <- lm(update(design, PC ~ .), df)
    # anova    <- Anova(model, type = 3)
    anova    <- aov(update(design, PC ~ .), df)
    anova_df <- if (type == "experiment") {
      summary(anova)[[1]] # anova |> data.frame()
    } else {
      bind_rows(
        summary(anova)[["Error: plot"]][[1]], 
        summary(anova)[["Error: Within"]][[1]]
      )
    }
    # total_ss <- sum(anova_df["Sum Sq"])
    total_ss <- sum((pc_dat - mean(pc_dat))^2)
    data.table(
      kingdom = kingdom,
      source  = source,
      pc      = pc_name,
      factor  = rownames(anova_df) |> trimws(),
      df      = anova_df[["Df"]],
      SS      = anova_df[["Sum Sq"]],
      MS      = anova_df[["Mean Sq"]],
      F       = anova_df[["F value"]],
      P       = anova_df[["Pr(>F)"]],
      pc_var  = variance[i],
      var     = anova_df[["Sum Sq"]] / total_ss * 100,
      var_adj = (anova_df[["Sum Sq"]] / total_ss * 100) * variance[i]
    )
  })

  # Bind rows and rename columns
  final_results <- bind_rows(results)
  colnames(final_results) <- sub("\\..*", "", colnames(final_results))

  return(final_results)
}


#' @title PCA linear mixed model anova
#'
pcaLMM <- function(data, kingdom, source, design, n_pcs = NULL) {
  pcs      <- 1:ifelse(!is.null(n_pcs), n_pcs, ncol(data$pca$x) - 1)
  pca_data <- as.data.frame(data$pca$x[, pcs])
  pc_names <- colnames(pca_data)
  variance <- data$pca$percentVar[pcs]

  results <- lapply(seq_along(pcs), function(i) {
    pc_name <- pc_names[i]
    df      <- cbind(PC = pca_data[[pc_name]], data$colData)
    message <- NULL
    model   <- withCallingHandlers(
      lmer(update(design, PC ~ .), df),
      message = function(m) {
        msg <- conditionMessage(m)
        message <<- sub(":.*", "", msg)
        invokeRestart("muffleMessage")
      }
    )
    anova   <- Anova(model, type = 3) |> data.frame()
    r2      <- r2_nakagawa(model, verbose = FALSE)
    # part_r2 <- partR2(model, partvars = c("treatment", "timepoint", "treatment:timepoint"))
    data.table(
      kingdom = kingdom,
      source  = source,
      pc      = pc_name,
      factor  = rownames(anova),
      Chisq   = anova$Chisq,
      df      = anova$Df,
      P       = anova$Pr..Chisq.,
      R2_marg = r2$R2_marginal,
      R2_cond = r2$R2_conditional,
      pc_var  = variance[i],
      message = message
    )
  })

  bind_rows(results)
}


#' @title PCA plot
#' @description Create a PCA plot for a given pair of PCs
#' @param data PCA data
#' @param shape A string representing factor to use for the shape aesthetic
#' @param filename A string representing the filename to save the plot
#' @return A ggplot object of the PCA plot
#'
pcaPlot <- function(data, pcs, shape, filename) {
  plot <- plotOrd(
    t(data.frame(t(data$pca$x) * data$pca$percentVar)),
    # data$pca, 
    data$colData, 
    design = "treatment", 
    shapes = shape, 
    axes = pcs,
    alpha = 0.75, 
    cbPalette = TRUE
  )
  ggsave(filename, path = FIGURES_DIR)
  return(plot)
}


#' @title Print table
#' 
printTable <- function(table, caption = NULL) {
  if (knitr::is_html_output()) {
    library(DT)
    DT::datatable(
      table,
      caption = caption,
      extensions = "Buttons",
      options = list(
        dom = "Bt",
        buttons = "copy",
        scrollX = TRUE,
        pageLength = 70,
        lengthMenu = c(10, 20, 50, 100, -1),
        autoWidth = TRUE
      ),
      autoHideNavigation = TRUE,
      rownames = FALSE
    )
  } else {
    table
  }
}
