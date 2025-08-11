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


#' @title Flatten and rowbind a nested list of tables
#' @param data_list A nested list of tables
#' @return A single table of data
#' 
combine_tables <- function(data_list) {
  unlist(data_list, recursive = FALSE) |> rbindlist()
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
# p_stars <- function(p) {
#   ifelse(is.na(p), "",
#     ifelse(p < 0.001, "***",
#       ifelse(p < 0.01, "**",
#         ifelse(p < 0.05, "*",
#           ifelse(p < 0.1, ".", "")
#         )
#       )
#     )
#   )
# }
