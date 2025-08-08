#' Generic iterator function
#' @description Applies a function to each subset of data
#' @param data_list A list of data subsets (e.g., FUN or BAC)
#' @param name A string representing the name of the data (e.g, "Fungi" or "Bacteria")
#' @param func A function to apply to each subset
#' @param ... Additional arguments to pass to the function
#' @return A list of results from applying the function to each subset
#'
apply_to_subsets <- function(data_list, func, name = NULL, ...) {
  # Iterate over sources (root, soil)
  result_list <- lapply(names(data_list), function(source_name) {
    data <- data_list[[source_name]]

    # Apply the provided function to the innermost list element and names
    # Pass the subset data, name, source_name, and any extra args
    if (!is.null(name)) {
      result <- func(data, name, source_name, ...)
    } else {
      result <- func(data, ...)
    }

    # Return the result (assuming func returns the modified list element or results)
    return(result)
  })
  names(result_list) <- names(data_list) # Set names for the result list
  return(result_list)
}
