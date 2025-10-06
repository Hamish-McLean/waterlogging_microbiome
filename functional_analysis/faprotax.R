# These are a set of functions to run faprotax and import the results into R

#inp <- data.table(taxonomy=apply(taxData[res.merge[padj<=ALPHA&log2FoldChange>0,OTU],-1],1,paste0,collapse = "; ")) # this is using the output from DESeq to grab lines from taxData - the format needs to be changed to a single column ";" seperated
# the above line might need adjusting to capture taxonomy (not confidence scores) info only

faprorun <- function(input, outfile, fapro_dir) {
  fwrite(input, "faprotaxa", sep = "\t", col.names = TRUE) # write to file
  system2(
    # "conda run -n FAPROTAX python",
    # "python --version"
    "conda",
    args = c(
      "run -n FAPROTAX python",
      here(fapro_dir, "collapse_table.py"),
      paste0("-i ", "faprotaxa"), # input taxonomy file
      paste0("-g ", fapro_dir, "FAPROTAX.txt"), # the FAPROTAX database location, this 
      paste0("-r ", outfile), # output file
      "-d 'taxonomy'", # column name containing taxonomy info
      # "--group_leftovers_as 'other'", # combine taxons which can't be assigned into their own group
      # "--verbose", # extra output to echo to the R console
      "--force"), # force over writing of the output file (otherwise it bombs out if the file exists)
    stdout = TRUE
  )
  outfile
}


parse_faprotax_raw <- function(file_path) {
  # Read file
  lines <- readLines(file_path)
  
  # Identify group header lines with "(x records)" and x > 0
  header_idx <- grep("\\([1-9][0-9]* records\\)", lines)  # matches "(1 records)" or more
  
  results <- lapply(seq_along(header_idx), function(i) {
    start <- header_idx[i]
    end <- ifelse(i < length(header_idx), header_idx[i + 1] - 1, length(lines))
    
    # Extract header line and clean it
    header_line <- str_trim(lines[start])
    
    # Function name (everything before " (x records)")
    function_group <- str_trim(str_remove(header_line, "\\s*\\([0-9]+ records\\)"))
    function_group <- str_remove(function_group, "^#\\s*|:$") # Remove leading '# '
    function_group <- str_remove(function_group, ":$")     # Remove trailing ':'
    
    # OTU/taxonomy lines are indented after header
    otus_raw <- lines[(start + 1):end]
    otus_raw <- otus_raw[grepl("^\\s", otus_raw)]  # keep only indented lines
    
    # Extract OTU names (last token after last ;)
    otu_names <- str_trim(otus_raw) %>%
      str_extract("[^;]+$")  # everything after last semicolon
    
    # Build result row
    data.frame(
      function_group = function_group,
      total_otus = length(otu_names),
    #   otu_names = paste(otu_names, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(results)
}


faprotaxit  <- function(taxData, outfile, fapro_dir){
  #. %>%
#   taxons <- subset(
#     taxData, 
#     row.names(taxData) %in% res.merge[padj <= ALPHA & filt(log2FoldChange), OTU]
#   ) %>% 
    # taxa_to_funguild(confidence = .8) %>%
  taxons <- taxData %>%
    mutate(taxonomy = paste(kingdom, phylum, class, order, family, genus, species, row.names(.), sep = ";") ) %>%
    dplyr::select(taxonomy) 
  if (nrow(taxons) > 0) {
    faprorun(taxons, outfile, fapro_dir) %>%
      parse_faprotax_raw()
  } else {
    print("No rows after select(taxonomy) - skipping faprotax run.") 
    #message("No rows after select(taxonomy) - skipping faprotax run.")
  }
}


#### RUN FAPROTAX ####

faprotax_main <- function() {

    library(data.table)
    library(here)
    library(tidyverse)

    fapro_dir <- here("/home/hmclean/apps/faprotax/FAPROTAX_1.2.12/")

    BAC  <- readRDS(here("objects", "BAC.rds"))
    taxData <- BAC$taxData

    outfile <- here("functional_analysis", "faprotax.txt")

    # set filter for increased abundance
    filt <- function(.).>0 # filt function will return T if > 0; filt(c(0,1,2,-3))# set to <0 for decreased abundence OTUs

    faproresults <- faprotaxit(taxData, outfile, fapro_dir)
    return(faproresults)
}
