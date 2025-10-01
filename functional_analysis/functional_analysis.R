# Libraries
library(data.table)
library(tidyverse)


## functions
# converts a taxonomy file to SINTAX output
taxa_to_funguild <- function(taxData,confidence=.6) {
  taxData <- phyloTaxaTidy(taxData[,-8],confidence)
  taxData$V <- gsub("(.*\\()([a-zA-Z])(\\))","\\2",taxData$rank)
  col_keep <- c(names(taxData)[1:7],"V")
  if(nrow(taxData) > 0) { 
    taxData <- as.data.frame(t(apply(taxData[,col_keep],1,function(x) {
      x[-1:-match(x[8],substr(col_keep[1:7],1,1))] <- "unidentified"
      x[-8]
    }))) 
  }
  return(taxData)
}



#~~~~~~~~~~~~~~~~~FUNGuild~~~~~~~~~~~~~~~~~~~~~# 


# Setup FUNGuild

#install.packages("remotes")
#remotes::install_github("brendanf/FUNGuildR")
library(FUNGuildR)

# Setup the funguild database
# Final table will retain OTU as a separate column
# Later steps will subset the database based on OTUs reported as differential
guilds <- taxData %>% # taxData is the fungal taxonomy file
  taxa_to_funguild(confidence = .8) %>% #reassign nodes with a confidence <= to "confidence" to unknown
  mutate(
    taxonomy = paste(kingdom, phylum, class, order, family, genus, species,row.names(.), sep = ",")
  ) %>%
  select(taxonomy) %>%
  unlist() %>%
  FUNGuildR::funguild_assign() %>%
  mutate(
    OTU=gsub(".*,","",Taxonomy), # add OTU as column
    Taxonomy=gsub(",OTU.*","",Taxonomy) # remove OTU from taxonomy
  )

# DESeq analysis stuff...
# dds is the DESeq object

# get DESeq results
res <- results(dds)

# merge results with taxonomy table - this uses data.table, but the same can be done in base or dplyr
res.merge <- as.data.table(res,keep.rownames="OTU")[as.data.table(taxData,keep.rownames="OTU"),on="OTU"]

#  subset the guilds database, print to screen and write to file
guilds %>%
  subset(OTU%in%res.merge[padj<=ALPHA&log2FoldChange>0,OTU]) %T>% # select OTUs with sig. padj and LFC > 0 (i.e. increased abundance)
  print() %>%
  fwrite("funguild.increased.txt")


#~~~~~~~~~~~~~~~~~~~Faprotax~~~~~~~~~~~~~~~~#


# These are a set of functions to run faprotax and import the results into R
# Requires wsl to be installed in Windows. 

# Path to run FAPROTAX from the current directory via wsl
PATH=gsub(".*OneDrive - NIAB","~/OneDrive",getwd()) # get the Windows path and convert to equivalent path under wsl

#inp <- data.table(taxonomy=apply(taxData[res.merge[padj<=ALPHA&log2FoldChange>0,OTU],-1],1,paste0,collapse = "; ")) # this is using the output from DESeq to grab lines from taxData - the format needs to be changed to a single column ";" seperated
# the above line might need adjusting to capture taxonomy (not confidence scores) info only

faprorun <- function(input,outfile) {
  fwrite(input,"faprotaxa",sep="\t",col.names = T) # write to file
  system2("wsl", args = c("collapse_table.py", # run "collapse_table.py" script under wsl - depending on setup/OS the script might be ble to be run directly
                          paste0("-i ","faprotaxa"), # input taxonomy file
                          "-g /usr/local/bin/FAPROTAX.txt", # the FAPROTAX database location, this 
                          paste0("-r ",outfile), # output file
                          "-d 'taxonomy'", # column name containing taxonomy info
                          # "--group_leftovers_as 'other'", # combine taxons which can't be assigned into their own group
                          # "--verbose", # extra output to echo to the R console
                          "--force") # force over writing of the output file (otherwise it bombs out if the file exists)
          , stdout = TRUE)
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
      otu_names = paste(otu_names, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(results)
}

faprotaxit  <- function(taxData){
  #. %>%
  taxons <- subset(taxData,row.names(taxData)%in%res.merge[padj<=ALPHA&filt(log2FoldChange),OTU]) %>%
    taxa_to_funguild(confidence = .8) %>%
    mutate(taxonomy = paste(kingdom, phylum, class, order, family, genus, species,row.names(.), sep = ";") ) %>%
    dplyr::select(taxonomy) 
  if (nrow(taxons) > 0) {
    faprorun(taxons,outfile = paste(RHB, "faprotax.increased", paste(contr, collapse = "_"), "txt", sep = ".")) %>%
      parse_faprotax_raw() %>%
      print()
  } else {
    print("No rows after select(taxonomy) - skipping faprotax run.") 
    #message("No rows after select(taxonomy) - skipping faprotax run.")
  }
} 

#### RUN FAPROTAX ####
# set filter for increased abundance
filt <- function(.).>0 # filt function will return T if > 0; filt(c(0,1,2,-3))# set to <0 for decreased abundence OTUs

faproresults <- faprotaxit(taxData)

