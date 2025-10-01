# Libraries
library(data.table)
library(here)
library(tidyverse)

FUNC_DIR <- here("functions")

source(here(FUNC_DIR, "metafuncs.R"))


DATA_DIR   <- normalizePath(here("..", "Data"))
OBJECT_DIR <- here("objects")

# taxData <- fread(here(DATA_DIR, "DADA2", "FUN.sintax.taxa"))

FUN  <- readRDS(here(OBJECT_DIR, "FUN.rds"))
taxData <- FUN$taxData


## functions
# converts a taxonomy file to SINTAX output
taxa_to_funguild <- function(taxData,confidence=.6) {
#   taxData <- phyloTaxaTidy(taxData[,-8],confidence)
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

# Setup FUNGuild

# renv::install("brendanf/FUNGuildR")
library(FUNGuildR)

# Setup the funguild database
# Final table will retain OTU as a separate column
# Later steps will subset the database based on OTUs reported as differential
guilds <- taxData %>% # taxData is the fungal taxonomy file
  taxa_to_funguild(confidence = .6) %>% #reassign nodes with a confidence <= to "confidence" to unknown
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


saveRDS(guilds, here("objects", "FUNGuilds.rds"))




# # DESeq analysis stuff...
# # dds is the DESeq object

# # get DESeq results
# res <- results(dds)

# # merge results with taxonomy table - this uses data.table, but the same can be done in base or dplyr
# res.merge <- as.data.table(res,keep.rownames="OTU")[as.data.table(taxData,keep.rownames="OTU"),on="OTU"]

# #  subset the guilds database, print to screen and write to file
# guilds %>%
#   subset(OTU%in%res.merge[padj<=ALPHA&log2FoldChange>0,OTU]) %T>% # select OTUs with sig. padj and LFC > 0 (i.e. increased abundance)
#   print() %>%
#   fwrite("funguild.increased.txt")