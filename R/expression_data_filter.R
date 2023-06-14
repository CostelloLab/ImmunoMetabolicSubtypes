#### This program preps the HTP data for downstream analysis
## Processing workflow:
##    - load count data
##    - log transform
##    - apply variance filter
##    - pivot to wider


library(tidyverse)
library(optparse
library(datawizard)

option_list <- list(
    make_option(c("-i", "--input_file"), type = "character", default = NULL, help = "counts input file, long format", metavar="character"),
    make_option(c("-o", "--output_directory"), type = "character", default = NULL, help = "output file location", metavar="character"),
    make_option(c("-c", "--clinical_metadata_file"), type = "character", default = NULL, help = "clinical_metadata_file", metavar="character"),
    make_option(c("-m", "--missingness"), type = "numeric", default = .2, help = "filter genes with missingness less than defined threshold, default is .2", metavar="numeric"),
    make_option(c("-v", "--variance"), type = "numeric", default = .2, help = "filter genes with log variance less than defined threshold, default is .2", metavar="numeric"))


opt_parse <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parse)

## load expression data
expression <- read.delim(opt$input_file)

## find the number of subjects
n <- length(unique(expression$LabID))

## filter for those with less than 20% missingness
missingness <- expression %>%
    group_by(EnsemblID) %>%
    summarise(missingness = sum(Value ==0)/n) 

miss_filtered <- expression %>%
    filter(EnsemblID %in% (missingness %>% filter(missingness < opt$missingness) %>% .$EnsemblID))

## logtranformation
miss_filtered$logValue <- log(miss_filtered$Value+1)


## filtering based on variance threshold
var_filtered <- miss_filtered %>%
    group_by(EnsemblID) %>%
    summarise(variance = var(logValue)) %>%
    filter(variance > opt$variance)
    

variance_filtered <- miss_filtered %>%
    filter(EnsemblID %in% var_filtered$EnsemblID) %>%
    as_tibble() 


### protein_filtered
prot_filtered <- variance_filtered %>%
    filter(Gene_type == "protein_coding")

## pivot wider. Have to include the unique identifier for row numbers to avoid any errors. 
htp_expr <- prot_filtered %>%
    dplyr::select(LabID, Gene_name,Value) %>%
    distinct(LabID, Gene_name, .keep_all = T) %>%
    pivot_wider(names_from = Gene_name, id_cols = LabID, values_from = Value)



### get the mappings of gene_name, ensembl, and chr from expression
gene_map <- expression %>%
    dplyr::select(Gene_name, EnsemblID,Chr) %>%
    filter(EnsemblID %in% variance_filtered$EnsemblID) %>%
    as_tibble()



## Adjust the transcriptomic data for age, sex, and sample source as the cytokine and metabolite data
metadata <- read.delim(opt$clinical_metadata_file)

expr_metadata <- htp_expr %>%
    inner_join(metadata, by = "LabID")

original_names <- names(htp_expr)
names(expr_metadata) <- make.names(names(expr_metadata))

to_exclude <- make.names(names(metadata))

mixed <- adjust(expr_metadata, effect = "Sample_source",
                exclude = to_exclude ,
                multilevel = T, keep_intercept = TRUE)

mixed <- adjust(mixed, effect = "Age_at_visit", exclude = to_exclude, keep_intercept = TRUE)
mixed <- adjust(mixed, effect = "Sex", exclude = to_exclude, keep_intercept = TRUE)

names(mixed) <- names(expr_metadata)
mixed <- mixed %>% dplyr::select(!(names(metadata)))

htp_expr <- mixed

expression_list<- list(expression = htp_expr,gene_map = gene_map)
save(expression_list, file = paste0(opt$output_directory,"/processed_transcriptome.RData"))



