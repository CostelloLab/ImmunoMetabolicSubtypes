#### This program preps the HTP data for downstream analysis
## Pipeline:
##    - load count data
##    - log transform
##    - apply variance filter
##    - pivot to wider

library(tidyverse)
library(caret)
library(sva)
## load expression data
expression <- read.delim("../data/HTP_WholeBlood_RNAseq_FPKMs_Synapse.txt")

## find the number of subjects
n <- length(unique(expression$LabID))

## filter for those with less than 20% missingness
missingness <- expression %>%
    group_by(EnsemblID) %>%
    summarise(missingness = sum(Value ==0)/n) 

miss_filtered <- expression %>%
    filter(EnsemblID %in% (missingness %>% filter(missingness < .2) %>% .$EnsemblID))

## logtranformation
expression$logValue <- log(expression$Value+1)


## filtering based on variance threshold
var_filtered <- miss_filtered %>%
    group_by(EnsemblID) %>%
    summarise(variance = var(logValue)) %>%
    filter(variance > .2)
    

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

expression_list<- list(expression = htp_expr,gene_map = gene_map)
save(expression_list, file = "../data/HTP_transcription_FPKMs_wide_missing_variance_filtered.Rdata")

