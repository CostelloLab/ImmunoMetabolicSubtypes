### This is a program that just combines the transcriptomic profiles with the processed omics data
### Is it a problem that the transcriptomic data is not processed like the cytokine and metabolite data?

library(tidyverse)
library(optparse)

## process the arguments
option_list <- list(
    make_option(c("-r", "--RNA"), type = "character", default = NULL, help = "processed transcriptome data ", metavar="character"),
    make_option(c("-t", "--T21_data"), type = "character", default = NULL, help = "processed T21 cytokine and metabolite data ", metavar="character"),
    make_option(c("-d", "--D21_data"), type = "character", default = NULL, help = "processed D21 cytokine and metabolite data ", metavar="character"),
     make_option(c("-c", "--clusterings"), type = "character", default = NULL, help = "clustering results ", metavar="character"),
    make_option(c("-o", "--output_directory"), type = "character", default = NULL, help = "output file location", metavar="character"))


opt_parse <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parse)

### Load the data
load(opt$RNA)
load(opt$T21_data)
load(opt$D21_data)
load(opt$clusterings)



## Further process the processed data
T21.omics.data <- rbind(t(T21.mixed[[1]]), t(T21.mixed[[2]]))
D21.omics.data <- rbind(t(D21.mixed[[1]]), t(D21.mixed[[2]]))
omics_data <- cbind(T21.omics.data, D21.omics.data)
data_names <- colnames(omics_data)

## standardize omics_data
omics_data <- as.data.frame(apply(omics_data, 1, scale))
rownames(omics_data) <- data_names
omics_data_clustering = as.data.frame(all_clust$clustering)
names(omics_data_clustering) <- "clustering"
omics_data <- merge(omics_data, omics_data_clustering, by = 0, all = TRUE)
omics_data$clustering[is.na(omics_data$clustering)] <-  "D21"
omics_data$clustering <- as.factor(omics_data$clustering)

## Add the RNAseq data                                 
exp <- expression_list$expression
## log transform the FPKMs
exp[, 2:ncol(exp)] <- log2(exp[, 2:ncol(exp)] + 1e-6)
## scale the log transformed FPKMs
exp[, 2:ncol(exp)] <- apply(exp[, 2:ncol(exp)], 2, scale)

exp <- exp %>%
    filter(LabID %in% omics_data$Row.names) %>%
    arrange(match(LabID, omics_data$Row.names))

omics_data <- omics_data %>%
    filter(Row.names %in% exp$LabID) %>%
    arrange(match( Row.names, exp$LabID))

## create dataset with all the omic measurements
omics_data <- exp %>%
    inner_join(omics_data, by = c("LabID" = "Row.names"))

## save the final dataset
save(omics_data, file = paste0(opt$output_directory,"/processed_omics_data.RData"))
