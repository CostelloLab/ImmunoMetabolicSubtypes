library(tidyverse)
library(ggpubr)

LongCorr <- function(psych_object) {

    psych_object[["long"]] <- psych_object[["r"]] %>%
        as.data.frame() %>%
        rownames_to_column("A") %>%
        pivot_longer(cols = -A, names_to = "B") %>%
        filter(A!=B)

    psych_object$p_long <- psych_object$p %>%
        as.data.frame() %>%
        rownames_to_column("A") %>%
        pivot_longer(cols = -A, names_to = "B") %>%
        filter(A!=B)



    psych_object$long <- cbind(psych_object$long, p = psych_object$p_long$value, fdr = psych_object$p.adj)

    psych_object$long <- psych_object$long %>%
        filter(fdr < .1) %>%
        arrange(-value)

    names(psych_object$long)[3] <- "r"

    return(psych_object)
}

corrScatterPlotsSingle <- function(psych_object, molecule_class, omics_data) {

    options(ggpubr.parse_aes = FALSE)

    sapply(1:10, function(i) {
        tmp.1 = psych_object[["long"]][["A"]][i]
        tmp.2 = psych_object[["long"]][["B"]][i]
        tmp.plot =  ggscatter(data = omics_data[[molecule_class]],
                              x =  tmp.1,
                              y = tmp.2,
                              add = "reg.line",  # Add regression line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
                              ) +
            stat_cor( aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "spearman",label.sep = "\n", size = 8)+
            theme(axis.text = element_text(size =18),
                  axis.title = element_text(size =20))
        pdf(paste0("~/OneDrive - The University of Colorado Denver/Projects/ImmunoMetabolicSubtypes/results/correlation_scatter_plots/",molecule_class,"/",make.names(tmp.1), " - ", make.names(tmp.2), ".pdf"))
        print(tmp.plot)
        dev.off()
    })

     n = nrow(psych_object$long) + 1
    sapply(1:10, function(i) {

        tmp.1 = psych_object[["long"]][["A"]][n-i]
        tmp.2 = psych_object[["long"]][["B"]][n-i]
        tmp.plot =  ggscatter(data = omics_data[[molecule_class]],
                              x =  tmp.1,
                              y = tmp.2,
                              add = "reg.line",  # Add regression line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
                              ) +
            stat_cor( aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "spearman",label.sep = "\n", size = 8)+
            theme(axis.text = element_text(size =18),
                  axis.title = element_text(size =20))
        pdf(paste0("~/OneDrive - The University of Colorado Denver/Projects/ImmunoMetabolicSubtypes/results/correlation_scatter_plots/",molecule_class,"/",make.names(tmp.1), " - ", make.names(tmp.2), ".pdf"))
        print(tmp.plot)
        dev.off()
    })
    
}


corrScatterPlotsMultiple <- function(psych_object, molecule_classes, omics_data) {

    options(ggpubr.parse_aes = FALSE)

    tmp_data <- cbind(omics_data[[molecule_classes[1]]], omics_data[[molecule_classes[2]]])
    
    sapply(1:10, function(i) {
        tmp.1 = psych_object[["long"]][["A"]][i]
        tmp.2 = psych_object[["long"]][["B"]][i]
        tmp.plot =  ggscatter(data = tmp_data,
                              x =  tmp.1,
                              y = tmp.2,
                              add = "reg.line",  # Add regression line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
                              ) +
            stat_cor( aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "spearman",label.sep = "\n", size = 8)+
            theme(axis.text = element_text(size =18),
                  axis.title = element_text(size =20))
        pdf(paste0("~/OneDrive - The University of Colorado Denver/Projects/ImmunoMetabolicSubtypes/results/correlation_scatter_plots/multi/",make.names(tmp.1), " - ", make.names(tmp.2), ".pdf"))
        print(tmp.plot)
        dev.off()
    })

    n = nrow(psych_object$long) + 1
    sapply(1:10, function(i) {

        tmp.1 = psych_object[["long"]][["A"]][n-i]
        tmp.2 = psych_object[["long"]][["B"]][n-i]
        tmp.plot =  ggscatter(data = tmp_data,
                              x =  tmp.1,
                              y = tmp.2,
                              add = "reg.line",  # Add regression line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
                              ) +
            stat_cor( aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "spearman",label.sep = "\n", size = 8)+
            theme(axis.text = element_text(size =18),
                  axis.title = element_text(size =20))
        pdf(paste0("~/OneDrive - The University of Colorado Denver/Projects/ImmunoMetabolicSubtypes/results/correlation_scatter_plots/multi/",make.names(tmp.1), " - ", make.names(tmp.2), ".pdf"))
        print(tmp.plot)
        dev.off()
    })
    
}
s
