library(tidyverse)
library(ggpubr)
library(psych)

LongCorr <- function(psych_object) {

    psych_object[["long"]] <- psych_object[["r"]] %>%
        as.data.frame() %>%
        rownames_to_column("A") %>%
        pivot_longer(cols = -A, names_to = "B") %>%
        filter(A!=B) %>%
        rowwise() %>%
        mutate(cyt_met = list(sort(c(A,B)))) %>%
        ungroup() %>%
        distinct(cyt_met, .keep_all = T) %>%
        dplyr::select(-cyt_met)

    psych_object[["p_long"]] <- psych_object$p %>%
        as.data.frame() %>%
        rownames_to_column("A") %>%
        pivot_longer(cols = -A, names_to = "B") %>%
        filter(A!=B) %>%
        rowwise() %>%
        mutate(cyt_met = list(sort(c(A,B)))) %>%
        ungroup() %>%
        distinct(cyt_met, .keep_all = T) %>%
        dplyr::select(-cyt_met)



    psych_object$long <- psych_object$long %>%
        full_join(psych_object$p_long, by = c("A", "B")) %>%
        arrange(value.y)

    psych_object[["long"]] <- psych_object[["long"]] %>%
        arrange()


    names(psych_object[["long"]])[3:4] <- c("r", "fdr")

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



# pathEnrich is a function for identifying enriched metabolite pathways
pathEnrich = function(correlations, met_path) {
	
	enriched = data.frame(pathway = character(), p_value = numeric(), num_hits = numeric(), perc = numeric(), avg_corr = numeric())

	metabolites = unique(met_path$name)
	num_met = length(metabolites)
	num_sig_met = length(metabolites[metabolites %in% correlations[["B"]]])

	for(path in unique(met_path$Pathway)){
		tmp = met_path[met_path$Pathway ==path ,]
		num_path = dim(tmp)[1]
		if(num_path > 1){
                    sig_in_path  = sum(tmp$name %in% correlations[["B"]])
                    avg_corr <- correlations %>%
                        filter(B %in% tmp$name) %>%
                        summarise(res = mean(r)) %>%
                        .$res
                    
			perc = sig_in_path/num_sig_met
			sig_out_path =  num_sig_met - sig_in_path
			un_in_path = num_path - sig_in_path
			un_out_path = num_met - (sig_in_path+sig_out_path+un_in_path)

			cntg_tab = matrix(ncol = 2, nrow = 2,
							c(sig_in_path, sig_out_path,un_in_path, un_out_path))
			res = fisher.test(cntg_tab, alternative = "greater")
			met_p = res$p.value
			#met_OR  = res$estimate
			enriched = rbind(enriched, c(path,met_p, sig_in_path,perc,avg_corr))
		}

	}


	names(enriched) = c("pathway", "p.value", "num_hits","percent", "avg_corr")
    enriched[2:5] = apply(enriched[2:5],2,as.numeric)
    enriched$fdr <- p.adjust(enriched$p.value, method = "fdr")

	return(enriched)

}


correlationsEnrichmentWrapper <- function(correlations_long, met_path) {

    cytokines <- correlations_long %>%
        filter(fdr < .1) %>%
        distinct(A) %>%
        .$A
    
    results <- lapply(cytokines, function(x) {

        correlations <- correlations_long  %>%
            filter(A == x) %>%
            filter(fdr < .1) %>%
            as.data.frame()

        pathEnrich(correlations,met_path)
    })

    names(results) <- cytokines
    return(results)
}


enrichmentPlotPrep <- function(enrichment_list) {

    sig_enriched <- lapply(enrichment_list, function(x) {
        x %>%
            filter(fdr < .1)
    })

    cytokines <- sapply(sig_enriched,nrow)
    cytokines <- names(cytokines)[cytokines >0]
    pathways <- unique(unlist(sapply(sig_enriched, function(x) x$pathway)))

    toplot <- lapply(cytokines, function(x) {
        tmp <- enrichment_list[[x]]
        tmp$cytokine <- x
        tmp <- tmp %>%
            filter(pathway %in% pathways)
        return(tmp)
    })

    toplot <- toplot %>%
        bind_rows()

    return(toplot)
}

    
                     
                     
