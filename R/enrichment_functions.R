clusterEnrichWrapper <- function(cluster_diff_sig, met_path) {

    num_clusters <- ncol(cluster_diff_sig) - 1

    results <- lapply(1:num_clusters, function(x) {
        col <- names(cluster_diff_sig)[grepl(x,names(cluster_diff_sig))]
        tmp <- cluster_diff_sig %>%
            select(c("feature",col)) %>%
            filter((!!sym(col)) < .1) %>%
            filter(feature %in% unique(met_path$name))
        sig_mets <- tmp$feature
        res <- metEnrich(sig_mets,met_path)
        return(res)
    })

    return(results)
}




metEnrich <- function(sig_mets, met_path) {
    metabolites = unique(met_path$name)
    num_met = length(metabolites)
    num_sig_met = length(sig_mets)
    results <- sapply(unique(met_path$Pathway), function(x) {
        tmp = met_path[met_path$Pathway ==x ,]
        num_path = dim(tmp)[1]
        if(num_path > 1){
            sig_in_path  = sum(tmp$name %in% sig_mets)
            sig_out_path =  num_sig_met - sig_in_path
            un_in_path = num_path - sig_in_path
            un_out_path = num_met - (sig_in_path+sig_out_path+un_in_path)
            cntg_tab = matrix(ncol = 2, nrow = 2,
                              c(sig_in_path, sig_out_path,un_in_path, un_out_path))
            res = fisher.test(cntg_tab, alternative = "greater")
            res_formatted <- as.data.frame(c(x,res$p.value, res$estimate))
            return(res_formatted)
        }
    })
    test <- suppressMessages(as.data.frame(t(bind_cols(results))))
    rownames(test) <- test[,1]
    names(test) <- c("pathway","p.value", "odd.ratio")
    test$p.value <- as.numeric(test$p.value)
    test <- test %>%
        arrange(p.value)
    return(test)
}
            
