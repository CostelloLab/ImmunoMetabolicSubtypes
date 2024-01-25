#### Cluster evaluation
library(ggplot2)
library(limma)
library(RColorBrewer)
library(ggbeeswarm)
library(cluster)
library(clusterSim)
library(gridExtra)
library(circlize)
library(ComplexHeatmap)
library(openxlsx)
library(tidyverse)
library(randomcoloR)

### adjustedMutualInformation is a function that uses the aricode function "AMI" to create a table of the mutual information indeces between the lists of clusterings. 

adjustedMutualInformation = function(clustering_list, plot = F){

	all_names = lapply(clustering_list, names)
	keep = Reduce(intersect, all_names)

	clustering_list = lapply(clustering_list, function(x) {
					x = x[keep]
	})


	MI = lapply(1:length(clustering_list), function(x){
				lapply(1:length(clustering_list), function(y){
					ARI(clustering_list[[x]],clustering_list[[y]])
				})
		})

	if(is.null(names(clustering_list))){
		dimnames = paste("solution ", 1:length(clustering_list))
	} else {
		dimnames = names(clustering_list)
	}

	MI_tab = matrix(
					unlist(MI), 
					unique(lengths(MI)), 
					dimnames = list(dimnames,dimnames) 
			)

	if(plot == T){
		toplot = suppressWarnings(melt(MI_tab))
		toplot$value = round(toplot$value,2)
		AMI_plot = ggplot(toplot, aes( x= Var1, y = Var2, fill = value)) +
						 geom_tile() +
						 scale_fill_gradient(low = "lightblue", high = "darkblue") +
						 geom_text(aes(label = value)) +
						 theme_minimal()+
						 xlab("") + ylab("") + theme(legend.position = "none", 
						 						     plot.title = element_text(hjust = .5), 
						 						     axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1) ) +
						 ggtitle("Adjusted Mutual information")
		return(list(AMI= MI_tab, plot = AMI_plot))
	} else {
		return(list(AMI = MI_tab))
	}

}




#### pairwise_mi is a function that finds the average pairwise mutual information over a series of clustering solutions. 

pairwise_mi = function(clustering_list){
	n = 1
	pmi_tab = data.frame(item1 = character(), item2 = character(), ave_pmi = numeric() )

	for(z in 1:(length(clustering_list)-1)){

		tab1 = as.data.frame(clustering_list[[z]])
		names(tab1)[1] = 'clustering'
		design1 = as.data.frame(model.matrix(~0+as.factor(clustering), data = tab1))
		design1$LabID = rownames(design1)
		design1 = tibble(design1)
		

		for(x in 1:(length(clustering_list)-z)){
			tab2 = as.data.frame(clustering_list[[z+x]])
			names(tab2)[1] = 'clustering'
			design2 = as.data.frame(model.matrix(~0+as.factor(clustering), data = tab2))
			design2$LabID = rownames(design2)
			design2 = tibble(design2)
			



			for(col1 in names(design1)[-ncol(design1)]) {
				tmp1 = design1[,c(col1,"LabID")]
				names(tmp1)[1] = 'clustering' 
				
				for(col2 in names(design2)[-ncol(design2)]) {
					tmp2 = design2[,c(col2,"LabID")]
					names(tmp2)[1] = 'clustering'
					tab = rbind(tmp1,tmp2)
					pmi = pairwise_pmi(tab, LabID, clustering)
					pmi_tab =merge(x = pmi_tab,  y = pmi, by = c("item1", "item2"), all = T)
					pmi_tab$ave_pmi = apply(pmi_tab[,3:4],1, function(x)mean(x,na.rm = T ))
					pmi_tab = pmi_tab[, 1:3]
					print(n)
					n = n+1
				}

			}
		}

	}

	return(pmi_tab)

}


robust_clust <- function(clustering_list) {
    n <- 1
    pmi_tab <- data.frame(item1 <- character(), item2 <- character(), ave_pmi <- numeric() )

    sim <- list()
    for(z in 1:(length(clustering_list)-1)){

        tab1 <- as.data.frame(clustering_list[[z]])
        names(tab1)[1] <- 'clustering'
        design1 <- as.data.frame(model.matrix(~0+as.factor(clustering), data <- tab1))
        design1$LabID <- rownames(design1)
        design1 <- tibble(design1)
        

        for(x in 1:(length(clustering_list)-z)){
            tab2 <- as.data.frame(clustering_list[[z+x]])
            names(tab2)[1] <- 'clustering'
            design2 <- as.data.frame(model.matrix(~0+as.factor(clustering), data <- tab2))
            design2$LabID <- rownames(design2)
            design2 <- tibble(design2)
            



            for(col1 in names(design1)[-ncol(design1)]) {
                tmp1 <- design1[,c(col1,"LabID")]
                names(tmp1)[1] <- 'clustering' 
                
                for(col2 in names(design2)[-ncol(design2)]) {
                    tmp2 <- design2[,c(col2,"LabID")]
                    names(tmp2)[1] <- 'clustering'
                    tab <- merge(tmp1,tmp2, by <- "LabID")
                    score <- AMI(tab$clustering.x, tab$clustering.y)
                                        #					sub_tab <- tab[tab$clustering !<-0,]
                                        #					tot <- nrow(sub_tab)
                                        #					dup <- 2*sum(duplicated(sub_tab$LabID))
                                        #					score <- dup/tot
                    sim[[n]] <- data.frame(clust1 <- sprintf("solution %s: %s",z ,col1), clust2 <- sprintf("solution %s: %s",x+z ,col2) , score <- score)
                    n <- n+1

                }

            }
        }

    }

    res <- rbindlist(sim)
    res <- res[order(-res$score),]

    return(res)

}

global_effect = function(omics.list, clustering){
	feat.list = list()
	n = 1
	for(i in 1:length(omics.list)){
		for(z in 1:nrow(omics.list[[i]])){
			print(n)
			feat = rownames(omics.list[[i]])[z]
			tmp.omics.list = lapply(omics.list, function(x) x[!(rownames(x) %in% feat), ])
			num.clusts = dim(table(clustering))
			tmp.dropout.clust = nemo.clustering(tmp.omics.list, num.clusters =num.clusts )
			sim = robust_clust(clustering_list = list(clustering, tmp.dropout.clust$clustering))
			sim = sim[!duplicated(sim$clust1),]
			sim$clust = substr(sim$clust1,nchar(sim$clust1),nchar(sim$clust1))
			sim = sim[order(sim$clust),]
			feat.list[[n]] = sim[, c("clust", "score") ]
			n = n+1
		}
	}

	names(feat.list) = unlist(sapply(omics.list,rownames))

	weights = as.numeric(table(clustering))/ sum(as.numeric(table(clustering)))

	global_effect  = lapply(feat.list, function(x) 1 - weighted.mean(x$score,weights))
	global_effect = as.data.frame(unlist(global_effect))
	global_effect$feature = rownames(global_effect)


	return(global_effect)


}


global_effect2 = function(omics.list, clustering){
	feat.list = data.frame(feature = character(), score = numeric())
	n = 1
	for(i in 1:length(omics.list)){
		for(z in 1:nrow(omics.list[[i]])){
			print(n)
			feat = rownames(omics.list[[i]])[z]
			tmp.omics.list = lapply(omics.list, function(x) x[!(rownames(x) %in% feat), ])
			num.clusts = dim(table(clustering))
			tmp.dropout.clust = nemo.clustering(tmp.omics.list, num.clusters =num.clusts )
			sim = as.data.frame(unlist(AMI(clustering, tmp.dropout.clust$clustering)))
			sim$feature = feat
			names(sim)[1] = "score"
			feat.list[n,] = sim[, c("feature", "score") ]
			n = n+1
		}
	}


	feat.list$score = 1-feat.list$score

	return(feat.list)


}

global_effect3 = function(omics.list, clustering){
	feat.list = list()
	n = 1
	for(i in 1:length(omics.list)){
		for(z in 1:nrow(omics.list[[i]])){
			print(n)
			feat = rownames(omics.list[[i]])[z]
			tmp.omics.list = lapply(omics.list, function(x) x[!(rownames(x) %in% feat), ])
			num.clusts = dim(table(clustering))
			tmp.dropout.clust = nemo.clustering(tmp.omics.list, num.clusters =num.clusts )
			sim = robust_clust(clustering_list = list(clustering, tmp.dropout.clust$clustering))
			sim = sim[!duplicated(sim$clust1),]
			sim$clust = substr(sim$clust1,nchar(sim$clust1),nchar(sim$clust1))
			sim = sim[order(sim$clust),]
			feat.list[[n]] = sim[, c("clust", "score") ]
			n = n+1
		}
	}

	
	feat.list = lapply(feat.list, as.data.frame)
	effect = Reduce(cbind, feat.list)
	effect = effect[, c(1,seq(2,ncol(effect),2))]
	names(effect)[2:ncol(effect)] = unlist(sapply(omics.list,rownames))
	t.effect = as.data.frame(t(effect))
	toplot = melt(effect)


	#weights = as.numeric(table(clustering))/ sum(as.numeric(table(clustering)))


	return(toplot)


}


#### within_cluster_pmi is a function for finding the within cluster pairwise mutual information.
#### Requires inputs from clustering algorithm and pairwise_mi
within_cluster_pmi = function(pmi_tab, clustering){
	cluster_pmi = list()
	for(i in 1:length(unique(clustering))){
		test.cluster = clustering[clustering == i]
		tmp.pmi = pmi_tab[pmi_tab$item1 %in% names(test.cluster) & pmi_tab$item2 %in% names(test.cluster),]
		tmp.pmi =  as.data.frame(t(apply(tmp.pmi,1,sort)))
		tmp.pmi = tmp.pmi[!duplicated(tmp.pmi),]
		tmp.pmi$V1 = as.numeric(tmp.pmi$V1)
		cluster_pmi[[i]] = aggregate(tmp.pmi$V1, by = list(tmp.pmi$V2), function(x) mean(x, na.rm = T))
		cluster_pmi[[i]]$flag = ifelse(cluster_pmi[[i]]$x < 0, 1, 0 )
	}

	cluster_pmi = as.data.frame(rbindlist(cluster_pmi))
	names(cluster_pmi) = c("LabID", "within_clust_pmi", "flag")
	return(cluster_pmi)
}



##### permute_associations is a function for performing permutations over cluster association tests


association_test = function(clustering, clinic, phenotypes) {

    clustering = as.data.frame(clustering)
    dat = merge(clinic, clustering, by.x = "LabID", by.y = 0)



    res <- matrix(0,nrow =  length(unique(dat$clustering)), ncol = 2*length(phenotypes))
    dir = matrix(0,nrow =  length(unique(dat$clustering)), ncol = length(phenotypes))
    lcl = matrix(0,nrow =  length(unique(dat$clustering)), ncol = length(phenotypes))
    ucl = matrix(0,nrow =  length(unique(dat$clustering)), ncol = length(phenotypes))
    ratio = matrix(0,nrow =  length(unique(dat$clustering)), ncol = length(phenotypes))

    for(i in 1:length(phenotypes)) {

                                        # counter for results
	n = 1

                                        # test variable type
	if(length(unique(dat[,phenotypes[i]])) < 7| class(dat[,phenotypes[i]]) != "numeric"){ # categorical
            
            tab <- table(dat[ ,phenotypes[i]],dat$clustering)
            if(nrow(tab) == 1){next}

            cntg_tab <- matrix(0,nrow = nrow(tab), ncol = 2)

            for(j in unique(dat$clustering)){
                for(z in 1:nrow(tab)){
                    cntg_tab[z,1] <- sum(tab[z,])- tab[z,paste(j)]
                    cntg_tab[z,2] <- tab[z,paste(j)]

                }
                ft<-  fisher.test(cntg_tab, simulate.p.value = T)
 		if(sum(!is.na(unique(dat[,phenotypes[i]]))) ==2){
                    
                    res[j,(i*2-1)] <- sprintf("Clust:%s%% | All:%s%%",round((tab[2,paste(j)]/sum(tab[,paste(j)]))*100,0), round((sum(tab[2,])/sum(tab))*100,0))
                    res[j,(i*2)] <- ft$p.value
                    ratio[j,i] = tab[2,paste(j)]/sum(tab[, paste(j)])
                    dir[j,i] = ft$estimate
                    lcl[j,i] = ft$conf.int[1]
                    ucl[j,i] = ft$conf.int[2]
                    
                    


                }
                else { 
                    

                    ft<-  fisher.test(cntg_tab, simulate.p.value = T)
                    rownames(cntg_tab) <- names(tab[,paste(j)])
                    tmp = rownames(cntg_tab)
                    perc = round((cntg_tab[,1]/sum(cntg_tab[,1]))*100,2)
                    est = paste0(tmp,":",perc, "%")

                    res[j,(i*2-1)] <- paste(est, collapse = " | ")
                    res[j,(i*2)] <- ft$p.value
                    ratio[j,i] = tab[2,paste(j)]/sum(tab[, paste(j)])
                    dir[j,i] = "not calculated"
                    lcl[j,i] = NA
                    ucl[j,i] = NA
                    
                    
                }
            }
            

	} else { # numeric
            for(j in unique(dat$clustering)){

                tt <- t.test(dat[dat$clustering == j, phenotypes[i]], dat[dat$clustering != j, phenotypes[i]])
                res[n,(i*2-1)] <- sprintf("Clust %s mean:%s | Other %s mean:%s", phenotypes[i], round(mean(dat[dat$clustering == j, phenotypes[i]],na.rm = T), 1), phenotypes[i], round(mean(dat[dat$clustering != j, phenotypes[i]], na.rm = T),1))
                
                res[j,(i*2)] <-tt$p.value
                ratio[j,i] = mean(dat[dat$clustering == j, phenotypes[i]],na.rm = T)
                dir[j,i] =round(mean(dat[dat$clustering == j, phenotypes[i]],na.rm = T), 1) - round(mean(dat[dat$clustering != j, phenotypes[i]],na.rm = T), 1)
                lcl[j,i] = NA
                ucl[j,i] = NA
                n = n+1		
            }
            

	}
        if(nrow(res > 2)){res[,(i*2)] = p.adjust(res[,(i*2)], method = "fdr")}
        
    }

    name = character()
    for(i in 1:length(phenotypes)){
	name[i*2-1] <- paste(phenotypes[i],"values")
	name[i*2] <- paste(phenotypes[i],"q.value")
    }
    colnames(res) <- name
    pvals= apply(res[,seq(2,ncol(res),2)],2, as.numeric)
                                        # dir = apply(dir, 2, as.numeric)
                                        # lcl = apply(lcl, 2, as.numeric)
                                        # ucl = apply(ucl, 2, as.numeric)


    rownames(res) <- paste("cluster", 1:length(unique(dat$clustering)))
    rownames(pvals) <- paste("cluster", 1:length(unique(dat$clustering)))
    rownames(dir) <- paste("cluster", 1:length(unique(dat$clustering)))
    rownames(ratio) = paste("cluster", 1:length(unique(dat$clustering)))
    colnames(dir) = paste(phenotypes, "OR")
    colnames(lcl) = paste(phenotypes, "lcl")
    colnames(ucl) = paste(phenotypes, "ucl")
    colnames(ratio) = phenotypes


                                        # comining pvals and dir (Odds ratio) into a single table
    dir.pvals = cbind(dir,pvals)
    dir.pvals = dir.pvals[order(rownames(dir.pvals)), order(colnames(dir.pvals))]


    dir.pvals.confint = cbind(dir.pvals,lcl,ucl)
    dir.pvals.confint = dir.pvals.confint[order(rownames(dir.pvals.confint)), order(colnames(dir.pvals.confint))]


    return(list(full.res = res, pvals = pvals, dir = dir, dir.pvals = dir.pvals, dir.pvals.confint = dir.pvals.confint, ratio = ratio ))

}


#### permute_associations is a function for performing permutations over cluster association tests



permute_associations = function(clustering, clinic, phenotypes, comorbidity_meta, num_permutations){
	pvals = association_test(clustering, clinic, phenotypes)$pvals

	tab = table(clustering)
	prob = tab/sum(tab)
	pdist = lapply(1:num_permutations,function(x) {
		random.clustering = sample(x = clustering )
		names(random.clustering) = names(clustering);

		association_test(random.clustering, clinic, phenotypes)$pvals

		})
	empiricalP = pvals
	for(y in 1:nrow(pvals)) {
		for(z in 1:ncol(pvals)) {
			if(pvals[y,z] < .3 ){
				empiricalP[y,z] = sum(unlist(lapply(pdist, function(x) x[y,z] < pvals[y,z])))/sum(length(pdist),1)
			} else{
			 	empiricalP[y,z] = 1
			}

		}
	}

	#empiricalP = apply(empiricalP,2, function(x) p.adjust(x, method = "fdr"))
	return(empiricalP)
}


# tmp = unlist(lapply(pdist, function(x) x[2,3]))
# plot(density(tmp), main = sprintf("cluster 2 Karyotype - empiricalP = %s", empiricalP[2,3]), 
# 	xlab = "p value")
# abline(v = pvals[2,3], col = "red")


# tmp = unlist(lapply(pdist, function(x) x[2,4]))
# plot(density(tmp), main = sprintf("cluster 2 Age - empiricalP = %s", empiricalP[2,4]), 
# 	xlab = "p value")
# abline(v = pvals[2,4], col = "red")

##### metap_condition is a function for finding the metap value by cluster across condtions
metap_cluster = function(cluster_assoc){
	meanp = suppressWarnings(as.data.frame(apply(cluster_assoc,1, function(x) sumz(x)$p)))
	names(meanp) = "cluster_metap"
	meanp[is.na(meanp)] = 0
	return(meanp)
}


##### metap_condition is a function for finding the metap value by cluster across condtions
metap_condition = function(cluster_assoc, weights){
	meanp = suppressWarnings(as.data.frame(apply(cluster_assoc,2, function(x) sumz(x, weights = weights)$p)))
	names(meanp) = "condition_metap"
	meanp[is.na(meanp)] = 1
	return(meanp)
}

###### sample_parameters is a function for sampling number of neighbors and clusters in nemo clustering/

sample_param = function(omics.list, num.clusters = NULL , num.neighbors = NA, plot = T, clinic, phenotypes, comorbidity_meta, num_permutations = 1e3 ){
	
	metap.cond = list()
	modularity = list()

	if(!is.null(num.clusters)) {
		for(i in num.clusters) {
			print(paste("Testing", i, "clusters"))
			clusts = nemo.clustering(omics.list, num.clusters = i, num.neighbors = NA)
			tmp.clusts = clusts$clustering			
			tmp.cluster_assoc = permute_associations(tmp.clusts,clinic, phenotypes, comorbidity_meta, num_permutations)
			weights = sqrt(table(tmp.clusts))
			metap.cond[[i]] = metap_condition(tmp.cluster_assoc, weights)
			g = graph_from_adjacency_matrix(clusts$graph, mode = "undirected", weighted = T)
			g = simplify(g)
			mod  = modularity(g, membership = tmp.clusts )
			modularity[[i]] = mod
		}

		metap.cond = metap.cond[num.clusters]
		metap.cond = Reduce(cbind, metap.cond)
		names(metap.cond) = paste(num.clusters,"clusters" )
	
		if(plot) {
			metap.cond$condition = rownames(metap.cond)
			toplot = melt(metap.cond, id.vars = 'condition')
			cols = distinctColorPalette(k = length(phenotypes))
			toplot$condition = gsub("p.value","", toplot$condition )
			toplot$value = ifelse(toplot$value == 0, -log10(min(toplot$value[toplot$value != 0])), -log10(toplot$value) )
			toplot$value = ifelse(toplot$value < -log10(.1), NA, toplot$value)


			p1 = ggplot(toplot, aes(x= variable, y = condition, fill = value)) +
				geom_tile() + theme_classic() + 
				theme(legend.position = "bottom" ) +
				ylab("meta p value") + xlab("") +
				ggtitle("") + labs(fill = "-log10(metap)")

		# 	p1 = ggplot(toplot, aes(x= variable, y = value, group = condition )) +
		# 		geom_line(aes(color = condition),size = 1.2) + 
		# 		geom_point(aes(color = condition), size = 2) + 
		# 		scale_color_manual(values = cols)+
		# 		theme_classic() + 
		# 		theme(legend.position = "bottom") +
		# 		ylab("meta p value") + xlab("") +
		# 		ggtitle("") 
		}
	}
	
	if(length(num.neighbors) > 1) {
		for(i in num.neighbors) {
			print(paste("Testing", i, "neighbors"))
			tmp.clusts = nemo.clustering(omics.list, num.clusters = NULL, num.neighbors = i)$clustering
			tmp.cluster_assoc = permute_associations(tmp.clusts,clinic, phenotypes,comorbidity_meta, num_permutations)$pvals
			metap.cond[[i]] = metap_condition(tmp.cluster_assoc)
		}

		metap.cond = metap.cond[num.neighbors]
		metap.cond = Reduce(cbind, metap.cond)
		names(metap.cond) = paste(num.neighbors,"neighbors" )
	
		if(plot) {
			metap.cond$condition = rownames(metap.cond)
			toplot = melt(metap.cond, id.vars = 'condition')
			cols = distinctColorPalette(k = length(phenotypes))
			p1 = ggplot(toplot, aes(x= variable, y = value, group = condition )) +
				geom_line(aes(color = condition),size = 1.2) + 
				geom_point(aes(color = condition), size = 2) + 
				scale_color_manual(values = cols)+
				theme_classic() + 
				theme(legend.position = "bottom") +
				ylab("meta p value") + xlab("") +
				ggtitle("") 
		}
	}
	



	return(list(metap.table = metap.cond, metap.plot = p1, modularity = unlist(modularity)))

}

#### diff_expr is a function for finding the deferentially expressed omics levels by cluster.
diff_expr <- function(omics.data,clustering){
	all.siglist = list()

	meta = as.data.frame(clustering)
	design = model.matrix(~0+as.factor(clustering), data = meta)
	colnames(design) = paste("c", 1:length(unique(clustering)), sep = "")

	for(i in 1:length(unique(clustering))){
		c0 = sprintf("contrast =  makeContrasts(contrast= c%s, levels = design)",i)
		eval(parse(text = c0))	
		fit = eBayes(contrasts.fit(lmFit(omics.data,design),contrast))
		all.siglist[[i]] = topTable(fit, adjust = 'f', number = nrow(omics.data)+1,
                          confint = T)
		all.siglist[[i]]$feature = rownames(all.siglist[[i]])
		names(all.siglist[[i]])[1:8] = paste0(i,".",names(all.siglist[[i]])[1:8])
	}

	res = Reduce(function(x,y) merge(x,y, by = "feature"), all.siglist)


	return(res)
}



#### diff_expr is a function for finding the differentially expressed omics levels by cluster.
diff_expr_cluster <- function(omics.data,clustering){
	cluster.siglist = list()
	meta = as.data.frame(clustering)
	design = model.matrix(~0+as.factor(clustering), data = meta)
	colnames(design) = paste("c", 1:length(unique(clustering)), sep = "")


	for(i in 1:length(unique(clustering))){
		sig.list = list()
		for(j in (1:(length(unique(clustering))))[-i]){
			c0 = sprintf("contrast =  makeContrasts(contrast= c%s - c%s, levels = design)",i,j)
			eval(parse(text = c0))	
			fit = eBayes(contrasts.fit(lmFit(omics.data,design),contrast))
			sig.list[[j]] = topTable(fit, adjust = 'fdr', number = nrow(omics.data)+1,
	                          confint = T)
			sig.list[[j]] = sig.list[[j]][sig.list[[j]]$adj.P.Val < .05, ]
		}
		sig.list = sig.list[(1:(length(unique(clustering))))[-i]]
		sig.rows = lapply(sig.list, rownames)
		sig = Reduce(intersect, sig.rows)
		cluster.siglist[[i]] = sig
	}

	return(cluster.siglist)

}


#### wilcoxon differential expression test

diff_expr_wilcoxon <- function(omics.data,clustering){
	cluster.siglist = list()
	meta = as.data.frame(clustering)
	design = model.matrix(~0+as.factor(clustering), data = meta)
	colnames(design) = 1:length(unique(clustering))

	results = list()
	for(i in colnames(design)){
		results[[i]] = diff_expr_wrapper(omics.data,design, i, 0, test = "wilcoxon")
		results[[i]]$feature = rownames(results[[i]])
	}

	res =  suppressWarnings(Reduce(function(x,y) merge(x,y, by = "feature"), results))
	sig = res[,  grepl("q.value|feature",names(res) )]
	sig[,2:ncol(sig)] = t(apply(sig[,2:ncol(sig)], 1, function(x) p.adjust(x,method = "fdr")))
	names(sig)[2:ncol(sig)] = paste0("q.value.",colnames(design))

	FC = res[, grepl("FC|feature",names(res) )]
	names(FC)[2:ncol(FC)] = paste0("FC.",colnames(design))


	return(list(sig = sig, FC = FC))

}





#### diff_expr is a function for finding the differentially expressed omics levels by cluster.
diff_expr_cluster_ranked <- function(omics.data,clustering) {
	cluster.siglist = list()
	meta = as.data.frame(clustering)
	design = model.matrix(~0+as.factor(clustering), data = meta)
	colnames(design) = paste("c", 1:length(unique(clustering)), sep = "")

	sig.list = list()
	for(i in 1:length(unique(clustering))){
		
	
			comp.group = paste0("c", (1:(length(unique(clustering))))[-i], collapse = "+")
			c0 = sprintf("contrast =  makeContrasts(contrast= c%s - ((%s)/%s), levels = design)",i,comp.group, (length(unique(clustering))-1))
			eval(parse(text = c0))	
			fit = eBayes(contrasts.fit(lmFit(omics.data,design),contrast))
			sig.list[[i]] = topTable(fit, adjust = 'fdr', number = nrow(omics.data)+1,
	                          confint = T)
			#sig.list[[j]] = sig.list[[j]][sig.list[[j]]$adj.P.Val < .05, ]
	}
	sig.list = lapply(sig.list, function(x) cbind(x, feature = rownames(x)))
	sig = Reduce(function(x,y) merge(x,y, by = "feature"), sig.list)
	sig = sig[, c(1, grep("dj.P.Val", names(sig)))]
	names(sig)[2:ncol(sig)] = paste("Cluster", 1:(ncol(sig)-1))
	

	return(sig)

}

uniqueFeatures = function(diff.expr.cluster){
	sig_features = list()
	for(i in 2:ncol(diff.expr.cluster)){
		sig_features[[i]] = diff.expr.cluster$feature[diff.expr.cluster[,i] < 1e-5]
	}

	sig_features = sig_features[-1]

	unique_features = list()
	for(i in 1:length(sig_features) ){
		common_features = unlist(sig_features[-i])
		unique_features[[i]] = sig_features[[i]][!sig_features[[i]] %in% common_features]
	}

	return(unique_features)

}

### plot_diff_expr
plot_diff_expr = function(diff.expr, output_dir){
	for(i in 1:length(diff.expr)){
		p = EnhancedVolcano(diff.expr[[i]], 
				lab = rownames(diff.expr[[i]]), 
				x = 'logFC', y = 'adj.P.Val', 
				FCcutoff = .25, 
				pCutoff = .05,
				labSize = 4, 
				title = paste("Cluster",i), 
				drawConnectors = F, 
				boxedLabels = F)
	png(sprintf("%s/Cluster%s-%s.png", output_dir, i, Sys.Date()), units = 'in', res = 100, height =8, width = 8 )
	print(p)
	dev.off()
	}
}



plot_diff_expr_cluster = function(diff.expr.cluster, omics.data, clustering, output_dir){

	tmp.omics.data = as.data.frame(t(omics.data))
	tmp.clustering = as.data.frame(clustering)
	tmp.omics.data = merge(tmp.omics.data, tmp.clustering, by = 0)
	tmp.omics.data$clustering = as.factor(tmp.omics.data$clustering)
	
	for(i in 1:length(diff.expr.cluster)){
		if(length(diff.expr.cluster[[i]]) > 0){
		tmp.cluster = diff.expr.cluster[[i]]
	} else (next)
		for (j in 1: length(tmp.cluster)){

			p1 = ggplot(tmp.omics.data) + 
				geom_boxplot(aes_string(x = "clustering", y = sprintf("`%s`",tmp.cluster[[j]])), fill = "lightblue", alpha = .5) + 
				geom_beeswarm(aes_string(x = "clustering", y = sprintf("`%s`",tmp.cluster[[j]])))+
				coord_flip() + 
				geom_boxplot(data = tmp.omics.data[tmp.omics.data$clustering == i,], aes_string(x = "clustering", y = sprintf("`%s`",tmp.cluster[[j]])), fill = "red", alpha = .5) + theme_minimal()
				feat = sprintf("%s",tmp.cluster[[j]])
				feat = gsub("/", "--", feat)
			png(sprintf("%sCluster_%s_%s.png", output_dir, i, feat), units = 'in', res = 100, height =4, width = 6 )
			print(p1)
			dev.off()


		}


	}



}





###### plot_tsne
plot_tsne = function(graph, clustering,clinic, plot_var, output_file){
	tsne = Rtsne(graph, theta = 0, perplexity = 3, max_iter = 5000)	
	rownames(tsne$Y) <- rownames(graph)
	tsne <- merge(tsne$Y, as.data.frame(clustering), by= 0)
	tsne <- merge(tsne, clinic, by.x = "Row.names", by.y = "LabID" )
	names(tsne)[2:3] <- c("tsne1", "tsne2")
	tsne$clustering <- as.factor(tsne$clustering)
	cols = distinctColorPalette(k = length(unique(tsne$clustering)))
	p1 <- ggplot(tsne, aes(x = tsne1, y = tsne2, col = 	tsne[, plot_var])) + 
			geom_point() + 
			theme_classic() + 
			scale_color_manual(values = cols) + 
			labs(col = plot_var)

	svg(output_file)
	print(p1)
	dev.off()


}


plot_umap = function(graph, clusterings,clinic, plot_var, title, output_file){
	umap = umap(graph, theta = 0)	
	umap <- merge(umap$layout, clinic, by.x = 0, by.y = "LabID" )
	names(umap)[2:3] <- c("umap1", "umap2")
        plot_list <- list()
        cols = distinctColorPalette(k = max(sapply(clusterings, function(x) length(table(x)))))
        
        
        plot_list <- sapply(1:length(clusterings), function(i) {
            umap_tmp <- umap
            clustering = clusterings[[i]]
            umap_tmp <- merge(umap_tmp, as.data.frame(clustering), by.x = "Row.names", by.y = 0)
            umap_tmp$clustering <- as.factor(umap_tmp$clustering)
            ggplot(umap_tmp, aes(x = umap1, y = umap2, color = 	umap_tmp[, plot_var])) + 
                geom_point(size =3) + 
                guides(color =guide_legend(override.aes = list(size = 5) ) )+
                theme_classic() + 
                scale_color_manual(values = cols) + 
                labs(col = plot_var)+ 
                theme(legend.position = "bottom",
                      plot.title = element_text(hjust = 0.5 )) +
                guides(size = "none")+
                ggtitle(names(clusterings[i]))
        },simplify  = F)
            

       p1 <-  grid.arrange(grobs = plot_list, ncol = 3)

        ggsave(p1, file = output_file,width = 12, height = 5)


}

### clusterEval  is a wrapper function for creating a table to evaluate number of clusters



clusterEval <- function(omics.list, clinic, phenotypes, iterations = 100, NUMC = 2:7,  diff_sig_threshold = .01, clinical_threshold = .1, num_neighbors = 100) {

    W <- nemo.affinity.graph(omics.list, k = num_neighbors)
    omics.data <- as.data.frame(rbindlist(lapply(omics.list,as.data.frame)))
    rownames(omics.data) <- unlist(sapply(omics.list,rownames))

    min_cluster <- .05 * ncol(omics.data)

    final_NUMC <- numeric()
    results <- list()
    feature_diff <- list()
    condition_diff <- list()
    
    for(i in NUMC) {
        print(paste("Testing" ,i, "Clusters"))

        clustering <- spectralClustering(W, i)
        names(clustering) <- colnames(W)

        if(min(table(clustering)) > min_cluster) {final_NUMC <- c(final_NUMC,i)}
        
        cluster_eval <- association_test(clustering <- clustering, clinic, phenotypes)
        condition_diff[[i]] <- cluster_eval
        num_sig <- sum(cluster_eval$pvals[,!(colnames(cluster_eval$pvals) %in% c("Sex q.value","Age_at_visit q.value"))]  < clinical_threshold)
        if(i ==2){ num_sig <- num_sig/2}                                         # Divide by two in the case of 2 clusters, because each cluster will have the same score
            
        diff.expr.cluster <- diff_expr_wilcoxon(omics.data,clustering)
        feature_diff[[i]] <- diff.expr.cluster
        num_diff <- diff.expr.cluster$sig
        num_diff <- sum(apply(num_diff[2:ncol(num_diff)],1, function(x) any(x < diff_sig_threshold)))
        DB <- index.DB(x = W, cl = clustering)$DB ## Davies-Bouldin
        CH <- index.G1(x = W, cl = clustering) ## Calinkski-Harabasz

        bootstrap <- boostrapCluster(omics.list = omics.list,
                                     iterations = iterations,
                                     num.clusters = i,
                                     num.neighbors = num_neighbors)
        
        results[[i]] <- c(i, num_sig,num_diff, min(table(clustering)), DB, CH, bootstrap)
    } 
    

    eigengap <- estimateNumberOfClustersGivenGraph(W,final_NUMC)
    nemo_num <- nemo.num.clusters(W = W, NUMC = final_NUMC)

    final.results <-  data.frame(do.call(rbind,results))
    names(final.results) <- c("# of clusters", "# of enriched conditions","# of differential features", "Minimum cluster size", "DB", "CH", "bootstrap")


    return(list(summary  = final.results,  eigengap = eigengap, nemo = nemo_num, feature = feature_diff, condition = condition_diff))

}



 

featureSelection = function(all_clust,omics.data, diff_expr) {

	library(caret)
	library(randomForest)


	# define the control using a random forest selection function

	# create dummy variables for clustering 
	clustering = as.data.frame(all_clust$clustering)
	design = model.matrix(~0+as.factor(all_clust$clustering), data = clustering)
	colnames(design) = paste0("Cluster", 1:ncol(design))

	features = base::t(omics.data)
	identical(rownames(features), rownames(design))

	# define the control using a random forest selection function
	control <- rfeControl(functions=rfFuncs, method="cv", number=10, repeats = 5)

	important_features = list()

	for(i in 1:ncol(design)) {
		tmp.features = diff_expr$sig$feature[diff_expr$sig[,i+1] < .05]
	#	tmp.features = tmp[[i]]


		important_features[[i]] = rfe(features[,tmp.features], as.factor(design[,i]), sizes=c(1:length(tmp.features)), rfeControl=control)	
		# feat_mod_select=	cv.glmnet(features, design[,i], standardize = FALSE, alpha = 1  )
		# test = as.matrix(coef(feat_mod_select, feat_mod_select$lambda.min))	
	}

	predictors = lapply(important_features, predictors)
	n = max(sapply(predictors, length))
	predictors_df = data.frame(cluster1 = rep(NA, n))

	for(x in 1:length(predictors)) {
		z = predictors[[x]]
		length(z) = n	
		predictors_df = cbind(predictors_df, z)
		#predictors_df[1:length(predictors[[x]]),x] = predictors[[x]]
		tmp = diff_expr$FC[diff_expr$FC$feature %in% predictors[[x]], c(1,x+1)] 
		tmp = tmp[order(match(tmp$feature, predictors[[x]])),2]
		length(tmp) = n
		predictors_df = cbind(predictors_df,tmp)
	}

	predictors_df = predictors_df[,-1]

	clusters =  paste0("cluster",1:length(predictors))
	FCs =  paste0("log2FC",1:length(predictors))
	names(predictors_df) = c(rbind(clusters,FCs))
	return(predictors_df)
}





### bootstrap clustering
boostrapCluster <- function(omics.list, iterations = 100, num.clusters, num.neighbors) {
    W <- nemo.affinity.graph(omics.list,num.neighbors)
    gold.clustering <- spectralClustering(W,num.clusters)

    bootstrap_results <- sapply(1:iterations, function(x) {
        rand <- sample(colnames(omics.list[[1]]), replace = T)
        rand.omics.list <- lapply(omics.list, function(x) x[,rand])
        rand.W <- nemo.affinity.graph(rand.omics.list, num.neighbors)
        rand.gold.clustering <- spectralClustering(rand.W, num.clusters)
        return(comparing.Partitions(gold.clustering, rand.gold.clustering, type = "rand"))
    })
    return(mean(bootstrap_results, na.rm = T))
}



## heatmap of clusters
clusterHeatmap = function(omics.data,clustering,diff_expr, cluster_assoc,threshold, phenotypes, split_num, clinical_threshold, title, fontsize, dendrogram = TRUE,enrichment = NULL, column_names = TRUE, legend_pos_auto = TRUE ){
	library(ComplexHeatmap)
	library(RColorBrewer)
	library(circlize)

	omics = t(omics.data)
	sig = apply(diff_expr$sig[,2:ncol(diff_expr$sig)],2, function(x) diff_expr$sig$feature[x<threshold])
	sig = unlist(sig)
	unique_features = sig[!duplicated(sig)]

	FC = diff_expr$FC
	rownames(FC) = FC$feature 
	FC = FC[,-1]
	FC = FC[colnames(omics)[colnames(omics) %in% unique_features],]


	padj_mat = diff_expr$sig
	rownames(padj_mat) = padj_mat$feature
	padj_mat = padj_mat[,-1]
	padj_mat = padj_mat[colnames(omics)[colnames(omics) %in% unique_features],]
	padj_mat = t(padj_mat)

	padj_fun <- function(j, i, x, y, width, height, fill) {

	  if(padj_mat[i, j] < .01)
	    grid.text(
	      # "\u066D", # asterisk but not very well centered
	      "*", # asterisk but not very well centered
	      x, y, gp = gpar(fontsize = 8))

	}


names(FC) = paste(1:length(unique(clustering)))

 	clustering = as.factor(clustering)
	ha = rowAnnotation(
		Cluster = clustering, 
		col = list(Cluster = c("1" = "#7FC97F", "2" = "#BEAED4", "3" = "#FDC086", "4" = "#FFFF99", "5" = "#386CB0", "6" = "#F0027F","7" = "#BF5B17", "8" ="#666666")))

	col_fun_main = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
	if(dendrogram){
		ht <- Heatmap(t(FC),
				name = "log2FC", 
				show_row_names = FALSE, 
				show_column_names = column_names, 
				show_row_dend = FALSE, 
				row_dend_reorder = FALSE,
				column_title = title, 
				row_title = "", 
				column_km = split_num, 
				row_split = rownames(t(FC)), 
			#	column_names_rot = 60, 
				column_names_gp = gpar(fontsize = fontsize),
				heatmap_legend_param = list(legend_direction = "horizontal"),
				border_gp = gpar(col = "black", lty = 1),
				#cell_fun = padj_fun,
				row_names_side = "left",
				row_order = rownames(t(FC)),
				show_heatmap_legend = ,
                              col = col_fun_main,
                              height = ncol(FC) * unit(5,"mm")

			)
	} else{
		ht <- Heatmap(t(FC),
				name = "log2FC", 
				show_row_names = TRUE, 
				show_column_names = column_names, 
				row_order = rownames(t(FC)),
				#right_annotation = ha, h_list
				 
				column_title = title, 
				row_title = "Clusters", 
				column_km = split_num, 
				row_split = rownames(t(FC)), 
			#	column_names_rot = 60, 
				column_names_gp = gpar(fontsize = fontsize),
				heatmap_legend_param = list(legend_direction = "horizontal"),
				show_row_dend = FALSE, 
				row_dend_reorder = FALSE,
				show_column_dend = dendrogram, 
				column_dend_reorder = dendrogram,
				border_gp = gpar(col = "black", lty = 1),
				#cell_fun = padj_fun
				show_heatmap_legend = legend_pos_auto,
                              col = col_fun_main,
                              height = ncol(FC) * unit(5,"mm")

			)
	}

	#	ht = draw(ht)
	cluster_assoc$dir[is.infinite(cluster_assoc$dir)] = 1.5

	assoc_toplot = as.matrix(cluster_assoc$pvals)
	tokeep = apply(assoc_toplot, 2, function(x) any(x< .5))
	assoc_toplot = assoc_toplot[, tokeep]
	assoc_toplot = -log10(assoc_toplot)
	colnames(assoc_toplot) = gsub(" q.value","", colnames(assoc_toplot) )
	
	cluster_assoc$dir = cluster_assoc$dir[, tokeep]
	colnames(cluster_assoc$dir) = colnames(assoc_toplot)

	assoc_toplot = assoc_toplot[, colnames(assoc_toplot) %in% phenotypes]
	cluster_assoc$dir = cluster_assoc$dir[, colnames(cluster_assoc$dir) %in% phenotypes]



	for(i in colnames(assoc_toplot)){
		if(min(cluster_assoc$dir[,i]) < 0){
			assoc_toplot[,i] = ifelse(cluster_assoc$dir[,i] < 0, assoc_toplot[,i] * -1, assoc_toplot[,i])
		} else{
			assoc_toplot[,i] = ifelse(cluster_assoc$dir[,i] < 1, assoc_toplot[,i] * -1, assoc_toplot[,i])
		}
	}


	assoc_mat = cluster_assoc$pvals[,grep(paste(colnames(assoc_toplot), collapse = "|"), colnames(cluster_assoc$pvals))]

	padj_fun <- function(j, i, x, y, width, height, fill) {

	  if(assoc_mat[i, j] < clinical_threshold)

	    grid.text(

	      # "\u066D", # asterisk but not very well centered

	      "*", # asterisk but not very well centered

	      x, y, gp = gpar(fontsize = 8))

	}


	col_fun_phen = colorRamp2(c(-2, 0, 2), c("purple", "white", "orange4"))
	#col_fun = colorRamp2(c(min(assoc_toplot), 0, max(assoc_toplot)), c("purple", "white", "darkgreen"))
	colnames(assoc_toplot) = gsub("_"," ", colnames(assoc_toplot))
	colnames(assoc_toplot) = gsub("Any","", colnames(assoc_toplot))


	hp= Heatmap(assoc_toplot,
			name = "signed -log10(q)",
			col = col_fun_phen,
			row_order = unlist(row_order(ht)), 
			column_order = colnames(assoc_toplot), 
			column_title = "Phenotypes", 
			column_names_gp = gpar(fontsize = fontsize),
			cell_fun = padj_fun,
			heatmap_legend_param = list(legend_direction = "horizontal"),
			#column_names_rot = 60 ),
			row_dend_reorder = dendrogram,
			border_gp = gpar(col = "black", lty = 1),
			width = ncol(assoc_toplot)*unit(4, "mm"),
			show_row_names = FALSE,
			show_heatmap_legend = legend_pos_auto
	)


	if(!is.null(enrichment)){
		comp = suppressWarnings(Reduce(function(x,y) merge(x,y, by = "pathway"), enrichment))
		comp.names = lapply(1:length(unique(clustering)), function(x)
								c(	paste0("p.value",x),
								 	paste0("num_hits",x),
									paste0("percent",x),
									paste0("cluster ",x)))
		comp.names = Reduce(function(x,y) c(x,y), comp.names)
		names(comp)[2:ncol(comp)] = comp.names

		comp = comp[, c("pathway", names(comp)[grep("cluster",names(comp))])]
		rownames(comp) = comp$pathway
		comp = comp[,-1]

		tokeep = lapply(1:ncol(comp), function(x) which(comp[,x] > 1))
		tokeep = unique(unlist(tokeep))
		comp = t(comp[tokeep,])

		col_fun = colorRamp2(c(-3, 0, 3), c("darkgreen", "white", "orange"))

		he= Heatmap(comp,
			name = "pathways/classes signed -log10(q.value)",
			col = col_fun,
			row_order = unlist(row_order(ht)), 
			column_order = colnames(comp), 
			column_title = "Pathways/Classes", 
			column_names_gp = gpar(fontsize = fontsize),
			#cell_fun = padj_fun,
			heatmap_legend_param = list(legend_direction = "horizontal"),
			#column_names_rot = 60 ),
			row_dend_reorder = dendrogram,
			border_gp = gpar(col = "black", lty = 1),
			show_heatmap_legend  =legend_pos_auto
		)
		h_list = ht+hp+he

	} else{h_list = ht+hp}

	#draw(h_list)

	if(legend_pos_auto){
		return(list(heatmap = h_list))
	}else{
		lgd = packLegend(
			Legend(title = "log2FC", col_fun = col_fun_main, direction = "horizontal" ),
			Legend(title = "signed -log10(q)", col_fun = col_fun_phen, direction = "horizontal"),
			direction = "horizontal"
		)
		return(list(heatmap = h_list, lgd = lgd))
	}


}


chordPlot <- function(data,cols, col_fun, title, diff_thresh, r_thresh) {

    tmp <- data
    ## tmp$flag <- ifelse(abs(tmp$diff) >= diff_thresh & abs(tmp$r) >= r_thresh, "yes","no")
    ## tmp <- tmp %>%
    ##     dplyr::select(-diff)
    circos.clear()
    chordDiagram(tmp,  grid.col = cols,  col = col_fun ,annotationTrack = "grid", preAllocateTracks = 1, scale = T, link.sort = T, link.visible = abs(tmp$r) > .2)
    title(title)
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
        circos.axis(h = "top", labels.cex = 0.2, major.tick.length = 0.2, sector.index = sector.name, track.index = 2)
    }, bg.border = NA)
    
    lgd <- Legend(at = c(-.6, 0, .6), col_fun = col_fun, title = "r2", direction = "horizontal")
    
    draw( lgd, just = "bottom", x = unit(6, "in"),y = unit(6, "mm"))

}


plotCorr_1plot <- function(tmp.data, x_var, y_var,   met_path) {

    path <- met_path %>%
        filter(name == y_var) %>%
        .$Pathway

    
    p1 =  ggplot(tmp.data, aes_string(tmp.data[,x_var], y = tmp.data[,y_var], color = "cluster")) +
        geom_point(size = 2)    +
        geom_smooth(method = "lm",  se = F) +
        scale_color_manual(values = c("red", "grey", "lightblue")) +       
        stat_cor(method = "spearman", label.x.npc = .8, label.y.npc = .2, size = 4 ) + 
        ggpubr::theme_pubr() +
        ggtitle(paste( x_var, "-", paste0(y_var, " (",path,")")))+
        theme(plot.title = element_text(hjust = .5))+
        xlim(c(min(tmp.data[,x_var]),max(tmp.data[,x_var]))) + 
        ylim(c(min(tmp.data[,y_var]),max(tmp.data[,y_var])))+
        xlab(x_var)+ ylab(y_var)

    return(p1)
    
}








