#### NEMO functions 

# these functions are for running the nemo algorithm and some functions for evaluating the number of clusters

#library(devtools)
#devtools::install_github('Shamir-Lab/NEMO/NEMO', force = T)
library(NEMO)
library(SNFtool)


normalize.matrix <- function(data.matrix) {
  temp = data.matrix - rowMeans(data.matrix)
  should.keep = (apply(temp, 1, sd) != 0)
  return ((temp / apply(temp, 1, sd))[should.keep, ])
}



#' @title NEMO num clusters
#' @name nemo.num.clusters
#' @description Estimates the number of clusters in an affinity graph.
#' @param W the affinity graph.
#' @param NUMC possible values for the number of clusters. Defaults to 2:15.
#' @return the estimated number of clusters in the graph.
#' @export
nemo.num.clusters <- function(W, NUMC=2:15) {
  if (min(NUMC) == 1) {
    warning("Note that we always assume there are more than one cluster.")
    NUMC = NUMC[NUMC > 1]
  }
  W = (W + t(W))/2
  diag(W) = 0
  if (length(NUMC) > 0) {
    degs = rowSums(W)
    degs[degs == 0] = .Machine$double.eps
    D = diag(degs)
    L = D - W
    Di = diag(1/sqrt(degs))
    L = Di %*% L %*% Di
    print(dim(L))
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return = T)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]
    eigengap = abs(diff(eigs$values))
    eigengap = (1:length(eigengap)) * eigengap

    t1 <- sort(eigengap[NUMC], decreasing = TRUE, index.return = T)$ix
    return(NUMC[t1[1]])
  }
}


#' @title NEMO num clusters
#' @name elbow.num.clusters
#' @description Estimates the number of clusters in an affinity graph.
#' @param W the affinity graph.
#' @param NUMC possible values for the number of clusters. Defaults to 2:15.
#' @return the estimated number of clusters in the graph.
#' @export
elbow.num.clusters <- function(W, NUMC=2:15) {
  if (min(NUMC) == 1) {
    warning("Note that we always assume there are more than one cluster.")
    NUMC = NUMC[NUMC > 1]
  }
  W = (W + t(W))/2
  diag(W) = 0
  if (length(NUMC) > 0) {
    degs = rowSums(W)
    degs[degs == 0] = .Machine$double.eps
    D = diag(degs)
    L = D - W
    Di = diag(1/sqrt(degs))
    L = Di %*% L %*% Di
#    print(dim(L))
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return = T)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]

    plot(NUMC, eigs$values[NUMC], xlab = "index", ylab = 'eigenvalues', type = "b")
    myplot = recordPlot()
    # find the elbow index

    # this is to help with testing
    x = NUMC
    y = eigs$values[NUMC]
    threshold =.01
    # runt the function
    elbow = get.elbow(NUMC, eigs$values[NUMC], NUMC )

    return(list(num.clusts = elbow, plot = myplot ))
  }
}


# function to identify the elbow from https://stackoverflow.com/questions/41518870/finding-the-elbow-knee-in-a-curve 
get.elbow <- function(x, y, NUMC) {
  d1 <- diff(y) / diff(x)
  d2 = diff(d1) / diff(x[-1]) # second derivative
  indices <- sort(abs(d2), decreasing = T, index.return = T)$ix  + min(NUMC)


  return(indices)
}


#' @title NEMO affinity graph
#' @name nemo.affinity.graph
#' @description Constructs a single affinity graph measuring similarity across different omics.
#' @param raw.data A list of the data to be clustered, where each an entry is a matrix of features x samples.
#' @param k The number of neighbors to use for each omic. It can either be a number, a list of numbers
#' or NA. If it is a number, this is the number of neighbors used for all omics. If this is a list,
#' the number of neighbors are taken for each omic from that list. If it is NA, each omic chooses the
#' number of neighbors to be the number of samples divided by NUM.NEIGHBORS.RATIO.
#' @return A single matrix measuring similarity between the samples across all omics.
#' @export

nemo.affinity.graph <- function(raw.data, k=NA, num.clusters) {
  if (is.na(k)) { k =ncol(raw.data[[1]])/num.clusters}

  sim.data = lapply(1:length(raw.data), function(i) {affinityMatrix(SNFtool::dist2(as.matrix(t(raw.data[[i]])),
                                                                as.matrix(t(raw.data[[i]]))), k, 0.5)})
  affinity.per.omic = lapply(1:length(raw.data), function(i) {
      sim.datum = sim.data[[i]]
      non.sym.knn = apply(sim.datum, 1, function(sim.row) {
      returned.row = sim.row
      threshold = sort(sim.row, decreasing = T)[k]
      returned.row[sim.row < threshold] = 0
      row.sum = sum(returned.row)
      returned.row[sim.row >= threshold] = returned.row[sim.row >= threshold] / row.sum
      return(returned.row)
    })

    sym.knn = non.sym.knn + t(non.sym.knn)
    return(sym.knn)
  })
  patient.names = Reduce(union, lapply(raw.data, colnames))
  num.patients = length(patient.names)
  returned.affinity.matrix = matrix(0, ncol = num.patients, nrow=num.patients)
  rownames(returned.affinity.matrix) = patient.names
  colnames(returned.affinity.matrix) = patient.names

  shared.omic.count = matrix(0, ncol = num.patients, nrow=num.patients)
  rownames(shared.omic.count) = patient.names
  colnames(shared.omic.count) = patient.names

  for (j in 1:length(raw.data)) {
    curr.omic.patients = colnames(raw.data[[j]])
    returned.affinity.matrix[curr.omic.patients, curr.omic.patients] = returned.affinity.matrix[curr.omic.patients, curr.omic.patients] + affinity.per.omic[[j]][curr.omic.patients, curr.omic.patients]
    shared.omic.count[curr.omic.patients, curr.omic.patients] = shared.omic.count[curr.omic.patients, curr.omic.patients] + 1
  }

  final.ret = returned.affinity.matrix / shared.omic.count
  lower.tri.ret = final.ret[lower.tri(final.ret)]
  final.ret[shared.omic.count == 0] = mean(lower.tri.ret[!is.na(lower.tri.ret)])

  return(final.ret)
}



#' @title NEMO clustering
#' @name nemo.clustering
#' @description Performs multi-omic clustering on a datset using the NEMO algorithm.
#' Uses nemo.num.clusters to estimate the number of clusters.
#' @param omics.list A list of the data to be clustered, where each an entry is a matrix of features x samples.
#' @param k The number of neighbors to use for each omic. It can either be a number, a list of numbers
#' or NA. If it is a number, this is the number of neighbors used for all omics. If this is a list,
#' the number of neighbors are taken for each omic from that list. If it is NA, each omic chooses the
#' number of neighbors to be the number of samples divided by NUM.NEIGHBORS.RATIO.
#' @return A single matrix measuring similarity between the samples across all omics.
#' @export
nemo.clustering <- function(omics.list, num.clusters=NA, num.neighbors=NA) {
  graph = nemo.affinity.graph(omics.list, k = num.neighbors, num.clusters)

  clustering = spectralClustering(graph,num.clusters) 
	names(clustering) = colnames(graph)

  return(list(clustering = clustering, graph = graph,  num.clusters = length(unique(clustering))))
}


#### Robustness is a function for determining the robustness of the clusters. 
### This method employs bootstrapping (sampling with replacement)

## should I impose a fixed number of clusters each time resampled? I don't think this makes sense if I am resampling with replacement, because it would be expected that the same individual would always cluster in the same cluster, thus decreasing the number of clusters. 

# not sure if the above statement is true

RobustClust = function(omics.list, iterations = 100, num.clusters = NA){
	
	# find the minimum cluster size
	min.cluster.size = floor(.05*ncol(omics.list[[1]]))
	if(is.na(num.clusters)){

		# find clustering with full data
		full_clusters = nemo.clustering(omics.list, min.cluster.size = min.cluster.size) 
		full_clustering = full_clusters$clustering

	} else{
		full_clusters = nemo.clustering(omics.list, num.clusters = num.clusters, min.cluster.size = min.cluster.size) 
		full_clustering = full_clusters$clustering

	}

	if(min(table(full_clustering)) <= min.cluster.size){
		print("Smallest cluster < Minimum cluster size ")
		return(list(scores = NA, cluster_means = NA, results = NA))
	} else {

	# resample the subjects
		subjects = colnames(omics.list[[1]])
		num_subjects = ncol(omics.list[[1]])

		scores = list()

		for(i in 1:iterations){

			lessThanMin = TRUE 
			
			while(lessThanMin){
				resampled_subjects = sample(subjects,num_subjects, replace = T) 
				resampled.list = lapply(omics.list, function(x) x[,resampled_subjects])
				resampled_clusters = nemo.clustering(resampled.list, min.cluster.size = min.cluster.size, NUMC = 2:10)
				resampled_clustering = resampled_clusters$clustering
				if(min(table(resampled_clustering)) > min.cluster.size){lessThanMin = FALSE }
			}

			unique_resampled = unique(resampled_subjects)
			num_unique = length(unique_resampled)

			tmp.scores = numeric()

			for(j  in 1: max(full_clustering)){
				tmp.clust = full_clustering[names(full_clustering[full_clustering ==j])]
				tmp.clust = tmp.clust[names(tmp.clust) %in% unique_resampled]
				tmp.resampled = resampled_clustering[names(resampled_clustering) %in% names(tmp.clust)]
				counts = table(tmp.resampled)
				tmp.scores = c(tmp.scores,max(counts)/sum(counts))
			}

			scores[[i]]= tmp.scores
		}


		scores = data.frame(do.call(rbind,scores))
		names(scores) = paste0("cluster",1:max(full_clustering))

		cluster_means = colMeans(scores)



		return(list(scores = scores, cluster_means = cluster_means, results = full_clusters))

	}
}


## subClust is a function that defines subsets of the data and clusters them. 

# The subsampling will have a parameter that defines what pecentage of the data to subsample. There is no replacement, so each individiual will only be sampled once. Also, the number of clusters will be fixed to the same number as the full clustered data. This 

subClust = function(omics.list, iterations = 10, sampling.perc = .8, num.clusters = NA) {
    
                                        # find the minimum cluster size
    min.cluster.size = floor(.05*ncol(omics.list[[1]]))
                                        #min.cluster.size = 0

    full_clusters = nemo.clustering(omics.list, num.clusters = num.clusters) 
    full_clustering = full_clusters$clustering

    if(min(table(full_clustering)) <= min.cluster.size){
        print("Smallest cluster < Minimum cluster size ")
        return(list(scores = NA, cluster_means = NA, results = NA))
    } else {

                                        # resample the subjects
        subjects = colnames(omics.list[[1]])
        num_subjects = ncol(omics.list[[1]])

        scores = list()

        for(i in 1:iterations){

            lessThanMin = TRUE 
            
            while(lessThanMin) {
                resampled_subjects = sample(subjects,floor(num_subjects * sampling.perc) , replace = F) 
                resampled.list = lapply(omics.list, function(x) x[,resampled_subjects])
                                        #resampled.min.cluster.size = floor(.05*ncol(resampled.list[[1]]))
                resampled_clusters = nemo.clustering(omics.list = resampled.list,  num.clusters =full_clusters$num.clusters )
                resampled_clustering = resampled_clusters$clustering
                if(min(table(resampled_clustering)) > min.cluster.size) {
                    lessThanMin = FALSE }
            }

                                        #unique_resampled = unique(resampled_subjects)
                                        #num_unique = length(unique_resampled)

            tmp.scores = numeric()
            intersects <- list()
            resampled_length <- length(resampled_clustering)
            for(j  in 1: max(full_clustering)){
                full_cluster_names <- names(full_clustering)[full_clustering ==j]
               
                    jaccard <- sapply(1: max(full_clustering), function(x) {
                        resampled_cluster_names <- names(resampled_clustering)[resampled_clustering == x ]
                        length(intersect(full_cluster_names,resampled_cluster_names))/length(union(full_cluster_names,resampled_cluster_names))
                        
                    })
                print(jaccard)
                }

                tmp.scores[[j]] <- max(jaccard)


            }
            

            scores[[i]]= tmp.scores
        }


        scores = data.frame(do.call(rbind,scores))
        names(scores) = paste0("cluster",1:max(full_clustering))

        cluster_means = colMeans(scores)



        return(list(scores = scores, cluster_means = cluster_means, results = full_clusters))

    }
}
