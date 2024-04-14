### Reading in the data

# * Cytokine and targeted metabolomic data was obtained from Joaquin Espinosa and the [human trisomy project](http://www.trisome.org/) working group at the University of Colorado Anschutz Medical campus. 
# * Here is some information of the data processing and cleaning:
# 	* MSD cytokine profiling data. Plasma concentration values (pg/mL) for each of the cytokines and related immune factors measured across multiple MSD assay plates was imported to R, combined, and analytes with >10% of values outside of detection or fit curve range flagged. For each analyte, missing values were replaced with either the minimum (if below fit curve range) or maximum (if above fit curve range) calculated concentration and means of duplicate wells used in all further analysis. Use "Adjusted_Concentration", not really adjusted but we were trying to have similar names.
# 	* MS-metabolite data. Relative abundance values. 0/missing values were replaced with a random value sampled from between 0 and 0.5x the minimum non-zero intensity value for that metabolite. Data was then normalized using a scaling factor derived from the global median intensity value across all metabolites / sample median intensity across all proteins. Use "adjusted_relative_abundance".

library(ggplot2)
library(openxlsx)
library(limma)
library(data.table)

 options(stringsAsFactors = FALSE)



### read_HTP is a function for reading in the HTP data

read_HTP <- function(cytokine_file, metabolite_file, metadata_file, commorbidity_file){

	cyt  = read.delim(cytokine_file)
	met  = read.delim(metabolite_file)
 	meta = read.delim(metadata_file)
 	com  = read.delim(comorbidity_file, skip = 1)

	##Create matrices from molecular data for downstream analysis.

	## Cytokines

    
	cyt <- suppressWarnings(dcast(cyt[, c("LabID", "Analyte", "Adjusted_Concentration")], Analyte ~ LabID ))
	rownames(cyt) <- cyt$Analyte
	cyt <- cyt[,-1]
	cyt <- as.data.frame(t(cyt))


	## Metabolites
	met <- suppressWarnings(dcast(met[, c("LabID", "Analyte", "adjusted_relative_abundance")], Analyte ~ LabID ))
	rownames(met) <- met$Analyte
	met <- met[,-1]
	met <- as.data.frame(t(met))

	# matching metadata
	# including metadata from either cohort
	meta <- meta[meta$LabID %in% c(rownames(cyt), rownames(met)), ]


	# comorbidity meta data
	comorbidity_meta <- com[!duplicated(com$Condition), c("Condition", "Age.group", "min_Age", "max_Age")]

	# merge meta and comorbidity data to get clinic data
	com <- suppressWarnings(dcast(com[, c("RecordID", "Condition", "HasCondition")], RecordID ~ Condition))
	clinic <- merge(meta,com, by = "RecordID", all.x = TRUE) 


	# assign LABID's as rownames
	rownames(clinic) = clinic$LabID

	# # merge BMI with clinic 
	# clinic$Age.bin <- ordered(ifelse(clinic$Age >= 40, ">40",
	# 			  ifelse(clinic$Age >= 30, "30-39",	
	# 				ifelse(clinic$Age >= 20, "20-29",
	# 				ifelse(clinic$Age >= 10, "10-19", "<10")))), 
	# 				levels = c("<10", "10-19", "20-29", "30-39", ">40"))

	return(list(cytokines = cyt, metabolites = met, clinic = clinic, comorbidity_meta = comorbidity_meta))


}

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}
#omic_data = as.matrix(nona)

# library('hydroGOF')
imputationTest = function(omic_data){

	missing = seq(.01,.2,.01)
	rmse = numeric()
	n = 1
	for(z in missing){
		tmp.omic_data = omic_data
		tmp.omic_data[sample(nrow(tmp.omic_data)*ncol(tmp.omic_data), floor((nrow(tmp.omic_data)*ncol(tmp.omic_data))*z))] = NA

		impute.omic_data = impute.knn(tmp.omic_data, k = 10 )

		rmse[n] = mean(sapply(1:ncol(omic_data), function(x) RMSE(impute.omic_data$data[,x], omic_data[,x] )))
		n = n+1
}

	toplot = cbind(missing, rmse)
	toplot = as.data.frame(toplot)
	return(toplot)
}





### get_overlap is a function for reducing metabolite, cytokine, and clinical data to overlapping samples.
get_overlap <- function(cytokine_data, metabolite_data, clinic_data){

	cyt = cytokine_data
	met = metabolite_data
	clinic = clinic_data

	overlap = intersect(rownames(cyt), rownames(met))

	cyt = cyt[overlap,]
	met = met[overlap,]
	clinic = clinic[overlap,]


	return(list(cytokines = cyt, metabolites = met,  clinic = clinic))
}


#### remove_subj
remove_subj = function(data_list, toRemove){
		subj = which(rownames(data_list[[1]]) %in% toRemove)
		data.list = lapply(data_list, function(x) x[-subj,])
		return(data.list)
}

## adjust_batch is a function for adjusting datasets based on sample source
### Not working currently
adjust_vars <- function(omics.list, clinic, covars){
	
	covar.mat = eval(parse(text(sprintf("model.matrix(~%s)", paste0("clinic$", covars)))))
	tmp.omics.list = lapply(omics.list, function(x) removeBatchEffect(x, batch = factor(clinic$Sample_source)))
	return(tmp.omics.list)

}

reduce_features = function(omics.list, features){
	omics.list= lapply(omics.list, function(x) x[rownames(x) %in% features,])
	return(omics.list)
}
### downsample_T21 is a function for downsampling the T21 subjects to be more similar in counts to D21

downsample_T21 = function(data_list, T21_D21_ratio = 1 ){

	met = data_list$metabolites
	cyt = data_list$cytokines
	clinic = data_list$clinic

	keep = c(which(clinic$Karyotype == "Control"), 
			   sample(which(clinic$Karyotype == "T21"), floor(length(which(clinic$Karyotype == "Control")) * T21_D21_ratio))
			)
	met = met[keep,]
	cyt = cyt[keep,]

	return(list(cytokines = cyt, metabolites = met))
}


only_T21 = function(data_list){
	data_list = data_list[[1]]
	cyt = as.data.frame(t(data_list[[1]]))
	met = as.data.frame(t(data_list[[2]]))
	clinic = data_list[[3]]

	keep = c(which(clinic$Karyotype == "T21"))

	met = met[keep,]
	cyt = cyt[keep,]

	return(list(cytokines = cyt, metabolites = met, clinic = clinic[clinic$Karyotype == "T21",]))		
}


only_D21 = function(data_list){
	data_list = data_list[[1]]
	cyt = as.data.frame(t(data_list[[1]]))
	met = as.data.frame(t(data_list[[2]]))
	clinic = data_list[[3]]

	keep = c(which(clinic$Karyotype == "Control"))

	met = met[keep,]
	cyt = cyt[keep,]

	return(list(cytokines = cyt, metabolites = met, clinic = clinic[clinic$Karyotype == "Control",]))		
}


replace_outlier <- function(x, coef = 3){
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  .IQR <- IQR(x, na.rm = TRUE)
  upper.limit <- Q3 + (coef*.IQR)
  lower.limit <- Q1 - (coef*.IQR)
  low.res <- ifelse(x < lower.limit, NA, x )
  
  count = sum(is.na(low.res))

  min  = min(low.res, na.rm = T)
  min_75 = ifelse(min < 0, 1.25*min, .75*min)
  
  x <- ifelse(is.na(low.res), sample(runif(100, min = min_75, max =  min),1), x)

  high.res <- ifelse(x > upper.limit ,NA,x)

  count = count + sum(is.na(high.res))

  max  = max(high.res, na.rm = T)
  max_125 = ifelse(max < 0, .75*max, 1.25*max)


  x <- ifelse(is.na(high.res), sample(runif(100, min = max, max =max_125), 1), x)
  return(list(x = x, count = count))

}


diff_expr_single = function(x,y, test ){
	if(test == "t.test"){
		test = t.test(x,y)
	} else if(test == "wilcoxon"){
		test = wilcox.test(x,y)
	}
	p = test$p.value
	stat = test$statistic
        FC = mean(x, na.rm = T) - mean(y, na.rm = T)
	res = c(p, FC)
	return(res)
}







diff_expr_wrapper = function(data, clinic, phenotype, control, q_cutoff = 0.05, test ){
    data <- as.data.frame(t(data))
	vals = unique(clinic[, phenotype])
	which_x = which(clinic[, phenotype] == vals[vals != control])
	which_y = which(clinic[, phenotype] == vals[vals == control])

	x_dat = data[which_x,]
    y_dat = data[which_y,]

       res = sapply(1:ncol(data), function(i){
				diff_expr_single(x =as.numeric(x_dat[,i]), as.numeric(y_dat[,i]), test)
	})
	res = as.data.frame(t(res))
	rownames(res) = names(data)
	colnames(res) = c("p.value", "FC")
	res$q.value = p.adjust(res$p.value, method = "fdr")
    res <- res[, c("FC", "p.value", "q.value")]
	return(res)
}


# SPlit Violin Plots
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    # Original function by Jan Gleixner (@jan-glx)
    # Adjustments by Wouter van der Bijl (@Axeman)
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
