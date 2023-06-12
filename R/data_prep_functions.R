### Reading in the data

# * Cytokine and targeted metabolomic data was obtained from Joaquin Espinosa and the [human trisomy project](http://www.trisome.org/) working group at the University of Colorado Anschutz Medical campus. 
# * Here is some information of the data processing and cleaning:
# 	* MSD cytokine profiling data. Plasma concentration values (pg/mL) for each of the cytokines and related immune factors measured across multiple MSD assay plates was imported to R, combined, and analytes with >10% of values outside of detection or fit curve range flagged. For each analyte, missing values were replaced with either the minimum (if below fit curve range) or maximum (if above fit curve range) calculated concentration and means of duplicate wells used in all further analysis. Use "Adjusted_Concentration", not really adjusted but we were trying to have similar names.
# 	* MS-metabolite data. Relative abundance values. 0/missing values were replaced with a random value sampled from between 0 and 0.5x the minimum non-zero intensity value for that metabolite. Data was then normalized using a scaling factor derived from the global median intensity value across all metabolites / sample median intensity across all proteins. Use "adjusted_relative_abundance".



library(data.table)
library(reshape2)
library(ggplot2)
# library(sva)
library(openxlsx)
library(limma)
library(impute)
# library(rstatix)
# library(lme4)

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

	### Removing phenotypes with low case counts(Cataracts, Autism, Regression, Pulmonary Hyptertension)
	clinic = clinic[, !(names(clinic) %in%c("Autism spectrum disorder", "Cataracts", "Regression", "Pulmonary hypertension"))]

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
	#parameter = test$parameter
	FC = mean(x)- mean(y)
	res = c(p,stat, FC)
	return(res)
}


diff_expr_wrapper = function(data, clinic, phenotype, control, q_cutoff = 0.05, test ){

	vals = unique(clinic[, phenotype])
	which_x = which(clinic[, phenotype] == vals[vals != control])
	which_y = which(clinic[, phenotype] == vals[vals == control])

	x_dat = data[,which_x]
	y_dat = data[, which_y]

	res = sapply(1:nrow(data), function(i){
				diff_expr_single(as.numeric(x_dat[i,]), as.numeric(y_dat[i,]), test)
	})
	res = as.data.frame(t(res))
	rownames(res) = rownames(data)
	colnames(res) = c("p.value", "statistic", "logFC")
	res$q.value = p.adjust(res$p.value, method = "fdr")
	res$color = ifelse(res$q.value < q_cutoff & res$logFC > 0, "navy", 
				ifelse(res$q.value < q_cutoff & res$logFC < 0, "blue", "gray" ))
	res$label = ifelse(res$q.value < q_cutoff, rownames(res), "")
	return(res)
}


# pathway_annotator is a function for properly annotating features based on their kegg ID's

# ##Create matrices from molecular data for downstream analysis.

# ## Cytokines

# cyt <- dcast(cyt[, c("LabID", "Analyte", "Adjusted_Concentration")], Analyte ~ LabID )
# rownames(cyt) <- cyt$Analyte
# cyt <- cyt[,-1]
# cyt <- as.data.frame(t(cyt))




# ## Metabolites

# met <- dcast(met[, c("LabID", "Analyte", "adjusted_relative_abundance")], Analyte ~ LabID )
# rownames(met) <- met$Analyte
# met <- met[,-1]
# met <- as.data.frame(t(met))

# # Find which samples have cytokine and metabolomic data
# ids = intersect(rownames(cyt), rownames(met))


# cyt <- cyt[ids,]
# met <- met[ids,]
# clinic <- clinic[clinic$LabID %in% rownames(met),]
# clinic <- clinic[order(match(clinic$LabID, rownames(met))),]

# save(cyt, met, clinic, file = sprintf('%s/data/raw_cyt_met.RData', DIR))

# ## Density plots

#pdf('~/CPBS/HTP/results/cytokine_density_plots.pdf',
 #   height = 5, width = 5)
#	for(i in 1:54){
#		plot(density(cyt[,i]),
#		main = sprintf("Density of %s", names(cyt)[i])) 
#	}
#dev.off()

#pdf('~/CPBS/HTP/results/cytokine_adj_density_plots.pdf',
#    height = 5, width = 5)
#	for(i in 1:54){
#		plot(density(log.norm.batch[[2]][i,]),
#		main = sprintf("Density of Adjusted %s", rownames(T21.log.norm.batch[[2]])[i])) 
#	}
#dev.off()


# pdf('~/labNotebook/content/post/Down-Syndrome-Omic-Integration/metabolite_density_plots.pdf',
#     height = 5, width = 5)
# 	for(i in 1:174){
# 		plot(density(met[,i]),
# 		main = sprintf("Density of %s", names(met)[i])) 
# 	}
# dev.off()


# # log2  transformation
# cyt <- log2(cyt)
# # cyt <- merge(cyt, clinic[, c("LabID","Sex","Karyotype", "Age" )], by.x = 0, by.y = "LabID")
# # rownames(cyt) <- cyt$Row.names
# # cyt <- cyt[,-1]


# met <- log2(met)
# # met <- merge(met,  clinic[, c("LabID","Sex","Karyotype", "Age" )], by.x = 0, by.y = "LabID")
# # rownames(met) <- met$Row.names
# # met <- met[,-1]


# save(cyt,met, clinic,file = sprintf("%sdata/log2_cytokine_metabolite_Data.RData", DIR))

# write.csv(cyt, file = sprintf("%s/data/log2_cytokine.csv", DIR))
# write.csv(met, file = sprintf("%s/data/log2_metabolite.csv", DIR))

# ## Density Plots

# pdf('~/labNotebook/content/post/Down-Syndrome-Omic-Integration/cytokine_log2_density_plots.pdf',
#     height = 5, width = 5)
# 	for(i in 1:54){
# 		plot(density(cyt[,i]),
# 		main = sprintf("Density of log2 %s", names(cyt)[i])) 

# 		print(ggplot(cyt, aes(x = cyt[, i], fill = Karyotype, alpha = .3)) + 
# 			geom_density() + 
# 			xlab(sprintf(names(cyt)[i])) + 
# 			ggtitle(sprintf("Density of log2 %s by Karyotype", names(cyt)[i]))+
# 			theme_classic() + 
# 			theme(legend.position = c(.85,.85)) + 
# 			guides(alpha = FALSE)
# 		)
# 	}
# dev.off()


# # Metabolites
# pdf('~/labNotebook/content/post/Down-Syndrome-Omic-Integration/metabolites_log2_density_plots.pdf',
#     height = 5, width = 5)
# 	for(i in 1:174){
# 		plot(density(met[,i]),
# 		main = sprintf("Density of log2 %s", names(met)[i])) 

# 		print(ggplot(met, aes(x = met[, i], fill = Karyotype, alpha = .3)) + 
# 			geom_density() + 
# 			xlab(sprintf(names(met)[i])) + 
# 			ggtitle(sprintf("Density of log2 %s by Karyotype", names(met)[i]))+
# 			theme_classic() + 
# 			theme(legend.position = c(.85,.85)) + 
# 			guides(alpha = FALSE)
# 		)
# 	}
# dev.off()


# ## z-score


# met[,1:174] <- apply(met[,1:174],2,scale)
# cyt[,1:54] <- apply(cyt[,1:54],2,scale)


# save(cyt,met, clinic,file = sprintf("%sdata/scale_cytokine_metabolite_Data.RData", DIR))


# write.csv(cyt, file = sprintf("%s/data/scale_cytokine.csv", DIR))
# write.csv(met, file = sprintf("%s/data/scale_metabolite.csv", DIR))

# ##### ComBat Normalization
# load(file = sprintf('%s/data/raw_cyt_met.RData', DIR))


# met.corr.model <- model.matrix(~1, data = clinic)
# met.age.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(met)), 
# 					  batch = clinic$Age, 
# 					  mod = met.corr.model)))

# cyt.corr.model <- model.matrix(~1, data = clinic)
# cyt.age.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(cyt)), 
# 					  batch = clinic$Age, 
# 					  mod = cyt.corr.model)))

# save(cyt.age.adj,met.age.adj, clinic,file = sprintf("%sdata/age_adj_cytokine_metabolite_Data.RData", DIR))



# met.corr.model <- model.matrix(~1, data = clinic)
# met.sex.model <- model.matrix(~0+Sex, data = clinic)

# met.sex.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(met)), 
# 					  batch = met.sex.model[,1], 
# 					  mod = met.corr.model)))




# cyt.corr.model <- model.matrix(~1, data = clinic)
# cyt.sex.model <- model.matrix(~0+Sex, data = clinic)

# cyt.sex.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(cyt)), 
# 					  batch = cyt.sex.model[,1], 
# 					  mod = cyt.corr.model)))


# save(cyt.sex.adj,met.sex.adj, clinic,file = sprintf("%sdata/sex_adj_cytokine_metabolite_Data.RData", DIR))



# met.sex.age.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(met.sex.adj)), 
# 					  batch = clinic$Age, 
# 					  mod = met.corr.model)))


# cyt.sex.age.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(cyt.sex.adj)), 
# 					  batch = clinic$Age, 
# 					  mod = cyt.corr.model)))

# save(cyt.sex.age.adj,met.sex.age.adj, clinic,file = sprintf("%sdata/sex_age_adj_cytokine_metabolite_Data.RData", DIR))



# ### Separate karyotypes
# load(file = sprintf('%s/data/raw_cyt_met.RData', DIR))

# met.T21 <- met[rownames(met) %in% clinic$LabID[clinic$Karyotype == "T21"],]
# met.D21 <- met[rownames(met) %in% clinic$LabID[clinic$Karyotype == "Control"],]

# cyt.T21 <- cyt[rownames(cyt) %in% clinic$LabID[clinic$Karyotype == "T21"],]
# cyt.D21 <- cyt[rownames(cyt) %in% clinic$LabID[clinic$Karyotype == "Control"],]


# save(met.T21, met.D21, cyt.T21, cyt.D21, clinic,file = sprintf("%sdata/Karyotype_cytokine_metabolite_Data.RData", DIR))

### Reading in the data

# * Cytokine and targeted metabolomic data was obtained from Joaquin Espinosa and the [human trisomy project](http://www.trisome.org/) working group at the University of Colorado Anschutz Medical campus. 
# * Here is some information of the data processing and cleaning:
# 	* MSD cytokine profiling data. Plasma concentration values (pg/mL) for each of the cytokines and related immune factors measured across multiple MSD assay plates was imported to R, combined, and analytes with >10% of values outside of detection or fit curve range flagged. For each analyte, missing values were replaced with either the minimum (if below fit curve range) or maximum (if above fit curve range) calculated concentration and means of duplicate wells used in all further analysis. Use "Adjusted_Concentration", not really adjusted but we were trying to have similar names.
# 	* MS-metabolite data. Relative abundance values. 0/missing values were replaced with a random value sampled from between 0 and 0.5x the minimum non-zero intensity value for that metabolite. Data was then normalized using a scaling factor derived from the global median intensity value across all metabolites / sample median intensity across all proteins. Use "adjusted_relative_abundance".



library(data.table)
library(reshape2)
library(ggplot2)
# library(sva)
library(openxlsx)
library(limma)
library(impute)
# library(rstatix)
# library(lme4)

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

	### Removing phenotypes with low case counts(Cataracts, Autism, Regression, Pulmonary Hyptertension)
	clinic = clinic[, !(names(clinic) %in%c("Autism spectrum disorder", "Cataracts", "Regression", "Pulmonary hypertension"))]

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
	FC = mean(x)- mean(y)
	res = c(p,FC)
	return(res)
}


diff_expr_wrapper = function(data, clinic, phenotype, control, q_cutoff = 0.05, test ){

	vals = unique(clinic[, phenotype])
	which_x = which(clinic[, phenotype] == vals[vals != control])
	which_y = which(clinic[, phenotype] == vals[vals == control])

	x_dat = data[,which_x]
	y_dat = data[, which_y]

	res = sapply(1:nrow(data), function(i){
				diff_expr_single(as.numeric(x_dat[i,]), as.numeric(y_dat[i,]), test)
	})
	res = as.data.frame(t(res))
	rownames(res) = rownames(data)
	colnames(res) = c("p.value", "logFC")
	res$q.value = p.adjust(res$p.value, method = "fdr")
	res$color = ifelse(res$q.value < q_cutoff & res$logFC > 0, "navy", 
				ifelse(res$q.value < q_cutoff & res$logFC < 0, "blue", "gray" ))
	res$label = ifelse(res$q.value < q_cutoff, rownames(res), "")
	return(res)
}


# pathway_annotator is a function for properly annotating features based on their kegg ID's

# ##Create matrices from molecular data for downstream analysis.

# ## Cytokines

# cyt <- dcast(cyt[, c("LabID", "Analyte", "Adjusted_Concentration")], Analyte ~ LabID )
# rownames(cyt) <- cyt$Analyte
# cyt <- cyt[,-1]
# cyt <- as.data.frame(t(cyt))




# ## Metabolites

# met <- dcast(met[, c("LabID", "Analyte", "adjusted_relative_abundance")], Analyte ~ LabID )
# rownames(met) <- met$Analyte
# met <- met[,-1]
# met <- as.data.frame(t(met))

# # Find which samples have cytokine and metabolomic data
# ids = intersect(rownames(cyt), rownames(met))


# cyt <- cyt[ids,]
# met <- met[ids,]
# clinic <- clinic[clinic$LabID %in% rownames(met),]
# clinic <- clinic[order(match(clinic$LabID, rownames(met))),]

# save(cyt, met, clinic, file = sprintf('%s/data/raw_cyt_met.RData', DIR))

# ## Density plots

#pdf('~/CPBS/HTP/results/cytokine_density_plots.pdf',
 #   height = 5, width = 5)
#	for(i in 1:54){
#		plot(density(cyt[,i]),
#		main = sprintf("Density of %s", names(cyt)[i])) 
#	}
#dev.off()

#pdf('~/CPBS/HTP/results/cytokine_adj_density_plots.pdf',
#    height = 5, width = 5)
#	for(i in 1:54){
#		plot(density(log.norm.batch[[2]][i,]),
#		main = sprintf("Density of Adjusted %s", rownames(T21.log.norm.batch[[2]])[i])) 
#	}
#dev.off()


# pdf('~/labNotebook/content/post/Down-Syndrome-Omic-Integration/metabolite_density_plots.pdf',
#     height = 5, width = 5)
# 	for(i in 1:174){
# 		plot(density(met[,i]),
# 		main = sprintf("Density of %s", names(met)[i])) 
# 	}
# dev.off()


# # log2  transformation
# cyt <- log2(cyt)
# # cyt <- merge(cyt, clinic[, c("LabID","Sex","Karyotype", "Age" )], by.x = 0, by.y = "LabID")
# # rownames(cyt) <- cyt$Row.names
# # cyt <- cyt[,-1]


# met <- log2(met)
# # met <- merge(met,  clinic[, c("LabID","Sex","Karyotype", "Age" )], by.x = 0, by.y = "LabID")
# # rownames(met) <- met$Row.names
# # met <- met[,-1]


# save(cyt,met, clinic,file = sprintf("%sdata/log2_cytokine_metabolite_Data.RData", DIR))

# write.csv(cyt, file = sprintf("%s/data/log2_cytokine.csv", DIR))
# write.csv(met, file = sprintf("%s/data/log2_metabolite.csv", DIR))

# ## Density Plots

# pdf('~/labNotebook/content/post/Down-Syndrome-Omic-Integration/cytokine_log2_density_plots.pdf',
#     height = 5, width = 5)
# 	for(i in 1:54){
# 		plot(density(cyt[,i]),
# 		main = sprintf("Density of log2 %s", names(cyt)[i])) 

# 		print(ggplot(cyt, aes(x = cyt[, i], fill = Karyotype, alpha = .3)) + 
# 			geom_density() + 
# 			xlab(sprintf(names(cyt)[i])) + 
# 			ggtitle(sprintf("Density of log2 %s by Karyotype", names(cyt)[i]))+
# 			theme_classic() + 
# 			theme(legend.position = c(.85,.85)) + 
# 			guides(alpha = FALSE)
# 		)
# 	}
# dev.off()


# # Metabolites
# pdf('~/labNotebook/content/post/Down-Syndrome-Omic-Integration/metabolites_log2_density_plots.pdf',
#     height = 5, width = 5)
# 	for(i in 1:174){
# 		plot(density(met[,i]),
# 		main = sprintf("Density of log2 %s", names(met)[i])) 

# 		print(ggplot(met, aes(x = met[, i], fill = Karyotype, alpha = .3)) + 
# 			geom_density() + 
# 			xlab(sprintf(names(met)[i])) + 
# 			ggtitle(sprintf("Density of log2 %s by Karyotype", names(met)[i]))+
# 			theme_classic() + 
# 			theme(legend.position = c(.85,.85)) + 
# 			guides(alpha = FALSE)
# 		)
# 	}
# dev.off()


# ## z-score


# met[,1:174] <- apply(met[,1:174],2,scale)
# cyt[,1:54] <- apply(cyt[,1:54],2,scale)


# save(cyt,met, clinic,file = sprintf("%sdata/scale_cytokine_metabolite_Data.RData", DIR))


# write.csv(cyt, file = sprintf("%s/data/scale_cytokine.csv", DIR))
# write.csv(met, file = sprintf("%s/data/scale_metabolite.csv", DIR))

# ##### ComBat Normalization
# load(file = sprintf('%s/data/raw_cyt_met.RData', DIR))


# met.corr.model <- model.matrix(~1, data = clinic)
# met.age.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(met)), 
# 					  batch = clinic$Age, 
# 					  mod = met.corr.model)))

# cyt.corr.model <- model.matrix(~1, data = clinic)
# cyt.age.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(cyt)), 
# 					  batch = clinic$Age, 
# 					  mod = cyt.corr.model)))

# save(cyt.age.adj,met.age.adj, clinic,file = sprintf("%sdata/age_adj_cytokine_metabolite_Data.RData", DIR))



# met.corr.model <- model.matrix(~1, data = clinic)
# met.sex.model <- model.matrix(~0+Sex, data = clinic)

# met.sex.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(met)), 
# 					  batch = met.sex.model[,1], 
# 					  mod = met.corr.model)))




# cyt.corr.model <- model.matrix(~1, data = clinic)
# cyt.sex.model <- model.matrix(~0+Sex, data = clinic)

# cyt.sex.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(cyt)), 
# 					  batch = cyt.sex.model[,1], 
# 					  mod = cyt.corr.model)))


# save(cyt.sex.adj,met.sex.adj, clinic,file = sprintf("%sdata/sex_adj_cytokine_metabolite_Data.RData", DIR))



# met.sex.age.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(met.sex.adj)), 
# 					  batch = clinic$Age, 
# 					  mod = met.corr.model)))


# cyt.sex.age.adj <- as.data.frame(t(ComBat(dat = as.matrix(t(cyt.sex.adj)), 
# 					  batch = clinic$Age, 
# 					  mod = cyt.corr.model)))

# save(cyt.sex.age.adj,met.sex.age.adj, clinic,file = sprintf("%sdata/sex_age_adj_cytokine_metabolite_Data.RData", DIR))



# ### Separate karyotypes
# load(file = sprintf('%s/data/raw_cyt_met.RData', DIR))

# met.T21 <- met[rownames(met) %in% clinic$LabID[clinic$Karyotype == "T21"],]
# met.D21 <- met[rownames(met) %in% clinic$LabID[clinic$Karyotype == "Control"],]

# cyt.T21 <- cyt[rownames(cyt) %in% clinic$LabID[clinic$Karyotype == "T21"],]
# cyt.D21 <- cyt[rownames(cyt) %in% clinic$LabID[clinic$Karyotype == "Control"],]


# save(met.T21, met.D21, cyt.T21, cyt.D21, clinic,file = sprintf("%sdata/Karyotype_cytokine_metabolite_Data.RData", DIR))

