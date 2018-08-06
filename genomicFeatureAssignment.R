#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

library(GenomicRanges)
library(ggplot2) 
library(reshape) 
library(reshape2)
library(made4)
library(cluster)
library(spatstat) # "im" function 

retrieveCores <- function(dir){
  return(read.table(dir, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
}

retrieveOrganoidSlices <- function(dir){
  return(read.table(dir, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
}

retrieveTrainingSet <- function(loaded_samples, gene_list, impute = FALSE, reference, res_dir){
  target_samples <- load_samples(classes = c("T", "F", "M"), sampleList = "./resources/sampleList.csv")
  
  setwd("~/Git-Projects/Git-Research-Projects/FACETS_write_files/")
  matrix_training_set <- data.frame(stringsAsFactors = FALSE)
  for(sample in target_samples){
    print(sample)
    facets_segment_data <- retrieveFacetsSegments(sample, sample_subdir = "/", reference = reference, dir = res_dir)
    facets_segment_data <- segmentsToBedFormat(facets_segment_data, median.clust = FALSE)
    gene_score <- data.frame(stringsAsFactors = FALSE)
    for(gene.i in seq(nrow(gene_list))){
      gene <- gene_list[gene.i, ]
      segments <- facets_segment_data[facets_segment_data$chrom == gene$chrom & ((facets_segment_data$start <= gene$start & gene$start <= facets_segment_data$end) | (facets_segment_data$start <= gene$end & gene$end <= facets_segment_data$end)), ]
      # TODO: genes may intersect across two different segments - will create logic to determine value
      value <- ifelse(nrow(segments) > 0, segments$value, NA)
      gene_score <- rbind(gene_score, cbind(gene, value))
    }
    gene_score <- gene_score[,c(1,2,3,5,4)]
    colnames(gene_score)[4] <- "value"
    feature.entry <- t(data.frame(gene_score[,"value"], row.names = gene_score$name))
    row.names(feature.entry) <- sample
    matrix_training_set <- rbind(matrix_training_set, feature.entry)  
  }
  
  # Impute NAs with column mean
  if(impute == TRUE){
    for(i in 1:ncol(matrix_training_set)){
      matrix_training_set[is.na(matrix_training_set[,i]), i] <- median(matrix_training_set[,i], na.rm = TRUE)
    }
  }
  
  melted_training_set <- do.call(rbind, lapply(seq(nrow(matrix_training_set)), function(index){
    return(do.call(rbind, lapply(colnames(matrix_training_set[index, ]), function(coreId, index){
      coreEntry <- data.frame(score = matrix_training_set[index, coreId], coreId = coreId, sampleId = rownames(matrix_training_set[index, ])[1])
      return(coreEntry)
    }, index)))
  }))
  
  return(list(melted=melted_training_set, matrix=matrix_training_set))
}

attachLabelsToSet <- function(matrix_training_set, labelData){
  sampleList <- rownames(matrix_training_set)
  labelLists <- lapply(names(labelData), function(label){
    aucList <- unlist(sapply(sampleList, function(sample, label){
      labelMatrix <- labelData[[label]]  
      auc <- c(labelMatrix[labelMatrix$SampleId == sample, ]$AUC, NA)[1]
      return(auc)
    }, label))
    return(aucList)
  })
  names(labelLists) <- names(labelData)
  labelDataframe <- do.call(cbind.data.frame, labelLists)
  labeled_matrix_training_set <- cbind(labelDataframe, matrix_training_set)
  return(labeled_matrix_training_set)
}

visualizeUnclusteredHeatmap <- function(training_set){
  ggplot(data = training_set, aes(x = coreId, y = sampleId)) + 
    geom_tile(aes(fill = score), color = "white", size = 1) + 
    scale_fill_gradient2(low = "blue", mid="white", high = "tomato") + 
    xlab("core ID") + 
    theme_grey(base_size = 10) + 
    ggtitle("Heatmap (ggplot)") + 
    theme(axis.ticks = element_blank(), 
          panel.background = element_blank(), 
          plot.title = element_text(size = 12, colour = "gray50")) 
}

clusterTrainingSet <- function(training_set, visualize = FALSE){
  # Unmelt training set for correlation analysis
  training_set_matrix <- dcast(data = training_set,formula = sampleId~coreId,fun.aggregate = sum,value.var = "score")
  sampleIds <- training_set_matrix$sampleId
  training_set_matrix <- training_set_matrix[,-c(1)]
  training_set_matrix <- t(training_set_matrix)
  colnames(training_set_matrix) <- sampleIds
  training_set_matrix <- as.data.frame(training_set_matrix)
  
  # Remove samples with 0 variance
  nonzero_variance_samples <- unlist(lapply(colnames(training_set_matrix), function(sample){
    if(var(na.omit(training_set_matrix[, sample])) != 0){
      return(sample)
    }
  }))
  training_set_matrix <- training_set_matrix[,c(nonzero_variance_samples)]
  
  training_set_matrix <- training_set_matrix[complete.cases(training_set_matrix), ]
  # Calculate distance matrix
  corRaw <- cor(training_set_matrix)
  dissimilarity <- 1 - corRaw
  distance.sample <- as.dist(dissimilarity)
  distance.core <- as.dist(t(dissimilarity))
  
  # Run hierarchical clustering
  hc.sample <- hclust(distance.sample)
  hc.core <- hclust(distance.core)
  
  if(visualize == TRUE){
    color.palette  <- colorRampPalette(c("blue4", "white", "red4"))(n=600)
#    heatmap.2(t(training_set_matrix),  trace="none", dendrogram="row", Rowv = as.dendrogram(hc.sample), Colv = as.dendrogram(hc.core), density.info = 'none', scale='none', col = color.palette)
    heatmap.2(t(training_set_matrix), trace = "none", density="none", scale = "row", hclust=function(x) hclust(x, method="complete"), distfun=function(x) as.dist((1-cor(t(x)))/2), col = color.palette)
  }
  
  return(hc.sample)  
}