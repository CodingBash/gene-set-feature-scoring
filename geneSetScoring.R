setwd("~/Git-Projects/Git-Research-Projects/gene-set-feature-scoring/")
source("helperFunctions.R")
source("genomicFeatureAssignment.R")

gene_list <- read.table("./resources/genes_hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gene_list <- gene_list[, c(7,8,9,3)]
colnames(gene_list) <- c("chrom", "start", "end", "name")

reference <- "hN30"
res_dir <- "output/FACETS_Reference_hN30_8_2_18_1/" # Determine FACETS reference to use _ hN30

target_samples <- load_samples(classes = c("N","T", "F", "M"), sampleList = "./resources/sampleList.csv")
target_samples <- target_samples[target_samples != reference]
target_samples <- target_samples[target_samples != "hN32"]
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

# TODO: Impute NAs with column mean
for(i in 1:ncol(matrix_training_set)){
  matrix_training_set[is.na(matrix_training_set[,i]), i] <- median(matrix_training_set[,i], na.rm = TRUE)
}

melted_training_set <- do.call(rbind, lapply(seq(nrow(matrix_training_set)), function(index){
  return(do.call(rbind, lapply(colnames(matrix_training_set[index, ]), function(coreId, index){
    coreEntry <- data.frame(score = matrix_training_set[index, coreId], coreId = coreId, sampleId = rownames(matrix_training_set[index, ])[1])
    return(coreEntry)
  }, index)))
}))

# TODO: distance matrix to hc may be incorrect (based on Euclidean instead of pearsons)
visualizeUnclusteredHeatmap(melted_training_set)
hc <- clusterTrainingSet(melted_training_set, visualize = TRUE)

