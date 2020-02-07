library(Ckmeans.1d.dp)
library(tidyverse)

# A script that attempts to retrieve clusters from a vector representation of a set of sounds.


# The path to and base filename of the .data, .sounds, and .contexts files to read.
DATA_ROOT <- "C:/Users/Sheng-Fu/Documents/GitHub/distributional_learning/vector_data/quechua-pbase_W2V_0"


# Whether or not we save the discovered classes to a file.
SAVE_CLASSES = TRUE

# Where to save the list of discovered classes
OUTPUT_DIRECTORY = "C:/Users/Sheng-Fu/Documents/GitHub/distributional_learning/output/"

# Optional suffix to add to the names of saved files
FILE_SUFFIX = ""

# A parameter that controls what amount of variance a principal component must account for to
# be used in clustering. The threshold is VARIABILITY_SCALAR * (average amount of variance).
VARIABILITY_SCALAR <- 1

# A parameter that, if TRUE, sets restrictions on the initial partition of the data set: namely,
# any partition of the full set of sounds must be into two classes (e.g. consonants vs. vowels,
# voiced vs. voiceless, etc.)
CONSTRAIN_INITIAL_PARTITION_NUMBER <- TRUE

# A parameter that, if TRUE, restricts the initial partition of the data set: namely, 
# only the first principal component is considered. Setting this to FALSE will result
# in the same classes being detected as when it is TRUE, but with additional partitions
# of the data set potentially discovered as well. Similar results can be gained by
# increasing VARIABILITY_SCALAR, but this will apply to all recursive calls to the clusterer
# rather than just the top level call.
CONSTRAIN_INITIAL_PARTITION_PCS <- TRUE

# Read files
values <- read.csv(paste(DATA_ROOT,".data",sep=""), sep=" ", header=FALSE)
values

sounds <- read.csv(paste(DATA_ROOT,".sounds",sep=""), sep=" ", header=FALSE, 
                   stringsAsFactors=FALSE, colClasses = c("character"))

dim(values)

sounds
contexts <- read.csv(paste(DATA_ROOT,".contexts",sep=""), sep=" ", header=FALSE)
contexts

values

rownames(values) <- t(sounds)
colnames(values) <- t(contexts)

colnames(values)
rownames(values)

values = values[!rownames(values) %in% c('', '#'),]
sounds = sounds[! sounds %in% c('', '#')]

#####value analysis here

# values1 = select(values, matches("#-"))
# values1
# 
# values2 = log(values + 0.00001) 
# values2[,1]
# 
# rownames(values)
# consonants = rownames(values)[c(F, T, T, F, T, 
#                               F, T, T, T, T, 
#                               T, T, F, T, T,
#                               T, T, F, T, T,
#                               T, F, F, T, T,
#                               T, T, T, T)]
# consonants
# vowels = rownames(values)[!c(F, T, T, F, T, 
#                                 F, T, T, T, T, 
#                                 T, T, F, T, T,
#                                 T, T, F, T, T,
#                                 T, F, F, T, T,
#                                 T, T, T, T)]
# vowels
# 
# consonants
# colnames(values)
# 
# initial_CCs = c()
# for (i in consonants){
#   expression = paste("#", "_", i, sep = "-")
#   initial_CCs = c(initial_CCs, expression)
# }
# initial_CCs = initial_CCs[initial_CCs %in% colnames(values)]
# 
# values_CCs = select(values, initial_CCs) %>% filter()
# values_CCs_Cs = subset(values_CCs, rownames(values_CCs) %in% consonants)
# 
# #View(values_CCs)
# View(values_CCs_Cs)
# 
# values_CCs_Cs = values_CCs_Cs[, colSums(values_CCs_Cs != 0) > 0]

#####


get_classes <- function(input_data, 
                        constrain_pcs=FALSE,
                        constrain_partition_number=FALSE,
                        variability_scalar=1,
                        visited_classes=list()) {
  # input_data: the matrix to factor and cluster
  # constrain_pcs: if TRUE, this partition of the entire set of sounds will be
  #                restricted to only the first PC. Otherwise all PCs up to the 
  #                variance threshold will be explored.
  # constraint_partition_number: if TRUE, the partition of the entire set of 
  #                              sounds will be forced to be into two classes
  #                              (e.g. consonants and vowels). If FALSE, any
  #                              number of classes between 1 and 3 may be
  #                              retrieved.
  # variability scalar: Scales the amount of variance that a PC must capture to be
  #                     used for clustering. For example, 1 indicates that a PC
  #                     must capture at least the average amount of variance of
  #                     all PCs to be used. 2 means it must capture twice the average
  #                     amount of variance, etc.
  # visited_classes: Keeps track of which classes we've already produced to avoid
  #                  unnecessary duplication of effort. You shouldn't need to use this.

  full_classes_list <- list()
  local_sounds <- row.names(input_data)
  
  # Do a PCA on the input data.
  pca_data <- prcomp(na.omit(input_data),  center=TRUE)
  #print(pca_data)
  
  if (constrain_pcs) {
    highest_dim <- 1
  }
  else {
    # If we're looking at all PCs, calculate which ones we will examine
    # based on Kaiser's stopping criterion. Cluster into between 1-3 clusters.
    mean_sdev <- mean(pca_data$sdev) * variability_scalar
    # print(mean_sdev)
    if (mean_sdev == 0){highest_dim <- 1}
    else {highest_dim <- max(1, min(which(pca_data$sdev <  mean_sdev)) - 1)}
  }
  if (constrain_partition_number) {
    # Only cluster into a maximum of two classes.
    cluster_range <- c(1,2)
  }
  else {
    # Cluster into a maximum of three classes.
    cluster_range <- c(1,3)
  }
  

  
  # Go through all the PCs we want to cluster over
  for (i in 1:highest_dim) {
    classes_list <- list()
    sub_classes <- list()
    col <- pca_data$x[,i]
  
    # Do 1D k-means clustering on this PC
    # Ckmeans.1d.dp issues a warning when the maximum number of clusters is smaller
    # than the length of the list of points it's clustering. Suppress it for sanity.
    result <- suppressWarnings(Ckmeans.1d.dp(col, cluster_range))
    
    # Add each of the discovered clusters to the list of discovered clusters
    num_clusters <- max(result$cluster)
    for (j in 1:num_clusters) {
      classes_list <- c(classes_list, list(local_sounds[result$cluster == j]))
    }
    
    # Go through all the discovered classes and see if we should cluster over them.
    for (j in 1:length(classes_list)) {
      class_vec <- classes_list[[j]]
      # Check that the class contains more than two objects and is a subset of the input sounds
      meaningful_partition <- length(class_vec) > 2 && length(class_vec) < length(local_sounds)
      # Check that we haven't already clustered this subset. This isn't strictly necessary,
      # but saves some cycles.
      visited <- list(class_vec) %in% visited_classes
      if (meaningful_partition && !visited) {
        # If this class looks worth partitioning, flag it as seen and apply the clustering
        # algorithm recursively.
        visited_classes <- c(visited_classes, list(class_vec))
        sub_classes <- c(
          sub_classes, 
          get_classes(
            input_data[class_vec,],
            variability_scalar=variability_scalar,
            visited_classes=visited_classes
          )
        )
      }
    }
    # Add the classes we found in this call and all recursive calls to the list of
    # discovered classes.
    full_classes_list <- c(full_classes_list, classes_list, sub_classes)
    print(length(unique(full_classes_list)))
  }
  # Return the classes we found, removing any duplicates.
  return(unique(full_classes_list))
}

VARIABILITY_SCALAR <- 1
CONSTRAIN_INITIAL_PARTITION_NUMBER <- T
CONSTRAIN_INITIAL_PARTITION_PCS <- T


values

# Main call to do clustering
classes <- get_classes(
  values,
  variability_scalar=VARIABILITY_SCALAR,
  constrain_pcs=CONSTRAIN_INITIAL_PARTITION_PCS,
  constrain_partition_number=CONSTRAIN_INITIAL_PARTITION_NUMBER
)



# Print the results to console
for (i in 1:length(classes)) {
  print(classes[[i]])
}

if (SAVE_CLASSES) {
  # Write the results to a file and print them to the console
  file_prefix <- basename(DATA_ROOT)
  filename <- paste(OUTPUT_DIRECTORY, file_prefix, FILE_SUFFIX, '_V1_T_T.classes', sep="")
  handle <-file(filename, 'w+')
  # This line writes up all the segments into the first line of the file
  writeLines(paste(as.character(sounds[1,]), collapse=' '), handle)
  for (i in 1:length(classes)) {
    writeLines(paste(classes[[i]], collapse=' '), handle)
  }
  
  close(handle)
}

VARIABILITY_SCALAR <- 0.75
CONSTRAIN_INITIAL_PARTITION_NUMBER <- T
CONSTRAIN_INITIAL_PARTITION_PCS <- T


values

# Main call to do clustering
classes <- get_classes(
  values,
  variability_scalar=VARIABILITY_SCALAR,
  constrain_pcs=CONSTRAIN_INITIAL_PARTITION_PCS,
  constrain_partition_number=CONSTRAIN_INITIAL_PARTITION_NUMBER
)



# Print the results to console
for (i in 1:length(classes)) {
  print(classes[[i]])
}

if (SAVE_CLASSES) {
  # Write the results to a file and print them to the console
  file_prefix <- basename(DATA_ROOT)
  filename <- paste(OUTPUT_DIRECTORY, file_prefix, FILE_SUFFIX, '_V0.75_T_T.classes', sep="")
  handle <-file(filename, 'w+')
  # This line writes up all the segments into the first line of the file
  writeLines(paste(as.character(sounds[1,]), collapse=' '), handle)
  for (i in 1:length(classes)) {
    writeLines(paste(classes[[i]], collapse=' '), handle)
  }
  
  close(handle)
}

VARIABILITY_SCALAR <- 0.5
CONSTRAIN_INITIAL_PARTITION_NUMBER <- T
CONSTRAIN_INITIAL_PARTITION_PCS <- T


values

# Main call to do clustering
classes <- get_classes(
  values,
  variability_scalar=VARIABILITY_SCALAR,
  constrain_pcs=CONSTRAIN_INITIAL_PARTITION_PCS,
  constrain_partition_number=CONSTRAIN_INITIAL_PARTITION_NUMBER
)



# Print the results to console
for (i in 1:length(classes)) {
  print(classes[[i]])
}

if (SAVE_CLASSES) {
  # Write the results to a file and print them to the console
  file_prefix <- basename(DATA_ROOT)
  filename <- paste(OUTPUT_DIRECTORY, file_prefix, FILE_SUFFIX, '_V0.5_T_T.classes', sep="")
  handle <-file(filename, 'w+')
  # This line writes up all the segments into the first line of the file
  writeLines(paste(as.character(sounds[1,]), collapse=' '), handle)
  for (i in 1:length(classes)) {
    writeLines(paste(classes[[i]], collapse=' '), handle)
  }
  
  close(handle)
}

VARIABILITY_SCALAR <- 0.25
CONSTRAIN_INITIAL_PARTITION_NUMBER <- T
CONSTRAIN_INITIAL_PARTITION_PCS <- T




values

# Main call to do clustering
classes <- get_classes(
  values,
  variability_scalar=VARIABILITY_SCALAR,
  constrain_pcs=CONSTRAIN_INITIAL_PARTITION_PCS,
  constrain_partition_number=CONSTRAIN_INITIAL_PARTITION_NUMBER
)



# Print the results to console
for (i in 1:length(classes)) {
  print(classes[[i]])
}

if (SAVE_CLASSES) {
  # Write the results to a file and print them to the console
  file_prefix <- basename(DATA_ROOT)
  filename <- paste(OUTPUT_DIRECTORY, file_prefix, FILE_SUFFIX, '_V0.25_T_T.classes', sep="")
  handle <-file(filename, 'w+')
  # This line writes up all the segments into the first line of the file
  writeLines(paste(as.character(sounds[1,]), collapse=' '), handle)
  for (i in 1:length(classes)) {
    writeLines(paste(classes[[i]], collapse=' '), handle)
  }
  
  close(handle)
}



 VARIABILITY_SCALAR <- 1
 CONSTRAIN_INITIAL_PARTITION_NUMBER <- T
 CONSTRAIN_INITIAL_PARTITION_PCS <- F


 values

# Main call to do clustering
 classes <- get_classes(
   values,
   variability_scalar=VARIABILITY_SCALAR,
   constrain_pcs=CONSTRAIN_INITIAL_PARTITION_PCS,
   constrain_partition_number=CONSTRAIN_INITIAL_PARTITION_NUMBER
 )



 # Print the results to console
 for (i in 1:length(classes)) {
   print(classes[[i]])
 }

 if (SAVE_CLASSES) {
   # Write the results to a file and print them to the console
   file_prefix <- basename(DATA_ROOT)
   filename <- paste(OUTPUT_DIRECTORY, file_prefix, FILE_SUFFIX, '_V1_T_F.classes', sep="")
   handle <-file(filename, 'w+')
   # This line writes up all the segments into the first line of the file
   writeLines(paste(as.character(sounds[1,]), collapse=' '), handle)
   for (i in 1:length(classes)) {
     writeLines(paste(classes[[i]], collapse=' '), handle)
   }

   close(handle)
 }
# 
#   VARIABILITY_SCALAR <- 0.25
#  CONSTRAIN_INITIAL_PARTITION_NUMBER <- F
# CONSTRAIN_INITIAL_PARTITION_PCS <- F
# 
# 
# # values
# 
# # Main call to do clustering
#  classes <- get_classes(
#    values,
#    variability_scalar=VARIABILITY_SCALAR,
#    constrain_pcs=CONSTRAIN_INITIAL_PARTITION_PCS,
#    constrain_partition_number=CONSTRAIN_INITIAL_PARTITION_NUMBER
#  )
# 
# 
# 
#  # Print the results to console
#  for (i in 1:length(classes)) {
#    print(classes[[i]])
#  }
# 
#  if (SAVE_CLASSES) {
#    # Write the results to a file and print them to the console
#    file_prefix <- basename(DATA_ROOT)
#    filename <- paste(OUTPUT_DIRECTORY, file_prefix, FILE_SUFFIX, '_V0.25_F_F.classes', sep="")
#    handle <-file(filename, 'w+')
#    # This line writes up all the segments into the first line of the file
#    writeLines(paste(as.character(sounds[1,]), collapse=' '), handle)
#    for (i in 1:length(classes)) {
#      writeLines(paste(classes[[i]], collapse=' '), handle)
#    }
# 
#    close(handle)
#  }

#
