## Inputs

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)

num_sub <- as.numeric(args[1])
dataset_file <- args[2]
pvalue_th <- as.numeric(args[3])
fc_th <- as.numeric(args[4])
column_seq <- as.matrix(as.numeric(args[5:length(args)]))


## Load data

library('Matching')

# Read miRNA dataset
X <- read.delim(dataset_file,
                header = TRUE,
                row.names = 1,
                sep = "\t")
#X <- as.matrix(X)
# Load Clusters
cluster <-
    read.delim(
        "../data/mirbase_cluster_g2.txt",
        header = FALSE,
        row.names = 1,
        sep = "\t",
        stringsAsFactors = FALSE
    )
# Load conversion table between IDs and pre-miRNA IDs
conversion_mi_miamt <-
    as.matrix(
        read.delim(
            "../data/conversion_IDpremirna_namemimat.txt",
            header = TRUE,
            sep = "\t",
            stringsAsFactors = FALSE
        )
    )

## Computation of fold change and KS on the miRNA matrix

# Initialise an empty table for the results
summary_new <- data.frame(matrix(ncol = 7, nrow = 0))
colnames_summ <- c("miRNA_cluster", "subtype_or_class", "sign", "members", "found", "pval", "fold_change")
names(summary_new) <- colnames_summ

#Loop through each cluster
for (i in 1:nrow(cluster)) {
    
    # Select miRNA IDs in current cluster (remove empty names)
    cluster_i <- cluster[i, ,drop=FALSE]
    cluster_i <- cluster_i[,cluster_i != ""]
    
    # Convert each miRNA ID in the cluster
    pos_converted <- conversion_mi_miamt[, 2] %in% cluster_i
    list_mirna <- conversion_mi_miamt[pos_converted, 1]
    list_mirna <- gsub("-5p", "", list_mirna, fixed=TRUE)
    list_mirna <- gsub("-3p", "", list_mirna, fixed=TRUE)
    list_mirna <- apply(expand.grid(unique(list_mirna), c("", "-3p", "-5p")), MARGIN = 1, paste0, collapse="")
    
    # Select matching lines in dataset
    # m <- match(row.names(X), list_mirna)
    # w <- which(!is.na(m))
    # X_mir <- X[w, , drop=FALSE]
    X_mir <- X[row.names(X) %in% list_mirna, , drop=FALSE]
    
    # If there was at least one match
    if (nrow(X_mir) != 0) {
        
        # For each class 
        for (class_idx in 1:num_sub) {
            
            # starting/ending index position of a class in the data table
            begin_class <- column_seq[(class_idx-1)*2+1, 1]
            end_class <- column_seq[(class_idx-1)*2+2, 1]

            # Get all data values for current class and current cluster
            values_for_class <- X_mir[, begin_class:end_class]
            # Get all data values for other classes and current cluster
            values_for_other_classes <- X_mir[, -1*(begin_class:end_class)]

            # For each miRNA in the cluster test significance of current class compared to the rest of samples
            p1 <- vapply(1:nrow(X_mir), FUN.VALUE = c(1),
                   FUN=function(current_mir) {
                       ks.boot(values_for_class[current_mir,], 
                               values_for_other_classes[current_mir,], 
                               alternative = "two.sided")$ks.boot.pvalue
                   })
            
            # Compute and store fold change
            # fc_keep <- paste0("[", rowMeans(values_for_class), "-", rowMeans(values_for_other_classes), "=", rowMeans(values_for_class) - rowMeans(values_for_other_classes), "]")
            fc1 <- rowMeans(values_for_class) - rowMeans(values_for_other_classes)
            
            #####control on significance

            # For each miRNA
            which_below_thresholds <- p1 < pvalue_th & abs(fc1) > fc_th
            at_least_one_dn <- sum(fc1[which_below_thresholds] < 0)
            at_least_one_up <- sum(fc1[which_below_thresholds] > 0)
            
            if (at_least_one_up >= 2) {
                summary_new <- rbind(summary_new, 
                                 cbind(
                                     row.names(cluster[i, ]), 
                                     class_idx, 
                                     'up', 
                                     paste0(cluster_i,collapse = "/"),
                                     paste0(row.names(X_mir),collapse = "/"),
                                     paste0(p1,collapse = "/"), 
                                     paste0(fc1,collapse = "/") 
                                    ) )
            }
            if (at_least_one_dn >= 2) {
                summary_new <- rbind(summary_new, 
                                 cbind(
                                     row.names(cluster[i, ]), 
                                     class_idx, 
                                     'down', 
                                     paste0(cluster_i,collapse = "/"),
                                     paste0(row.names(X_mir),collapse = "/"),
                                     paste0(p1,collapse = "/"), 
                                     paste0(fc1,collapse = "/") 
                                 ) )
            }
            
            # if (at_least_one_up < 2 && at_least_one_dn <2) { # CEH: Added to have all results
            #     summary_new <- rbind(summary_new, 
            #                      cbind(
            #                          row.names(cluster[i,]), 
            #                          class_idx, 
            #                          'none', 
            #                          paste0(cluster_i,collapse = "/"),
            #                          paste0(row.names(X_mir),collapse = "/"),
            #                          paste0(p1,collapse = "/"), 
            #                          paste0(fc1,collapse = "/") 
            #                      ) )
            # } # CEH: Added to have all results

        }
    }
}

write.table(
    summary_new,
    file = "../results/step1_results.txt",
    sep = "\t",
    row.names = FALSE,
    col.names = c("miRNA_cluster", "subype/class", "sign(up/down)", "members", "found", "pval", "fold_change")
)
