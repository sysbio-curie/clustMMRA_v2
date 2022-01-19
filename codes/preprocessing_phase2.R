options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
dataset_mirna_file <- args[1]
dataset_mrna_file <- args[2]

# Load miRNA and RNA data
TCGA_mirna <-
    read.table(
        dataset_mirna_file,
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE,
        sep = '\t'
    )
TCGA_mrna <-
    read.table(
        dataset_mrna_file,
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE,
        sep = '\t'
    )
# Check that they are similarly ordered
if (all(colnames(TCGA_mirna) == colnames(TCGA_mrna))) {
    print("columns of the datasets correctly ordered")
} else{
    print("ATTENTION! the columns of the mRNA and microRNA dataset are not in the same order. Reorder them and rerun")
}

# Load results from Step 2
cluster <-
    as.matrix(read.delim(
        "../results/step2_results.txt",
        header = T,
        sep = "\t"
    ))

# Load miRNA IDs conversion table
conversion_mi_miamt <-
    as.matrix(
        read.delim(
            "../data/conversion_IDpremirna_namemimat.txt",
            header = TRUE,
            sep = "\t"
        )
    )
# Load miRNA cluster members
conv_clust_singlemir <-
    as.matrix(read.delim(
        "../data/mirbase_cluster_g2.txt",
        header = FALSE,
        sep = "\t"
    ))

library('preprocessCore')

# For each significantly DE cluster of miRNA with enriched targets
for (i in 1:nrow(cluster)) {
    
    print(i)
    
    # Find miRNA members of the cluster and (remove empty strings)
    cluster_i <- conv_clust_singlemir[ conv_clust_singlemir[, 1]==cluster[i, 1], ]
    cluster_i <- cluster_i[ cluster_i != ""]
    cluster_i <- matrix(cluster_i) # to remove
    
    # Convert identifiers    
    pos_converted <- conversion_mi_miamt[, 2] %in% cluster_i[-1,]
    list_mirna <- conversion_mi_miamt[pos_converted, 1]

    # list_mirna <- gsub("-5p", "", list_mirna)
    # list_mirna <- gsub("-3p", "", list_mirna)
    # list_mirna2 <- paste(list_mirna, '-1', sep = '')
    # list_mirna3 <- paste(list_mirna, '-2', sep = '')
    # list_mirna <-
    #     as.matrix(rbind(
    #         as.matrix(list_mirna),
    #         as.matrix(list_mirna2),
    #         as.matrix(list_mirna3)
    #     ))
    
    # miRNA expression
    mat <- TCGA_mirna[list_mirna,]
    names <- row.names(mat)
    dim_cluster <- nrow(mat)
    
    # Concatenates to RNA data
    matrice <- rbind(TCGA_mrna, mat)

    # Matrix normalization 
    # NB data from RNA-seq is in RPKM (and it's ok :) 
    matrice2 <- normalize.quantiles(as.matrix(exp(matrice)))
    matrice2 <- log2(matrice2)
    rownames(matrice2) <- rownames(matrice)
    colnames(matrice2) <- colnames(matrice)
    # Round very small neg values written as exponential
    matrice2[matrice2 < 0] <- 0
    
    # Filter to keep only RNA-seq lines with a certain variation (sd above 1.2)
    indexes_wo_mirna <- -1 * (nrow(matrice2)-dim_cluster+1):nrow(matrice2)
    matrice_new <- matrice2[indexes_wo_mirna, ][apply(matrice2[indexes_wo_mirna, ], MARGIN=1, sd) > 1.2, ]
    # Bind again miRNA cluster data 
    matrice_new <- rbind(matrice_new, matrice2[-1*indexes_wo_mirna, ])
    # Add 2 columns of annotations
    tmp <- cbind(rownames(matrice_new), matrice_new)
    # Add one line with column names
    matrice_new <- rbind(c('gene', colnames(matrice_new)), 
                         tmp)
    
    output_folder <- paste0('../results/boot/', cluster[i,1], "/" )
    dir.create(output_folder, showWarnings = FALSE, recursive = TRUE) # Don't warn if folder already exists
    write.table(
        matrice_new,
        paste0(output_folder, cluster[i,1], '_expr_new.txt'),
        sep = '\t', quote=FALSE,
        col.names = FALSE,
        row.names = FALSE
    )
    write.table(
        as.matrix(names),
        paste0(output_folder , cluster[i,1], '_mirna_new.txt'),
        sep = '\t', quote = FALSE,
        col.names = FALSE,
        row.names = FALSE
    )
}
