
## Load data

# Load results of step 1
ttest <-
    read.delim(
        "../results/step1_results.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
# Load clusters and miRNA members
cluster <-
    read.delim(
        "../data/mirbase_cluster_g2.txt",
        header = FALSE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
# Load miRNA ID conversion table
conversion_mi_miamt <-
    read.delim(
        "../data/conversion_IDpremirna_namemimat.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
# Load miRNA target genes (2 sources)
X <-
    read.delim(
        "../data/all_no_mTB.txt",
        header = FALSE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
Y <- read.delim("../data/sperimentale.txt",
                header = FALSE,
                sep = "\t")
# List of genes that can be targeted by miRNA
genes <- unique(X[, 2])


#####
## I Classify miRNA-target interaction according to number of DBs

DB_2 <- X[X[,3] >= 2, ]
DB_3 <- X[X[,3] >= 3, ]
DB_4 <- X[X[,3] == 4, ]

# Load all signatures
signatures_up <-
    lapply(1:max(ttest$subype.class), 
           FUN = function(class_nb) {
               as.vector(read.delim(
                   paste0("../data/signatures/up", class_nb, ".txt"),
                   header = FALSE,
                   stringsAsFactors = FALSE,
                   sep = "\t"
               ))
           })
signatures_down <-
    lapply(1:max(ttest$subype.class), 
           FUN = function(class_nb) {
               as.vector(read.delim(
                   paste0("../data/signatures/down", class_nb, ".txt"),
                   header = FALSE,
                   stringsAsFactors = FALSE,
                   sep = "\t"
               ))
           })



#####
## II miRNA target enrichment test

# Results
riassunto <- data.frame(matrix(ncol = 7, nrow = 0))
colnames_summ <- c("cluster", "Expression_class", "miRNAcluster_sign", "signature", "minimal_number_DB", "Cluster_members", "miRNA_targets")
names(riassunto) <- colnames_summ

# For each differential miRNA cluster found in step 1
#### count microR's targets
for (j in 1:nrow(ttest)) {
    print(j)
    
    # Class index where cluster was found significantly DE
    class_DE_cluster <- as.numeric(ttest[j, 2])
    
    # Name of the cluster
    cl <- ttest[j, 1]
    # Find miRNA members of the cluster and (remove empty strings)
    cluster_i <- cluster[ cluster[, 1]==cl, ]
    cluster_i <- cluster_i[ , cluster_i != ""]
    # Convert identifiers    
    pos_converted <- conversion_mi_miamt[, 2] %in% cluster_i[,-1]
    list_mirna <- conversion_mi_miamt[pos_converted, 1]
    
    for (segno in c('up', 'down')) {
        
        # Load one of the two signatures of the class
        if(segno=='up') {
            sign <- signatures_up[[class_DE_cluster]][[1]]
        } else {
            sign <- signatures_down[[class_DE_cluster]][[1]]
        }

        # Filter DBs to keep only the miRNA targets
        DB2_mirna <- DB_2[DB_2[, 1] %in% list_mirna, ]
        DB3_mirna <- DB_3[DB_3[, 1] %in% list_mirna, ]
        DB4_mirna <- DB_4[DB_4[, 1] %in% list_mirna, ]
        Y_mirna <- Y[Y[, 1] %in% list_mirna, ]
        
        # Compute number of targets of cluster members that are present in the signature
        t_num_2 <- sum(DB2_mirna[,2] %in% sign)
        t_num_3 <- sum(DB3_mirna[,2] %in% sign)
        t_num_4 <- sum(DB4_mirna[,2] %in% sign)
        t_num_s <- sum(Y_mirna[,2] %in% sign)
        
        ## Compute null model for random genes of same dimension of targets of each cluster

        # 1000 random samples of gene targets (= 1000 random signatures)
        # Caution: replicate binds samples by column
        length_sign <- length(sign)
        random_samples <- replicate(1000, sample(genes, length_sign, replace = FALSE, prob = NULL))
        
        # Filter DBs to keep only the random signature genes
        target_r <- 
            t(apply(random_samples, MARGIN = 2, 
                    function(sign_r) {
                        #### compute targets of the cluster in sign_r
                        # Compute number of targets of cluster members that are present in the random signature, for each database
                        return(c(
                            DB2=sum(DB2_mirna[, 2] %in% sign_r), 
                            DB3=sum(DB3_mirna[, 2] %in% sign_r), 
                            DB4=sum(DB4_mirna[, 2] %in% sign_r), 
                            S=sum(Y_mirna[, 2] %in% sign_r)
                        ))
                    })
            )
        
        # Compute 5% thresholds under the null distribution
        soglia <- apply(target_r, MARGIN=2, quantile, 0.95)
        
        # Save result if cluster value is above one of the thresholds
        # Note that being significant for db 2 is restrictive, as it's the biggest database.
        if (t_num_2 > soglia["DB2"]) {
            riassunto <-
                rbind(riassunto, cbind(cl, class_DE_cluster, ttest[j, 3], segno, 2, paste0(list_mirna, collapse = "/"),
                                       paste(DB2_mirna[DB2_mirna[,2] %in% sign, 1], DB2_mirna[DB2_mirna[,2] %in% sign, 2], collapse="/", sep=":")))

            if (t_num_3 > soglia["DB3"]) {
                riassunto <-
                    rbind(riassunto, cbind(cl, class_DE_cluster, ttest[j, 3], segno, 3, paste0(list_mirna, collapse = "/"),
                                           paste(DB3_mirna[DB3_mirna[,2] %in% sign, 1], DB3_mirna[DB3_mirna[,2] %in% sign, 2], collapse="/", sep=":")))
            }
            if (t_num_4 > soglia["DB4"]) {
                riassunto <-
                    rbind(riassunto, cbind(cl, class_DE_cluster, ttest[j, 3], segno, 4, paste0(list_mirna, collapse = "/"),
                                           paste(DB4_mirna[DB4_mirna[,2] %in% sign, 1], DB4_mirna[DB4_mirna[,2] %in% sign, 2], collapse="/", sep=":")))
            }
            if (t_num_s > soglia["S"]) {
                riassunto <-
                    rbind(riassunto, cbind(cl, class_DE_cluster, ttest[j, 3], segno, 's', paste0(list_mirna, collapse = "/"),
                                           paste(Y_mirna[Y_mirna[,2] %in% sign, 1], Y_mirna[Y_mirna[,2] %in% sign, 2], collapse="/", sep=":")))
            }
        }
    }
}

## Export results
write.table(
    riassunto,
    "../results/step2_results.txt",
    col.names = c("cluster", "Expression_class", "miRNAcluster_sign", "signature", "minimal_number_DB", "Cluster_members", "miRNA_targets"),
    row.names = FALSE,
    sep = "\t"
)
