options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
th_STR <- as.numeric(args[1])

# Load results of step2
step2 <-
    read.delim(
        "../results/step2_results.txt",
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
    )

# Container for results
res <- matrix(ncol = 6, nrow = 0)

# For each DE cluster
for (j in 1:nrow(step2)) {
    
    # Cluster name
    mirna_DB <- step2[j, 1]
    
    # Read consolidated network
    my.data <-
        read.delim(
            paste0('../results/boot/', mirna_DB, '/network.txt'),
            header = TRUE,
            sep = '\t',
            stringsAsFactors = FALSE
        )
    # Remove miRNA-miRNA interactions
    my.data <- my.data[!grepl('hsa-', my.data[, 2]), ]
    # Targets of this cluster (MI)
    net_targets <- unique(my.data[, 2])
    
    # Read expression data used for this ARACNe analysis
    TCGA <-
        read.delim(
            paste0('../results/boot/', mirna_DB, "/", mirna_DB, '_expr_new.txt'),
            header = TRUE,
            sep = '\t',
            stringsAsFactors = FALSE
        )[,1]
    #TCGA <- TCGA[!is.na(TCGA)]
    
    # Read the signature 
    class <-
        read.delim(
            paste0('../data/signatures/', step2[j, 4], step2[j, 2], '.txt'),
            header = FALSE,
            sep = '\t',
            stringsAsFactors = FALSE
        )[, 1]
    
    # Expr genes not in targets
    is_gene_in_network <- TCGA %in% net_targets
    non_net <- TCGA[!is_gene_in_network]
    
    # Expr genes not in signature
    is_gene_in_signature <- TCGA %in% class
    non_class <- TCGA[!is_gene_in_signature]
    
    # Expr genes not in targets, in signature 
    is_sign_and_nonTargetExpr <- class %in% non_net
    nonnet_class <- length(class[is_sign_and_nonTargetExpr])
    
    # Expr genes not in signature, in targets
    is_notSign_target <- non_class %in% net_targets
    noclass_net <- length(non_class[is_notSign_target])

    # Expr genes not in signature not in target
    is_notSign_notTarget <- non_class %in% non_net
    nonet_noclass <- length(non_class[is_notSign_notTarget])

    # Genes in signature and target
    is_signature_and_target <- class %in% net_targets
    is_signature_and_target_inexpr <- class[is_signature_and_target] %in% TCGA
    class_net <- length(class[is_signature_and_target][is_signature_and_target_inexpr])

    # Prepare matrix for Fisher test
    mat <- matrix(c(class_net, noclass_net, nonnet_class, nonet_noclass),
                  nrow=2, byrow = TRUE)
    # Fisher test
    FET <- fisher.test(mat, alternative = "two.sided")
    
    #if (FET$p.value <= th_STR) {
        res <-
            rbind(res, c(mirna_DB, step2[j, 2], step2[j, 3], step2[j, 4], FET$p.value, 
                         paste(c(class_net, nonnet_class), c(noclass_net, nonet_noclass), sep=":", collapse="/")))
    #}
    
}

res <- cbind(res, p.adjust(res[,5], "fdr"))
colnames(res) <- c("mirna", "class", "mirna_sig", "signature_sign", "p-value", "Fisher_test_values", "FDR")

# Write results
write.table(
    res,
    '../clustMMRA_output.txt',
    sep = '\t',
    row.names = FALSE,
    col.names = c("mirna", "class", "mirna_sig", "signature_sign", "p-value", "Fisher_test_values", "FDR")
)
