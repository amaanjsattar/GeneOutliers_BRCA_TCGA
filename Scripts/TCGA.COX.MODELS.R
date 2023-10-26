output_dir <- "tcga_cox_modeling_output";
dir.create(output_dir, showWarnings = FALSE);

library(BoutrosLab.plotting.survival);
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);
library(BoutrosLab.statistics.survival);
library(survival);
load('/Users/amaanjsattar/Desktop/2023-07-07_TCGA_BRCA_Outlier.rda')

tcga.fpkm <- fpkm.tumor.symbol.filter
tcga.outlier.tags <- outlier.patient.tag.01.t.p.order
tcga.fpkm <- tcga.fpkm[rownames(tcga.outlier.tags), ]
rownames(tcga.fpkm) <- tcga.fpkm$Symbol
tcga.fpkm$Symbol <- NULL
# Subset the clinical DataFrame to keep only the desired columns
tcga.subset <- brca.clinic[, c('Overall.Survival..Months.', 'Overall.Survival.Status')]

# transpose tcga.fpkm so that patients are rows and genes are columns
tcga.fpkm <- t(tcga.fpkm)

extract_first_12_chars <- function(row_names) {
    substr(row_names, 1, 12)
}

tcga.fpkm.substring <- extract_first_12_chars(rownames(tcga.fpkm))
sum(duplicated(tcga.fpkm.substring))

duplicate_indices <- which(duplicated(tcga.fpkm.substring) | duplicated(tcga.fpkm.substring, fromLast = TRUE))

rows_with_duplicates <- tcga.fpkm[duplicate_indices, ]
indices_to_remove <- duplicate_indices[!grepl("01A$", rownames(tcga.fpkm[duplicate_indices, ]))]
tcga.fpkm.filtered <- tcga.fpkm[-indices_to_remove, ]
tcga.fpkm.substring <- extract_first_12_chars(rownames(tcga.fpkm.filtered))

tcga.subset.filtered <- tcga.subset[rownames(tcga.subset) %in% tcga.fpkm.substring, ]
### FUNCTIONS (TAKEN FROM TCGA.FUNCTIONS.R) ########################################################

binary.survival <- function(brca.data, old.col.name, new.col.name) {
    brca.data[[new.col.name]] <- ifelse(
        grepl('^1:', brca.data[[old.col.name]]), 1, 0)
    return(brca.data)
};


create.surv <- function(brca.data, survival.time, survival.status) {
    surv.obj <- Surv(as.numeric(brca.data[[survival.time]]),
                     brca.data[[survival.status]])
    return(surv.obj)
};

####################################################################################################
tcga.subset.filtered <- binary.survival(
    tcga.subset.filtered, 'Overall.Survival.Status', 'OS.Status'
);

tcga.surv.cox <- create.surv(
    tcga.subset.filtered, 'Overall.Survival..Months.', 'OS.Status'
);

tcga.cox.results <- list()
tcga.genes <- colnames(tcga.fpkm.filtered)

for (gene in tcga.genes) {
    # Extract the FPKM values for the current gene (as a vector)
    fpkm.values <- log2(as.numeric(tcga.fpkm.filtered[, gene]) + 1)
    
    # Fit the Cox proportional hazards model using the 'coxph' function from the 'survival' package
    cox.model <- coxph(tcga.surv.cox ~ fpkm.values)
    
    # Store the model and its summary in the cox.model.results list
    tcga.cox.results[[gene]] <- list(model = cox.model, summary = summary(cox.model))
}

gene_names_vector <- c()
p_values_vector <- c()

# Loop through each gene in gene.names and extract the p-value from the summary
for (gene in tcga.genes) {
    # Access the summary object for the specific gene
    summary_for_gene <- tcga.cox.results[[gene]]$summary
    
    # Access the p-value for the gene
    p_value <- summary_for_gene$coefficients["fpkm.values", "Pr(>|z|)"]
    
    # Store the gene name and p-value in the corresponding vectors
    gene_names_vector <- c(gene_names_vector, gene)
    p_values_vector <- c(p_values_vector, p_value)
}

# Create a dataframe with p-value results
tcga_p_values <- data.frame(GeneName = gene_names_vector, PValue = p_values_vector)

# Sort p-values in ascending order
tcga_p_values <- tcga_p_values[order(tcga_p_values$PValue), ]

fdr_corrected_p_values <- p.adjust(tcga_p_values$PValue, method = 'fdr')

tcga_p_values$FDR.Corrected <- fdr_corrected_p_values
save(tcga.cox.results, tcga_p_values, file = 'tcga.cox.RData')

