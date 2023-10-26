### UNIVARIATE COX PROPORTIONAL HAZARDS MODELING FOR OUTLIER GENES AND OVERALL SURVIVAL ############
output_dir <- "cox_modeling_output";
dir.create(output_dir, showWarnings = FALSE);

library(BoutrosLab.plotting.survival);
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);
library(BoutrosLab.statistics.survival);
library(survival);
load("/Users/amaanjsattar/Desktop/2023-07-07_Metabric_Outlier.rda");




# Remove metadata
meta.clinic.patient <- meta.clinic.patient[-c(1:4), ];

# Remove rows lacking survival status
meta.clinic.patient <- subset(
    meta.clinic.patient, 
    meta.clinic.patient$Overall.Survival.Status != '')

names(meta.clinic.patient)[names(meta.clinic.patient) == 'Overall.Survival..Months.'] <- 'OS.Months'
names(meta.clinic.patient)[names(meta.clinic.patient) == 'Overall.Survival.Status'] <- 'OS.Status'


meta.fpkm <- fpkm.tumor.symbol.filter.symbol
meta.symbols <- as.data.frame(meta.fpkm$Symbol)
meta.fpkm <- meta.fpkm[rownames(outlier.patient.tag.01), ]
# Remove rows with NA values in the 'Symbol' column
meta.fpkm <- meta.fpkm[complete.cases(meta.fpkm$Symbol), ]
rownames(meta.fpkm) <- meta.fpkm$Symbol
meta.fpkm$Symbol <- NULL

# Convert the patient IDs in 'meta.clinic.patient' to the desired format (with dashes)
colnames(meta.fpkm) <- gsub("\\.", "-", colnames(meta.fpkm))

# Transpose the 'meta.fpkm' dataframe to have patients as rows and genes as columns
meta.fpkm <- t(meta.fpkm)

# Get the common patient IDs between 'meta.clinic.patient' and 'meta.fpkm'
common.patients <- intersect(rownames(meta.clinic.patient), rownames(meta.fpkm))

# Subset 'meta.fpkm' to include only the common patients
meta.fpkm <- meta.fpkm[common.patients, ]
meta.clinic.patient <- meta.clinic.patient[common.patients, ]

# Reorder meta.fpkm so that it matches the order of meta.clinic.patient
meta.fpkm <- meta.fpkm[match(rownames(meta.clinic.patient), rownames(meta.fpkm)), ]

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
meta.clinic.patient <- binary.survival(
    meta.clinic.patient, 'OS.Status', 'OS.Status'
);

meta.surv.cox <- create.surv(
    meta.clinic.patient, 'OS.Months', 'OS.Status'
);


# Create a list to store the results of the Cox models and their summaries
cox.model.results <- list()
gene.names <- colnames(meta.fpkm)

# Loop through each gene and fit univariate Cox models
for (gene in gene.names) {
    # Extract the FPKM values for the current gene (as a vector)
    fpkm.values <- log2(as.numeric(meta.fpkm[, gene]) + 1)
    
    # Fit the Cox proportional hazards model using the 'coxph' function from the 'survival' package
    cox.model <- coxph(meta.surv.cox ~ fpkm.values)
    
    # Store the model and its summary in the cox.model.results list
    cox.model.results[[gene]] <- list(model = cox.model, summary = summary(cox.model))
}


# Initialize empty vectors to store gene names and p-values
gene_names_vector <- c()
p_values_vector <- c()

# Loop through each gene in gene.names and extract the p-value from the summary
for (gene_name in gene.names) {
    # Access the summary object for the specific gene
    summary_for_gene <- cox.model.results[[gene_name]]$summary
    
    # Access the p-value for the gene
    p_value <- summary_for_gene$coefficients["fpkm.values", "Pr(>|z|)"]
    
    # Store the gene name and p-value in the corresponding vectors
    gene_names_vector <- c(gene_names_vector, gene_name)
    p_values_vector <- c(p_values_vector, p_value)
}

# Create a data frame with gene names and p-values
p_values_data <- data.frame(GeneName = gene_names_vector, PValue = p_values_vector)

# Print the data frame
print(p_values_data)
# Assuming p_values_data is a dataframe with columns "Gene" and "P_Value"

# Sort the dataframe in ascending order of p-values
meta_p_values <- p_values_data[order(p_values_data$PValue), ]

# View the sorted dataframe
View(meta_p_values)


meta_fdr_corrected_p_values <- p.adjust(meta_p_values$PValue, method = 'fdr')

meta_p_values$FDR.Corrected <- meta_fdr_corrected_p_values

save(cox.model.results, meta_p_values, file = 'meta.cox.RData')


