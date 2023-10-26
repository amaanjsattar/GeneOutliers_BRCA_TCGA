library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)
library(BoutrosLab.plotting.survival)
load('meta.cox.RData')

cox.gene.names <- names(cox.model.results)

result_df <- data.frame(GeneName = character(0), PValue = numeric(0), HR = numeric(0))

for (gene_name in cox.gene.names) {
    # Extract p-value from cox_model_results
    p_value <- cox.model.results[[gene_name]]$summary$coefficients["fpkm.values", "Pr(>|z|)"]
    
    # Extract hazard ratio from cox_model_results
    hr <- exp(coef(cox.model.results[[gene_name]]$model)["fpkm.values"])
    
    # Create a new row for the gene in the result dataframe
    new_row <- data.frame(GeneName = gene_name, PValue = p_value, HR = hr)
    
    # Append the new row to the result dataframe
    result_df <- rbind(result_df, new_row)
}

# SCATTER PLOTS: x = LOG2(HR) , y = -LOG10(PVAL)

colors <- ifelse(-log10(result_df$PValue) < -log10(0.05), "gray",
                 ifelse(-log10(result_df$PValue) >= -log10(0.05) & log2(result_df$HR) < 0, "blue", 
                        "red"))
result_df$Color <- colors
create.scatterplot(
    height = 15,
    width = 15,
    formula = -log10(PValue) ~ log2(HR),
    data = result_df,
    main = 'METABRIC: Univariate Cox Models',
    main.cex = 1.2,
    ylab.label = expression(-log[10]("P-Value")),
    xlab.label = expression(log[2]("Hazard Ratio")),
    ylab.cex = 1,
    xlab.cex = 1,
    xaxis.cex = 1,
    yaxis.cex = 1,
    alpha = 0.5,
    col = colors,
    abline.h = -log10(0.05),
    abline.col = 'red',
    abline.lwd = 1.2,
    top.padding = 1
    
)

result_df$FDR.Corrected <- p.adjust(result_df$PValue, method = 'fdr')
colors_meta_fdr <- ifelse(-log10(result_df$FDR.Corrected) < -log10(0.05), "gray",
                 ifelse(-log10(result_df$FDR.Corrected) >= -log10(0.05) & log2(result_df$HR) < 0, "blue", 
                        "red"))
create.scatterplot(
    height = 15,
    width = 15,
    formula = -log10(FDR.Corrected) ~ log2(HR),
    data = result_df,
    main = 'METABRIC: Univariate Cox Models',
    main.cex = 1.2,
    ylab.label = expression(-log[10]("FDR-Adjusted P-Value")),
    xlab.label = expression(log[2]("Hazard Ratio")),
    ylab.cex = 1,
    xlab.cex = 1,
    xaxis.cex = 1,
    yaxis.cex = 1,
    alpha = 0.5,
    col = colors_meta_fdr,
    top.padding = 1,
    xaxis.tck = c(1, 0),
    yaxis.tck = c(1, 0)
    
)



create.histogram(x = result_df$PValue,
                 data = result_df,
                 main = 'METABRIC Univariate Cox Models, Unadjusted P-Values',
                 main.cex = 1.2,
                 xlimits = c(0.0, 1.0),
                 breaks = seq(0, 1, 0.05),
                 xat = seq(0, 1, 0.1),
                 xaxis.cex = 0.9,
                 yat = seq(0, 40, 5),
                 yaxis.cex = 0.9,
                 ylab.label = 'Frequency',
                 xlab.cex = 1.25,
                 ylab.cex = 1.25,
                 xlab.label = 'P-Value',
                 ylimits = c(0, 35),
                 col = 'skyblue2',
                 abline.v = 0.05,
                 abline.col = 'red',
                 abline.lwd = 2
)



#TCGA SCATTER PLOTS
load('tcga.cox.RData')
tcga.gene.names <- names(tcga.cox.results)
tcga_result_df <- data.frame(GeneName = character(0), PValue = numeric(0), HR = numeric(0))


for (gene_name in tcga.gene.names) {
    # Extract p-value from cox_model_results
    p_value <- tcga.cox.results[[gene_name]]$summary$coefficients["fpkm.values", "Pr(>|z|)"]
    
    # Extract hazard ratio from cox_model_results
    hr <- exp(coef(tcga.cox.results[[gene_name]]$model)["fpkm.values"])
    
    # Create a new row for the gene in the result dataframe
    new_row <- data.frame(GeneName = gene_name, PValue = p_value, HR = hr)
    
    # Append the new row to the result dataframe
    tcga_result_df <- rbind(tcga_result_df, new_row)
}

# SCATTER PLOTS: x = LOG2(HR) , y = -LOG10(PVAL)
colors_tcga <- ifelse(-log10(tcga_result_df$PValue) < -log10(0.05), "gray",
                          ifelse(-log10(tcga_result_df$PValue) >= -log10(0.05) & log2(tcga_result_df$HR) < 0, "blue", 
                                 "red"))
create.scatterplot(
    height = 15,
    width = 15,
    formula = -log10(PValue) ~ log2(HR),
    data = tcga_result_df,
    main = 'TCGA: Univariate Cox Models',
    main.cex = 1.2,
    ylab.label = expression(-log[10]("P-Value")),
    xlab.label = expression(log[2]("Hazard Ratio")),
    ylab.cex = 1,
    xlab.cex = 1,
    xaxis.cex = 1,
    yaxis.cex = 1,
    alpha = 0.5,
    col = colors_tcga,
    abline.h = -log10(0.05),
    abline.col = 'red',
    abline.lwd = 1.2,
    top.padding = 1,
    xaxis.tck = c(1, 0),
    yaxis.tck = c(1, 0)
    
)

tcga_result_df$FDR.Corrected <- p.adjust(tcga_result_df$PValue, method = 'fdr')
create.scatterplot(
    height = 15,
    width = 15,
    formula = -log10(FDR.Corrected) ~ log2(HR),
    data = tcga_result_df,
    main = 'TCGA: Univariate Cox Models',
    main.cex = 1.2,
    ylab.label = expression(-log[10]("FDR-Adjusted P-Value")),
    xlab.label = expression(log[2]("Hazard Ratio")),
    ylab.cex = 1,
    xlab.cex = 1,
    xaxis.cex = 1,
    yaxis.cex = 1,
    alpha = 0.5,
    col = 'gray',
    top.padding = 1,
    xaxis.tck = c(1, 0),
    yaxis.tck = c(1, 0)
    
)

create.scatterplot(
    height = 10,
    width = 10,
    formula = -log10(FDR.Corrected) ~ log2(HR),
    data = tcga_result_df,
    main = 'TCGA: Univariate Cox Models',
    main.cex = 1.2,
    ylab.label = expression(-log[10]("FDR-Adjusted P-Value")),
    xlab.label = expression(log[2]("Hazard Ratio")),
    ylab.cex = 1,
    xlab.cex = 1,
    xaxis.cex = 1,
    yaxis.cex = 1,
    alpha = 0.5,
    col = 'gray',
    top.padding = 1,
    xaxis.tck = c(0, 1),
    yaxis.tck = c(0, 1)
    
)
create.histogram(x = tcga_result_df$PValue,
                 data = tcga_result_df,
                 main = 'TCGA Univariate Cox Models, P-Values',
                 main.cex = 1.2,
                 xlimits = c(0.0, 1.0),
                 breaks = seq(0, 1, 0.05),
                 xat = seq(0, 1, 0.1),
                 xaxis.cex = 0.9,
                 yaxis.cex = 0.9,
                 ylab.label = 'Frequency',
                 xlab.cex = 1.25,
                 ylab.cex = 1.25,
                 xlab.label = 'P-Value',
                 col = 'skyblue2',
                 abline.v = 0.05,
                 abline.col = 'red',
                 abline.lwd = 2
)
colnames(meta_p_values) <- c('GeneName', 'Meta.P.Value', 'Meta.FDR')
colnames(tcga_p_values) <- c('GeneName', 'TCGA.P.Value', 'TCGA.FDR')

# Get the common 'GeneName' values
common_genes <- intersect(result_df$GeneName, tcga_result_df$GeneName)

# Subset both dataframes to keep only the rows with common 'GeneName'
cumulative_result_df <- result_df[result_df$GeneName %in% common_genes, ]
cumulative_result_df <- merge(cumulative_result_df, tcga_result_df[tcga_result_df$GeneName %in% common_genes, ], by = 'GeneName')

colnames(cumulative_result_df) <- c('GeneName', 'Meta.P', 'Meta.HR', 'Meta.Color', 'Meta.FDR',
                                    'TCGA.P', 'TCGA.HR', 'TCGA.FDR')
create.scatterplot(
    height = 15,
    width = 15,
    formula = TCGA.P ~ Meta.P,
    data = cumulative_result_df,
    main = 'METABRIC-TCGA Comparison: P-Values',
    ylab.label = expression('TCGA P-Value'),
    xlab.label = expression('METABRIC P-Value'),
    main.cex = 1.2,
    ylab.cex = 1,
    xlab.cex = 1,
    xaxis.cex = 1,
    yaxis.cex = 1,
    alpha = 0.5,
    col = 'black',
    top.padding = 1,
    add.xyline = TRUE,
    xaxis.tck = c(1, 0),
    yaxis.tck = c(1, 0),
    
    
)

create.scatterplot(
    height = 15,
    width = 15,
    formula = log(TCGA.HR) ~ log(Meta.HR),
    data = cumulative_result_df,
    main = 'METABRIC-TCGA Comparison: Hazard Ratios',
    ylab.label = expression(ln("TCGA Hazard Ratio")),
    xlab.label = expression(ln('METABRIC Hazard Ratio')),
    main.cex = 1.2,
    ylab.cex = 1,
    xlab.cex = 1,
    xaxis.cex = 1,
    yaxis.cex = 1,
    alpha = 0.5,
    col = 'black',
    add.xyline = TRUE,
    top.padding = 1,
    xaxis.tck = c(1, 0),
    yaxis.tck = c(1, 0)
    
)

