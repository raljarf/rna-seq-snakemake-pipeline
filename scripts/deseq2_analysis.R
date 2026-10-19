library(DESeq2)

# Read the gene count matrix
counts <- read.delim(
    "counts/gene_count_matrix.tsv",
    row.names = 1,
    check.names = FALSE
)

# Read sample metadata
samples <- read.delim(
    "samples.tsv"
)

# Reorder count matrix columns to match sample metadata
counts <- counts[, samples$sample]

# Verify that the sample order matches
stopifnot(identical(colnames(counts), samples$sample))

# Set control as the reference condition
samples$condition <- factor(
    samples$condition,
    levels = c("control", "treatment")
)

# Create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = samples,
    design = ~ condition
)


# Display the DESeq2 dataset
print(dds)

# Display the sample information
print(colData(dds))

# Run differential expression analysis
dds <- DESeq(dds)

res <- results(
    dds,
    contrast = c("condition", "treatment", "control")
)

summary(res)

# Sort genes by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Display the first 10 genes
head(res_ordered, 10)

# Create a directory for differential expression results
dir.create("deseq2_results", showWarnings = FALSE)

# Save the results, including gene IDs
write.table(
    as.data.frame(res_ordered),
    file = "deseq2_results/differential_expression.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = NA
)


# Apply variance stabilizing transformation
vsd <- vst(dds, blind = FALSE, nsub = 1000)

# Get PCA coordinates
pca_data <- plotPCA(
    vsd,
    intgroup = "condition",
    returnData = TRUE
)

# Extract percentage of variance explained
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Load ggplot2
library(ggplot2)

# Create and save the PCA plot
pdf("deseq2_results/PCA_plot.pdf")

print(
    ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 3) +
        geom_text(aes(label = name), vjust = -1, size = 3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle("PCA: Control vs Treatment") +
        theme_bw()
)

# Close the PDF file
dev.off()

# Create a volcano plot
volcano_data <- as.data.frame(res)

volcano_data$significance <- "Not significant"

volcano_data$significance[
    !is.na(volcano_data$padj) &
    volcano_data$padj < 0.1 &
    volcano_data$log2FoldChange > 0
] <- "Upregulated"

volcano_data$significance[
    !is.na(volcano_data$padj) &
    volcano_data$padj < 0.1 &
    volcano_data$log2FoldChange < 0
] <- "Downregulated"

pdf("deseq2_results/volcano_plot.pdf")

print(
    ggplot(
        volcano_data,
        aes(
            x = log2FoldChange,
            y = -log10(padj),
            color = significance
        )
    ) +
        geom_point(alpha = 0.6, na.rm = TRUE) +
        scale_color_manual(
            values = c(
                "Upregulated" = "red",
                "Downregulated" = "blue",
                "Not significant" = "grey"
            )
        ) +
        xlab("Log2 fold change (Treatment vs Control)") +
        ylab("-Log10 adjusted p-value") +
        ggtitle("Volcano Plot: Control vs Treatment") +
        theme_bw()
)

dev.off()
