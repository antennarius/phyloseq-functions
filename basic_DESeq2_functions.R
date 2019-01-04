library(phyloseq)
library(DESeq2)
library(ggplot2)

sigdds <- function(physeq, by, test, fitType, sig.level = .05) {
  physeqdds <- phyloseq_to_deseq2(physeq, as.formula(paste(" ~", by)))
  physeqdds <- DESeq(physeqdds, test = test, fitType = fitType)
  res <- results(physeqdds, cooksCutoff = FALSE)
  alpha = sig.level
  sigtab <- res[which(res$padj < alpha), ]
  sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
  sigtab
}

resdds <- function(physeq, by, test, fitType, sig.level = .05) {
  physeqdds <- phyloseq_to_deseq2(physeq, as.formula(paste(" ~", by)))
  physeqdds <- DESeq(physeqdds, test = test, fitType = fitType)
  results(physeqdds, cooksCutoff = FALSE)
}

plotdds <- function(table) {
  x = tapply(table$log2FoldChange, table$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  table$Phylum = factor(as.character(table$Phylum), levels=names(x))
  # Genus order
  x = tapply(table$log2FoldChange, table$Genus, function(x) max(x))
  x = sort(x, TRUE)
  table$Genus = factor(as.character(table$Genus), levels=names(x))
  ggplot(table, aes(x = Genus, y = log2FoldChange, color = Phylum, fill = Phylum)) + 
    geom_bar(stat = "identity") +
    
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
    scale_color_manual(values = cbPalette) +
    scale_fill_manual(values = cbPalette) +
    coord_flip()
}
