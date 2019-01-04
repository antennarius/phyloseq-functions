correlate_2_genera <- function(physeq, taxon1, taxon2, logtrans = TRUE) {
  row1 <- which(tax_table(physeq)[,'Genus'] == taxon1)
  row2 <- which(tax_table(physeq)[,'Genus'] == taxon2)
  names1 <- taxa_names(tax_table(physeq))
  names2 <- taxa_names(tax_table(physeq))
  taxon1pseq <- prune_taxa(names1[row1], physeq)
  taxon2pseq <- prune_taxa(names2[row2], physeq)
  sum1 <- sample_sums(taxon1pseq)
  sum2 <- sample_sums(taxon2pseq)
  if(logtrans == TRUE) {
    sum1 <- log(sum1 + 1)
    sum2 <- log(sum2 + 2)
  }
  model <- lm(sum1 ~ sum2)
  summary(model)
}