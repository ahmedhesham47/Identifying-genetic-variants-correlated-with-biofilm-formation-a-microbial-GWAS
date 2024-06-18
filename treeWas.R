# Written by Ahmed Hesham Saadawy

library(devtools)
library(treeWAS)
library(rlist)

## Load CA_MRSA data:
## import files manually by 'Import Dataset' feature

# Data Preparation Steps

names(gene_presence_absence) = gsub(pattern = "_polish*", replacement = "", x = names(gene_presence_absence))
names(gene_presence_absence) = gsub(pattern = "_scaffolds*", replacement = "", x = names(gene_presence_absence))
names(gene_presence_absence) = gsub(pattern = "_pilosh*", replacement = "", x = names(gene_presence_absence))
names(gene_presence_absence) = gsub(pattern = "MRSA", replacement = "", x = names(gene_presence_absence))
# To be run twice
rownames(biofilm_phenotypes) = gsub(pattern = "_0", replacement = "_", x = rownames(biofilm_phenotypes))

# Saving the files in case we want them for further usage
write.csv(gene_presence_absence, "gene_presence_absence_updated.csv")
write.csv(biofilm_phenotypes, "Phenotypes_updated.csv")

# After loading the new, updated files, we do data cleaning by selecting only the columns we are interesting in
# And by transposing the gene presence absence dataframe (treeWAS format)

Phenotypes_updated = subset(Phenotypes_updated, select = X)

## Convert phenotypes to a vector:
traits <- as.vector(unlist(Phenotypes_updated))
names(traits) <- rownames(Phenotypes_updated)

# Removing unneeded columns
gene_presence_absence_updated = gene_presence_absence_updated[-c(3:13)]
gene_presence_absence_updated = gene_presence_absence_updated[-1]

# Converting gene presence absence dataframe to a matrix and transposing it
gene_presence_absence_updated <- as.matrix(gene_presence_absence_updated)

gene_presence_absence_updated = t(gene_presence_absence_updated)


# Removing hypothetical proteins
vector = c()

for (j in 1:ncol(gene_presence_absence_updated)) {
  if (gene_presence_absence_updated[1, j] == "hypothetical protein" | gene_presence_absence_updated[1, j] == "Hypothetical protein") {
    vector = append(vector, j, after=length(vector))
  }
}

gene_presence_absence_updated = gene_presence_absence_updated[,-vector]

# Convert gene presence absence matrix to binary (1s and 0s)
for (i in 1:nrow(gene_presence_absence_updated)) {
  # loop over columns
  for (j in 1:ncol(gene_presence_absence_updated)) {
    # check if value is not null
    if (startsWith(gene_presence_absence_updated[i, j], "CA")) {
      # replace value with 1
      gene_presence_absence_updated[i, j] <- 1
    } else {
      # replace value with 0
      gene_presence_absence_updated[i, j] <- 0
    }
  }
}

# Re-naming the phylogenetic tree tips to be matching with the SNP matrix
tree <- read.tree(file = "phylo_tree.newick")
tree$tip.label = gsub(pattern = "_polish*", replacement = "", x = tree$tip.label)
tree$tip.label = gsub(pattern = "_scaffolds*", replacement = "", x = tree$tip.label)
tree$tip.label = gsub(pattern = "_pilosh*", replacement = "", x = tree$tip.label)
tree$tip.label = gsub(pattern = "MRSA", replacement = "", x = tree$tip.label)

# Checking no valeus are null and cross-checking labels
is.null(tree$tip.label)
is.null(rownames(gene_presence_absence_updated))
is.null(names(traits))

all(tree$tip.label %in% rownames(gene_presence_absence_updated))
all(rownames(gene_presence_absence_updated) %in% tree$tip.label)
all(tree$tip.label %in% names(traits))
all(names(traits) %in% tree$tip.label)
all(names(traits) %in% rownames(gene_presence_absence_updated))
all(rownames(gene_presence_absence_updated) %in% names(traits))


# The main output variable that contains all the results
# Adjusting the parameters here
out <- treeWAS(snps = gene_presence_absence_updated,
               phen = traits,
               tree = tree,
               filename.plot = "Plots.pdf",
               plot.dist = TRUE,
               p.value = 0.05,
               correct.prop = TRUE,
)

#write.csv(out$simultaneous$p.vals, "p-vals_simultaneous.csv")

write.csv(out2$simultaneous$sig.snps, "significant_SNPs_simultaneous_without_hypotheticals.csv")

#write.csv(out$treeWAS.combined, "Output.csv")
