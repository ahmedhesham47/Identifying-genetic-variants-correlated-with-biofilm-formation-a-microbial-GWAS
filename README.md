# Identifying-genetic-variants-correlated-with-biofilm-formation-a-microbial-GWAS
This project aims to identify genetic variants correlated with biofilm formation using treeWAS R package.

## What is treeWAS?
TreeWAS is an R package that features a phylogenetic, pangenome-based approach (though supporting other variants’ types).  It implements simulated “null” genetic datasets to control for multiple confounders. These simulated loci make a null distribution of three different association scores (terminal, subsequent, and simultaneous), which are then used as a cutoff to identify the real (unsimulated; experimental) significant associations [6].

## What are its inputs?
First off, you need to install devtools, treeWAS, and rlist packages.
The package takes two required files as inputs: a gene variants matrix and a phenotype vector (discrete or continuous). In the gene variants matrix, which can be SNP data, gene presence/absence, INDELs, or anything else, the variants must be in columns and samples, which must have matching names with the phenotype file, must be in rows. Each variant value should be either a 1 or 0, indicating the presence of the variant or its absence. Moreover, the package takes an optional, though recommended, phylogenetic tree, based on which it calculates the association scores.
All data cleaning and preparation must be within R, according to the source of data and how it looks like. Typically, you would have to ensure that the names are matching, no nulls are found, remove unnecessary columns/rows, remove hypothetical-protein-coding genes, and so on, depending on the nature of the variant file. Remember to convert your variants file to a matrix using as.matrix() function and your phenotypes file to a vector using as.vector().
After finalizing all the preparations, you can call the treeWAS function to run the association tests based on the parameters you set. An example is shown below. Out is an object that stores all the analysis results and can be used to retrieve whichever association score of the three. “Plots.pdf” is a pdf file called “Plots” that will contain all the plots.

![image](https://github.com/ahmedhesham47/Identifying-genetic-variants-correlated-with-biofilm-formation-a-microbial-GWAS/assets/44484663/84850816-fb97-4fab-bdd2-e28f58ac80eb)


## What are its outputs?
The output of treeWAS is an object of type “treeWAS”, called “out” in the picture above. This object stores all the statistical analyses of the GWAS as well as all the necessary plots that are needed for interpretation. If you put a dollar sign after the word “out”, as follows “out$”, you will find five options that you can choose from. These different options give you different outputs and results. Most importantly, the plots file will be saved in the same directory of the R file.

### Sample output
After running treeWAS on our sample dataset, the following three plots were generated, one for each association score. Only the simultaneous association score shows significant genes.

![image](https://github.com/ahmedhesham47/Identifying-genetic-variants-correlated-with-biofilm-formation-a-microbial-GWAS/assets/44484663/89d0e29b-c186-4751-afdd-f25d27c31e9e)
![image](https://github.com/ahmedhesham47/Identifying-genetic-variants-correlated-with-biofilm-formation-a-microbial-GWAS/assets/44484663/324465ba-4f5f-4b3c-9f8a-74c41df4e987)
![image](https://github.com/ahmedhesham47/Identifying-genetic-variants-correlated-with-biofilm-formation-a-microbial-GWAS/assets/44484663/8a89c9a1-e802-47c5-9caa-8bee949ddd86)

Moreover, the table below can be retrieved from the treeWAS object and saved to a .csv file. It illustrates the p-adjusted value of each gene, its association score test, and its abundance in each group of the four, G0P0, G0P1, G1P0, G1P1, where G0P0 means the gene was not found and the phenotype was non-biofilm former, and so on for the rest.

![image](https://github.com/ahmedhesham47/Identifying-genetic-variants-correlated-with-biofilm-formation-a-microbial-GWAS/assets/44484663/7a55abf5-712d-49f3-b86a-5b79d57976cf)

## Author:
Ahmed H. Saadawy
