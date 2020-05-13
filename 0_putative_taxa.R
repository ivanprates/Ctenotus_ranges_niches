####
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### May 2020.

### The goals of this script are:
### To try to match clusters inferred by sNMF and Ctenotus taxa.
### These assigments will then be used to label plots.

## Clearing working space:
rm(list = ls())

## Number of clusters considered:
minK = 13
maxK = 15

## Creating empty lists to populate with matrices of qscores:
putative_taxa_list <- vector("list", length(minK:maxK))
names(putative_taxa_list) <- c(minK:maxK)

## For K = 13:
putative_taxa_list[[1]] <- c("coggeri", 
                             "arnhemensis",
                             "vertebralis",
                             "essingtonii 1",
                             "essingtonii 2",
                             "decaneurus",
                             "joanae",
                             "robustus 1",
                             "spaldingi",
                             "robustus 2",
                             "superciliaris",
                             "lateralis",
                             "inornatus")

## For K = 14:
putative_taxa_list[[2]] <- c("coggeri", 
                             "arnhemensis",
                             "vertebralis",
                             "essingtonii",
                              "decaneurus",
                             "joanae",
                             "robustus 1",
                             "robustus 2",
                             "spaldingi 1",
                             "spaldingi 2",
                             "robustus 3",
                             "superciliaris",
                              "lateralis",
                             "inornatus")

## For K = 15:
putative_taxa_list[[3]] <- c("coggeri", 
                              "arnhemensis",
                              "vertebralis",
                              "essingtonii 1",
                              "essingtonii 2",
                              "decaneurus",
                              "joanae",
                              "robustus 1",
                              "robustus 2",
                              "spaldingi 1",
                              "spaldingi 2",
                              "robustus 3",
                              "superciliaris 2",
                              "superciliaris 1",
                              "inornatus")

## End of script.
