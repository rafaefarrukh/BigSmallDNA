###################################################################################################
# Script 5
# Clustering sequences using h-clustering
###################################################################################################

# working directory
wd <- "data"

# libraries
library(tidyverse) # general data manipulation and visualization
library(ape) # read, store (DNAbin), and write fasta
library(kmer) # kdistance()
library(dendextend) # sub-clustering

# import ORFs
ORF <- read.FASTA(paste0(wd, "/raw/ORF.fasta"))

# trim labels for ease
## check your data before trimming
## current: IDs usually are 10 char long
labels(ORF) <- substr(labels(ORF), 1, 10)

# trim results table
res <- res[which(res$Accession %in% labels(ORF)),]

# kmean clustering (time consuming)
kmean <- cluster(ORF)

# hierarchical clustering
## set k accordingly
## current: k=5 
cluster <- cutree(as.dendrogram(kmean), k = 5)
View(table(cluster))

# convert to table
cluster <- data.frame(id = names(cluster), cluster = cluster)

# number of clusters (will be used as index)
var <- 1:length(unique(temp.cluster$cluster))

# write
for (i in 1:length(var)) {
    
    print(i) # progress
    
    # write
    write.FASTA(
        x = ORF[filter(cluster, cluster == i)$id],
        file = paste0(wd, "/fasta/cluster.", temp[1], ".fasta"),
    )
    
    # update index
    c.var <- c(c.var, paste0(temp, ".cluster.", i))
    
}

# remove objects and clear cache
rm(list = c("kmean", "cluster"))
gc()

###################################################################################################