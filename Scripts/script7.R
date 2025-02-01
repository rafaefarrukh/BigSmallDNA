###################################################################################################
# Script 7
# Improving genetic diversity per cluster
###################################################################################################

# working directory
wd <- "data"

# libraries
library(tidyverse) # general data manipulation and visualization
library(ape) # read, store (DNAbin), and write fasta
library(pegas) # nuc.div()

# read MSA
for (i in 1:length(var)) {
    assign(
        paste0("cluster.", var[i]),
        read.FASTA(paste0(wd, "/msa/cluster.", var[i], ".msa.fasta"))
    )
}

# dataframe to summarize genetic diversity improvements
nd <- data.frame(
    cluster = paste0("cluster.", var), # cluster
    old.nd = NA, new.nd = NA, change.nd = NA, # old, new, and change in nuc div
    old.seq = NA, new.seq = NA, change.seq = NA, # old, new, and change in num of seq
    th = NA, # threshold
    time = NA # time taken
)

# improve diversity
for (i in 1:length(.var)) { # check index before running
    
    print(i) # progress 
    
    print(Sys.time()) # time of start
    
    try({ # ignore errors
        
        time.old <- Sys.time() # time for time complexity
        
        # get cluster
        temp <- get(paste0("cluster.", c.var[i])) # run for i>=3
        
        # convert to distance matrix
        ## choose model according to data
        ## current: K80 (default)
        dist <- dist.dna(temp, model = "K80", as.matrix = TRUE)
        # convert diagonal entries to NA for ease
        diag(dist) <- NA
        
        # temporary dataframe to evaluate optimal threshold
        ## chose threshold range according to data
        ## current: 0.9 - 1
        ## if unsure use c(0.5, 0.6, 0.7, 0.8, 0.9, 1) to get an idea
        threshold <- data.frame(th = c(0.9, 0.95, 0.99, 0.999, 0.9999, 1), nd = NA)
        
        # nd of temp (calc once to prevent recalculations)
        temp.nd <- nuc.div(temp)
        
        # evaluate optimal threshold
        # new = nd of seqs with identity >= threshold * max score
        # nd = new - old
        for (j in 1:nrow(threshold)) {
            threshold$nd[j] <-
                nuc.div(
                    temp[
                        unique(c(which(
                            dist >= threshold$th[j]*max(dist, na.rm = TRUE),
                            arr.ind = TRUE)))
                    ]
                ) 
            - temp.nd
        }
        
        # determine optimal threshold and save seqs
        assign(
            paste0("cluster", var[i], ".i"),
            temp[
                unique(c(which(
                    dist >= threshold$th[
                        order(threshold$nd, decreasing = TRUE)[1]
                    ]*max(dist, na.rm = TRUE),
                    arr.ind = TRUE)))
            ]
        )
        
        # nd
        nd$old.nd[i] <- temp.nd
        nd$new.nd[i] <- max(threshold$nd) + temp.nd
        nd$change.nd[i] <- (nd$new.nd[i] - nd$old.nd[i]) / nd$old.nd[i]
        nd$th[i] <- threshold$th[order(threshold$nd, decreasing = TRUE)[1]]
        nd$time[i] <- as.numeric(Sys.time() - time.old)
        nd$old.seq[i] <- length(temp)
        nd$new.seq[i] <- length(get(paste0(var[i], ".i")))
        nd$change.seq[i] <- (nd$new.seq[i] - nd$old.seq[i]) / nd$old.seq[i]
        
    })
    
}

# check if nd was decreased for any cluster
View(nd[which(nd$change.nd <= 0),])

# For negative cases, lower the threshold range and repeat the loop above

# merge clusters
temp <- c()
for (i in 1:length(var)) {
    temp <- c(temp, get(paste0("cluster", var[i], ".i")))
}

# trim results table
res <- res[res$Accession %in% labels(temp),]

# get un-aligned seqs
dataset <- c(ORF[labels(ORF) %in% res$Accession])

# write
write.FASTA(c12, file = paste0(wd, "/fasta/dataset.fasta"))

# remove objects and clear cache
rm(list=c("temp", paste0("cluster.", var), paste0("cluster", var[i], ".i"), "var"))
gc()

###################################################################################################