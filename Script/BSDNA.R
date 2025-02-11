###################################################################################################

# BigSmallDNA - Improving Genetic Diversity

###################################################################################################

BSDNA <- function(data, change = T, change.name = "changes", dist.model = "K80", t.coeff = c(0.5, 0.75, 1.0), prefix = data, suffix = "reduced") {
    
    # parameters
    # data         # vector containing name(s) of DNAbin objects
    # change        # if true, then it makes a dataframe which tracks changes in dataset
    # change.name   # name of dataframe tracking changes
    # t.coeff       # vector with threshold coefficients
    # prefix        # prefix for name of reduced dataset
    # suffix        # suffix for name of reduced dataset
    
    # checks
    require(ape); require(pegas)
    if(!(is.vector(data))) {stop("data must be the names of the DNAobjects", .call = F)}
    if(sum(unlist(lapply(lapply(data, get), class)) == "DNAbin") != length(data)) {stop("the data must be DNAbin object")}
    if(!(dist.model %in% c("raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", "indelblock"))) {stop("dist.model does not exist", call. = F)}
    if(class(t.coeff) != "numeric") {stop("t.coeff is not numeric", call. = F)}
    if(sum(1 < range(t.coeff) | range(t.coeff) <= 0) > 0) {stop("range of t.coeff must be betweeen 0 and 1", call. = F)}
    if(sum(length(prefix) == length(data) | length(suffix) == length(data)) == 0) {stop("atleast one of prefix or suffix must have same length as data", .call = F)}
    
    # if change = T, then construct a dataframe for logging changes.
    if (change == T) {
        nd <- data.frame(
            cluster = data, # cluster name
            old.nd = NA, new.nd = NA, change.nd = NA, # old, new, and change in nuc div
            old.seq = NA, new.seq = NA, change.seq = NA, # old, new, and change in num of seq
            threshold = NA, # threshold
            time = NA, # time taken
            length = NA # length of sequence
        )
    }
    
    # if multiple datasets, then for loop is required
    for (i in 1:length(data)) {
        
        # get DNAbin object
        temp <- get(data[i])
        
        # progress report
        print(paste(i, "out of", length(data), "with", length(temp), "sequences")); print(Sys.time())
        
        try({ # ignore errors
            
            # if change = T, then track time
            if (change == T) {time.old <- Sys.time()}
            
            # convert to distance matrix
            dist <- dist.dna(temp, model = dist.model, as.matrix = TRUE)
            
            # convert diagonal to NA for ease
            diag(dist) <- NA
            
            # dataframe to evaluate optimal threshold
            threshold <- data.frame(
                t = t.coeff * max(dist, na.rm = TRUE), # threshold
                ind = NA, # sequences used
                nd = NA # nd
                )
            
            # evaluate nd
            for (j in 1:nrow(threshold)) {
                threshold$ind[j] <- list(unique(c(which(dist >= threshold$t[j], arr.ind = TRUE))))
                threshold$nd[j] <- nuc.div(temp[unlist(threshold$ind[j])])
                }
            
            # get sequences of improved cluster
            if (length(suffix) == 1) {
                assign(paste0(prefix[i], ".", suffix),
                       temp[unlist(threshold$ind[order(threshold$nd, decreasing = TRUE)[1]])], envir = .GlobalEnv
                       )}
            if (length(prefix) == 1) {
                assign(paste0(prefix, ".",
                              suffix[i]), temp[unlist(threshold$ind[order(threshold$nd, decreasing = TRUE)[1]])], envir = .GlobalEnv
                       )}
            
            # if change = T, then track changes
            if (change == T){
                nd$time[i] <- as.numeric(Sys.time() - time.old)
                nd$old.nd[i] <- nuc.div(temp)
                nd$new.nd[i] <- nuc.div(get(paste0(data[i], ".reduced")))
                nd$change.nd[i] <- ((nd$new.nd[i] - nd$old.nd[i]) / nd$old.nd[i])
                nd$old.seq[i] <- length(temp)
                nd$new.seq[i] <- length(get(paste0(data[i], ".reduced")))
                nd$change.seq[i] <- ((nd$old.seq[i] - nd$new.seq[i]) / nd$old.seq[i])
                nd$threshold[i] <- threshold$t[order(threshold$nd, decreasing = TRUE)[1]]
                nd$length[i] <- length(temp[[1]])
                }
            
            })
        
        }
    
    # if changes = T, then save changes
    if (change == T) {assign(change.name, nd, envir = .GlobalEnv)}
    
}

###################################################################################################
