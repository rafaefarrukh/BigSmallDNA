###################################################################################################
# Script 2
# Extract IDs from csv results table and write as batches of files
###################################################################################################

# wd
wd <- "data"

# read csv results table
res <- read.csv(paste0(wd, "/raw/res.csv"))

# make index for writing IDS
# current: 10000IDs per file
temp <- c(seq(1, nrow(res), 10000), nrow(res))

# write IDs
for (i in 1:length(temp)) {
    write.table(
        res$Accession[temp[i-1]:temp[i]],
        paste0(wd, "/raw/ID.", i-1, ".txt"),
        row.names = FALSE, col.names = FALSE, quote = FALSE
    )
}

###################################################################################################