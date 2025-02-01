###################################################################################################
# Script 4
# Extract ORFs from dataset using EMBOSS getORF
###################################################################################################

# working directories
wd="data" # default

# wd
cd $wd/raw

# check EMBOSS manual for optimal command
## current: all ORFs are extracted (not ideal)
getorf -sequence dataset.fasta -outseq ORF

# remove duplicates
seqkit rmdup -s < ORF > ORF.fasta

# remove file with duplicates
rm ORFs

###################################################################################################