###################################################################################################
# Script 6
# Rough MSA using loose parameters
###################################################################################################

# working directories
wd="data" # default

# wd
cd $wd

# MSA with loose parameters
for i in "fasta"/*".fasta"
do
mafft --retree 1 --thread 12 $i > msa/$(basename $i .fasta).msa.fasta
done

###################################################################################################