###################################################################################################
# Script 3
# Download large number of sequences using ncbi datasets and IDs
###################################################################################################

# working directories
wd="data" # default
wdDatasets="Downloads" # containing datasets bin file

# wd
cd $wdDatasets

# assuming 100 files of IDs, so change index accordingly

# download sequences
## check datasets manual for optimal command
## it may fail for some files, in which case rerun command for those files

for i in {1..100}
do
./datasets download genome accession --inputfile $wd/raw/res.$i.txt --filename $wd/raw/res.$i.zip
done


for i in {1..100}
do
unzip $wd/raw/res.$i.zip -d $wd/raw/res.$i
done

# empty file to merge all seqs
> $wd/raw/dataset

# add all seqs to empty file
for i in {1..100}
do 
cat "$wd/raw/res.$i/ncbi_dataset/data/genomic.fna" >> "$wd/raw/dataset"
done

# remove duplicates
seqkit rmdup -s < $wd/raw/dataset > $wd/raw/dataset.fasta

# remove unwanted files
for i in {1..100}
do 
rm $wd/raw/res.$i.zip # downloaded zip files
rm $wd/raw/res.$i # extracted folders
rm $wd/raw/dataset # merged seqs with duplicates
done

###################################################################################################