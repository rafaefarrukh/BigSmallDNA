# BigSmallDNA
Reduce large nucleotide sequence dataset without sacrificing data of significance by clustering and improving nucleotide diversity of each cluster.

See Overview for more details.

# General Workflow
1. Filter genomes on NCBI using the following criteria and download the csv results table:
2. Extract accession IDs from the results tables.
3. Download genomes using NCBI datasets and remove duplicates using SeqKit.
4. Extract ORFs using EMBOSS and remove duplicates using SeqKit.
5. Cluster the sequences using h-clustering.
6. Align clusters roughly using MAFFT and remove duplicates using SeqKit.
7. Improve genetic diversity for each lineage.

# Scripts
CompletePipeline is a Rmd file with code chunks so that the entire pipeline can be run from a single file. But since this may not be desirable depending on the usecase, each step is provided as a seperate script labelled as mentioned in the general workflow above.

# Example
ExamplePipeline is made using Severe Acute Respiratory Syndrome Coronavirus-2 (SARS-CoV) data available on NCBI. It reduced the dataset size by 143,190% and improved average nucleotide diversity per cluster by 26%.

|Original Dataset|Downloaded from NCBI|Unique ORFs|Reduced Dataset|Reduction (%)|
|:--:|:----:|:----:|:----:|:----:|
|9,064,523|1,365,419|181,940|6,326|143,190|

|Dataset|Average nucleotide diversity per cluster|Minimum nucleotide diversity per cluster|Maximum nucleotide diversity per cluster|
|:--:|:--:|:--:|:--:|
|Original|0.0012775|0.0000000|0.4823382|
|Reduced|0.0082952|0.0000000|0.6081363|
|Change (%)|549.3|0.0|26.0|
