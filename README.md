# CSU_Capstone
Project data for CSU capstone


## Intersect
This folder contains subsamples of Illumina's bam files from their public dataset. They were subsampled due to the sheer size of the original files. Also, the python script for the intersect tool is available here as well. Calling this script while passing in a bam file and the reference genome will result in desired output of indexed bam files. These indexed files were converted to JSON's so they could be passed through Cromwell's workflow and the output of Cromwell's analysis can be seen in the Comparison folder.

## Comparison
This folder contains the output from the Cromwell analysis. The output was put into tsv files and used for comparisons between masked and unmasked data to ensure no data was lost in the process. The bar plot to visually see any differences is in the jupyter notebook in this folder.

## PCA
Contains the python script for running PCA along with the necessary tsv files that need to be imported to retrieve the metrics.
