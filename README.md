# seq-targets

Extract the interval list for variant calling from the gene table of FoundationOneCDx. Specifically, the list of genes are avaiable here:

* FoundationOneCDx (page 4-5): https://www.accessdata.fda.gov/cdrh_docs/pdf17/P170019S006C.pdf
* FoundationOne Liquid CDx (page 3-4): https://www.accessdata.fda.gov/cdrh_docs/pdf19/P190032C.pdf


`gtf` folder contains the code to get gtf file that contains the annotation of genes from EBI. `src` folder contains all source code for interval list extraction. all required dependencies are listed in `./requirements.txt`

# Usage:

All the code for the job are in the `src` folder. The file `generate_interval_list.sh` conduct the full process for interval list generation. Interval list based on hg37 and hg38 will be generated. 



