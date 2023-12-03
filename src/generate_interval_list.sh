# /bin/bash

set -euo pipefail

# ! Step 1: Download the pdf file
if [[ ! -f "../output/gene_list_meta.pdf" ]]; then
    wget -O ../output/gene_list_meta.pdf https://www.accessdata.fda.gov/cdrh_docs/pdf17/P170019S006C.pdf
fi

# ! Step 2: Clean and generated the gene list + code for 
python extract_table.py

# ! Step 3: Generate raw interval list
(cd ../output && bash ../output/subset_gff3.sh)

# ! Step 4: Refine the interval list
python refine_interval.py

# ! Step 5: Merge all things together and sort bed file:
cd ../output
cat *.interval_list > target_regions.bed
bedtools sort -i target_regions.bed > target_regions.sorted.bed
#rm *.interval_list
#rm *.interval_list.raw