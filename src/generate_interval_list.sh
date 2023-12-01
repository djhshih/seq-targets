# /bin/bash

set -euo pipefail

# ! Step 1: Download the pdf file

wget -O ../output/gene_list_meta.pdf https://www.accessdata.fda.gov/cdrh_docs/pdf17/P170019S006C.pdf

# ! Step 2: Clean and generated the gene list
rm ../data/*.interval_list
python extract_table.py


# ! Step 3: Generate interval list

