# /bin/bash

set -euo pipefail

gff38="../gtf/gencode.v44.annotation.gff3"
gff37="../gtf/gencode.v44lift37.annotation.gff3"

# ! Step 1: Download the pdf file
if [[ ! -f "../interval_list/gene_list_meta_cdx.pdf" ]]; then
    wget -O ../interval_list/gene_list_meta_cdx.pdf https://www.accessdata.fda.gov/cdrh_docs/pdf17/P170019S006C.pdf
fi

if [[ ! -f "../interval_list/gene_list_meta_liquid_cdx.pdf" ]]; then
    wget -O ../interval_list/gene_list_meta_liquid_cdx.pdf https://www.accessdata.fda.gov/cdrh_docs/pdf19/P190032C.pdf
fi

# ! Step 2: Clean and generated the gene list + code for 
python extract_table.py --target_pdf ../interval_list/gene_list_meta_cdx.pdf --output ../interval_list/cdx_target.tsv --target_page 4 5
python extract_table.py --target_pdf ../interval_list/gene_list_meta_liquid_cdx.pdf --output ../interval_list/liquid_cdx_target.tsv --target_page 3 4

# ! Step 3: Generate raw interval list
rawpath="../interval_list/raw"
mkdir -p ${rawpath}
python subset_gff.py --gff ${gff38} --output_path ${rawpath} --output_prefix cdx_hg38 --target_table ../interval_list/cdx_target.tsv
python subset_gff.py --gff ${gff37} --output_path ${rawpath} --output_prefix cdx_hg37 --target_table ../interval_list/cdx_target.tsv
python subset_gff.py --gff ${gff38} --output_path ${rawpath} --output_prefix liquid_cdx_hg38 --target_table ../interval_list/liquid_cdx_target.tsv
python subset_gff.py --gff ${gff37} --output_path ${rawpath} --output_prefix liquid_cdx_hg37 --target_table ../interval_list/liquid_cdx_target.tsv

# ! Step 4: Refine the interval list
python refine_interval.py --input_path ${rawpath} --input_prefix cdx_hg38 --target_table ../interval_list/cdx_target.tsv --output_path ../interval_list
python refine_interval.py --input_path ${rawpath} --input_prefix cdx_hg37 --target_table ../interval_list/cdx_target.tsv --output_path ../interval_list
python refine_interval.py --input_path ${rawpath} --input_prefix liquid_cdx_hg38 --target_table ../interval_list/liquid_cdx_target.tsv --output_path ../interval_list
python refine_interval.py --input_path ${rawpath} --input_prefix liquid_cdx_hg37 --target_table ../interval_list/liquid_cdx_target.tsv --output_path ../interval_list

# ! Step 5: Merge all things + sort & merge bed file:
cd ../interval_list
cat cdx_hg38.*.interval_list > cdx_hg38.bed
cat cdx_hg37.*.interval_list > cdx_hg37.bed
cat liquid_cdx_hg38.*.interval_list > liquid_cdx_hg38.bed
cat liquid_cdx_hg37.*.interval_list > liquid_cdx_hg37.bed

rm *.interval_list

bedtools sort -i cdx_hg38.bed | bedtools merge -i - > cdx_hg38.sorted.bed
bedtools sort -i cdx_hg37.bed | bedtools merge -i - > cdx_hg37.sorted.bed
bedtools sort -i liquid_cdx_hg38.bed | bedtools merge -i - > liquid_cdx_hg38.sorted.bed
bedtools sort -i liquid_cdx_hg37.bed | bedtools merge -i - > liquid_cdx_hg37.sorted.bed
