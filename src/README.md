# Usage

`generate_interval_list.sh` warup up all relevant code together. A one time run of this script will finish the interval list generation from the pdf file. 

Following codes are wrapped up: 
* The script `extract_table.py` extract names and target region of gene from the FoundationOne pdf file:

```
python extract_table.py \
  --target_pdf </path/to/foundation1_pdf> \
  --target_page <list of page> \
  --output </path/to/output>
```
* The script `subset_gff.py` extract related line from gff file based on the output of `extract_table.py`. Records that related to exon, intron, promoter, utr and ncrna will be extracted. They will be saved in 5 files. 

```
python subset_gff.py \
  --gff <gff_file> \
  --output_path </path/to/output> \
  --output_prefix <output_file_prefix> \
  --target_table <output_table_from_extract_table>
```

* The script `refine_interval.py` use the 5 files generated in last step to extract the correct region 

```
python refine_interval.py \
  --input_path <output_path_from_subset_gff> \
  --input_prefix <output_prefix_from_subset_gff> \
  --target_table .<output_table_from_extract_table> \
  --output_path </path/to/output>
```
