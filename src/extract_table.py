import pdfplumber
import re

target_genome = "../gtf/gencode.v44lift37.annotation.gff3"

## ~ Input:
### ~ Target input
#### ~ P4: Full coding extrons
#### ~ P5: Selected introns
pdf_in = pdfplumber.open("../output/gene_list_meta.pdf")
extrons_raw = pdf_in.pages[3].extract_table()
introns_raw = pdf_in.pages[4].extract_table()

## ~ Target output
### ~ 1: a human-readable table for all extracted information
#### ~ gene_name_1 \t alt_name  \t region \t info 
#### ~ gene_name_1: name of gene
#### ~ alt_name: alternative gene name (if have)
#### ~ region: intron/extron/3utr/promoter
#### ~ info: additional information

### ~ 2: A bash code
#### ~ Generate subset of gff3 files contains all target genes
#### ~ for further refinement.

## ~ Functions
# @ intron entry:
# @@ Input: intron(s) xx,xx-xx,xx 
# @@ Output: xx,xx,xx,xx, include all intervals
def clean_intron_entry(input_entry):
    target_intervals = input_entry.replace("intron ","").\
        replace("introns ","").split(",")

    full_interval = []
    for sub_interval in target_intervals:
        if re.search("-", sub_interval):
            start, end = sub_interval.split("-")
            start = int(start)
            while(start <= int(end)):
                full_interval.append(str(start))
                start = start + 1
        else:
            full_interval.append(sub_interval)

    return(",".join(full_interval))

# @ Clean the full input in various cases: 
def clean_input_string(input_string, extron = True):
    alt_name = ""
    # @ Detect possible alternative name of gene:
    if re.search("\\(",input_string):
        # @ Extract alternative name
        alt_name = re.findall("\\((.*?)\\)", input_string)[0]
        alt_name = alt_name.strip()
        # @ Replace alternative name
        input_string = \
            input_string.replace("(","").\
            replace(alt_name,"").replace(")","")
    
    # @ Split by the first /n, get the gene name
    input_split = input_string.split("\n",1)
    gene_name = input_split[0].strip()
    
    # @ Extron:
    if extron:
        return(gene_name, alt_name, "extron", "")


    info = input_split[1].lower().\
            replace("\n","").replace("*","").\
            replace("- ", "-"). replace(", ",",").\
            strip()
    # @ utr, ncrna, promoter
    if re.search("(utr)|(ncrna)|(promoter)",info):
        return(gene_name, alt_name, info, "")

    # @ introns:
    else:
        info = clean_intron_entry(info)
        return(gene_name, alt_name, "intron", info)

# @ Generate awk commands to extract target gene from gtf
def match_gene_name(gene_name, alt_name):
    if alt_name == "":
        gene_name_cmd = f'/;gene_name={gene_name};/'
    else:
        gene_name_cmd = f'(/;gene_name={gene_name};/||/;gene_name={alt_name};/)'
    return(gene_name_cmd)

## ! -1 for the coordinate to generate 0-based interval list

def process_extron(gene_name, alt_name, genome):
    gene_name_cmd = match_gene_name(gene_name, alt_name)
    ## ! Directly extracted:
    full_command = f'awk \'BEGIN {{FS=\"[\\t;]\"; OFS=\"\\t\"}} {gene_name_cmd} && $3==\"exon\" {{print $1,$4-1,$5-1}}\' {genome} >> extron.interval_list\n' 
    return(full_command)

def process_intron(gene_name, alt_name, genome):
    # @ Just extract, calculate will be in another file
    gene_name_cmd = match_gene_name(gene_name, alt_name)
    ## ! Keep the numbering of extron for intron numbering:
    full_command = f'awk \'BEGIN {{FS=\"[\\t;]\"; OFS=\"\\t\"}} {gene_name_cmd} && $3==\"exon\" {{print $1,$4-1,$5-1,$7, $14,$16,$17}}\' {genome} >> intron.interval_list.raw\n' 
    return(full_command)

def process_lnrna(gene_name, alt_name, genome):
    # ? Which entry to include (i.e. gene/transcript/exon)
    gene_name_cmd = match_gene_name(gene_name, alt_name)
    ## ! Keep the numbering of extron for intron numbering:
    full_command = f'awk \'BEGIN {{FS=\"[\\t;]\"; OFS=\"\\t\"}} {gene_name_cmd} && $3==\"gene\" && /;gene_type=lncRNA;/ {{print $1,$4-1,$5-1}}\' {genome} >> lnrna.interval_list\n' 
    return(full_command)

def process_3utr(gene_name, alt_name, genome):
    gene_name_cmd = match_gene_name(gene_name, alt_name)
    ## ! Direct extract:
    full_command = f'awk \'BEGIN {{FS=\"[\\t;]\"; OFS=\"\\t\"}} {gene_name_cmd} && $3==\"three_prime_UTR\" {{print $1,$4-1,$5-1}}\' {genome} >> 3utr.interval_list\n' 
    return(full_command)

def process_promoter(gene_name, alt_name, genome):
    gene_name_cmd = match_gene_name(gene_name, alt_name)
    ## ! Extract with orientation: 
    full_command = f'awk \'BEGIN {{FS=\"[\\t;]\"; OFS=\"\\t\"}} {gene_name_cmd} && $3==\"five_prime_UTR\" {{print $1,$4-1,$5-1,$7}}\' {genome} >> promoter.interval_list.raw\n' 
    return(full_command)

## ~ Running:
tsv_out = open("../output/gene_list.tsv", "w+")
bash_out= open("../output/subset_gff3.sh","w+")
tsv_out.write("gene_name\talt_name\tregion\tinfo\n")
bash_out.write("#/bin/bash\n")
bash_out.write("rm *.interval_list\n")
bash_out.write("rm *.interval_list.raw\n")

### ~ Tidy coding extrons
for sub_list in extrons_raw:
    for sub_gene in sub_list:
        if sub_gene == "":
            continue
        gene_name, alt_name, region, info = clean_input_string(sub_gene, extron = True)
        tsv_out.write(f'{gene_name}\t{alt_name}\t{region}\t{info}\n')
        bash_out.write(process_extron(gene_name, alt_name, target_genome))

### ~ Tidy selected introns
for sub_list in introns_raw:
    for sub_gene in sub_list:
        if sub_gene == "":
            continue
        gene_name, alt_name, region, info = clean_input_string(sub_gene, extron = False)
        tsv_out.write(f'{gene_name}\t{alt_name}\t{region}\t{info}\n')
        
        if re.search("utr",region):
            bash_out.write(process_3utr(gene_name, alt_name, target_genome))
        elif re.search("ncrna",region):
            bash_out.write(process_lnrna(gene_name, alt_name, target_genome))
        elif re.search("promoter",region):
            bash_out.write(process_promoter(gene_name, alt_name, target_genome))
        else:
            bash_out.write(process_intron(gene_name, alt_name, target_genome))