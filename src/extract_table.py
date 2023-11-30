import pdfplumber
import re

pdf_in = pdfplumber.open("gene_list_meta.pdf")

## ~ Target input
#### ~ P4: Full coding extrons
#### ~ P5: Selected introns
extrons_raw = pdf_in.pages[3].extract_table()
introns_raw = pdf_in.pages[4].extract_table()
# print(introns_raw)


## ~ Target output
#### ~ Gene_name_1 \t Gene_name_2  \t Group \t INFO 
#### ~ Gene_name_1: name of gene
#### ~ Gene name_2: alternative gene name (if have)
#### ~ Group: intron/extron
#### ~ INFO: additional information


def clean_input_string(input_string, group_name =""):
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
    
    # @ Tidy rest piece
    if len(input_split) == 2:
        info = \
            input_split[1].replace("\n","").\
            replace("*","").replace("- ", "-").\
            replace(", ",",").lower().strip()
    else:
        info = ""
    return(f'{gene_name}\t{alt_name}\t{group_name}\t{info}\n')


with open("gene_list.tsv", "w+") as f:

    ## ~ Tidy coding extrons
    for sub_list in extrons_raw:
        for sub_gene in sub_list:
            if sub_gene == "":
                continue
            f.write(clean_input_string(sub_gene, "all_extron")) 

    ## ~ Tidy selected introns
    for sub_list in introns_raw:
        for sub_gene in sub_list:
            if sub_gene == "":
                continue
            f.write(clean_input_string(sub_gene, "selected_intron"))
