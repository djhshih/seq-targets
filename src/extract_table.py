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
#### ~ region: intron/extron/3utr/promoter
#### ~ INFO: additional information


# @ intron xx,xx-xx,xx => numbers
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
    




def clean_input_string(input_string, target_region):
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
    if target_region == "extron":
        return(f'{gene_name}\t{alt_name}\textron\t\n')


    info = input_split[1].lower().\
            replace("\n","").replace("*","").\
            replace("- ", "-"). replace(", ",",").\
            strip()
    # @ utr, ncrna, promoter
    if re.search("(utr)|(ncrna)|(promoter)",info):
        return(f'{gene_name}\t{alt_name}\t{info}\t\n')

    # @ introns:
    else:
        info = clean_intron_entry(info)
        return(f'{gene_name}\t{alt_name}\tintron\t{info}\n')

    


with open("gene_list.tsv", "w+") as f:
    f.write("gene_name_1\tgene_name_2\tregion\tinfo\n")

    ## ~ Tidy coding extrons
    for sub_list in extrons_raw:
        for sub_gene in sub_list:
            if sub_gene == "":
                continue
            f.write(clean_input_string(sub_gene, "extron")) 

    ## ~ Tidy selected introns
    for sub_list in introns_raw:
        for sub_gene in sub_list:
            if sub_gene == "":
                continue
            f.write(clean_input_string(sub_gene, "selected_intron"))
