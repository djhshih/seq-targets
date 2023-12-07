import pdfplumber
import re
import argparse

# ! Args:
parser = argparse.ArgumentParser(description="Process input pdf file, extract table and info for each gene")

### ~ Add the arguments
parser.add_argument("--target_pdf", type=str, required=True, help="Path to the target PDF file.")
parser.add_argument("--target_page", type=int, nargs='+', required=True, help="List of target pages.")
parser.add_argument("--output", type=str, required=True, help="Path to the output table.")


### ~ Parse the arguments
args = parser.parse_args()

### ~ Example usage:
# * python extract_table.py --target_pdf /path/to/pdf --target_page 1 2 3 --output /path/to/output

### ! Expected patter in the annotation:
# * Information in <> are optional
# * GENE_NAME <(ALTERNATIVE_NAME)> <[ or { or (> <indicator of region> INFORMATION <] or } or )> 
# * Also some \n are expected. 


# ! ~ Target output
### ~  A human-readable table for all extracted information
#### ~ gene_name_1: name of gene
#### ~ alt_name: alternative gene name (if have)
#### ~ target_intron: intron number to be extract
#### ~ other_info: additional information

# ! Functions
#### ~ Transform input interval like xx, xx-xx, xx 
#### ~ to xx,xx,xx,xx,xx that include all numbers
def transform_interval_to_list(input_interval):
    full_interval = []
    input_interval = input_interval.split(",")
    for sub_interval in input_interval:
        if re.search("-", sub_interval):
            start, end = sub_interval.split("-")
            start = int(start)
            while(start <= int(end)):
                full_interval.append(str(start))
                start = start + 1
        else:
            full_interval.append(sub_interval)

    return(",".join(full_interval))


#### ~ Clean the string and get all information we need:
#### ~ Expected return: gene_name, alt_name, target_extrons, other_info
def clean_input_raw_gene_annotation(input_gene_string):
    # @ Remove non-sense items:
    input_gene_string = input_gene_string.replace("\n", " ").replace("*","").strip()

    # @ Bug fix for ALOX12B
    input_gene_string = input_gene_string.replace("A LOX12B", "ALOX12B")
    
    # @ Refine following information:
    ## @ (Coding Exons ...)
    ## @ (alternative designation exon ...)
    input_gene_string = re.sub(r"\(Coding Exons[^)]*\)", "", input_gene_string)
    matched_extron = re.search(r"exon (\d+)", input_gene_string)
    #print(matched_extron)
    if matched_extron:
        input_gene_string = re.sub(r"\(alternative designation exon (\d+)\)", f',{matched_extron.group(1)}', input_gene_string)
    
    # @ Refine bracket:
    input_gene_string = re.sub(r"\(|\{|\[","(", input_gene_string)
    input_gene_string = re.sub(r"\)|\}|\]",")", input_gene_string)

    # @ Split by the first whitespace and get the gene name:
    input_split = input_gene_string.split(" ", 1)
    target_gene = input_split[0]

    if len(input_split) == 1:
        return(f'{target_gene}\t\t\t\n')
    
    # @ Add proper brackets:
    ### @ After this step, all information are separated by ()
    new_annotation_string = ""
    ### @ Flag indicate whether we need to add ( or )
    input_split[1] = input_split[1].replace(" ","").strip()
    left_flag = input_split[1][0] != "("
    for i in input_split[1]:
        if i != "(" and left_flag:
            left_flag = False
            new_annotation_string += "("
        if i == ")":
            left_flag = True
        new_annotation_string += i

    if i != ")":
        new_annotation_string += ")"


    # @ Breakdown the bracket
    bracket_pattern = r"\(.*?\)"
    all_item_in_bracket = re.findall(bracket_pattern, new_annotation_string)
    all_item_in_bracket = [re.sub(r"\(|\)", "", item) for item in all_item_in_bracket]

    #print(f'{target_gene} \t {all_item_in_bracket}')
    # @ Further clean things in the bracket:

    target_intron = ""
    other_info = ""
    target_gene_alt_name = ""
    for description in all_item_in_bracket:
        # @ intron or extron
        if re.match(r"[eE]xtron|[eE]xon", description):
            continue
        elif re.match(r"[iI]ntron", description):
            target_intron = transform_interval_to_list(re.sub(r"[a-zA-Z]","",description)).strip()
        # @ utr, ncrna, promoter
        elif re.match(r"(.*UTR)|(ncRNA)|([Pp]romoter)",description):
            other_info = description.strip()
        # @ alternative_name
        else:
            target_gene_alt_name = description.strip()

    return(f'{target_gene}\t{target_gene_alt_name}\t{target_intron}\t{other_info}\n')



# ! Workflow:
pdf_in = pdfplumber.open(args.target_pdf)
table_out = open(args.output, "w+")
table_out.write(f'gene_name\talt_name\ttarget_introns\tother_info\n')
for page in args.target_page:
    temp_table = pdf_in.pages[page-1].extract_table()
    for item in temp_table:
        for sub_item in item:
            if sub_item == "":
                continue
            else: 
                sub_item = sub_item.replace("\n", " ").replace("*","").replace("- ","-").strip()
                table_out.write(clean_input_raw_gene_annotation(sub_item))
