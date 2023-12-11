import re
import argparse

# ! Args:
parser = argparse.ArgumentParser(description="Process 5 subset of gff file and generate interval list in .bed format")

### ~ Add the arguments
parser.add_argument("--target_table", type=str, required=True, help="Path to the gene information table.")
parser.add_argument("--input_path", type=str, required=True, help="Path to the input file")
parser.add_argument("--input_prefix", type=str, required=True, help="Prefix for 5 raw interval list")
parser.add_argument("--output_path", type=str, required=True, help="Path to the output table.")

### ~ Parse the arguments
args = parser.parse_args()

input_path = args.input_path
if input_path[-1] != "/":
    input_path += "/"


output_path = args.output_path
if output_path[-1] != "/":
    output_path += "/"

# ! Read in the gene list and grab all targets:
target_gene_dict = {}

with open(args.target_table) as f:
    next(f)
    for line in f:
        current_line = line.split("\t")
        if current_line[2] != "":
            intron_vec = [int(i) for i in current_line[2].split(',')]
            target_gene_dict[current_line[0]] = intron_vec

# ! Clean the intron.interval_list.raw
## @ Read line by line of interval list, calculate the exon and exon number:
## * Pointers initialization
target_transcript = ""
target_gene = ""
transcript_start_coordinate = 0
transcript_end_coordinate = 0

## * A flag to indicate whether we have record the end of previous
record_flag = False

## * Keep a set to remove duplicate line for output:
all_output_lines = set()

with open(f'{input_path}{args.input_prefix}.intron.interval_list.raw') as f:
    for line in f:
        # @ current_line: CHROM, START, END, STRAND, GENE, TRANSCRIPT, EXON_NUMBER
        current_line = line.split()
        current_line = [re.sub("^.+?=","",i) for i in current_line]
        
        # @ 0. Include this extron first:
        all_output_lines.add(f'\t{current_line[1]}\t{current_line[2]}\n')
        
        # @ 1. Check whether changed the gene:
        if current_line[4] != target_gene:
            record_flag = False
            target_gene = current_line[4]
            target_intron=target_gene_dict[target_gene]
              
        # @ 2. Check whether changed the transcript:
        if current_line[5] != target_transcript:
            record_flag = False
            target_transcript = current_line[5]
            target_intron = target_gene_dict[target_gene]
            
        # @ 3. Check whether we need to record:
        if record_flag:
            record_flag = False
            # * Get the end coordinate
            if current_line[3] == "+":
                transcript_end_coordinate = int(current_line[1]) 
            else:
                transcript_start_coordinate = int(current_line[2]) 
            # * Record target exon:
            all_output_lines.add(f'{current_line[0]}\t{transcript_start_coordinate}\t{transcript_end_coordinate}\n')

        # @ Check whether we need to calculate this exon:
        # @@ Correct transcript + extron number:
        current_exon_number = int((current_line[6].split("="))[0])

        if len(target_intron) > 0 and current_exon_number == target_intron[0]:
            # * Record the start coordinate:
            record_flag = True
            if current_line[3] == "+":
                transcript_start_coordinate = int(current_line[2]) 
            else:
                transcript_end_coordinate = int(current_line[1]) 
            # * Remove current intron index
            target_intron = target_intron[1:]

# @ Write the final output
all_output_lines = list(all_output_lines)
with open(f'{output_path}{args.input_prefix}.intron.interval_list', "w+") as f:
    f.writelines(all_output_lines)

# ! Clean the promoter.interval_list.raw
all_output_lines = set()
with open(f'{input_path}{args.input_prefix}.promoter.interval_list.raw') as f:
    for line in f:
        current_line = line.split()
        #current_line[0] = current_line[0].replace("chr", "")
        # @ Determine by strand:
        if current_line[3] == "+":
            all_output_lines.add(f'{current_line[0]}\t{int(current_line[1]) - 500}\t{int(current_line[2])}\n')
        else:
            all_output_lines.add(f'{current_line[0]}\t{int(current_line[1])}\t{int(current_line[2]) + 500}\n')

all_output_lines = list(all_output_lines)
with open(f'{output_path}{args.input_prefix}.promoter.interval_list', "w+") as f:
    f.writelines(all_output_lines)


# ! Clean rest interval_list.raw
## @ Remove duplication + remove dummy line
def clean_interval_list(input_path, interval_type):
    all_output_lines = set()
    with open(input_path) as f:
        for line in f:
            current_line = line.split()
            #current_line[0] = current_line[0].replace("chr", "")
            if current_line[1] == current_line[2]:
                continue
            all_output_lines.add(f'{current_line[0]}\t{current_line[1]}\t{current_line[2]}\n')

    all_output_lines = list(all_output_lines)
    with open(f'{output_path}{args.input_prefix}.{interval_type}.interval_list', "w+") as f:
        f.writelines(all_output_lines)

clean_interval_list(f'{input_path}{args.input_prefix}.exon.interval_list.raw', "exon")
clean_interval_list(f'{input_path}{args.input_prefix}.3utr.interval_list.raw', "3utr")
clean_interval_list(f'{input_path}{args.input_prefix}.ncrna.interval_list.raw', "ncrna")