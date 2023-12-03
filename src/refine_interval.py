import re

# @ Clean the intron.interval_list
## @ Read in the gene list and grab all targets:
target_gene_list = [[],[]]
with open("../output/gene_list.tsv") as f:
    for line in f:
        current_line = line.split()
        if 'intron' in current_line:
            target_gene_list[0].append(current_line[0])
            intron_vec = [int(i) for i in current_line[-1].split(',')]
            target_gene_list[1].append(intron_vec)

## @ Read line by line of interval list, calculate the exon and exon number:

## * Pointers initialization
target_transcript = ""
target_gene_pointer = 0
target_gene = target_gene_list[0][target_gene_pointer]
target_intron=target_gene_list[1][target_gene_pointer]
transcript_start_coordinate = 0
transcript_end_coordinate = 0

## * A flag to indicate whether we have record the end of previous
record_flag = False

## * Keep a set to remove duplicate line for output:
all_output_lines = set()

with open("../output/intron.interval_list.raw") as f:
    for line in f:
        # @ Current_line: CHROM, START, END, STRAND, GENE, TRANSCRIPT, EXON_NUMBER
        current_line = line.split()
        current_line = [re.sub("^.+?=","",i) for i in current_line]
        
        # @ 1. Check whether changed the gene:
        if current_line[4] != target_gene:
            record_flag = False
            target_gene_pointer = target_gene_pointer + 1
            target_gene = target_gene_list[0][target_gene_pointer]
            target_intron=target_gene_list[1][target_gene_pointer]
              
        # @ 2. Check whether changed the transcript:
        if current_line[5] != target_transcript:
            record_flag = False
            target_transcript = current_line[5]
            target_intron = target_gene_list[1][target_gene_pointer]
            
        # @ 3. Check whether we need to record:
        if record_flag:
            record_flag = False
            # * Get the end coordinate
            if current_line[3] == "+":
                transcript_end_coordinate = int(current_line[1]) - 1
            else:
                transcript_start_coordinate = int(current_line[2]) + 1
            # * Record target exon:
            all_output_lines.add(f'{current_line[0]}\t{transcript_start_coordinate}\t{transcript_end_coordinate}\n')

        # @ Check whether we need to calculate this exon:
        # @@ Correct transcript + extron number:
        current_exon_number = int((current_line[6].split("="))[0])

        if len(target_intron) > 0 and current_exon_number == target_intron[0]:
            # * Record the start coordinate:
            record_flag = True
            if current_line[3] == "+":
                transcript_start_coordinate = int(current_line[2]) + 1
            else:
                transcript_end_coordinate = int(current_line[1]) - 1
            # * Remove current intron index
            target_intron = target_intron[1:]

# @ Write the final output
all_output_lines = list(all_output_lines)
with open("../output/intron.interval_list", "w+") as f:
    f.writelines(all_output_lines)

# @ Clean the promoter.interval_list
promoter_output = open("../output/promoter.interval_list", "w+")
with open("../output/promoter.interval_list.raw") as f:
    for line in f:
        current_line = line.split()
        # @ Determine by strand:
        if current_line[3] == "+":
            promoter_output.write(f'{current_line[0]}\t{int(current_line[1]) - 500}\t{int(current_line[2])}\n')
        else:
            promoter_output.write(f'{current_line[0]}\t{int(current_line[1])}\t{int(current_line[2]) + 500}\n')
