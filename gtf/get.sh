# /bin/bash

# hg19:

fversion=44

if [[ ! -f "gencode.v${fversion}lift37.annotation.gff3" ]]; then
    if [[ ! -f "gencode.v${fversion}lift37.annotation.gff3.gz" ]]; then
        wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${fversion}/GRCh37_mapping/gencode.v${fversion}lift37.annotation.gff3.gz"
    fi
    gunzip "gencode.v${fversion}lift37.annotation.gff3.gz"
fi

# hg38

if [[ ! -f "gencode.v${fversion}.annotation.gff3" ]]; then
    if [[ ! -f "gencode.v${fversion}lift37.annotation.gtf.gz" ]]; then
        wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${fversion}/gencode.v${fversion}.annotation.gff3.gz"
    fi
    gunzip "gencode.v${fversion}.annotation.gff3.gz"
fi
