#/bin/bash
rm *.interval_list
rm *.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ABL1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BRAF;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CDKN1A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=EPHA3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGFR4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=IKZF1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MCL1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NKX2-1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PMS2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RNF43;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TET2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ACVR1B;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BRCA1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CDKN1B;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=EPHB1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FH;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=INPP4B;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MDM2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NOTCH1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=POLD1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ROS1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TGFBR2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=AKT1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BRCA2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CDKN2A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=EPHB4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FLCN;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=IRF2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MDM4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NOTCH2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=POLE;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RPTOR;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TIPARP;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=AKT2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BRD4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CDKN2B;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ERBB2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FLT1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=IRF4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MED12;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NOTCH3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PPARG;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SDHA;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TNFAIP3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=AKT3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BRIP1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CDKN2C;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ERBB3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FLT3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=IRS2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MEF2B;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NPM1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PPP2R1A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SDHB;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TNFRSF14;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ALK;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BTG1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CEBPA;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ERBB4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FOXL2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=JAK1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MEN1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NRAS;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PPP2R2A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SDHC;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TP53;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ALOX12B;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BTG2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CHEK1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ERCC4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FUBP1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=JAK2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MERTK;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NT5C2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PRDM1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SDHD;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TSC1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=AMER1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BTK;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CHEK2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ERG;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=GABRA6;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=JAK3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MET;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NTRK1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PRKAR1A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SETD2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TSC2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=APC;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=C11orf30;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CIC;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ERRFI1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=GATA3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=JUN;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MITF;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NTRK2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PRKCI;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SF3B1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TYRO3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=AR;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CALR;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CREBBP;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ESR1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=GATA4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=KDM5A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MKNK1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NTRK3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PTCH1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SGK1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=U2AF1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ARAF;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CARD11;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CRKL;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=EZH2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=GATA6;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=KDM5C;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MLH1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=P2RY8;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PTEN;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SMAD2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=VEGFA;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ARFRP1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CASP8;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CSF1R;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FAM46C;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} (/;gene_name=GID4;/||/;gene_name=C17orf39;/) && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=KDM6A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MPL;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PALB2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PTPN11;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SMAD4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=VHL;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ARID1A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CBFB;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CSF3R;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FANCA;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=GNA11;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=KDR;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MRE11A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PARK2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PTPRO;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SMARCA4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=WHSC1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ASXL1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CBL;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CTCF;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FANCC;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=GNA13;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=KEAP1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MSH2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PARP1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=QKI;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SMARCB1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=WHSC1L1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ATM;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CCND1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CTNNA1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FANCG;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=GNAQ;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=KEL;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MSH3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PARP2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RAC1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SMO;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=WT1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ATR;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CCND2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CTNNB1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FANCL;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=GNAS;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=KIT;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MSH6;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PARP3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RAD21;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SNCAIP;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=XPO1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ATRX;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CCND3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CUL3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FAS;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=GRM3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=KLHL6;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MST1R;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PAX5;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RAD51;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SOCS1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=XRCC2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=AURKA;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CCNE1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CUL4A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FBXW7;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=GSK3B;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} (/;gene_name=KMT2A;/||/;gene_name=MLL;/) && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MTAP;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PBRM1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RAD51B;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SOX2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ZNF217;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=AURKB;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CD22;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CXCR4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGF10;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=H3F3A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} (/;gene_name=KMT2D;/||/;gene_name=MLL2;/) && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MTOR;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PDCD1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RAD51C;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SOX9;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ZNF703;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=AXIN1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CD274;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CYP17A1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGF12;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=HDAC1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=KRAS;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MUTYH;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PDCD1LG2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RAD51D;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SPEN;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=AXL;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CD70;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=DAXX;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGF14;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=HGF;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=LTK;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MYC;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PDGFRA;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RAD52;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SPOP;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BAP1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CD79A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=DDR1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGF19;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=HNF1A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=LYN;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MYCL;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PDGFRB;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RAD54L;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SRC;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BARD1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CD79B;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=DDR2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGF23;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=HRAS;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MAF;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MYCN;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PDK1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RAF1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=STAG2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BCL2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CDC73;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=DIS3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGF3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=HSD3B1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MAP2K1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MYD88;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PIK3C2B;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RARA;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=STAT3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BCL2L1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CDH1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=DNMT3A;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGF4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ID3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MAP2K2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NBN;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PIK3C2G;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RB1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=STK11;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BCL2L2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CDK12;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=DOT1L;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGF6;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=IDH1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MAP2K4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NF1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PIK3CA;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RBM10;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SUFU;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BCL6;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CDK4;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=EED;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGFR1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=IDH2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MAP3K1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NF2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PIK3CB;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=REL;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SYK;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BCOR;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CDK6;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=EGFR;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGFR2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=IGF1R;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MAP3K13;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NFE2L2;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PIK3R1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RET;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TBX3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BCORL1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CDK8;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=EP300;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGFR3;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=IKBKE;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MAPK1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NFKBIA;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PIM1;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RICTOR;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TEK;/ && $3=="exon" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> extron.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ALK;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BRCA1;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ETV4;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=EZR;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=KIT;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MYC;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NUTM1;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RET;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SLC34A2;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BCL2;/ && $3=="three_prime_UTR" {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> 3utr.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BRCA2;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ETV5;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGFR1;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} (/;gene_name=KMT2A;/||/;gene_name=MLL;/) && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NOTCH2;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=PDGFRA;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ROS1;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TERC;/ && $3=="gene" && /;gene_type=lncRNA;/ {print $1,$4-1,$5-1}' ../gtf/gencode.v44lift37.annotation.gff3 >> lnrna.interval_list
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BCR;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=CD74;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=ETV6;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGFR2;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MSH2;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NTRK1;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RAF1;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RSPO2;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TERT;/ && $3=="five_prime_UTR" {print $1,$4-1,$5-1,$7}' ../gtf/gencode.v44lift37.annotation.gff3 >> promoter.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=BRAF;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=EGFR;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=EWSR1;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=FGFR3;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=MYB;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=NTRK2;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=RARA;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=SDC4;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
awk 'BEGIN {FS="[\t;]"; OFS="\t"} /;gene_name=TMPRSS2;/ && $3=="exon" {print $1,$4-1,$5-1,$7, $14,$16,$17}' ../gtf/gencode.v44lift37.annotation.gff3 >> intron.interval_list.raw
