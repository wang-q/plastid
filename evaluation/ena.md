# Download fastq files from ENA

[TOC levels=1-3]: # ""

- [Download fastq files from ENA](#download-fastq-files-from-ena)
  - [evaluation](#evaluation)


## evaluation

```shell script
mkdir -p ~/data/plastid/ena
cd ~/data/plastid/ena

cat << EOF > source.csv
SRX202246,Atha_Col_0_1,HiSeq 2000 PE100
SRX2527206,Atha_Col_0_2,MiSeq 2000 PE300
SRR616965,Atha_Ler_0,HiSeq 2000 PE100
SRX179254,Osat_Nip,HiSeq 2000 PE100
SRX673852,Mtru_A17,HiSeq 2000 PE150
SRX150254,Pper_Lovell,Illumina Genome Analyzer IIx PE100
SRX698770,Slyc_H1706,Illumina HiSeq 2000 PE100
SRX7009428,Gmax_W82,HiSeq X Ten
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml --ascp

mlr --icsv --omd cat ena_info.csv

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 1 "{}"

```

| name         | srx        | platform | layout | ilength | srr         | spots     | bases  |
|:-------------|:-----------|:---------|:-------|:--------|:------------|:----------|:-------|
| Atha_Col_0_1 | SRX202246  | ILLUMINA | PAIRED | 450     | SRR611086   | 49891349  | 9.29G  |
| Atha_Col_0_1 | SRX202246  | ILLUMINA | PAIRED | 450     | SRR616966   | 24851796  | 4.63G  |
| Atha_Col_0_2 | SRX2527206 | ILLUMINA | PAIRED |         | SRR5216995  | 26893065  | 14.46G |
| Atha_Ler_0   | SRX202247  | ILLUMINA | PAIRED | 450     | SRR611087   | 50791450  | 9.46G  |
| Atha_Ler_0   | SRX202247  | ILLUMINA | PAIRED | 450     | SRR616965   | 25436255  | 4.74G  |
| Gmax_W82     | SRX7009428 | ILLUMINA | PAIRED |         | SRR10296600 | 162110355 | 45.29G |
| Mtru_A17     | SRX673852  | ILLUMINA | PAIRED | 360     | SRR1542422  | 99418334  | 16.67G |
| Mtru_A17     | SRX673852  | ILLUMINA | PAIRED | 360     | SRR1542423  | 29663436  | 8.34G  |
| Osat_Nip     | SRX179254  | ILLUMINA | PAIRED | 300     | SRR545059   | 85148124  | 7.93G  |
| Osat_Nip     | SRX179254  | ILLUMINA | PAIRED | 300     | SRR545231   | 85251097  | 16.04G |
| Pper_Lovell  | SRX150254  | ILLUMINA | PAIRED | 400     | SRR502985   | 123590441 | 23.25G |
| Slyc_H1706   | SRX698770  | ILLUMINA | PAIRED |         | SRR1572628  | 24198345  | 4.51G  |

