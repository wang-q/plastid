# 

## evaluation

```shell script
mkdir -p ~/data/organelles/ena
cd ~/data/organelles/ena

cat << EOF > source.csv
SRX202246,Atha_Col_0,HiSeq 2000 PE100
SRX179254,Osat_nip,HiSeq 2000 PE100
SRX673852,Mtru_A17,HiSeq 2000 PE150
SAMN02946808,Gmax_Soybean_100,HiSeq 2000 PE100
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 4 -s 2 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name             | srx       | platform | layout | ilength | srr        | spot      | base   |
|:-----------------|:----------|:---------|:-------|:--------|:-----------|:----------|:-------|
| Atha_Col_0       | SRX202246 | ILLUMINA | PAIRED | 450     | SRR611086  | 49891349  | 9.29G  |
| Atha_Col_0       | SRX202246 | ILLUMINA | PAIRED | 450     | SRR616966  | 24851796  | 4.63G  |
| Gmax_Soybean_100 | SRX666906 | ILLUMINA | PAIRED |         | SRR1533313 | 148394234 | 27.64G |
| Mtru_A17         | SRX673852 | ILLUMINA | PAIRED | 360     | SRR1542422 | 99418334  | 16.67G |
| Mtru_A17         | SRX673852 | ILLUMINA | PAIRED | 360     | SRR1542423 | 29663436  | 8.34G  |
| Osat_nip         | SRX179254 | ILLUMINA | PAIRED | 300     | SRR545059  | 85148124  | 7.93G  |
| Osat_nip         | SRX179254 | ILLUMINA | PAIRED | 300     | SRR545231  | 85251097  | 16.04G |

