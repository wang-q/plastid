# 

## evaluation

```shell script
mkdir -p ~/data/organelles/ena
cd ~/data/organelles/ena

cat << EOF > source.csv
SRX202246,Atha_Col_0,HiSeq 2000 PE100
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 4 -s 2 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name       | srx       | platform | layout | ilength | srr       | spot     | base  |
|:-----------|:----------|:---------|:-------|:--------|:----------|:---------|:------|
| Atha_Col_0 | SRX202246 | ILLUMINA | PAIRED | 450     | SRR611086 | 49891349 | 9.29G |
| Atha_Col_0 | SRX202246 | ILLUMINA | PAIRED | 450     | SRR616966 | 24851796 | 4.63G |

