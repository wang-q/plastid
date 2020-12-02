# Download fastq files from ENA

[TOC levels=1-3]: # ""

- [Download fastq files from ENA](#download-fastq-files-from-ena)
  - [evaluation](#evaluation)
  - [*Arabidopsis thaliana* 1001 Genomes Project](#arabidopsis-thaliana-1001-genomes-project)



## evaluation

```shell script
mkdir -p ~/data/plastid/ena
cd ~/data/plastid/ena

cat << EOF > source.csv
SRX202246,Atha_Col_0_1,HiSeq 2000 PE100
SRX2527206,Atha_Col_0_2,MiSeq 2000 PE300
SRX179254,Osat_Nip,HiSeq 2000 PE100
SRX673852,Mtru_A17,HiSeq 2000 PE150
SRX7009428,Gmax_W82,HiSeq X Ten
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 4 -s 2 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name         | srx        | platform | layout | ilength | srr         | spot      | base   |
|:-------------|:-----------|:---------|:-------|:--------|:------------|:----------|:-------|
| Atha_Col_0_1 | SRX202246  | ILLUMINA | PAIRED | 450     | SRR611086   | 49891349  | 9.29G  |
| Atha_Col_0_1 | SRX202246  | ILLUMINA | PAIRED | 450     | SRR616966   | 24851796  | 4.63G  |
| Atha_Col_0_2 | SRX2527206 | ILLUMINA | PAIRED |         | SRR5216995  | 26893065  | 14.46G |
| Gmax_W82     | SRX7009428 | ILLUMINA | PAIRED |         | SRR10296600 | 162110355 | 45.29G |
| Mtru_A17     | SRX673852  | ILLUMINA | PAIRED | 360     | SRR1542422  | 99418334  | 16.67G |
| Mtru_A17     | SRX673852  | ILLUMINA | PAIRED | 360     | SRR1542423  | 29663436  | 8.34G  |
| Osat_nip     | SRX179254  | ILLUMINA | PAIRED | 300     | SRR545059   | 85148124  | 7.93G  |
| Osat_nip     | SRX179254  | ILLUMINA | PAIRED | 300     | SRR545231   | 85251097  | 16.04G |

## *Arabidopsis thaliana* 1001 Genomes Project

* Download `Metadata` from NCBI SRA Run Selector via a web browser
  * https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA273563
  * Save it to `SraRunTable.txt`

* Download Accessions info
  * http://1001genomes.org/accessions.html
  * There's no header line
  * Rename `Accession ID` to `Ecotype`

```shell script
mkdir -p ~/data/plastid/1001/ena
cd ~/data/plastid/1001/ena

cat SraRunTable.txt |
    mlr --icsv --otsv cat \
    > SraRunTable.tsv

cat SraRunTable.tsv |
    tsv-filter -H \
        --str-in-fld TISSUE:leaf \
        --ge AvgSpotLen:90

curl -o accessions.csv \
    https://tools.1001genomes.org/api/accessions.csv?query=SELECT%20*%20FROM%20tg_accessions%20ORDER%20BY%20id

echo -e 'Ecotype\tName\tCS Number\tCountry\tLat\tLong\tCollector\tAdmixture Group\tSequenced by' \
    > accessions.tsv

cat accessions.csv |
    mlr --icsv --implicit-csv-header --otsv cat |
    sed 1d |
    tsv-select -f 1,3,10,4,6,7,8,11,2 \
    >> accessions.tsv

cat accessions.tsv |
    tsv-join -H --key-fields "Ecotype" \
        -f SraRunTable.tsv \
        --append-fields TISSUE,AvgSpotLen,Instrument,Bases,Experiment |
    tsv-filter -H \
        --str-in-fld TISSUE:leaf \
        --str-not-in-fld Instrument:"Genome Analyzer" \
        --ge Bases:2000000000 \
        --le Bases:10000000000 \
        --regex Name:'^[\w\d-]+$' |
    tsv-uniq -H \
        -f Country --max 2 |
    tsv-select -H -f Experiment,Name,Country |
    mlr --itsv --ocsv cat |
    sed 1d \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 4 -s 2 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name         | srx       | platform | layout | ilength | srr        | spot     | base  |
|:-------------|:----------|:---------|:-------|:--------|:-----------|:---------|:------|
| ANH-1        | SRX972422 | ILLUMINA | PAIRED |         | SRR1945739 | 10783643 | 2.03G |
| Aitba-1      | SRX972878 | ILLUMINA | PAIRED |         | SRR1946195 | 11730031 | 2.18G |
| Amel-1       | SRX972487 | ILLUMINA | PAIRED |         | SRR1945804 | 12458066 | 2.34G |
| Ang-0        | SRX972488 | ILLUMINA | PAIRED |         | SRR1945805 | 10482055 | 1.97G |
| BRR12        | SRX972131 | ILLUMINA | PAIRED |         | SRR1945448 | 63237282 | 5.71G |
| BRR4         | SRX972130 | ILLUMINA | PAIRED |         | SRR1945447 | 70270246 | 6.25G |
| Bela-1       | SRX972989 | ILLUMINA | PAIRED |         | SRR1946306 | 11510427 | 2.14G |
| Bik-1        | SRX973011 | ILLUMINA | PAIRED |         | SRR1946328 | 14726330 | 2.77G |
| Bl-1         | SRX972496 | ILLUMINA | PAIRED |         | SRR1945813 | 10355956 | 1.95G |
| Bs-1         | SRX972492 | ILLUMINA | PAIRED |         | SRR1945809 | 10547692 | 1.98G |
| CYR          | SRX972118 | ILLUMINA | PAIRED |         | SRR1945435 | 15022786 | 2.8G  |
| Co           | SRX972512 | ILLUMINA | PAIRED |         | SRR1945829 | 31073234 | 5.85G |
| Co-1         | SRX972511 | ILLUMINA | PAIRED |         | SRR1945828 | 9909746  | 1.86G |
| Cvi-0        | SRX972441 | ILLUMINA | PAIRED |         | SRR1945758 | 36851486 | 6.93G |
| Dja-1        | SRX972148 | ILLUMINA | PAIRED |         | SRR1945465 | 10193882 | 1.92G |
| Dolen-1      | SRX972959 | ILLUMINA | PAIRED |         | SRR1946276 | 21745054 | 4.05G |
| Doubravnik7  | SRX972126 | ILLUMINA | PAIRED |         | SRR1945443 | 18807933 | 3.5G  |
| Ei-2         | SRX972443 | ILLUMINA | PAIRED |         | SRR1945760 | 11999881 | 2.26G |
| Es-0         | SRX972526 | ILLUMINA | PAIRED |         | SRR1945843 | 12184312 | 2.29G |
| Est          | SRX972527 | ILLUMINA | PAIRED |         | SRR1945844 | 10419012 | 1.96G |
| Faneronemi-3 | SRX972985 | ILLUMINA | PAIRED |         | SRR1946302 | 11164743 | 2.08G |
| Geg-14       | SRX972730 | ILLUMINA | PAIRED |         | SRR1946047 | 28986447 | 5.4G  |
| Goced-1      | SRX972960 | ILLUMINA | PAIRED |         | SRR1946277 | 14331777 | 2.67G |
| Gr-1         | SRX972129 | ILLUMINA | PAIRED |         | SRR1945446 | 9957191  | 1.87G |
| Gradi-1      | SRX972914 | ILLUMINA | PAIRED |         | SRR1946231 | 10677561 | 1.99G |
| Kar-1        | SRX972146 | ILLUMINA | PAIRED |         | SRR1945463 | 18918780 | 3.56G |
| Kas-1        | SRX972543 | ILLUMINA | PAIRED |         | SRR1945860 | 36472683 | 6.86G |
| Ko-2         | SRX972617 | ILLUMINA | PAIRED |         | SRR1945934 | 10379035 | 1.95G |
| Kondara      | SRX972453 | ILLUMINA | PAIRED |         | SRR1945770 | 11713120 | 2.2G  |
| Kz-9         | SRX972454 | ILLUMINA | PAIRED |         | SRR1945771 | 12026979 | 2.26G |
| LDV-18       | SRX972119 | ILLUMINA | PAIRED |         | SRR1945436 | 25591692 | 4.77G |
| La-0         | SRX972551 | ILLUMINA | PAIRED |         | SRR1945868 | 20083405 | 3.78G |
| Lag1-2       | SRX972719 | ILLUMINA | PAIRED |         | SRR1946036 | 34465657 | 6.42G |
| Lag1-4       | SRX972720 | ILLUMINA | PAIRED |         | SRR1946037 | 27090130 | 5.05G |
| Ms-0         | SRX972457 | ILLUMINA | PAIRED |         | SRR1945774 | 16166537 | 3.04G |
| Neo-6        | SRX972150 | ILLUMINA | PAIRED |         | SRR1945467 | 13690770 | 2.58G |
| Nok-3        | SRX972461 | ILLUMINA | PAIRED |         | SRR1945778 | 20380029 | 3.83G |
| Olympia-2    | SRX972986 | ILLUMINA | PAIRED |         | SRR1946303 | 11570327 | 2.16G |
| PHW-2        | SRX972661 | ILLUMINA | PAIRED |         | SRR1945978 | 18690981 | 3.52G |
| Rubezhnoe-1  | SRX972580 | ILLUMINA | PAIRED |         | SRR1945897 | 16257115 | 3.06G |
| Se-0         | SRX972468 | ILLUMINA | PAIRED |         | SRR1945785 | 13228869 | 2.49G |
| Stiav-2      | SRX972988 | ILLUMINA | PAIRED |         | SRR1946305 | 20683246 | 3.85G |
| Strand-1     | SRX973236 | ILLUMINA | PAIRED |         | SRR1946553 | 29763531 | 4.44G |
| Tamm-27      | SRX972473 | ILLUMINA | PAIRED |         | SRR1945790 | 24125976 | 4.54G |
| Teiu-2       | SRX972994 | ILLUMINA | PAIRED |         | SRR1946311 | 10261930 | 1.91G |
| Ts-5         | SRX972475 | ILLUMINA | PAIRED |         | SRR1945792 | 32959846 | 6.2G  |
| UKSW06-179   | SRX972227 | ILLUMINA | PAIRED |         | SRR1945544 | 17662847 | 3.29G |
| UKSW06-207   | SRX972228 | ILLUMINA | PAIRED |         | SRR1945545 | 15152715 | 2.82G |
| Ulies-1      | SRX972995 | ILLUMINA | PAIRED |         | SRR1946312 | 10003075 | 1.86G |
| Uod-1        | SRX972478 | ILLUMINA | PAIRED |         | SRR1945795 | 11518606 | 2.17G |
| Van-0        | SRX972603 | ILLUMINA | PAIRED |         | SRR1945920 | 16055225 | 3.02G |
| Wa-1         | SRX972606 | ILLUMINA | PAIRED |         | SRR1945923 | 11537790 | 2.17G |
| Wei-0        | SRX972480 | ILLUMINA | PAIRED |         | SRR1945797 | 14332334 | 2.7G  |
| Wil-1        | SRX972693 | ILLUMINA | PAIRED |         | SRR1946010 | 22149228 | 4.17G |
| Ws-2         | SRX972481 | ILLUMINA | PAIRED |         | SRR1945798 | 11488247 | 2.16G |
| Xan-3        | SRX972706 | ILLUMINA | PAIRED |         | SRR1946023 | 14306256 | 2.66G |
| Xan-5        | SRX972707 | ILLUMINA | PAIRED |         | SRR1946024 | 39523387 | 7.29G |
| Yeg-2        | SRX972731 | ILLUMINA | PAIRED |         | SRR1946048 | 20491320 | 3.82G |
| Zabar-1      | SRX973002 | ILLUMINA | PAIRED |         | SRR1946319 | 13427571 | 2.5G  |
| Zagub-1      | SRX973003 | ILLUMINA | PAIRED |         | SRR1946320 | 10484877 | 1.95G |
| Zdarec3      | SRX972125 | ILLUMINA | PAIRED |         | SRR1945442 | 24356076 | 4.54G |
| Zupan-1      | SRX972913 | ILLUMINA | PAIRED |         | SRR1946230 | 13300921 | 2.48G |

