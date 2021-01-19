# *Arabidopsis thaliana* Columbia X Ler Genome sequencing

[TOC levels=1-3]: # ""

- [*Arabidopsis thaliana* Columbia X Ler Genome sequencing](#arabidopsis-thaliana-columbia-x-ler-genome-sequencing)
  - [基本信息](#基本信息)
  - [项目信息](#项目信息)
  - [数据下载](#数据下载)
    - [Reference](#reference)
    - [Illumina](#illumina)
  - [Symlink](#symlink)
  - [Run](#run)
  - [Pack and clean](#pack-and-clean)
  - [VCF](#vcf)


## 基本信息

* Genome: GCF_000001735.3, TAIR10, 119.668 Mb
* Chloroplast: [NC_000932](https://www.ncbi.nlm.nih.gov/nuccore/NC_000932), **Columbia**, 154478 bp
* Chloroplast: [KX551970](https://www.ncbi.nlm.nih.gov/nuccore/KX551970), **Landsberg erecta**,
  154515 bp
* Mitochondrion: [Y08501](https://www.ncbi.nlm.nih.gov/nuccore/Y08501), 366924 bp


## 项目信息

* <https://www.pnas.org/content/109/51/20992.long>
* PRJNA178613

* <https://www.nature.com/articles/nature14649>
* PRJNA243018
* PRJNA232554 - rice
* PRJNA252997 - bee


## 数据下载

### Reference

```shell script
mkdir -p ~/data/plastid/Atha_cross/genome
cd ~/data/plastid/Atha_cross/genome

for ACCESSION in "NC_000932" "Y08501"; do
    URL=$(printf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=%s&id=%s&retmode=text" "fasta" "${ACCESSION}");
    curl $URL -o ${ACCESSION}.fa
done

TAB=$'\t'
cat <<EOF > replace.tsv
NC_000932${TAB}Pt
Y08501${TAB}Mt
EOF

cat NC_000932.fa Y08501.fa |
    faops filter -s stdin stdout |
    faops replace stdin replace.tsv stdout |
    faops order stdin <(echo Pt; echo Mt) genome.fa

```

### Illumina

* Download `Metadata` from NCBI SRA Run Selector via a web browser
  * <https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA178613>
  * Save it to `SraRunTable.txt`
  * <https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA243018>
  * Save it to `SraRunTable.2.txt`

```shell script
mkdir -p ~/data/plastid/Atha_cross/ena
cd ~/data/plastid/Atha_cross/ena

cat SraRunTable.txt |
    mlr --icsv --otsv cat |
    tsv-select -H -f Experiment,"Sample\ Name",Bases \
    > SraRunTable.tsv

cat SraRunTable.2.txt |
    mlr --icsv --otsv cat |
    tsv-filter -H --str-eq Organism:"Arabidopsis thaliana" |
    tsv-select -H -f Experiment,"Sample\ Name",Bases |
    sed '1d' \
    >> SraRunTable.tsv

cat SraRunTable.tsv |
    sed '1 s/^/#/' |
    keep-header -- tsv-sort -k2,2 -k3,3nr |
    tsv-uniq -H -f "Sample\ Name" --max 1 |
    mlr --itsv --ocsv cat \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml --ascp

mlr --icsv --omd cat ena_info.csv

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 1 "{}"

```

| name            | srx        | platform | layout | ilength | srr        | spots    | bases  |
|:----------------|:-----------|:---------|:-------|:--------|:-----------|:---------|:-------|
| Col3-1          | SRX1003302 | ILLUMINA | PAIRED | 500     | SRR1984953 | 27101014 | 5.05G  |
| Col3-10         | SRX1003311 | ILLUMINA | PAIRED | 500     | SRR1984962 | 26353507 | 4.91G  |
| Col3-11         | SRX1003312 | ILLUMINA | PAIRED | 500     | SRR1984963 | 26796420 | 4.99G  |
| Col3-12         | SRX1003313 | ILLUMINA | PAIRED | 500     | SRR1984964 | 26830987 | 5G     |
| Col3-2          | SRX1003303 | ILLUMINA | PAIRED | 500     | SRR1984954 | 27243562 | 5.07G  |
| Col3-3          | SRX1003304 | ILLUMINA | PAIRED | 500     | SRR1984955 | 27243165 | 5.07G  |
| Col3-4          | SRX1003305 | ILLUMINA | PAIRED | 500     | SRR1984956 | 27162752 | 5.06G  |
| Col3-5          | SRX1003306 | ILLUMINA | PAIRED | 500     | SRR1984957 | 27206235 | 5.07G  |
| Col3-6          | SRX1003307 | ILLUMINA | PAIRED | 500     | SRR1984958 | 26372134 | 4.91G  |
| Col3-7          | SRX1003308 | ILLUMINA | PAIRED | 500     | SRR1984959 | 27073410 | 5.04G  |
| Col3-8          | SRX1003309 | ILLUMINA | PAIRED | 500     | SRR1984960 | 27270554 | 5.08G  |
| Col3-9          | SRX1003310 | ILLUMINA | PAIRED | 500     | SRR1984961 | 26514685 | 4.94G  |
| Col3_L1         | SRX1003300 | ILLUMINA | PAIRED | 500     | SRR1984951 | 27259770 | 5.08G  |
| Col3_L2         | SRX1003301 | ILLUMINA | PAIRED | 500     | SRR1984952 | 13615431 | 2.54G  |
| Col4_new        | SRX1003315 | ILLUMINA | PAIRED | 500     | SRR1984965 | 27355896 | 5.1G   |
| Col5            | SRX1003316 | ILLUMINA | PAIRED | 500     | SRR1984966 | 27290718 | 5.08G  |
| Ler1-1          | SRX1003319 | ILLUMINA | PAIRED | 500     | SRR1984970 | 31619686 | 5.89G  |
| Ler1-2          | SRX1003320 | ILLUMINA | PAIRED | 500     | SRR1984971 | 31642314 | 5.89G  |
| Ler1-3          | SRX1003321 | ILLUMINA | PAIRED | 500     | SRR1984972 | 31952689 | 5.95G  |
| Ler1-4          | SRX1003322 | ILLUMINA | PAIRED | 500     | SRR1984973 | 31974318 | 5.96G  |
| Ler1-5          | SRX1003323 | ILLUMINA | PAIRED | 500     | SRR1984974 | 31950307 | 5.95G  |
| Ler1_L1         | SRX1003317 | ILLUMINA | PAIRED | 500     | SRR1984967 | 27889098 | 5.19G  |
| Ler1_L2         | SRX1003318 | ILLUMINA | PAIRED | 500     | SRR1984968 | 13452095 | 2.51G  |
| Sample_14       | SRX202214  | ILLUMINA | PAIRED | 450     | SRR611088  | 25212780 | 4.7G   |
| Sample_18       | SRX202197  | ILLUMINA | PAIRED | 450     | SRR611072  | 12000000 | 2.23G  |
| Sample_19       | SRX202198  | ILLUMINA | PAIRED | 450     | SRR611073  | 12000000 | 2.23G  |
| Sample_20       | SRX202201  | ILLUMINA | PAIRED | 450     | SRR611074  | 12329694 | 2.3G   |
| Sample_21       | SRX202202  | ILLUMINA | PAIRED | 450     | SRR611075  | 12000000 | 2.23G  |
| Sample_4        | SRX202195  | ILLUMINA | PAIRED | 450     | SRR611076  | 12000000 | 2.23G  |
| Sample_5        | SRX202211  | ILLUMINA | PAIRED | 450     | SRR611089  | 12000000 | 2.23G  |
| Sample_5        | SRX202211  | ILLUMINA | PAIRED | 450     | SRR616963  | 13935870 | 2.6G   |
| Sample_6        | SRX202212  | ILLUMINA | PAIRED | 450     | SRR611090  | 23590718 | 4.39G  |
| Sample_7        | SRX202213  | ILLUMINA | PAIRED | 450     | SRR611091  | 25939674 | 4.83G  |
| Sample_8        | SRX202196  | ILLUMINA | PAIRED | 450     | SRR611077  | 12000000 | 2.23G  |
| Sample_Col_G    | SRX202246  | ILLUMINA | PAIRED | 450     | SRR611086  | 49891349 | 9.29G  |
| Sample_Col_G    | SRX202246  | ILLUMINA | PAIRED | 450     | SRR616966  | 24851796 | 4.63G  |
| Sample_Ler_XL_4 | SRX202247  | ILLUMINA | PAIRED | 450     | SRR611087  | 50791450 | 9.46G  |
| Sample_Ler_XL_4 | SRX202247  | ILLUMINA | PAIRED | 450     | SRR616965  | 25436255 | 4.74G  |
| Sample_c1c2     | SRX202204  | ILLUMINA | PAIRED | 450     | SRR611078  | 13983882 | 2.6G   |
| Sample_c41      | SRX202215  | ILLUMINA | PAIRED | 450     | SRR611092  | 25604873 | 4.77G  |
| Sample_c42      | SRX202216  | ILLUMINA | PAIRED | 450     | SRR611093  | 25505766 | 4.75G  |
| Sample_c45      | SRX202217  | ILLUMINA | PAIRED | 450     | SRR611094  | 25572903 | 4.76G  |
| Sample_c47      | SRX202203  | ILLUMINA | PAIRED | 450     | SRR611079  | 11642813 | 2.17G  |
| Sample_c48      | SRX202218  | ILLUMINA | PAIRED | 450     | SRR611095  | 25722535 | 4.79G  |
| Sample_c51      | SRX202219  | ILLUMINA | PAIRED | 450     | SRR611096  | 25759432 | 4.8G   |
| Sample_c52      | SRX202220  | ILLUMINA | PAIRED | 450     | SRR611097  | 25769891 | 4.8G   |
| Sample_c54      | SRX202221  | ILLUMINA | PAIRED | 450     | SRR611098  | 25762016 | 4.8G   |
| Sample_c57      | SRX202222  | ILLUMINA | PAIRED | 450     | SRR611099  | 25534424 | 4.76G  |
| Sample_c61      | SRX202223  | ILLUMINA | PAIRED | 450     | SRR611100  | 25601360 | 4.77G  |
| Sample_c62      | SRX202224  | ILLUMINA | PAIRED | 450     | SRR611101  | 25682651 | 4.78G  |
| Sample_c63      | SRX202225  | ILLUMINA | PAIRED | 450     | SRR611102  | 25535887 | 4.76G  |
| Sample_c64      | SRX202226  | ILLUMINA | PAIRED | 450     | SRR611103  | 25536189 | 4.76G  |
| Sample_c65      | SRX202227  | ILLUMINA | PAIRED | 450     | SRR611104  | 25493386 | 4.75G  |
| Sample_c66      | SRX202228  | ILLUMINA | PAIRED | 450     | SRR611105  | 25586286 | 4.77G  |
| Sample_c73      | SRX202229  | ILLUMINA | PAIRED | 450     | SRR611106  | 23348564 | 4.35G  |
| Sample_c81      | SRX202230  | ILLUMINA | PAIRED | 450     | SRR611107  | 25760277 | 4.8G   |
| Sample_c82      | SRX202231  | ILLUMINA | PAIRED | 450     | SRR611108  | 25625314 | 4.77G  |
| Sample_c83      | SRX202232  | ILLUMINA | PAIRED | 450     | SRR611109  | 25800250 | 4.81G  |
| Sample_c84      | SRX202233  | ILLUMINA | PAIRED | 450     | SRR611110  | 25880633 | 4.82G  |
| Sample_c85      | SRX202234  | ILLUMINA | PAIRED | 450     | SRR611111  | 25794536 | 4.8G   |
| Sample_c87      | SRX202235  | ILLUMINA | PAIRED | 450     | SRR611112  | 24917455 | 4.64G  |
| Sample_c88      | SRX202236  | ILLUMINA | PAIRED | 450     | SRR611113  | 24818024 | 4.62G  |
| Sample_c89      | SRX202237  | ILLUMINA | PAIRED | 450     | SRR611114  | 25551606 | 4.76G  |
| Sample_c90      | SRX202238  | ILLUMINA | PAIRED | 450     | SRR611115  | 25161363 | 4.69G  |
| Sample_c91      | SRX202240  | ILLUMINA | PAIRED | 450     | SRR611116  | 25624037 | 4.77G  |
| Sample_c92      | SRX202241  | ILLUMINA | PAIRED | 450     | SRR611117  | 25688992 | 4.78G  |
| Sample_c93      | SRX202242  | ILLUMINA | PAIRED | 450     | SRR611118  | 13983461 | 2.6G   |
| Sample_c93      | SRX202242  | ILLUMINA | PAIRED | 450     | SRR616982  | 11739658 | 2.19G  |
| Sample_c94      | SRX202243  | ILLUMINA | PAIRED | 450     | SRR611084  | 75625099 | 14.09G |
| Sample_c95      | SRX202245  | ILLUMINA | PAIRED | 450     | SRR611085  | 74530836 | 13.88G |
| Sample_l2c2     | SRX202205  | ILLUMINA | PAIRED | 450     | SRR611080  | 13984182 | 2.6G   |
| Sample_l2l3     | SRX202206  | ILLUMINA | PAIRED | 450     | SRR611081  | 13980212 | 2.6G   |
| Sample_l4c1     | SRX202207  | ILLUMINA | PAIRED | 450     | SRR611082  | 13984884 | 2.6G   |
| Sample_l4l3     | SRX202210  | ILLUMINA | PAIRED | 450     | SRR611083  | 13964774 | 2.6G   |
| c32             | SRX1003324 | ILLUMINA | PAIRED | 500     | SRR1984975 | 26484443 | 4.93G  |
| c33             | SRX1003325 | ILLUMINA | PAIRED | 500     | SRR1984976 | 27162437 | 5.06G  |
| c34             | SRX1003326 | ILLUMINA | PAIRED | 500     | SRR1984977 | 27224391 | 5.07G  |
| c35             | SRX1003327 | ILLUMINA | PAIRED | 500     | SRR1984978 | 27211820 | 5.07G  |
| c36             | SRX1003328 | ILLUMINA | PAIRED | 500     | SRR1984979 | 27309477 | 5.09G  |
| c37             | SRX1003329 | ILLUMINA | PAIRED | 500     | SRR1984980 | 27299393 | 5.08G  |
| c38             | SRX1003330 | ILLUMINA | PAIRED | 500     | SRR1984981 | 27184544 | 5.06G  |
| c39             | SRX1003331 | ILLUMINA | PAIRED |         | SRR1984982 | 27161575 | 5.06G  |
| c40             | SRX1003332 | ILLUMINA | PAIRED |         | SRR1984986 | 27166879 | 5.06G  |
| c43             | SRX1003333 | ILLUMINA | PAIRED |         | SRR1984988 | 26234630 | 4.89G  |
| c44             | SRX1003334 | ILLUMINA | PAIRED |         | SRR1984991 | 27315463 | 5.09G  |
| c46             | SRX1003335 | ILLUMINA | PAIRED |         | SRR1984993 | 27201317 | 5.07G  |
| c49             | SRX1003336 | ILLUMINA | PAIRED |         | SRR1984996 | 27048405 | 5.04G  |
| c50             | SRX1003337 | ILLUMINA | PAIRED |         | SRR1984997 | 27279919 | 5.08G  |
| c52-F4-1-3      | SRX1008009 | ILLUMINA | PAIRED |         | SRR1994588 | 26648552 | 4.96G  |
| c52-F4-1-4      | SRX1008011 | ILLUMINA | PAIRED |         | SRR1994590 | 26607356 | 4.96G  |
| c52-F4-11-2     | SRX1008003 | ILLUMINA | PAIRED |         | SRR1994584 | 26254971 | 4.89G  |
| c52-F4-11-3     | SRX1008006 | ILLUMINA | PAIRED |         | SRR1994585 | 26442401 | 4.93G  |
| c52-F4-12-2     | SRX1008007 | ILLUMINA | PAIRED |         | SRR1994586 | 27295717 | 5.08G  |
| c52-F4-12-3     | SRX1008008 | ILLUMINA | PAIRED |         | SRR1994587 | 27062342 | 5.04G  |
| c52-F4-13-1     | SRX1008010 | ILLUMINA | PAIRED |         | SRR1994589 | 27294669 | 5.08G  |
| c52-F4-14-2     | SRX1008012 | ILLUMINA | PAIRED |         | SRR1994591 | 26422833 | 4.92G  |
| c52-F4-14-3     | SRX1008013 | ILLUMINA | PAIRED |         | SRR1994592 | 26469653 | 4.93G  |
| c52-F4-15-1     | SRX1008014 | ILLUMINA | PAIRED |         | SRR1994593 | 26340857 | 4.91G  |
| c52-F4-16-2     | SRX1008015 | ILLUMINA | PAIRED |         | SRR1994594 | 27342095 | 5.09G  |
| c52-F4-16-3     | SRX1008016 | ILLUMINA | PAIRED |         | SRR1994595 | 26372590 | 4.91G  |
| c52-F4-20-1     | SRX1008017 | ILLUMINA | PAIRED |         | SRR1994596 | 26356850 | 4.91G  |
| c52-F4-20-2     | SRX1008018 | ILLUMINA | PAIRED |         | SRR1994597 | 26289494 | 4.9G   |
| c52-F4-4-2      | SRX1008019 | ILLUMINA | PAIRED |         | SRR1994598 | 27276319 | 5.08G  |
| c52-F4-4-3      | SRX1008020 | ILLUMINA | PAIRED |         | SRR1994599 | 26248628 | 4.89G  |
| c52-F4-6-1      | SRX1008021 | ILLUMINA | PAIRED |         | SRR1994600 | 27295867 | 5.08G  |
| c52-F4-6-2      | SRX1008022 | ILLUMINA | PAIRED |         | SRR1994601 | 27096012 | 5.05G  |
| c52-F4-8-2      | SRX1008023 | ILLUMINA | PAIRED |         | SRR1994602 | 26744977 | 4.98G  |
| c52-F4-8-3      | SRX1008024 | ILLUMINA | PAIRED |         | SRR1994603 | 26607908 | 4.96G  |
| c52-F4-9-1      | SRX1008025 | ILLUMINA | PAIRED |         | SRR1994604 | 26404489 | 4.92G  |
| c52-F4-9-3      | SRX1008026 | ILLUMINA | PAIRED |         | SRR1994605 | 27108095 | 5.05G  |
| c53             | SRX1003338 | ILLUMINA | PAIRED |         | SRR1984998 | 27307127 | 5.09G  |
| c55             | SRX1003339 | ILLUMINA | PAIRED |         | SRR1984999 | 27236774 | 5.07G  |
| c56             | SRX1003340 | ILLUMINA | PAIRED |         | SRR1985000 | 27332277 | 5.09G  |
| c58             | SRX1003342 | ILLUMINA | PAIRED |         | SRR1985001 | 26679600 | 4.97G  |
| c59             | SRX1003343 | ILLUMINA | PAIRED |         | SRR1985002 | 26214984 | 4.88G  |
| c60             | SRX1003344 | ILLUMINA | PAIRED |         | SRR1985004 | 26280677 | 4.9G   |
| c64-F4-18-2     | SRX1008027 | ILLUMINA | PAIRED |         | SRR1994606 | 26267733 | 4.89G  |
| c64-F4-19-2     | SRX1008028 | ILLUMINA | PAIRED |         | SRR1994607 | 27322502 | 5.09G  |
| c64-F4-3-3      | SRX1008029 | ILLUMINA | PAIRED |         | SRR1994608 | 27217905 | 5.07G  |
| c64-F4-5-1      | SRX1008030 | ILLUMINA | PAIRED |         | SRR1994609 | 27324609 | 5.09G  |
| c64-F4-5-2      | SRX1008031 | ILLUMINA | PAIRED |         | SRR1994610 | 26443978 | 4.93G  |
| c64-F4-5-5      | SRX1008032 | ILLUMINA | PAIRED |         | SRR1994611 | 27342895 | 5.09G  |
| c64-F4-8-2      | SRX1008033 | ILLUMINA | PAIRED |         | SRR1994612 | 26279734 | 4.89G  |
| c64-F4-9-1      | SRX1008034 | ILLUMINA | PAIRED |         | SRR1994613 | 27177988 | 5.06G  |
| c64-F4-9-4      | SRX1008035 | ILLUMINA | PAIRED |         | SRR1994614 | 27128176 | 5.05G  |
| c64-F4-9-5      | SRX1008036 | ILLUMINA | PAIRED |         | SRR1994615 | 26338648 | 4.91G  |
| c67             | SRX1003345 | ILLUMINA | PAIRED |         | SRR1985005 | 26223688 | 4.88G  |
| c68             | SRX1003346 | ILLUMINA | PAIRED |         | SRR1985006 | 26291494 | 4.9G   |
| c69             | SRX1003348 | ILLUMINA | PAIRED |         | SRR1985007 | 26222454 | 4.88G  |
| c70             | SRX1003349 | ILLUMINA | PAIRED |         | SRR1985009 | 27221518 | 5.07G  |
| c71             | SRX1003350 | ILLUMINA | PAIRED |         | SRR1985010 | 26278280 | 4.89G  |
| c72             | SRX1003351 | ILLUMINA | PAIRED |         | SRR1985013 | 27271550 | 5.08G  |
| c74             | SRX1003353 | ILLUMINA | PAIRED |         | SRR1985014 | 26244328 | 4.89G  |
| c75             | SRX1003354 | ILLUMINA | PAIRED |         | SRR1985016 | 26291054 | 4.9G   |
| c76             | SRX1003355 | ILLUMINA | PAIRED |         | SRR1985018 | 26693773 | 4.97G  |
| c77             | SRX1003356 | ILLUMINA | PAIRED |         | SRR1985022 | 26282852 | 4.9G   |
| c78             | SRX1003357 | ILLUMINA | PAIRED |         | SRR1985025 | 27113642 | 5.05G  |
| c79             | SRX1007997 | ILLUMINA | PAIRED |         | SRR1994579 | 27209386 | 5.07G  |
| c80             | SRX1008001 | ILLUMINA | PAIRED |         | SRR1994580 | 26977625 | 5.02G  |
| c86             | SRX1008002 | ILLUMINA | PAIRED |         | SRR1994581 | 27166830 | 5.06G  |

## Symlink

* 采用的倍数因子值: `2`

* Rsync to hpcc

```shell script
rsync -avP \
    ~/data/plastid/Atha_cross/ \
    wangq@202.119.37.251:data/plastid/Atha_cross

# rsync -avP wangq@202.119.37.251:data/plastid/Atha_cross/ ~/data/plastid/Atha_cross

```

```shell script
cd ~/data/plastid/Atha_cross/

export FOLD=2
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/col_0/chr.sizes |
        tsv-summarize --sum 2
)

cat ena/ena_info.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f name,srr,bases |
    perl -nla -F'\t' -e '
        BEGIN { our %seen }
        /^name/ and next;
        $seen{$F[0]} and next;
        my $bases = $F[2];
        $bases =~ s/G$//;
        my $cutoff = $bases * 1000 * 1000 * 1000 / $ENV{GENOME_SIZE} * $ENV{FOLD};
        $cutoff = int $cutoff;
        print join qq(\t), ($F[0], $F[1], $cutoff);
        $seen{$F[0]}++;
    ' \
    > opts.tsv

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ ! -f ena/{2}_1.fastq.gz ]; then
            exit;
        fi
        if [ ! -f ena/{2}_2.fastq.gz ]; then
            exit;
        fi

        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        mkdir -p {1}/1_genome
        pushd {1}/1_genome

        cp ../../genome/genome.fa genome.fa
        popd

        mkdir -p {1}/2_illumina
        pushd {1}/2_illumina

        ln -fs ../../ena/{2}_1.fastq.gz R1.fq.gz
        ln -fs ../../ena/{2}_2.fastq.gz R2.fq.gz
        popd
    '

```

* Rsync non-processed files to hpcc

```shell script
cd ~/data/plastid/Atha_cross/

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        find ena -type f -name "*{2}*"
    ' \
    > rsync.lst

rsync -avP \
    --files-from=rsync.lst \
    ~/data/plastid/Atha_cross/ \
    wangq@202.119.37.251:data/plastid/Atha_cross

```

## Run

```shell script
cd ~/data/plastid/Atha_cross/

cat opts.tsv | # head -n 170 | #tail -n 10 |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        if [ ! -d {1} ]; then
            exit;
        fi

        if [ ! -e {1}/2_illumina/R1.fq.gz ]; then
            exit;
        fi
        if [ ! -e {1}/2_illumina/R2.fq.gz ]; then
            exit;
        fi

        if bjobs -w | tr -s " " | cut -d " " -f 7 | grep -w "^{1}$"; then
            echo Job {1} exists
            exit;
        fi

        cd {1}

        echo {1}

        rm *.sh
        anchr template \
            --genome 1000000 \
            --parallel 24 \
            --xmx 80g \
            \
            --fastqc \
            --insertsize \
            --kat \
            \
            --trim "--dedupe --cutoff {3} --cutk 31" \
            --qual "25" \
            --len "60" \
            --filter "adapter artifact" \
            \
            --bwa Q25L60 \
            --gatk

        bsub -q mpi -n 24 -J "{1}" "
            bash 2_fastqc.sh
            bash 2_insert_size.sh
            bash 2_kat.sh
            bash 2_trim.sh
            bash 9_stat_reads.sh
            bash 3_bwa.sh
            bash 3_gatk.sh
        "
    '

```

## Pack and clean

```shell script
cd ~/data/plastid/Atha_cross/

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            echo "==> {1} .tar.gz"
            exit;
        fi

        if [ ! -f {1}/3_gatk/R.filtered.vcf ]; then
            echo "==> {1} 3_gatk"
            exit;
        fi

#        if [ ! -f {1}/7_merge_anchors/anchor.merge.fasta ]; then
#            echo "==> {1} 7_merge_anchors"
#            exit;
#        fi
#
#        if [ ! -d {1}/9_quast ]; then
#            echo "==> {1} 9_quast"
#            exit;
#        fi

        echo "==> Clean {1}"
        bash {1}/0_cleanup.sh

        echo "==> Create {1}.tar.gz"

        tar -czvf {1}.tar.gz \
            {1}/1_genome/genome.fa \
            {1}/2_illumina/fastqc \
            {1}/2_illumina/insert_size \
            {1}/2_illumina/kat \
            {1}/3_bwa/join.tsv \
            {1}/3_bwa/R.dedup.metrics \
            {1}/3_bwa/R.wgs.metrics \
            {1}/3_bwa/R.sort.bam \
            {1}/3_bwa/R.sort.bai \
            {1}/3_gatk \
            {1}/7_merge_anchors/anchor.merge.fasta \
            {1}/7_merge_anchors/others.non-contained.fasta \
            {1}/8_megahit/anchor/anchor.fasta \
            {1}/8_megahit/megahit.non-contained.fasta \
            {1}/8_mr_megahit/anchor/anchor.fasta \
            {1}/8_mr_megahit/megahit.non-contained.fasta \
            {1}/8_spades/anchor/anchor.fasta \
            {1}/8_spades/spades.non-contained.fasta \
            {1}/8_mr_spades/anchor/anchor.fasta \
            {1}/8_mr_spades/spades.non-contained.fasta \
            {1}/8_platanus/anchor/anchor.fasta \
            {1}/8_platanus/platanus.non-contained.fasta \
            {1}/9_quast \
            {1}/*.md

        echo
    '

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            if [ -d {1} ]; then
                echo "==> Remove {1}/"
                rm -fr {1}
            fi
        fi
    '

```

* Unpack

```shell script
cd ~/data/plastid/Atha_cross/

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -d {1} ]; then
            echo "==> {1} exists"
            exit;
        fi

        if [ ! -f {1}.tar.gz ]; then
            echo "==> {1}.tar.gz not exists"
            exit;
        fi

        tar -xzvf {1}.tar.gz
        rm {1}.tar.gz
    '

```

## VCF

```shell script
cd ~/data/plastid/Atha_cross/

mkdir -p vcf

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -f {1}.tar.gz ]; then
            echo "==> {1}.tar.gz not exists"
            exit;
        fi

        tar -xOzvf {1}.tar.gz {1}/3_gatk/R.filtered.vcf |
            bcftools reheader --samples <(echo {1}) |
            bcftools view \
                --apply-filters PASS --types snps --max-alleles 2 --targets Pt -Oz |
            bcftools view --include "AF>0.01" -Oz -o vcf/{1}.vcf.gz

        bcftools index -f vcf/{1}.vcf.gz
    '

bcftools merge --merge all -l <(
        cat opts.tsv |
            cut -f 1 |
            parallel -k -j 1 ' [ -f vcf/{}.vcf.gz ] && echo "vcf/{}.vcf.gz" '
    ) \
    > Atha_cross.vcf

```

