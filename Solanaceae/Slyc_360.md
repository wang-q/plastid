# 🍅 *Solanum lycopersicum* 360 accessions

[TOC levels=1-3]: # ""

- [🍅 *Solanum lycopersicum* 360 accessions](#-solanum-lycopersicum-360-accessions)
  - [基本信息](#基本信息)
  - [项目信息](#项目信息)
  - [其他可能可用的项目](#其他可能可用的项目)
  - [数据下载](#数据下载)
    - [Reference](#reference)
    - [Illumina](#illumina)
  - [采用的倍数因子值](#采用的倍数因子值)
  - [Symlink](#symlink)
  - [Run](#run)
  - [Pack and clean](#pack-and-clean)


## 基本信息

* Genome: GCF_000188115.4, SL3.0, 828.349 Mb
* Chloroplast: [NC_007898](https://www.ncbi.nlm.nih.gov/nuccore/NC_007898), 155461 bp
* Mitochondrion: [NC_035963](https://www.ncbi.nlm.nih.gov/nuccore/NC_035963), 446257 bp


## 项目信息

* [PRJNA259308](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA259308)

* <https://www.nature.com/articles/ng.3117>


## 其他可能可用的项目


https://solgenomics.net/projects/varitome

https://www.ncbi.nlm.nih.gov/bioproject/PRJNA454805


## 数据下载

### Reference

```shell script
mkdir -p ~/data/plastid/360/genome
cd ~/data/plastid/360/genome

for ACCESSION in "NC_007898" "NC_035963"; do
    URL=$(printf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=%s&id=%s&retmode=text" "fasta" "${ACCESSION}");
    curl $URL -o ${ACCESSION}.fa
done

TAB=$'\t'
cat <<EOF > replace.tsv
NC_007898${TAB}Pt
NC_035963${TAB}Mt
EOF

cat NC_007898.fa NC_035963.fa |
    faops filter -s stdin stdout |
    faops replace stdin replace.tsv stdout |
    faops order stdin <(echo Pt; echo Mt) genome.fa

# bowtie2 index
bowtie2-build --threads 20 genome.fa genome.fa

# chr.sizes
faops size genome.fa > chr.sizes

```

### Illumina


* Download `Metadata` from NCBI SRA Run Selector via a web browser
  * https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA259308
  * Save it to `SraRunTable.txt`

* "Supplementary Table 1" of <https://www.nature.com/articles/ng.3117> provides detailed info

```shell script
mkdir -p ~/data/plastid/360/ena
cd ~/data/plastid/360/ena

cat SraRunTable.txt |
    mlr --icsv --otsv cat \
    > SraRunTable.tsv

cat SraRunTable.tsv |
    tsv-select -H -f Experiment,"Sample\ Name",Population,collected_by,Bases |
    tsv-filter -H \
        --not-blank collected_by \
        --istr-ne collected_by:missing \
        --ge Bases:4000000000 \
        --le Bases:10000000000 |
    sed '1 s/^/#/' |
    keep-header -- tsv-sort -k2,2 -k5,5nr |
    tsv-uniq -H -f "Sample\ Name" --max 1 |
    tsv-uniq -H -f Population,collected_by --max 5 | #wc -l
    cut -f 1-4 |
    mlr --itsv --ocsv cat \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

sed -i".bak" "s/ftp:/http:/" ena_info.ftp.txt # ftp server busy
aria2c -j 4 -x 4 -s 2 --file-allocation=none -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name   | srx       | platform | layout | ilength | srr        | spot     | base  |
|:-------|:----------|:---------|:-------|:--------|:-----------|:---------|:------|
| TS-1   | SRX698594 | ILLUMINA | PAIRED |         | SRR1572452 | 34236150 | 4.78G |
| TS-10  | SRX698603 | ILLUMINA | PAIRED |         | SRR1572461 | 38680247 | 7.2G  |
| TS-100 | SRX698641 | ILLUMINA | PAIRED |         | SRR1572499 | 47200811 | 8.79G |
| TS-102 | SRX698644 | ILLUMINA | PAIRED |         | SRR1572502 | 30800219 | 5.74G |
| TS-104 | SRX698647 | ILLUMINA | PAIRED |         | SRR1572505 | 37228653 | 6.93G |
| TS-105 | SRX698503 | ILLUMINA | PAIRED |         | SRR1572361 | 42135898 | 7.85G |
| TS-106 | SRX698504 | ILLUMINA | PAIRED |         | SRR1572362 | 35746785 | 6.66G |
| TS-107 | SRX698505 | ILLUMINA | PAIRED |         | SRR1572363 | 40467328 | 7.54G |
| TS-108 | SRX698648 | ILLUMINA | PAIRED |         | SRR1572506 | 23184310 | 4.32G |
| TS-110 | SRX698649 | ILLUMINA | PAIRED |         | SRR1572507 | 33544567 | 6.25G |
| TS-112 | SRX698652 | ILLUMINA | PAIRED |         | SRR1572510 | 47680880 | 8.88G |
| TS-114 | SRX698654 | ILLUMINA | PAIRED |         | SRR1572512 | 27504680 | 5.12G |
| TS-116 | SRX698507 | ILLUMINA | PAIRED |         | SRR1572365 | 33299923 | 6.2G  |
| TS-117 | SRX698656 | ILLUMINA | PAIRED |         | SRR1572514 | 28610182 | 5.33G |
| TS-118 | SRX698508 | ILLUMINA | PAIRED |         | SRR1572366 | 31511933 | 5.87G |
| TS-119 | SRX698511 | ILLUMINA | PAIRED |         | SRR1572369 | 27822589 | 5.18G |
| TS-121 | SRX698658 | ILLUMINA | PAIRED |         | SRR1572516 | 31601728 | 5.89G |
| TS-122 | SRX698660 | ILLUMINA | PAIRED |         | SRR1572518 | 35012707 | 6.52G |
| TS-123 | SRX698389 | ILLUMINA | PAIRED |         | SRR1572247 | 27018914 | 5.03G |
| TS-124 | SRX698390 | ILLUMINA | PAIRED |         | SRR1572248 | 35784650 | 6.67G |
| TS-126 | SRX698664 | ILLUMINA | PAIRED |         | SRR1572522 | 33475492 | 6.24G |
| TS-127 | SRX698666 | ILLUMINA | PAIRED |         | SRR1572524 | 26260820 | 4.89G |
| TS-129 | SRX698514 | ILLUMINA | PAIRED |         | SRR1572372 | 27916425 | 5.2G  |
| TS-130 | SRX698668 | ILLUMINA | PAIRED |         | SRR1572526 | 35282813 | 6.57G |
| TS-131 | SRX698515 | ILLUMINA | PAIRED |         | SRR1572373 | 25573376 | 4.76G |
| TS-132 | SRX698669 | ILLUMINA | PAIRED |         | SRR1572527 | 30633839 | 5.71G |
| TS-134 | SRX698516 | ILLUMINA | PAIRED |         | SRR1572374 | 39931892 | 7.44G |
| TS-135 | SRX698672 | ILLUMINA | PAIRED |         | SRR1572530 | 31948053 | 5.95G |
| TS-136 | SRX698673 | ILLUMINA | PAIRED |         | SRR1572531 | 27129484 | 5.05G |
| TS-137 | SRX698674 | ILLUMINA | PAIRED |         | SRR1572532 | 46773573 | 8.71G |
| TS-138 | SRX698517 | ILLUMINA | PAIRED |         | SRR1572375 | 36597411 | 6.82G |
| TS-139 | SRX698675 | ILLUMINA | PAIRED |         | SRR1572533 | 42090663 | 7.84G |
| TS-14  | SRX697461 | ILLUMINA | PAIRED | 500     | SRR1571027 | 32740751 | 4.57G |
| TS-144 | SRX698391 | ILLUMINA | PAIRED |         | SRR1572249 | 35211925 | 6.56G |
| TS-145 | SRX698392 | ILLUMINA | PAIRED |         | SRR1572250 | 27499487 | 5.12G |
| TS-146 | SRX698827 | ILLUMINA | PAIRED |         | SRR1572685 | 36211186 | 6.74G |
| TS-148 | SRX698518 | ILLUMINA | PAIRED |         | SRR1572376 | 28454608 | 5.3G  |
| TS-149 | SRX698519 | ILLUMINA | PAIRED |         | SRR1572377 | 35612333 | 6.63G |
| TS-15  | SRX698375 | ILLUMINA | PAIRED |         | SRR1572233 | 26058243 | 4.85G |
| TS-150 | SRX698683 | ILLUMINA | PAIRED |         | SRR1572541 | 35195197 | 6.56G |
| TS-152 | SRX698687 | ILLUMINA | PAIRED |         | SRR1572545 | 27080451 | 5.04G |
| TS-153 | SRX698688 | ILLUMINA | PAIRED |         | SRR1572546 | 31807728 | 5.92G |
| TS-154 | SRX698520 | ILLUMINA | PAIRED |         | SRR1572378 | 33629478 | 6.26G |
| TS-155 | SRX698690 | ILLUMINA | PAIRED |         | SRR1572548 | 35891255 | 6.69G |
| TS-156 | SRX698393 | ILLUMINA | PAIRED |         | SRR1572251 | 37512030 | 6.99G |
| TS-157 | SRX698691 | ILLUMINA | PAIRED |         | SRR1572549 | 36780978 | 6.85G |
| TS-158 | SRX698521 | ILLUMINA | PAIRED |         | SRR1572379 | 49693204 | 9.26G |
| TS-16  | SRX698376 | ILLUMINA | PAIRED |         | SRR1572234 | 32503044 | 4.54G |
| TS-162 | SRX698695 | ILLUMINA | PAIRED |         | SRR1572553 | 32461066 | 6.05G |
| TS-163 | SRX698696 | ILLUMINA | PAIRED |         | SRR1572554 | 23197653 | 4.32G |
| TS-164 | SRX698395 | ILLUMINA | PAIRED |         | SRR1572253 | 45371126 | 8.45G |
| TS-165 | SRX698522 | ILLUMINA | PAIRED |         | SRR1572380 | 33131989 | 6.17G |
| TS-166 | SRX698698 | ILLUMINA | PAIRED |         | SRR1572556 | 25066492 | 4.67G |
| TS-167 | SRX698700 | ILLUMINA | PAIRED |         | SRR1572558 | 31831754 | 5.93G |
| TS-17  | SRX698377 | ILLUMINA | PAIRED |         | SRR1572235 | 31988044 | 4.47G |
| TS-172 | SRX698706 | ILLUMINA | PAIRED |         | SRR1572564 | 38529623 | 7.18G |
| TS-176 | SRX698710 | ILLUMINA | PAIRED |         | SRR1572568 | 44434358 | 8.28G |
| TS-178 | SRX698712 | ILLUMINA | PAIRED |         | SRR1572570 | 41113580 | 7.66G |
| TS-181 | SRX698524 | ILLUMINA | PAIRED |         | SRR1572382 | 42926816 | 8G    |
| TS-183 | SRX698716 | ILLUMINA | PAIRED |         | SRR1572574 | 43651247 | 8.13G |
| TS-184 | SRX698717 | ILLUMINA | PAIRED |         | SRR1572575 | 27650944 | 5.15G |
| TS-185 | SRX698718 | ILLUMINA | PAIRED |         | SRR1572576 | 65237939 | 6.08G |
| TS-189 | SRX698526 | ILLUMINA | PAIRED |         | SRR1572384 | 33643947 | 6.27G |
| TS-197 | SRX698731 | ILLUMINA | PAIRED |         | SRR1572589 | 38228724 | 7.12G |
| TS-199 | SRX698830 | ILLUMINA | PAIRED |         | SRR1572688 | 44574679 | 8.3G  |
| TS-2   | SRX698595 | ILLUMINA | PAIRED |         | SRR1572453 | 31251247 | 4.37G |
| TS-201 | SRX698735 | ILLUMINA | PAIRED |         | SRR1572593 | 29086923 | 5.42G |
| TS-202 | SRX698527 | ILLUMINA | PAIRED |         | SRR1572385 | 35716706 | 6.65G |
| TS-207 | SRX698831 | ILLUMINA | PAIRED |         | SRR1572689 | 36657289 | 6.83G |
| TS-208 | SRX698828 | ILLUMINA | PAIRED |         | SRR1572686 | 22964915 | 4.28G |
| TS-209 | SRX698529 | ILLUMINA | PAIRED |         | SRR1572387 | 31169855 | 5.81G |
| TS-213 | SRX698531 | ILLUMINA | PAIRED |         | SRR1572389 | 32071916 | 5.97G |
| TS-214 | SRX698743 | ILLUMINA | PAIRED |         | SRR1572601 | 35822838 | 6.67G |
| TS-215 | SRX698744 | ILLUMINA | PAIRED |         | SRR1572602 | 46335168 | 8.63G |
| TS-217 | SRX698832 | ILLUMINA | PAIRED |         | SRR1572690 | 20336072 | 3.79G |
| TS-218 | SRX698747 | ILLUMINA | PAIRED |         | SRR1572605 | 22500992 | 4.19G |
| TS-219 | SRX698533 | ILLUMINA | PAIRED |         | SRR1572391 | 30753811 | 5.73G |
| TS-221 | SRX698534 | ILLUMINA | PAIRED |         | SRR1572392 | 28897252 | 5.38G |
| TS-223 | SRX698535 | ILLUMINA | PAIRED |         | SRR1572393 | 38337682 | 7.14G |
| TS-224 | SRX698750 | ILLUMINA | PAIRED |         | SRR1572608 | 24313414 | 4.53G |
| TS-225 | SRX698752 | ILLUMINA | PAIRED |         | SRR1572610 | 43668194 | 8.13G |
| TS-226 | SRX698753 | ILLUMINA | PAIRED |         | SRR1572611 | 31274973 | 5.83G |
| TS-229 | SRX698539 | ILLUMINA | PAIRED |         | SRR1572397 | 33539226 | 6.25G |
| TS-233 | SRX698543 | ILLUMINA | PAIRED |         | SRR1572401 | 30354597 | 5.65G |
| TS-235 | SRX698759 | ILLUMINA | PAIRED |         | SRR1572617 | 31856196 | 5.93G |
| TS-236 | SRX698760 | ILLUMINA | PAIRED |         | SRR1572618 | 23104524 | 4.3G  |
| TS-237 | SRX698761 | ILLUMINA | PAIRED |         | SRR1572619 | 42184062 | 7.86G |
| TS-238 | SRX698544 | ILLUMINA | PAIRED |         | SRR1572402 | 34242485 | 6.38G |
| TS-240 | SRX698546 | ILLUMINA | PAIRED |         | SRR1572404 | 38918481 | 7.25G |
| TS-241 | SRX698764 | ILLUMINA | PAIRED |         | SRR1572622 | 33138357 | 6.17G |
| TS-245 | SRX698766 | ILLUMINA | PAIRED |         | SRR1572624 | 35122043 | 6.54G |
| TS-249 | SRX698768 | ILLUMINA | PAIRED |         | SRR1572626 | 29652192 | 5.52G |
| TS-250 | SRX698550 | ILLUMINA | PAIRED |         | SRR1572408 | 38807052 | 7.23G |
| TS-251 | SRX698769 | ILLUMINA | PAIRED |         | SRR1572627 | 31732515 | 5.91G |
| TS-252 | SRX698551 | ILLUMINA | PAIRED |         | SRR1572409 | 32841436 | 6.12G |
| TS-255 | SRX698771 | ILLUMINA | PAIRED |         | SRR1572629 | 38031941 | 7.08G |
| TS-257 | SRX698554 | ILLUMINA | PAIRED |         | SRR1572412 | 32934461 | 6.13G |
| TS-26  | SRX698440 | ILLUMINA | PAIRED |         | SRR1572298 | 30246593 | 5.62G |
| TS-260 | SRX698556 | ILLUMINA | PAIRED |         | SRR1572414 | 29988136 | 5.59G |
| TS-261 | SRX698774 | ILLUMINA | PAIRED |         | SRR1572632 | 32007578 | 5.96G |
| TS-262 | SRX698558 | ILLUMINA | PAIRED |         | SRR1572416 | 30435032 | 5.67G |
| TS-267 | SRX698401 | ILLUMINA | PAIRED |         | SRR1572259 | 41415815 | 7.71G |
| TS-27  | SRX698442 | ILLUMINA | PAIRED |         | SRR1572300 | 36209787 | 6.73G |
| TS-270 | SRX698783 | ILLUMINA | PAIRED |         | SRR1572641 | 37507041 | 6.99G |
| TS-271 | SRX698559 | ILLUMINA | PAIRED |         | SRR1572417 | 31598011 | 5.89G |
| TS-274 | SRX698787 | ILLUMINA | PAIRED |         | SRR1572645 | 28891564 | 5.38G |
| TS-276 | SRX698791 | ILLUMINA | PAIRED |         | SRR1572649 | 44237422 | 8.24G |
| TS-282 | SRX698796 | ILLUMINA | PAIRED |         | SRR1572654 | 29847761 | 5.56G |
| TS-284 | SRX698566 | ILLUMINA | PAIRED |         | SRR1572424 | 33911423 | 6.32G |
| TS-286 | SRX698567 | ILLUMINA | PAIRED |         | SRR1572425 | 35742687 | 6.66G |
| TS-288 | SRX698798 | ILLUMINA | PAIRED |         | SRR1572656 | 30201506 | 5.63G |
| TS-289 | SRX698569 | ILLUMINA | PAIRED |         | SRR1572427 | 40452727 | 7.53G |
| TS-29  | SRX698447 | ILLUMINA | PAIRED |         | SRR1572305 | 31850843 | 5.92G |
| TS-294 | SRX698572 | ILLUMINA | PAIRED |         | SRR1572430 | 45100513 | 8.4G  |
| TS-297 | SRX698802 | ILLUMINA | PAIRED |         | SRR1572660 | 31789721 | 5.92G |
| TS-3   | SRX698596 | ILLUMINA | PAIRED |         | SRR1572454 | 34860050 | 4.87G |
| TS-30  | SRX698448 | ILLUMINA | PAIRED |         | SRR1572306 | 28175185 | 5.24G |
| TS-302 | SRX698581 | ILLUMINA | PAIRED |         | SRR1572439 | 35797300 | 6.67G |
| TS-303 | SRX698582 | ILLUMINA | PAIRED |         | SRR1572440 | 29000798 | 5.4G  |
| TS-31  | SRX698451 | ILLUMINA | PAIRED |         | SRR1572309 | 32839664 | 6.11G |
| TS-32  | SRX698453 | ILLUMINA | PAIRED |         | SRR1572311 | 35676098 | 6.63G |
| TS-36  | SRX698458 | ILLUMINA | PAIRED |         | SRR1572316 | 35922243 | 6.68G |
| TS-38  | SRX698463 | ILLUMINA | PAIRED |         | SRR1572321 | 33704928 | 6.27G |
| TS-400 | SRX698804 | ILLUMINA | PAIRED |         | SRR1572662 | 38635347 | 7.2G  |
| TS-401 | SRX698807 | ILLUMINA | PAIRED |         | SRR1572665 | 32941871 | 6.14G |
| TS-402 | SRX698835 | ILLUMINA | PAIRED |         | SRR1572693 | 31645378 | 5.89G |
| TS-403 | SRX698836 | ILLUMINA | PAIRED |         | SRR1572694 | 28149381 | 5.24G |
| TS-404 | SRX698837 | ILLUMINA | PAIRED |         | SRR1572695 | 34358030 | 6.4G  |
| TS-407 | SRX698839 | ILLUMINA | PAIRED |         | SRR1572697 | 34551270 | 6.44G |
| TS-408 | SRX698838 | ILLUMINA | PAIRED |         | SRR1572696 | 31750650 | 5.91G |
| TS-409 | SRX698808 | ILLUMINA | PAIRED |         | SRR1572666 | 39508515 | 7.36G |
| TS-41  | SRX698606 | ILLUMINA | PAIRED |         | SRR1572464 | 32535528 | 6.06G |
| TS-410 | SRX698405 | ILLUMINA | PAIRED |         | SRR1572263 | 28829002 | 5.37G |
| TS-432 | SRX698426 | ILLUMINA | PAIRED |         | SRR1572284 | 23514458 | 4.38G |
| TS-433 | SRX698427 | ILLUMINA | PAIRED |         | SRR1572285 | 23982839 | 4.47G |
| TS-48  | SRX698613 | ILLUMINA | PAIRED |         | SRR1572471 | 30333463 | 5.65G |
| TS-5   | SRX698598 | ILLUMINA | PAIRED |         | SRR1572456 | 34802765 | 4.86G |
| TS-59  | SRX698620 | ILLUMINA | PAIRED |         | SRR1572478 | 36746826 | 6.84G |
| TS-62  | SRX698475 | ILLUMINA | PAIRED |         | SRR1572333 | 33018690 | 6.15G |
| TS-64  | SRX698477 | ILLUMINA | PAIRED |         | SRR1572335 | 28718548 | 5.35G |
| TS-65  | SRX698478 | ILLUMINA | PAIRED |         | SRR1572336 | 26108405 | 4.86G |
| TS-67  | SRX698481 | ILLUMINA | PAIRED |         | SRR1572339 | 37273229 | 6.94G |
| TS-69  | SRX698624 | ILLUMINA | PAIRED |         | SRR1572482 | 28443142 | 5.3G  |
| TS-82  | SRX698631 | ILLUMINA | PAIRED |         | SRR1572489 | 29675046 | 5.53G |
| TS-83  | SRX698488 | ILLUMINA | PAIRED |         | SRR1572346 | 46431966 | 8.65G |
| TS-9   | SRX698602 | ILLUMINA | PAIRED |         | SRR1572460 | 31998972 | 4.47G |
| TS-91  | SRX698491 | ILLUMINA | PAIRED |         | SRR1572349 | 26571451 | 4.95G |
| TS-94  | SRX698493 | ILLUMINA | PAIRED |         | SRR1572351 | 48573204 | 9.05G |
| TS-97  | SRX698498 | ILLUMINA | PAIRED |         | SRR1572356 | 31539812 | 5.87G |
| TS-98  | SRX698500 | ILLUMINA | PAIRED |         | SRR1572358 | 31695055 | 5.9G  |


## Symlink

```shell script
cd ~/data/plastid/360/

export FOLD=2
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/h1706/chr.sizes |
        tsv-summarize --sum 2
)

cat ena/ena_info.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f name,srr,base |
    perl -nla -F'\t' -e '
        /^name/ and next;
        my $bases = $F[2];
        $bases =~ s/G$//;
        my $cutoff = $bases * 1000 * 1000 * 1000 / $ENV{GENOME_SIZE} * $ENV{FOLD};
        $cutoff = int $cutoff;
        print join qq(\t), ($F[0], $F[1], $cutoff);
    ' \
    > opts.tsv

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        mkdir -p {1}/1_genome
        pushd {1}/1_genome
    
        ln -fs ../../genome/genome.fa genome.fa
        popd
        
        mkdir -p {1}/2_illumina
        pushd {1}/2_illumina
        
        ln -fs ../../ena/{2}_1.fastq.gz R1.fq.gz
        ln -fs ../../ena/{2}_2.fastq.gz R2.fq.gz
        popd
    '

```

## Run

* Rsync to hpcc

```shell script
rsync -avP \
    ~/data/plastid/360/ \
    wangq@202.119.37.251:data/plastid/360

# rsync -avP wangq@202.119.37.251:data/plastid/360/ ~/data/plastid/360

```

```shell script
cd ~/data/plastid/360/

cat opts.tsv | head -n 10 | tail -n 10 |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        cd {1}
        
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
            --quorum \
            --merge \
            --ecphase "1 2 3" \
            \
            --cov "40 80 120 160 240 320" \
            --unitigger "superreads bcalm tadpole" \
            --splitp 100 \
            --statp 1 \
            --readl 100 \
            --uscale 50 \
            --lscale 5 \
            --redo \
            \
            --extend
            
        bsub -q mpi -n 24 -J "{1}" "
            bash 0_master.sh
        "
    '

```

## Pack and clean

```shell script
cd ~/data/plastid/360/

cat opts.tsv | head -n 30 |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        cd {1}
        
        if [ -d 4_down_sampling ]; then
            rm -fr 4_down_sampling
            rm -fr 6_down_sampling
        fi
    '

cat opts.tsv | head -n 30 | tail -n 10 |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '        
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        if [ ! -f {1}/7_merge_anchors/anchor.merge.fasta ]; then
            exit;
        fi

        if [ ! -d {1}/9_quast ]; then
            exit;
        fi

        echo "==> Clean {1}"
        bash {1}/0_cleanup.sh
        
        echo "==> Create {1}.tar.gz"

        tar -czvf {1}.tar.gz \
            {1}/2_illumina/fastqc \
            {1}/2_illumina/insert_size \
            {1}/2_illumina/kat \
            {1}/2_illumina/trim/Q25L60/pe.cor.fa.gz \
            {1}/2_illumina/trim/Q25L60/env.json \
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

cat opts.tsv | head -n 150 |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            if [ -d {1} ]; then
                echo "==> Remove {1}/"
                rm -fr {1}
            fi
        fi
    '

```

