# *Medicago truncatula* Hapmap Project

[TOC levels=1-3]: # ""

- [*Medicago truncatula* Hapmap Project](#medicago-truncatula-hapmap-project)
  - [基本信息](#基本信息)
  - [项目信息](#项目信息)
  - [其他可能可用的项目](#其他可能可用的项目)
  - [数据下载](#数据下载)
    - [Reference](#reference)
    - [Illumina](#illumina)


## 基本信息

* Genome: GCF_000219495.3, MedtrA17_4.0, 412.924 Mb
* Chloroplast: [NC_003119](https://www.ncbi.nlm.nih.gov/nuccore/NC_003119), 124033 bp
* Mitochondrion: [NC_029641](https://www.ncbi.nlm.nih.gov/nuccore/NC_029641), 271618 bp


## 项目信息

* [PRJNA256006](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA256006)

<http://www.medicagohapmap.org/>

> Briefly, 384 inbred lines spanning the range of Medicago diversity are being resequenced using
> Illumina next generation technology. This provides a foundation for discovering single nucleotide
> polymorphisms (SNPs), insertions/deletions (INDELs) and copy number variants (CNV) at very high
> resolution among the Medicago lines. Thirty of these lines have been deeply resequenced (20X
> coverage or more), while the remainder are sequenced at least 5X coverage. The resulting database
> of sequence variants establishes a basis for describing population structure and identifying
> genome segments with shared ancestry (haplotypes) - and thereby creating a long-term,
> community-accessible genome-wide association (GWA) mapping resource.


## 其他可能可用的项目

PRJNA170333


## 数据下载

### Reference


```shell script
mkdir -p ~/data/plastid/384/genome
cd ~/data/plastid/384/genome

for ACCESSION in "NC_003119" "NC_029641"; do
    URL=$(printf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=%s&id=%s&retmode=text" "fasta" "${ACCESSION}");
    curl $URL -o ${ACCESSION}.fa
done

TAB=$'\t'
cat <<EOF > replace.tsv
NC_003119${TAB}Pt
NC_029641${TAB}Mt
EOF

cat NC_003119.fa NC_029641.fa |
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
  * https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA256006
  * Save it to `SraRunTable.txt`

* http://www.medicagohapmap.org/hapmap/germplasm
  * Extract table via `pup`
  * Convert xls to csv via `excel`

```shell script
mkdir -p ~/data/plastid/384/ena
cd ~/data/plastid/384/ena

cat SraRunTable.txt |
    mlr --icsv --otsv cat \
    > SraRunTable.tsv

curl http://www.medicagohapmap.org/hapmap/germplasm |
    pup 'table#germplasmData' \
    > germplasm.xls

# Convert xls to csv via `excel`
cat germplasm.csv | wc -l
#338

cat germplasm.csv |
    head -n 11 |
    mlr --icsv --omd cat

cat SraRunTable.tsv |
    tsv-filter -H \
        --istr-in-fld "Instrument":'Illumina' \
        --istr-in-fld "DATASTORE\ filetype":'fastq' \
        --regex "Library\ Name":'^HM' \
        --ge Bases:4000000000 \
        --le Bases:10000000000 \
        --ge AvgSpotLen:90 |
    tsv-filter -H --istr-not-in-fld "Library\ Name":'Mate' |
    tsv-filter -H --istr-not-in-fld "Library\ Name":'Nex' |
    tsv-filter -H --istr-not-in-fld "Library\ Name":'MP' |
    perl -nla -F'\t' -e '
        $F[14] =~ s/^HM_(\d+)/HM$1/;
        $F[14] =~ s/^(HM\d+).*/$1/;
        print join qq(\t), @F;
    ' \
    > corrected.tsv

cat germplasm.csv |
    perl -p -e 's/\r\n/\n/g; s/ ,/,/g' |
    mlr --icsv --otsv cat |
    tsv-join -H --data-fields "ID" --key-fields "Library\ Name" \
        -f corrected.tsv \
        --append-fields AvgSpotLen,Instrument,Bases,Experiment |
    tsv-select -H -f Experiment,ID,"Country\ of\ Origin" |
    mlr --itsv --ocsv cat |
    sed 1d \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -j 4 -x 4 -s 2 --file-allocation=none -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| ID    | Line    | Population of Origin | Country of Origin | Category | Seeds From       | Status    |
|:------|:--------|:---------------------|:------------------|:---------|:-----------------|:----------|
| HM001 | L000163 | SA22322              | Syria             | CC8      | INRA-Montpellier | Processed |
| HM002 | L000174 | SA28064              | Cyprus            | CC8      | INRA-Montpellier | Processed |
| HM003 | L000544 | ESP105-L             | Spain             | CC8      | INRA-Montpellier | Processed |
| HM004 | L000736 | DZA045-6             | Algeria           | CC8      | INRA-Montpellier | Processed |
| HM005 | L000734 | DZA315-16            | Algeria           | CC8      | INRA-Montpellier | Processed |
| HM006 | L000530 | F83005-5             | France            | CC8      | INRA-Montpellier | Processed |
| HM007 | L000651 | Salses71B            | France            | CC8      | INRA-Montpellier | Processed |
| HM008 | L000368 | DZA012-J             | Algeria           | CC8      | INRA-Montpellier | Processed |
| HM009 | L000555 | GRC020-B             | Greece            | CC16     | INRA-Montpellier | Processed |
| HM010 | L000154 | SA24714              | Italy             | CC16     | INRA-Montpellier | Processed |


| name  | srx       | platform | layout | ilength | srr        | spot     | base  |
|:------|:----------|:---------|:-------|:--------|:-----------|:---------|:------|
| HM014 | SRX375910 | ILLUMINA | PAIRED | 293     | SRR1034070 | 29010895 | 4.86G |
| HM016 | SRX375953 | ILLUMINA | PAIRED | 299     | SRR1034113 | 22543648 | 3.78G |
| HM021 | SRX376029 | ILLUMINA | PAIRED | 358     | SRR1034189 | 27466998 | 4.6G  |
| HM023 | SRX376011 | ILLUMINA | PAIRED | 405     | SRR1034171 | 25507163 | 4.28G |
| HM024 | SRX376094 | ILLUMINA | PAIRED | 268     | SRR1034254 | 22320687 | 3.74G |
| HM029 | SRX375983 | ILLUMINA | PAIRED | 312     | SRR1034143 | 23073809 | 3.87G |
| HM102 | SRX375988 | ILLUMINA | PAIRED | 361     | SRR1034148 | 23768646 | 3.98G |
| HM106 | SRX376158 | ILLUMINA | PAIRED | 229     | SRR1034318 | 23504842 | 3.94G |
| HM107 | SRX376163 | ILLUMINA | PAIRED | 246     | SRR1034323 | 22413221 | 3.76G |
| HM146 | SRX376196 | ILLUMINA | PAIRED | 258     | SRR1034356 | 35653183 | 5.98G |
| HM150 | SRX376197 | ILLUMINA | PAIRED | 241     | SRR1034357 | 34880743 | 5.85G |
| HM155 | SRX376187 | ILLUMINA | PAIRED | 232     | SRR1034347 | 31231373 | 5.24G |
| HM156 | SRX376199 | ILLUMINA | PAIRED | 243     | SRR1034359 | 35047508 | 5.88G |
| HM157 | SRX376188 | ILLUMINA | PAIRED | 251     | SRR1034348 | 32434856 | 5.44G |
| HM160 | SRX376204 | ILLUMINA | PAIRED | 239     | SRR1034364 | 33712431 | 5.65G |
| HM161 | SRX376202 | ILLUMINA | PAIRED | 235     | SRR1034362 | 31891151 | 5.35G |
| HM162 | SRX376203 | ILLUMINA | PAIRED | 224     | SRR1034363 | 29198476 | 4.89G |
| HM163 | SRX376205 | ILLUMINA | PAIRED | 169     | SRR1034365 | 32748569 | 5.49G |
| HM164 | SRX376206 | ILLUMINA | PAIRED | 175     | SRR1034366 | 31184827 | 5.23G |
| HM165 | SRX376207 | ILLUMINA | PAIRED | 180     | SRR1034367 | 33172082 | 5.56G |
| HM166 | SRX376208 | ILLUMINA | PAIRED | 198     | SRR1034368 | 29551522 | 4.95G |
| HM167 | SRX376209 | ILLUMINA | PAIRED | 209     | SRR1034369 | 28502982 | 4.78G |
| HM168 | SRX376210 | ILLUMINA | PAIRED | 221     | SRR1034370 | 30190885 | 5.06G |
| HM169 | SRX376211 | ILLUMINA | PAIRED | 179     | SRR1034371 | 29657558 | 4.97G |
| HM170 | SRX376212 | ILLUMINA | PAIRED | 193     | SRR1034372 | 29076479 | 4.87G |
| HM172 | SRX376213 | ILLUMINA | PAIRED | 207     | SRR1034373 | 27983365 | 4.69G |
| HM175 | SRX376215 | ILLUMINA | PAIRED | 211     | SRR1034375 | 27182645 | 4.56G |
| HM176 | SRX376216 | ILLUMINA | PAIRED | 228     | SRR1034376 | 27582074 | 4.62G |
| HM178 | SRX376218 | ILLUMINA | PAIRED | 216     | SRR1034378 | 27543515 | 4.62G |
| HM179 | SRX376219 | ILLUMINA | PAIRED | 227     | SRR1034379 | 29050894 | 4.87G |
| HM180 | SRX376200 | ILLUMINA | PAIRED | 224     | SRR1034360 | 33094981 | 5.55G |
| HM181 | SRX376201 | ILLUMINA | PAIRED | 202     | SRR1034361 | 33340119 | 5.59G |
| HM183 | SRX376220 | ILLUMINA | PAIRED | 256     | SRR1034380 | 28752670 | 4.82G |
| HM184 | SRX376221 | ILLUMINA | PAIRED | 273     | SRR1034381 | 25530259 | 4.28G |
| HM185 | SRX376222 | ILLUMINA | PAIRED | 246     | SRR1034382 | 31875679 | 5.34G |
| HM186 | SRX376223 | ILLUMINA | PAIRED | 246     | SRR1034383 | 31032893 | 5.2G  |
| HM187 | SRX376224 | ILLUMINA | PAIRED | 207     | SRR1034384 | 30913811 | 5.18G |

