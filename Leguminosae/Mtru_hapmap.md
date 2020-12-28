# *Medicago truncatula* Hapmap Project

[TOC levels=1-3]: # ""

- [*Medicago truncatula* Hapmap Project](#medicago-truncatula-hapmap-project)
  - [基本信息](#基本信息)
  - [项目信息](#项目信息)
  - [其他可能可用的项目](#其他可能可用的项目)
  - [数据下载](#数据下载)
    - [Reference](#reference)
    - [Illumina](#illumina)
  - [Symlink](#symlink)
  - [Run](#run)
  - [Pack and clean](#pack-and-clean)


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

# Convert xls(html) to csv via `excel`
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
        --ge Bases:2000000000 \
        --le Bases:10000000000 \
        --ge AvgSpotLen:90 |
    tsv-filter -H --istr-not-in-fld "Library\ Name":'Mate' |
    tsv-filter -H --istr-not-in-fld "Library\ Name":'Nex' |
    tsv-filter -H --istr-not-in-fld "Library\ Name":'MP' |
    perl -nla -F'\t' -e '
        $F[14] =~ s/^HM_(\d+)/HM$1/;
        $F[14] =~ s/^(HM\d+).*/$1/;
        print join qq(\t), @F;
    ' |
    keep-header -- tsv-sort -k4,4nr -k14,14 | # grep HM001 # sort by bases
    tsv-uniq -H \
        -f "Library\ Name" --max 1 \
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

#sed -i".bak" "s/ftp:/http:/" ena_info.ftp.txt # when ftp server is busy
aria2c -j 4 -x 4 -s 1 --file-allocation=none -c -i ena_info.ftp.txt

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
| HM001 | SRX375894 | ILLUMINA | PAIRED | 257     | SRR1034054 | 18738482 | 3.14G |
| HM002 | SRX375896 | ILLUMINA | PAIRED | 219     | SRR1034056 | 20086982 | 3.37G |
| HM003 | SRX375917 | ILLUMINA | PAIRED | 283     | SRR1034077 | 18599526 | 3.12G |
| HM004 | SRX375905 | ILLUMINA | PAIRED | 273     | SRR1034065 | 19031676 | 3.19G |
| HM005 | SRX375923 | ILLUMINA | PAIRED | 247     | SRR1034083 | 18903805 | 3.17G |
| HM006 | SRX376009 | ILLUMINA | PAIRED | 359     | SRR1034169 | 18531405 | 3.11G |
| HM007 | SRX375930 | ILLUMINA | PAIRED | 229     | SRR1034090 | 22213827 | 3.72G |
| HM008 | SRX375937 | ILLUMINA | PAIRED | 214     | SRR1034097 | 21766454 | 3.65G |
| HM009 | SRX375970 | ILLUMINA | PAIRED | 262     | SRR1034130 | 16090912 | 2.7G  |
| HM010 | SRX376080 | ILLUMINA | PAIRED | 241     | SRR1034240 | 14071825 | 2.36G |
| HM011 | SRX375943 | ILLUMINA | PAIRED | 228     | SRR1034103 | 15624226 | 2.62G |
| HM012 | SRX375962 | ILLUMINA | PAIRED | 351     | SRR1034122 | 16190747 | 2.71G |
| HM013 | SRX375998 | ILLUMINA | PAIRED | 365     | SRR1034158 | 18067454 | 3.03G |
| HM014 | SRX375910 | ILLUMINA | PAIRED | 293     | SRR1034070 | 29010895 | 4.86G |
| HM015 | SRX375948 | ILLUMINA | PAIRED | 298     | SRR1034108 | 20171851 | 3.38G |
| HM016 | SRX375953 | ILLUMINA | PAIRED | 299     | SRR1034113 | 22543648 | 3.78G |
| HM017 | SRX376069 | ILLUMINA | PAIRED | 256     | SRR1034229 | 18409555 | 3.09G |
| HM018 | SRX376025 | ILLUMINA | PAIRED | 339     | SRR1034185 | 21545311 | 3.61G |
| HM019 | SRX375978 | ILLUMINA | PAIRED | 341     | SRR1034138 | 17367340 | 2.91G |
| HM020 | SRX376071 | ILLUMINA | PAIRED | 241     | SRR1034231 | 17181113 | 2.88G |
| HM021 | SRX376029 | ILLUMINA | PAIRED | 358     | SRR1034189 | 27466998 | 4.6G  |
| HM022 | SRX376075 | ILLUMINA | PAIRED | 266     | SRR1034235 | 18096959 | 3.03G |
| HM023 | SRX376011 | ILLUMINA | PAIRED | 405     | SRR1034171 | 25507163 | 4.28G |
| HM024 | SRX376094 | ILLUMINA | PAIRED | 268     | SRR1034254 | 22320687 | 3.74G |
| HM025 | SRX376095 | ILLUMINA | PAIRED | 248     | SRR1034255 | 20793012 | 3.49G |
| HM026 | SRX375990 | ILLUMINA | PAIRED | 315     | SRR1034150 | 12887454 | 2.16G |
| HM027 | SRX376004 | ILLUMINA | PAIRED | 292     | SRR1034164 | 20996668 | 3.52G |
| HM028 | SRX375968 | ILLUMINA | PAIRED | 441     | SRR1034128 | 13428840 | 2.25G |
| HM029 | SRX375983 | ILLUMINA | PAIRED | 312     | SRR1034143 | 23073809 | 3.87G |
| HM030 | SRX376036 | ILLUMINA | PAIRED | 314     | SRR1034196 | 19120675 | 3.21G |
| HM031 | SRX376039 | ILLUMINA | PAIRED | 435     | SRR1034199 | 19618533 | 3.29G |
| HM033 | SRX376125 | ILLUMINA | PAIRED | 387     | SRR1034285 | 20847474 | 3.49G |
| HM034 | SRX376131 | ILLUMINA | PAIRED | 377     | SRR1034291 | 16893530 | 2.83G |
| HM035 | SRX376130 | ILLUMINA | PAIRED | 398     | SRR1034290 | 20384946 | 3.42G |
| HM036 | SRX376056 | ILLUMINA | PAIRED | 456     | SRR1034216 | 13603172 | 2.28G |
| HM037 | SRX376129 | ILLUMINA | PAIRED | 362     | SRR1034289 | 18805564 | 3.15G |
| HM038 | SRX376040 | ILLUMINA | PAIRED | 480     | SRR1034200 | 18289823 | 3.07G |
| HM039 | SRX376057 | ILLUMINA | PAIRED | 472     | SRR1034217 | 16305059 | 2.73G |
| HM040 | SRX376042 | ILLUMINA | PAIRED | 419     | SRR1034202 | 20313648 | 3.41G |
| HM041 | SRX376043 | ILLUMINA | PAIRED | 476     | SRR1034203 | 19764690 | 3.31G |
| HM043 | SRX376142 | ILLUMINA | PAIRED | 363     | SRR1034302 | 16587057 | 2.78G |
| HM044 | SRX376059 | ILLUMINA | PAIRED | 369     | SRR1034219 | 19249373 | 3.23G |
| HM045 | SRX376061 | ILLUMINA | PAIRED | 362     | SRR1034221 | 21367867 | 3.58G |
| HM046 | SRX376127 | ILLUMINA | PAIRED | 353     | SRR1034287 | 16205869 | 2.72G |
| HM047 | SRX376137 | ILLUMINA | PAIRED | 340     | SRR1034297 | 19864858 | 3.33G |
| HM048 | SRX376126 | ILLUMINA | PAIRED | 338     | SRR1034286 | 19427996 | 3.26G |
| HM049 | SRX376136 | ILLUMINA | PAIRED | 369     | SRR1034296 | 21152110 | 3.55G |
| HM050 | SRX376132 | ILLUMINA | PAIRED | 373     | SRR1034292 | 20307663 | 3.4G  |
| HM051 | SRX376143 | ILLUMINA | PAIRED | 379     | SRR1034303 | 19059567 | 3.2G  |
| HM052 | SRX376144 | ILLUMINA | PAIRED | 235     | SRR1034304 | 17538344 | 2.94G |
| HM053 | SRX376148 | ILLUMINA | PAIRED | 236     | SRR1034308 | 17262697 | 2.89G |
| HM054 | SRX376152 | ILLUMINA | PAIRED | 252     | SRR1034312 | 14403878 | 2.41G |
| HM057 | SRX376156 | ILLUMINA | PAIRED | 271     | SRR1034316 | 19172878 | 3.21G |
| HM058 | SRX376044 | ILLUMINA | PAIRED | 522     | SRR1034204 | 16647508 | 2.79G |
| HM059 | SRX376164 | ILLUMINA | PAIRED | 268     | SRR1034324 | 17072204 | 2.86G |
| HM060 | SRX376159 | ILLUMINA | PAIRED | 245     | SRR1034319 | 15942892 | 2.67G |
| HM061 | SRX376145 | ILLUMINA | PAIRED | 241     | SRR1034305 | 13649083 | 2.29G |
| HM062 | SRX376062 | ILLUMINA | PAIRED | 376     | SRR1034222 | 17025579 | 2.85G |
| HM063 | SRX376045 | ILLUMINA | PAIRED | 410     | SRR1034205 | 14997697 | 2.51G |
| HM065 | SRX376063 | ILLUMINA | PAIRED | 494     | SRR1034223 | 12409127 | 2.08G |
| HM068 | SRX376134 | ILLUMINA | PAIRED | 386     | SRR1034294 | 17315467 | 2.9G  |
| HM070 | SRX376149 | ILLUMINA | PAIRED | 277     | SRR1034309 | 17546142 | 2.94G |
| HM071 | SRX376065 | ILLUMINA | PAIRED | 492     | SRR1034225 | 15538352 | 2.6G  |
| HM072 | SRX376128 | ILLUMINA | PAIRED | 377     | SRR1034288 | 16598239 | 2.78G |
| HM073 | SRX376124 | ILLUMINA | PAIRED | 366     | SRR1034284 | 16515292 | 2.77G |
| HM075 | SRX376048 | ILLUMINA | PAIRED | 425     | SRR1034208 | 17207913 | 2.88G |
| HM076 | SRX376141 | ILLUMINA | PAIRED | 374     | SRR1034301 | 17627304 | 2.96G |
| HM077 | SRX376123 | ILLUMINA | PAIRED | 387     | SRR1034283 | 14680655 | 2.46G |
| HM080 | SRX376168 | ILLUMINA | PAIRED | 234     | SRR1034328 | 16582882 | 2.78G |
| HM081 | SRX376157 | ILLUMINA | PAIRED | 251     | SRR1034317 | 16027834 | 2.69G |
| HM082 | SRX376162 | ILLUMINA | PAIRED | 282     | SRR1034322 | 19933264 | 3.34G |
| HM085 | SRX376147 | ILLUMINA | PAIRED | 207     | SRR1034307 | 14042233 | 2.35G |
| HM087 | SRX376155 | ILLUMINA | PAIRED | 247     | SRR1034315 | 16113807 | 2.7G  |
| HM088 | SRX376150 | ILLUMINA | PAIRED | 251     | SRR1034310 | 17881903 | 3G    |
| HM093 | SRX376154 | ILLUMINA | PAIRED | 244     | SRR1034314 | 20167574 | 3.38G |
| HM102 | SRX375988 | ILLUMINA | PAIRED | 361     | SRR1034148 | 23768646 | 3.98G |
| HM105 | SRX376140 | ILLUMINA | PAIRED | 317     | SRR1034300 | 18920661 | 3.17G |
| HM106 | SRX376158 | ILLUMINA | PAIRED | 229     | SRR1034318 | 23504842 | 3.94G |
| HM107 | SRX376163 | ILLUMINA | PAIRED | 246     | SRR1034323 | 22413221 | 3.76G |
| HM108 | SRX376165 | ILLUMINA | PAIRED | 226     | SRR1034325 | 21834727 | 3.66G |
| HM109 | SRX376167 | ILLUMINA | PAIRED | 277     | SRR1034327 | 16374331 | 2.74G |
| HM111 | SRX376135 | ILLUMINA | PAIRED | 344     | SRR1034295 | 16638222 | 2.79G |
| HM112 | SRX376146 | ILLUMINA | PAIRED | 276     | SRR1034306 | 16506517 | 2.77G |
| HM114 | SRX376194 | ILLUMINA | PAIRED | 257     | SRR1034354 | 17922954 | 3G    |
| HM118 | SRX376183 | ILLUMINA | PAIRED | 262     | SRR1034343 | 14819910 | 2.48G |
| HM119 | SRX376182 | ILLUMINA | PAIRED | 221     | SRR1034342 | 16233468 | 2.72G |
| HM126 | SRX376181 | ILLUMINA | PAIRED | 168     | SRR1034341 | 16939866 | 2.84G |
| HM130 | SRX376173 | ILLUMINA | PAIRED | 210     | SRR1034333 | 14593193 | 2.45G |
| HM131 | SRX376189 | ILLUMINA | PAIRED | 257     | SRR1034349 | 17856111 | 2.99G |
| HM133 | SRX376190 | ILLUMINA | PAIRED | 234     | SRR1034350 | 17087795 | 2.86G |
| HM138 | SRX376193 | ILLUMINA | PAIRED | 269     | SRR1034353 | 16094171 | 2.7G  |
| HM141 | SRX376177 | ILLUMINA | PAIRED | 213     | SRR1034337 | 14365982 | 2.41G |
| HM146 | SRX376196 | ILLUMINA | PAIRED | 258     | SRR1034356 | 35653183 | 5.98G |
| HM147 | SRX376180 | ILLUMINA | PAIRED | 202     | SRR1034340 | 18938657 | 3.17G |
| HM148 | SRX376172 | ILLUMINA | PAIRED | 233     | SRR1034332 | 20532423 | 3.44G |
| HM149 | SRX376169 | ILLUMINA | PAIRED | 207     | SRR1034329 | 20952891 | 3.51G |
| HM150 | SRX376197 | ILLUMINA | PAIRED | 241     | SRR1034357 | 34880743 | 5.85G |
| HM151 | SRX376184 | ILLUMINA | PAIRED | 204     | SRR1034344 | 18812797 | 3.15G |
| HM152 | SRX376179 | ILLUMINA | PAIRED | 208     | SRR1034339 | 19916982 | 3.34G |
| HM153 | SRX376185 | ILLUMINA | PAIRED | 228     | SRR1034345 | 21067557 | 3.53G |
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
| HM207 | SRX376052 | ILLUMINA | PAIRED | 408     | SRR1034212 | 16226937 | 2.72G |
| HM208 | SRX376054 | ILLUMINA | PAIRED | 429     | SRR1034214 | 17933612 | 3.01G |
| HM211 | SRX376085 | ILLUMINA | PAIRED | 260     | SRR1034245 | 19598873 | 3.23G |
| HM212 | SRX376086 | ILLUMINA | PAIRED | 321     | SRR1034246 | 19255256 | 3.17G |
| HM213 | SRX376087 | ILLUMINA | PAIRED | 278     | SRR1034247 | 18669194 | 3.08G |
| HM214 | SRX376088 | ILLUMINA | PAIRED | 288     | SRR1034248 | 20977571 | 3.52G |
| HM215 | SRX376089 | ILLUMINA | PAIRED | 273     | SRR1034249 | 20823028 | 3.49G |
| HM216 | SRX376090 | ILLUMINA | PAIRED | 328     | SRR1034250 | 19759770 | 3.2G  |
| HM217 | SRX376091 | ILLUMINA | PAIRED | 338     | SRR1034251 | 18976294 | 3.13G |
| HM218 | SRX376092 | ILLUMINA | PAIRED | 248     | SRR1034252 | 19098942 | 3.2G  |
| HM219 | SRX376093 | ILLUMINA | PAIRED | 326     | SRR1034253 | 19638730 | 3.29G |
| HM222 | SRX376098 | ILLUMINA | PAIRED | 303     | SRR1034258 | 19295507 | 3.23G |
| HM223 | SRX376100 | ILLUMINA | PAIRED | 278     | SRR1034260 | 18786715 | 3.15G |
| HM224 | SRX376101 | ILLUMINA | PAIRED | 288     | SRR1034261 | 18273944 | 3.06G |
| HM226 | SRX376102 | ILLUMINA | PAIRED | 373     | SRR1034262 | 13789065 | 2.31G |
| HM227 | SRX376103 | ILLUMINA | PAIRED | 298     | SRR1034263 | 21434460 | 3.59G |
| HM228 | SRX376104 | ILLUMINA | PAIRED | 348     | SRR1034264 | 13074843 | 2.19G |
| HM230 | SRX376106 | ILLUMINA | PAIRED | 268     | SRR1034266 | 19539756 | 3.28G |
| HM232 | SRX376108 | ILLUMINA | PAIRED | 348     | SRR1034268 | 16031153 | 2.69G |
| HM234 | SRX376110 | ILLUMINA | PAIRED | 293     | SRR1034270 | 13200102 | 2.21G |
| HM235 | SRX376099 | ILLUMINA | PAIRED | 268     | SRR1034259 | 20231756 | 3.39G |
| HM237 | SRX376113 | ILLUMINA | PAIRED | 328     | SRR1034273 | 15241718 | 2.56G |
| HM238 | SRX376114 | ILLUMINA | PAIRED | 278     | SRR1034274 | 21531986 | 3.61G |
| HM239 | SRX376115 | ILLUMINA | PAIRED | 328     | SRR1034275 | 20990330 | 3.5G  |
| HM240 | SRX376116 | ILLUMINA | PAIRED | 328     | SRR1034276 | 17143013 | 2.87G |
| HM241 | SRX376117 | ILLUMINA | PAIRED | 328     | SRR1034277 | 18603094 | 3.12G |
| HM242 | SRX376118 | ILLUMINA | PAIRED | 328     | SRR1034278 | 20703978 | 3.45G |
| HM244 | SRX376120 | ILLUMINA | PAIRED | 328     | SRR1034280 | 19585639 | 3.28G |
| HM245 | SRX376121 | ILLUMINA | PAIRED | 328     | SRR1034281 | 20632655 | 3.42G |
| HM247 | SRX376160 | ILLUMINA | PAIRED | 257     | SRR1034320 | 14603330 | 2.45G |

* Failed to assemble
    * HM016 - Bad quality of reads

## Symlink

* 采用的倍数因子值: `2`

* Rsync to hpcc

```shell script
rsync -avP \
    ~/data/plastid/384/ \
    wangq@202.119.37.251:data/plastid/384

# rsync -avP wangq@202.119.37.251:data/plastid/384/ ~/data/plastid/384

```

```shell script
cd ~/data/plastid/384/

export FOLD=2
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/a17/chr.sizes |
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

## Run

```shell script
cd ~/data/plastid/384/

cat opts.tsv | head -n 100 |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi
        
        if [ ! -d {1} ]; then
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
            --quorum \
            --merge \
            --ecphase "1 2 3" \
            \
            --bowtie "Q25L60" \
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
            rm -fr 4_down_sampling
            rm -fr 6_down_sampling
        "
    '

```

## Pack and clean

```shell script
cd ~/data/plastid/384/

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '        
        if [ -f {1}.tar.gz ]; then
            echo "==> {1} .tar.gz"
            exit;
        fi

        if [ ! -f {1}/7_merge_anchors/anchor.merge.fasta ]; then
            echo "==> {1} 7_merge_anchors"
            exit;
        fi

        if [ ! -d {1}/9_quast ]; then
            echo "==> {1} 9_quast"
            exit;
        fi

        echo "==> Clean {1}"
        bash {1}/0_cleanup.sh
        
        echo "==> Create {1}.tar.gz"

        tar -czvf {1}.tar.gz \
            {1}/1_genome/genome.fa \
            {1}/2_illumina/fastqc \
            {1}/2_illumina/insert_size \
            {1}/2_illumina/kat \
            {1}/2_illumina/trim/Q25L60/pe.cor.fa.gz \
            {1}/2_illumina/trim/Q25L60/env.json \
            {1}/3_bowtie \
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
cd ~/data/plastid/384/

cat opts.tsv | head -n 60 |
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
