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
  - [VCF](#vcf)


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
mkdir -p ~/data/plastid/Mtru_384/genome
cd ~/data/plastid/Mtru_384/genome

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

```

### Illumina

* Download `Metadata` from NCBI SRA Run Selector via a web browser
  * https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA256006
  * Save it to `SraRunTable.txt`

* http://www.medicagohapmap.org/hapmap/germplasm
  * Extract table via `pup`
  * Convert xls to csv via `excel`

```shell script
mkdir -p ~/data/plastid/Mtru_384/ena
cd ~/data/plastid/Mtru_384/ena

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
anchr ena prep | perl - ena_info.yml --ascp

mlr --icsv --omd cat ena_info.csv | head -n 20

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 2 "{}"

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


| name  | srx       | platform | layout | ilength | srr        | spots    | bases |
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

* Failed to assemble
  * HM016 -

## Symlink

* 采用的倍数因子值: `2`

* Rsync to hpcc

```shell script
rsync -avP \
    ~/data/plastid/Mtru_384/ \
    wangq@202.119.37.251:data/plastid/Mtru_384

# rsync -avP wangq@202.119.37.251:data/plastid/Mtru_384/ ~/data/plastid/Mtru_384

```

```shell script
cd ~/data/plastid/Mtru_384/

export FOLD=2
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/a17/chr.sizes |
        tsv-summarize --sum 2
)

cat ena/ena_info.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f name,srr,bases |
    grep -v -w 'HM016' | # Bad quality of reads
    grep -v -w 'HM207' |
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

## Run

```shell script
cd ~/data/plastid/Mtru_384/

cat opts.tsv |
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
cd ~/data/plastid/Mtru_384/

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

* Remove processed files

```shell script
cd ~/data/plastid/Mtru_384/

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ ! -f {1}.tar.gz ]; then
            exit;
        fi

        find ena -type f -name "*{2}*"
    ' |
    xargs rm

```

* Unpack

```shell script
cd ~/data/plastid/Mtru_384/

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
cd ~/data/plastid/Mtru_384/

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
    > Mtru_384.vcf

```

