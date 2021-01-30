# 🍅 *Solanum lycopersicum* 360 accessions

[TOC levels=1-3]: # ""

- [🍅 *Solanum lycopersicum* 360 accessions](#-solanum-lycopersicum-360-accessions)
  - [Basic info](#basic-info)
  - [Project](#project)
  - [Other Projects](#other-projects)
  - [Download](#download)
    - [Reference](#reference)
    - [Illumina](#illumina)
  - [Symlink](#symlink)
  - [Run](#run)
  - [Pack and clean](#pack-and-clean)
  - [VCF](#vcf)


## Basic info

* Genome: GCF_000188115.4, SL3.0, 828.349 Mb
* Chloroplast: [NC_007898](https://www.ncbi.nlm.nih.gov/nuccore/NC_007898), 155461 bp
* Mitochondrion: [NC_035963](https://www.ncbi.nlm.nih.gov/nuccore/NC_035963), 446257 bp


## Project

* [PRJNA259308](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA259308)

* <https://www.nature.com/articles/ng.3117>


## Other Projects


https://solgenomics.net/projects/varitome

https://www.ncbi.nlm.nih.gov/bioproject/PRJNA454805


## Download

### Reference

```shell script
mkdir -p ~/data/plastid/Slyc_360/genome
cd ~/data/plastid/Slyc_360/genome

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

```

### Illumina

* Download `Metadata` from NCBI SRA Run Selector via a web browser
  * https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA259308
  * Save it to `SraRunTable.txt`

* "Supplementary Table 1" of <https://www.nature.com/articles/ng.3117> provides detailed info

```shell script
mkdir -p ~/data/plastid/Slyc_360/ena
cd ~/data/plastid/Slyc_360/ena

cat SraRunTable.txt |
    mlr --icsv --otsv cat \
    > SraRunTable.tsv

cat SraRunTable.tsv |
    tsv-select -H -f Experiment,"Sample\ Name",Population,collected_by,Bases,AvgSpotLen |
    tsv-filter -H \
        --not-blank collected_by \
        --istr-ne collected_by:missing \
        --ge Bases:4000000000 \
        --le Bases:10000000000 \
        --ge AvgSpotLen:90 |
    sed '1 s/^/#/' |
    keep-header -- tsv-sort -k2,2 -k5,5nr |
    tsv-uniq -H -f "Sample\ Name" --max 1 |
    cut -f 1-4 |
    mlr --itsv --ocsv cat \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml --ascp

mlr --icsv --omd cat ena_info.csv | head -n 20

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 1 "{}"

# Skips
cat ena_info.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f name,srr,bases |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -f {2}_1.fastq.gz ]; then
            echo 1>&2 {1}
            echo {1}
            exit;
        fi
        if [ ! -f {2}_2.fastq.gz ]; then
            echo 1>&2 {1}
            echo {1}
            exit;
        fi

        LENGTH=$(gzip -dcf {2}_1.fastq.gz |
            head -n 100 |
            faops n50 -H stdin
        )
        [[ $LENGTH -le 90 ]] && ( echo 1>&2 {1}; echo {1}; exit; )

        LENGTH=$(gzip -dcf {2}_2.fastq.gz |
            head -n 100 |
            faops n50 -H stdin
        )
        [[ $LENGTH -le 90 ]] && ( echo 1>&2 {1}; echo {1}; exit;  )
    ' |
    sort -r |
    uniq \
    > skip.lst

```

| name   | srx       | platform | layout | ilength | srr        | spots    | bases |
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


## Symlink

* 采用的倍数因子值: `2`

* Rsync to hpcc

```shell script
rsync -avP \
    ~/data/plastid/Slyc_360/ \
    wangq@202.119.37.251:data/plastid/Slyc_360

# rsync -avP wangq@202.119.37.251:data/plastid/Slyc_360/ ~/data/plastid/Slyc_360

```

```shell script
cd ~/data/plastid/Slyc_360/

export FOLD=2
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/h1706/chr.sizes |
        tsv-summarize --sum 2
)

cat ena/ena_info.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f name,srr,bases |
    tsv-join -H --filter-file ena/skip.lst --key-fields name --exclude |
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

wc -l opts.tsv
# 250

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
        popd > /dev/null

        mkdir -p {1}/2_illumina
        pushd {1}/2_illumina

        ln -fs ../../ena/{2}_1.fastq.gz R1.fq.gz
        ln -fs ../../ena/{2}_2.fastq.gz R2.fq.gz
        popd > /dev/null
    '

```

* Rsync non-processed files to hpcc

```shell script
cd ~/data/plastid/Slyc_360/

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        find ena -type f -name "*{2}_*"
    ' \
    > rsync.lst

rsync -avP \
    --files-from=rsync.lst \
    ~/data/plastid/Slyc_360/ \
    wangq@202.119.37.251:data/plastid/Slyc_360

```

## Run

```shell script
cd ~/data/plastid/Slyc_360/

cat opts.tsv | #head -n 100 | #tail -n 10 |
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
cd ~/data/plastid/Slyc_360/

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
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
cd ~/data/plastid/Slyc_360/

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
cd ~/data/plastid/Slyc_360/

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
cd ~/data/plastid/Slyc_360/

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
    > Slyc_360.vcf

rm -fr vcf

```

