# 🌾 *Oryza sativa* 9311 X PA64s Genome sequencing

[TOC levels=1-3]: # ""

- [🌾 *Oryza sativa* 9311 X PA64s Genome sequencing](#-oryza-sativa-9311-x-pa64s-genome-sequencing)
  - [Basic info](#basic-info)
  - [Project](#project)
  - [Download](#download)
    - [Reference](#reference)
    - [Illumina](#illumina)
  - [Symlink](#symlink)
  - [Run](#run)
  - [Pack and clean](#pack-and-clean)
  - [VCF](#vcf)


## Basic info

* Genome: GCA_000004655.2, Cultivar: 93-11, 426.337 Mb
* Chloroplast: [NC_008155](https://www.ncbi.nlm.nih.gov/nuccore/NC_008155), **Indica**, 134496 bp
* Mitochondrion: [NC_007886](https://www.ncbi.nlm.nih.gov/nuccore/NC_007886), **Indica**, 491515 bp


## Project

* <https://www.nature.com/articles/nature14649>
* PRJNA232554 - rice LYP9
* PRJNA243018 - 9311, PA64s

## Download

### Reference

```shell script
mkdir -p ~/data/plastid/Osat_cross/genome
cd ~/data/plastid/Osat_cross/genome

for ACCESSION in "NC_008155" "NC_007886"; do
    URL=$(printf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=%s&id=%s&retmode=text" "fasta" "${ACCESSION}");
    curl $URL -o ${ACCESSION}.fa
done

TAB=$'\t'
cat <<EOF > replace.tsv
NC_008155${TAB}Pt
NC_007886${TAB}Mt
EOF

cat NC_008155.fa NC_007886.fa |
    faops filter -s stdin stdout |
    faops replace stdin replace.tsv stdout |
    faops order stdin <(echo Pt; echo Mt) genome.fa

```

### Illumina

* Download `Metadata` from NCBI SRA Run Selector via a web browser
  * <https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA232554>
  * Save it to `SraRunTable.txt`
  * <https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA243018>
  * Save it to `SraRunTable.2.txt`


```shell script
mkdir -p ~/data/plastid/Osat_cross/ena
cd ~/data/plastid/Osat_cross/ena

cat SraRunTable.txt |
    mlr --icsv --otsv cat |
    tsv-select -H -f Experiment,"Sample\ Name",Bases \
    > SraRunTable.tsv

cat SraRunTable.2.txt |
    mlr --icsv --otsv cat |
    tsv-filter -H --str-eq Organism:"Oryza sativa" |
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
    parallel --no-run-if-empty -j 2 "{}"

```

| name              | srx        | platform | layout | ilength | srr         | spots    | bases |
|:------------------|:-----------|:---------|:-------|:--------|:------------|:---------|:------|
| Aihuangmi         | SRX9014014 | ILLUMINA | PAIRED |         | SRR12524330 | 26502344 | 7.4G  |
| Aizizhan          | SRX9014081 | ILLUMINA | PAIRED |         | SRR12524263 | 24847302 | 6.94G |
| Ce64              | SRX9013320 | ILLUMINA | PAIRED |         | SRR12523666 | 29513684 | 8.25G |
| Chenghui448       | SRX9013350 | ILLUMINA | PAIRED |         | SRR12523636 | 26196738 | 7.32G |
| Gang46A           | SRX9013449 | ILLUMINA | PAIRED |         | SRR12523537 | 24368210 | 6.81G |
| Guangchang13      | SRX9013477 | ILLUMINA | PAIRED |         | SRR12523509 | 23726716 | 6.63G |
| Guangchangai      | SRX9013478 | ILLUMINA | PAIRED |         | SRR12523508 | 21711376 | 6.07G |
| Guanger104        | SRX9013480 | ILLUMINA | PAIRED |         | SRR12523506 | 23926052 | 6.68G |
| Guanghui3550      | SRX9013489 | ILLUMINA | PAIRED |         | SRR12523497 | 24548440 | 6.86G |
| Guangqiuai        | SRX9013497 | ILLUMINA | PAIRED |         | SRR12523489 | 20491477 | 5.73G |
| Gui630            | SRX9013509 | ILLUMINA | PAIRED |         | SRR12523477 | 24729656 | 6.91G |
| Guichao2hao       | SRX9013512 | ILLUMINA | PAIRED |         | SRR12523474 | 25728794 | 7.19G |
| IR127-80-1-10     | SRX9013598 | ILLUMINA | PAIRED |         | SRR12523388 | 20664447 | 5.77G |
| IR24-1            | SRX9013585 | ILLUMINA | PAIRED |         | SRR12523401 | 31908719 | 8.92G |
| IR24-2            | SRX9013586 | ILLUMINA | PAIRED |         | SRR12523400 | 21187347 | 5.92G |
| IR24609-134-3-3   | SRX9013652 | ILLUMINA | PAIRED |         | SRR12523334 | 27609578 | 7.71G |
| IR26-1            | SRX9013587 | ILLUMINA | PAIRED |         | SRR12523399 | 27253723 | 7.61G |
| IR26-2            | SRX9013588 | ILLUMINA | PAIRED |         | SRR12523398 | 25374053 | 7.09G |
| IR28-1            | SRX9013589 | ILLUMINA | PAIRED |         | SRR12523397 | 23789070 | 6.65G |
| IR28-2            | SRX9013590 | ILLUMINA | PAIRED |         | SRR12523396 | 22598944 | 6.31G |
| IR28212-71-4-2-3  | SRX9013663 | ILLUMINA | PAIRED |         | SRR12524032 | 22454542 | 6.27G |
| IR2863-6-3        | SRX9013615 | ILLUMINA | PAIRED |         | SRR12523371 | 25135346 | 7.02G |
| IR30              | SRX9013591 | ILLUMINA | PAIRED |         | SRR12523395 | 24713580 | 6.9G  |
| IR36              | SRX9013593 | ILLUMINA | PAIRED |         | SRR12523393 | 23394866 | 6.54G |
| IR56              | SRX9013595 | ILLUMINA | PAIRED |         | SRR12523391 | 23020110 | 6.43G |
| IR5657-33-2       | SRX9013618 | ILLUMINA | PAIRED |         | SRR12523368 | 22033730 | 6.16G |
| IR58              | SRX9013596 | ILLUMINA | PAIRED |         | SRR12523390 | 26952170 | 7.53G |
| IR64              | SRX9013597 | ILLUMINA | PAIRED |         | SRR12523389 | 25322564 | 7.08G |
| IR661             | SRX9013599 | ILLUMINA | PAIRED |         | SRR12523387 | 24898618 | 6.96G |
| IR8               | SRX9013584 | ILLUMINA | PAIRED |         | SRR12523402 | 23525870 | 6.57G |
| IR841             | SRX9013602 | ILLUMINA | PAIRED |         | SRR12523384 | 29063595 | 8.12G |
| Mianhui501        | SRX9013848 | ILLUMINA | PAIRED |         | SRR12523847 | 27237888 | 7.61G |
| Mianhui725        | SRX9013850 | ILLUMINA | PAIRED |         | SRR12523845 | 21948987 | 6.13G |
| Minghui86         | SRX9013867 | ILLUMINA | PAIRED |         | SRR12523828 | 26851702 | 7.5G  |
| Muquanzhong       | SRX9013884 | ILLUMINA | PAIRED |         | SRR12523811 | 30996097 | 8.66G |
| Peiai64           | SRX9013930 | ILLUMINA | PAIRED |         | SRR12523765 | 27027905 | 7.55G |
| Peiai64S          | SRX9013931 | ILLUMINA | PAIRED |         | SRR12523764 | 26948568 | 7.53G |
| Peidi             | SRX9013932 | ILLUMINA | PAIRED |         | SRR12523763 | 27727480 | 7.75G |
| Qiguizao25        | SRX9013953 | ILLUMINA | PAIRED |         | SRR12523742 | 28807308 | 8.05G |
| Qingsiai16A       | SRX9013961 | ILLUMINA | PAIRED |         | SRR12523734 | 21280102 | 5.95G |
| Taiyin1hao        | SRX9014167 | ILLUMINA | PAIRED |         | SRR12524177 | 31109455 | 8.69G |
| Yuejingsimiao2hao | SRX8946546 | ILLUMINA | PAIRED |         | SRR12452051 | 22244199 | 6.21G |
| Yuexiangzhan      | SRX8946551 | ILLUMINA | PAIRED |         | SRR12452046 | 20644825 | 5.77G |
| Zhaiyeqing8hao    | SRX8946574 | ILLUMINA | PAIRED |         | SRR12452023 | 23501708 | 6.57G |
| Zhenguiai1hao     | SRX8946581 | ILLUMINA | PAIRED |         | SRR12452016 | 24285367 | 6.79G |


## Symlink

* 采用的倍数因子值: `2`

* Rsync to hpcc

```shell script
rsync -avP \
    ~/data/plastid/Osat_cross/ \
    wangq@202.119.37.251:data/plastid/Osat_cross

# rsync -avP wangq@202.119.37.251:data/plastid/Osat_cross/ ~/data/plastid/Osat_cross

```

```shell script
cd ~/data/plastid/Osat_cross/

export FOLD=2
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/nip/chr.sizes |
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
        popd > /dev/null

        mkdir -p {1}/2_illumina
        pushd {1}/2_illumina

        ln -fs ../../ena/{2}_1.fastq.gz R1.fq.gz
        ln -fs ../../ena/{2}_2.fastq.gz R2.fq.gz
        popd > /dev/null
    '

```

## Run

```shell script
cd ~/data/plastid/Osat_cross/

cat opts.tsv | #head -n 20 | #tail -n 10 |
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
cd ~/data/plastid/Osat_cross/

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
cd ~/data/plastid/Osat_cross/

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
cd ~/data/plastid/Osat_cross/

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
cd ~/data/plastid/Osat_cross/

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
    > Osat_cross.vcf

rm -fr vcf

```

