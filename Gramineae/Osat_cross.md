# 🌾 *Oryza sativa* 9311 X PA64s Genome sequencing

<!-- TOC -->
* [🌾 *Oryza sativa* 9311 X PA64s Genome sequencing](#-oryza-sativa-9311-x-pa64s-genome-sequencing)
  * [Basic info](#basic-info)
  * [Project](#project)
  * [Download](#download)
    * [Reference](#reference)
    * [Illumina](#illumina)
  * [Symlink](#symlink)
  * [Run](#run)
  * [Pack and clean](#pack-and-clean)
  * [VCF](#vcf)
<!-- TOC -->

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
    hnsm filter -s stdin |
    hnsm replace stdin replace.tsv |
    hnsm order stdin <(echo Pt; echo Mt) -o genome.fa

```

### Illumina

* Download `Metadata` from NCBI SRA Run Selector via a web browser
    * <https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA232554>
    * Save it to `SraRunTable.csv`
    * <https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA243018>
    * Save it to `SraRunTable.2.csv`

```shell script
mkdir -p ~/data/plastid/Osat_cross/ena
cd ~/data/plastid/Osat_cross/ena

cat SraRunTable.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f Experiment,"Sample\ Name",Bases \
    > SraRunTable.tsv

cat SraRunTable.2.csv |
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

rgr md ena_info.tsv --fmt

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 2 "{}"

```

| name         | srx       | platform | layout | ilength | srr        |      spots | bases  |
|--------------|-----------|----------|--------|---------|------------|-----------:|--------|
| 9311         | SRX399032 | ILLUMINA | PAIRED | 450     | SRR1060330 | 47,893,596 | 8.92G  |
| 9311-1       | SRX505067 | ILLUMINA | PAIRED | 470     | SRR1210600 | 43,391,116 | 8.08G  |
| 9311-2       | SRX505069 | ILLUMINA | PAIRED | 470     | SRR1210602 | 44,561,271 | 8.3G   |
| 9311-3       | SRX505071 | ILLUMINA | PAIRED | 470     | SRR1210605 | 44,344,722 | 8.26G  |
| DZ-11        | SRX831999 | ILLUMINA | PAIRED |         | SRR1924249 | 85,323,671 | 15.89G |
| DZ-12        | SRX832000 | ILLUMINA | PAIRED |         | SRR1942959 | 85,923,172 | 16G    |
| DZ-19        | SRX832001 | ILLUMINA | PAIRED |         | SRR2002789 | 86,018,191 | 16.02G |
| F2-37        | SRX832002 | ILLUMINA | PAIRED |         | SRR2002790 | 85,313,518 | 15.89G |
| F2-41        | SRX832003 | ILLUMINA | PAIRED |         | SRR2002791 | 85,525,890 | 15.93G |
| LYP9_F1.1    | SRX399035 | ILLUMINA | PAIRED | 450     | SRR1060365 | 46,418,102 | 8.65G  |
| LYP9_F1.2    | SRX399036 | ILLUMINA | PAIRED | 450     | SRR1060366 | 44,878,715 | 8.36G  |
| LYP9_F1.3    | SRX399037 | ILLUMINA | PAIRED | 450     | SRR1060367 | 46,551,221 | 8.67G  |
| LYP9_F2_22.1 | SRX399038 | ILLUMINA | PAIRED | 450     | SRR1060368 | 45,132,915 | 8.41G  |
| LYP9_F2_22.2 | SRX399039 | ILLUMINA | PAIRED | 450     | SRR1060369 | 46,583,397 | 8.68G  |
| LYP9_F2_23.1 | SRX399040 | ILLUMINA | PAIRED | 450     | SRR1060370 | 45,202,202 | 8.42G  |
| LYP9_F2_23.2 | SRX399042 | ILLUMINA | PAIRED | 450     | SRR1060371 | 44,702,616 | 8.33G  |
| LYP9_F2_24.1 | SRX399043 | ILLUMINA | PAIRED | 450     | SRR1060372 | 45,181,206 | 8.42G  |
| LYP9_F2_24.2 | SRX399045 | ILLUMINA | PAIRED | 450     | SRR1060373 | 47,042,532 | 8.76G  |
| LYP9_F2_25.1 | SRX399048 | ILLUMINA | PAIRED | 450     | SRR1060374 | 44,529,836 | 8.29G  |
| LYP9_F2_25.2 | SRX399049 | ILLUMINA | PAIRED | 450     | SRR1060375 | 44,183,328 | 8.23G  |
| LYP9_F2_26.1 | SRX399051 | ILLUMINA | PAIRED | 450     | SRR1060377 | 44,734,088 | 8.33G  |
| LYP9_F2_27.1 | SRX399052 | ILLUMINA | PAIRED | 450     | SRR1060376 | 41,678,053 | 7.76G  |
| LYP9_F2_27.2 | SRX399053 | ILLUMINA | PAIRED | 450     | SRR1060378 | 46,305,884 | 8.63G  |
| LYP9_F2_30.1 | SRX399054 | ILLUMINA | PAIRED | 450     | SRR1060379 | 44,930,781 | 8.37G  |
| LYP9_F2_30.2 | SRX399056 | ILLUMINA | PAIRED | 450     | SRR1060380 | 44,574,102 | 8.3G   |
| LYP9_F2_31.1 | SRX399057 | ILLUMINA | PAIRED | 450     | SRR1060381 | 44,874,811 | 8.36G  |
| LYP9_F2_32.1 | SRX399059 | ILLUMINA | PAIRED | 450     | SRR1060382 | 45,394,532 | 8.46G  |
| LYP9_F2_32.2 | SRX399060 | ILLUMINA | PAIRED | 450     | SRR1060383 | 46,983,743 | 8.75G  |
| LYP9_F2_33.1 | SRX399061 | ILLUMINA | PAIRED | 450     | SRR1060384 | 41,897,760 | 7.8G   |
| LYP9_F2_46.1 | SRX399062 | ILLUMINA | PAIRED | 450     | SRR1060385 | 45,691,393 | 8.51G  |
| LYP9_F2_47.1 | SRX399091 | ILLUMINA | PAIRED | 450     | SRR1060387 | 46,555,725 | 8.67G  |
| LYP9_F2_56.1 | SRX399092 | ILLUMINA | PAIRED | 450     | SRR1060386 | 46,914,710 | 8.74G  |
| LYP9_F2_56.2 | SRX399094 | ILLUMINA | PAIRED | 450     | SRR1060388 | 42,629,618 | 7.94G  |
| LYP9_F2_87.1 | SRX399095 | ILLUMINA | PAIRED | 450     | SRR1060389 | 46,902,296 | 8.74G  |
| LYP9_F2_88.1 | SRX399097 | ILLUMINA | PAIRED | 450     | SRR1060390 | 46,073,642 | 8.58G  |
| LYP9_F2_88.2 | SRX399098 | ILLUMINA | PAIRED | 450     | SRR1060391 | 48,851,306 | 9.1G   |
| LYP9_F2_89.1 | SRX399099 | ILLUMINA | PAIRED | 450     | SRR1060392 | 45,153,709 | 8.41G  |
| LYP9_F2_89.2 | SRX399100 | ILLUMINA | PAIRED | 450     | SRR1060393 | 43,269,777 | 8.06G  |
| LYP9_F2_90.1 | SRX399101 | ILLUMINA | PAIRED | 450     | SRR1060394 | 43,099,482 | 8.03G  |
| LYP9_F2_90.2 | SRX399102 | ILLUMINA | PAIRED | 450     | SRR1060395 | 46,821,210 | 8.72G  |
| LYP9_F2_93.1 | SRX399103 | ILLUMINA | PAIRED | 450     | SRR1060396 | 43,214,069 | 8.05G  |
| LYP9_F2_94.1 | SRX399105 | ILLUMINA | PAIRED | 450     | SRR1060397 | 43,419,871 | 8.09G  |
| PA64         | SRX399034 | ILLUMINA | PAIRED | 450     | SRR1060331 | 48,774,641 | 9.08G  |
| PA64s-1      | SRX505073 | ILLUMINA | PAIRED | 470     | SRR1210618 | 43,688,101 | 8.14G  |
| PA64s-2      | SRX505092 | ILLUMINA | PAIRED | 470     | SRR1210626 | 44,430,332 | 8.28G  |
| PA64s-3      | SRX505097 | ILLUMINA | PAIRED | 470     | SRR1210632 | 39,615,687 | 7.38G  |

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

