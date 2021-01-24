# *Arabidopsis thaliana* Mutation Accumulation Lines

[TOC levels=1-3]: # ""

- [*Arabidopsis thaliana* Mutation Accumulation Lines](#arabidopsis-thaliana-mutation-accumulation-lines)
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

* Genome: GCF_000001735.3, TAIR10, 119.668 Mb
* Chloroplast: [NC_000932](https://www.ncbi.nlm.nih.gov/nuccore/NC_000932), **Columbia**, 154478 bp
* Mitochondrion: [Y08501](https://www.ncbi.nlm.nih.gov/nuccore/Y08501), 366924 bp


## Project

* <https://www.genetics.org/content/211/2/703>
* PRJNA434660


## Download

### Reference

```shell script
mkdir -p ~/data/plastid/Atha_ma/genome
cd ~/data/plastid/Atha_ma/genome

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
  * <https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA434660>
  * Save it to `SraRunTable.txt`

```shell script
mkdir -p ~/data/plastid/Atha_ma/ena
cd ~/data/plastid/Atha_ma/ena

cat SraRunTable.txt |
    mlr --icsv --otsv cat |
    tsv-select -H -f Experiment,"Library\ Name",Bases,AssemblyName \
    > SraRunTable.tsv

cat SraRunTable.tsv |
    sed '1 s/^/#/' |
    tsv-filter -H --empty "AssemblyName" | # filter out aligned reads
    keep-header -- tsv-sort -k2,2n -k3,3nr |
    tsv-uniq -H -f "Library\ Name" --max 1 |
    mlr --itsv --ocsv cat \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml --ascp

mlr --icsv --omd cat ena_info.csv

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 2 "{}"

```

| name | srx        | platform | layout | ilength | srr        | spots    | bases |
|:-----|:-----------|:---------|:-------|:--------|:-----------|:---------|:------|
| 1    | SRX3722757 | ILLUMINA | PAIRED |         | SRR6750188 | 19373010 | 5.45G |
| 100  | SRX3722807 | ILLUMINA | PAIRED |         | SRR6750138 | 5117910  | 1.43G |
| 101  | SRX3722769 | ILLUMINA | PAIRED |         | SRR6750176 | 6270781  | 1.76G |
| 102  | SRX3722768 | ILLUMINA | PAIRED |         | SRR6750177 | 7347849  | 2.06G |
| 103  | SRX3722767 | ILLUMINA | PAIRED |         | SRR6750178 | 7756818  | 2.17G |
| 105  | SRX3722766 | ILLUMINA | PAIRED |         | SRR6750179 | 25373012 | 7.1G  |
| 106  | SRX3722765 | ILLUMINA | PAIRED |         | SRR6750180 | 15679744 | 4.39G |
| 108  | SRX3722764 | ILLUMINA | PAIRED |         | SRR6750181 | 24420035 | 6.83G |
| 109  | SRX3722763 | ILLUMINA | PAIRED |         | SRR6750182 | 29576243 | 8.28G |
| 11   | SRX3722752 | ILLUMINA | PAIRED |         | SRR6750193 | 9678123  | 2.72G |
| 110  | SRX3722762 | ILLUMINA | PAIRED |         | SRR6750183 | 18479365 | 5.17G |
| 111  | SRX3722777 | ILLUMINA | PAIRED |         | SRR6750168 | 21107376 | 5.91G |
| 112  | SRX3722751 | ILLUMINA | PAIRED |         | SRR6750194 | 22462991 | 6.29G |
| 113  | SRX3722770 | ILLUMINA | PAIRED |         | SRR6750175 | 16933984 | 4.74G |
| 114  | SRX3722771 | ILLUMINA | PAIRED |         | SRR6750174 | 14382676 | 4.03G |
| 115  | SRX3722772 | ILLUMINA | PAIRED |         | SRR6750173 | 14386338 | 4.03G |
| 116  | SRX3722773 | ILLUMINA | PAIRED |         | SRR6750172 | 13879317 | 3.88G |
| 117  | SRX3722774 | ILLUMINA | PAIRED |         | SRR6750171 | 15671158 | 4.38G |
| 118  | SRX3722775 | ILLUMINA | PAIRED |         | SRR6750170 | 10573153 | 2.96G |
| 119  | SRX3722776 | ILLUMINA | PAIRED |         | SRR6750169 | 12457109 | 3.49G |
| 13   | SRX3722787 | ILLUMINA | PAIRED |         | SRR6750158 | 9761456  | 2.75G |
| 14   | SRX3722786 | ILLUMINA | PAIRED |         | SRR6750159 | 20246568 | 5.69G |
| 15   | SRX3722785 | ILLUMINA | PAIRED |         | SRR6750160 | 14300936 | 4.02G |
| 16   | SRX3722784 | ILLUMINA | PAIRED |         | SRR6750161 | 7079241  | 1.99G |
| 17   | SRX3722783 | ILLUMINA | PAIRED |         | SRR6750162 | 12964429 | 3.65G |
| 18   | SRX3722782 | ILLUMINA | PAIRED |         | SRR6750163 | 11057558 | 3.11G |
| 19   | SRX3722781 | ILLUMINA | PAIRED |         | SRR6750164 | 13794328 | 3.88G |
| 2    | SRX3722756 | ILLUMINA | PAIRED |         | SRR6750189 | 22898416 | 6.44G |
| 20   | SRX3722780 | ILLUMINA | PAIRED |         | SRR6750165 | 17273335 | 4.86G |
| 22   | SRX3722779 | ILLUMINA | PAIRED |         | SRR6750166 | 11222359 | 3.16G |
| 23   | SRX3722778 | ILLUMINA | PAIRED |         | SRR6750167 | 12132917 | 3.41G |
| 24   | SRX3722743 | ILLUMINA | PAIRED |         | SRR6750202 | 13142115 | 3.7G  |
| 25   | SRX3722744 | ILLUMINA | PAIRED |         | SRR6750201 | 16667704 | 4.69G |
| 26   | SRX3722741 | ILLUMINA | PAIRED |         | SRR6750204 | 12391990 | 3.49G |
| 27   | SRX3722742 | ILLUMINA | PAIRED |         | SRR6750203 | 11980358 | 3.37G |
| 28   | SRX3722747 | ILLUMINA | PAIRED |         | SRR6750198 | 7271193  | 2.05G |
| 29   | SRX3722748 | ILLUMINA | PAIRED |         | SRR6750197 | 7566356  | 2.13G |
| 3    | SRX3722755 | ILLUMINA | PAIRED |         | SRR6750190 | 19073071 | 5.36G |
| 31   | SRX3722745 | ILLUMINA | PAIRED |         | SRR6750200 | 6307628  | 1.77G |
| 32   | SRX3722746 | ILLUMINA | PAIRED |         | SRR6750199 | 8743376  | 2.46G |
| 33   | SRX3722749 | ILLUMINA | PAIRED |         | SRR6750196 | 8438537  | 2.37G |
| 34   | SRX3722750 | ILLUMINA | PAIRED |         | SRR6750195 | 15970940 | 4.49G |
| 35   | SRX3722839 | ILLUMINA | PAIRED |         | SRR6750106 | 16842383 | 4.74G |
| 36   | SRX3722838 | ILLUMINA | PAIRED |         | SRR6750107 | 4730674  | 1.33G |
| 37   | SRX3722841 | ILLUMINA | PAIRED |         | SRR6750104 | 6673071  | 1.88G |
| 38   | SRX3722840 | ILLUMINA | PAIRED |         | SRR6750105 | 6874418  | 1.93G |
| 39   | SRX3722843 | ILLUMINA | PAIRED |         | SRR6750102 | 6178608  | 1.73G |
| 4    | SRX3722754 | ILLUMINA | PAIRED |         | SRR6750191 | 19782170 | 5.56G |
| 40   | SRX3722842 | ILLUMINA | PAIRED |         | SRR6750103 | 30213719 | 8.46G |
| 41   | SRX3722845 | ILLUMINA | PAIRED |         | SRR6750100 | 34624502 | 9.74G |
| 42   | SRX3722844 | ILLUMINA | PAIRED |         | SRR6750101 | 32655703 | 9.18G |
| 43   | SRX3722847 | ILLUMINA | PAIRED |         | SRR6750098 | 21296310 | 5.99G |
| 44   | SRX3722846 | ILLUMINA | PAIRED |         | SRR6750099 | 20284562 | 5.71G |
| 45   | SRX3722832 | ILLUMINA | PAIRED |         | SRR6750113 | 19715321 | 5.55G |
| 46   | SRX3722833 | ILLUMINA | PAIRED |         | SRR6750112 | 17587441 | 4.95G |
| 47   | SRX3722834 | ILLUMINA | PAIRED |         | SRR6750111 | 16884938 | 4.75G |
| 48   | SRX3722835 | ILLUMINA | PAIRED |         | SRR6750110 | 31038257 | 8.73G |
| 49   | SRX3722828 | ILLUMINA | PAIRED |         | SRR6750117 | 27632916 | 7.77G |
| 5    | SRX3722761 | ILLUMINA | PAIRED |         | SRR6750184 | 24584146 | 6.91G |
| 50   | SRX3722829 | ILLUMINA | PAIRED |         | SRR6750116 | 32484942 | 9.14G |
| 51   | SRX3722830 | ILLUMINA | PAIRED |         | SRR6750115 | 16144296 | 4.54G |
| 52   | SRX3722831 | ILLUMINA | PAIRED |         | SRR6750114 | 10249750 | 2.88G |
| 53   | SRX3722836 | ILLUMINA | PAIRED |         | SRR6750109 | 17853999 | 5.02G |
| 54   | SRX3722837 | ILLUMINA | PAIRED |         | SRR6750108 | 14895309 | 4.19G |
| 55   | SRX3722823 | ILLUMINA | PAIRED |         | SRR6750122 | 16336521 | 4.59G |
| 56   | SRX3722822 | ILLUMINA | PAIRED |         | SRR6750123 | 12467043 | 3.51G |
| 57   | SRX3722821 | ILLUMINA | PAIRED |         | SRR6750124 | 11727847 | 3.3G  |
| 58   | SRX3722820 | ILLUMINA | PAIRED |         | SRR6750125 | 15062121 | 4.24G |
| 59   | SRX3722827 | ILLUMINA | PAIRED |         | SRR6750118 | 16604189 | 4.67G |
| 6    | SRX3722760 | ILLUMINA | PAIRED |         | SRR6750185 | 17843413 | 5.02G |
| 60   | SRX3722826 | ILLUMINA | PAIRED |         | SRR6750119 | 16073531 | 4.52G |
| 61   | SRX3722825 | ILLUMINA | PAIRED |         | SRR6750120 | 7287471  | 2.05G |
| 62   | SRX3722824 | ILLUMINA | PAIRED |         | SRR6750121 | 8883773  | 2.5G  |
| 63   | SRX3722819 | ILLUMINA | PAIRED |         | SRR6750126 | 7464946  | 2.1G  |
| 64   | SRX3722818 | ILLUMINA | PAIRED |         | SRR6750127 | 13087075 | 3.68G |
| 65   | SRX3722816 | ILLUMINA | PAIRED |         | SRR6750129 | 10847788 | 3.05G |
| 66   | SRX3722817 | ILLUMINA | PAIRED |         | SRR6750128 | 23217469 | 6.53G |
| 67   | SRX3722814 | ILLUMINA | PAIRED |         | SRR6750131 | 8709682  | 2.45G |
| 68   | SRX3722815 | ILLUMINA | PAIRED |         | SRR6750130 | 25060688 | 7.05G |
| 69   | SRX3722812 | ILLUMINA | PAIRED |         | SRR6750133 | 24155909 | 6.79G |
| 7    | SRX3722759 | ILLUMINA | PAIRED |         | SRR6750186 | 22660053 | 6.37G |
| 71   | SRX3722813 | ILLUMINA | PAIRED |         | SRR6750132 | 11256329 | 3.17G |
| 72   | SRX3722810 | ILLUMINA | PAIRED |         | SRR6750135 | 18843823 | 5.3G  |
| 73   | SRX3722811 | ILLUMINA | PAIRED |         | SRR6750134 | 20992892 | 5.9G  |
| 74   | SRX3722808 | ILLUMINA | PAIRED |         | SRR6750137 | 20394780 | 5.74G |
| 75   | SRX3722809 | ILLUMINA | PAIRED |         | SRR6750136 | 17185651 | 4.83G |
| 76   | SRX3722795 | ILLUMINA | PAIRED |         | SRR6750150 | 31594750 | 8.89G |
| 77   | SRX3722794 | ILLUMINA | PAIRED |         | SRR6750151 | 21260653 | 5.98G |
| 78   | SRX3722797 | ILLUMINA | PAIRED |         | SRR6750148 | 20129173 | 5.66G |
| 79   | SRX3722796 | ILLUMINA | PAIRED |         | SRR6750149 | 19283317 | 5.42G |
| 8    | SRX3722758 | ILLUMINA | PAIRED |         | SRR6750187 | 21114102 | 5.94G |
| 80   | SRX3722791 | ILLUMINA | PAIRED |         | SRR6750154 | 25663976 | 7.22G |
| 81   | SRX3722790 | ILLUMINA | PAIRED |         | SRR6750155 | 29672464 | 8.35G |
| 82   | SRX3722793 | ILLUMINA | PAIRED |         | SRR6750152 | 22156150 | 6.23G |
| 83   | SRX3722792 | ILLUMINA | PAIRED |         | SRR6750153 | 31121423 | 8.75G |
| 84   | SRX3722789 | ILLUMINA | PAIRED |         | SRR6750156 | 22660815 | 6.37G |
| 85   | SRX3722788 | ILLUMINA | PAIRED |         | SRR6750157 | 17865412 | 5.02G |
| 86   | SRX3722798 | ILLUMINA | PAIRED |         | SRR6750147 | 12898203 | 3.63G |
| 88   | SRX3722799 | ILLUMINA | PAIRED |         | SRR6750146 | 11027799 | 3.1G  |
| 89   | SRX3722800 | ILLUMINA | PAIRED |         | SRR6750145 | 9923887  | 2.79G |
| 9    | SRX3722753 | ILLUMINA | PAIRED |         | SRR6750192 | 24966582 | 7.02G |
| 91   | SRX3722801 | ILLUMINA | PAIRED |         | SRR6750144 | 16907458 | 4.76G |
| 92   | SRX3722802 | ILLUMINA | PAIRED |         | SRR6750143 | 18640836 | 5.24G |
| 94   | SRX3722803 | ILLUMINA | PAIRED |         | SRR6750142 | 10663211 | 3G    |
| 96   | SRX3722804 | ILLUMINA | PAIRED |         | SRR6750141 | 15697461 | 4.42G |
| 98   | SRX3722805 | ILLUMINA | PAIRED |         | SRR6750140 | 11142151 | 3.13G |
| 99   | SRX3722806 | ILLUMINA | PAIRED |         | SRR6750139 | 14290149 | 4.02G |


## Symlink

* 采用的倍数因子值: `2`

* Rsync to hpcc

```shell script
rsync -avP \
    ~/data/plastid/Atha_ma/ \
    wangq@202.119.37.251:data/plastid/Atha_ma

# rsync -avP wangq@202.119.37.251:data/plastid/Atha_ma/ ~/data/plastid/Atha_ma

```

```shell script
cd ~/data/plastid/Atha_ma/

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
cd ~/data/plastid/Atha_ma/

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
    ~/data/plastid/Atha_ma/ \
    wangq@202.119.37.251:data/plastid/Atha_ma

```

## Run

```shell script
cd ~/data/plastid/Atha_ma/

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
cd ~/data/plastid/Atha_ma/

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
cd ~/data/plastid/Atha_ma/

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
cd ~/data/plastid/Atha_ma/

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
cd ~/data/plastid/Atha_ma/

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
    > Atha_ma.vcf

```

