# Peaches (*Prunus* spp.)

<!-- TOC -->
* [Peaches (*Prunus* spp.)](#peaches-prunus-spp)
  * [基本信息](#基本信息)
  * [项目信息](#项目信息)
  * [其他可能可用的项目](#其他可能可用的项目)
  * [数据下载](#数据下载)
    * [Reference](#reference)
    * [Illumina](#illumina)
  * [采用的倍数因子值](#采用的倍数因子值)
  * [Symlink](#symlink)
  * [Run](#run)
  * [Pack and clean](#pack-and-clean)
<!-- TOC -->

## 基本信息

* Genome: GCF_000346465.2, Prunus_persica_NCBIv2, 227.569 Mb
* Chloroplast: [NC_014697](https://www.ncbi.nlm.nih.gov/nuccore/NC_014697), **Lovell**, 157790 bp

## 项目信息

* [PRJNA34817](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA34817)
    * [SRA Run Selector](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA34817)
    * https://www.nature.com/articles/ng.2586 Supplementary Table 18.

* [PRJNA310042](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA310042)
    * [SRA Run Selector](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA310042)
    * https://www.nature.com/articles/s41467-018-07744-3 Supplementary Data 2

## 其他可能可用的项目

* [PRJNA497989](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA497989)

* [PRJNA312014](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA312014)
    * Peach meiosis mutation and recombination analyisis

* [PRJNA663114](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA663114)
    * Genomic Structure Variation Analyses of Peach

## 数据下载

### Reference

```shell script
mkdir -p ~/data/plastid/peach/genome
cd ~/data/plastid/peach/genome

for ACCESSION in "NC_014697"; do
    URL=$(printf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=%s&id=%s&retmode=text" "fasta" "${ACCESSION}");
    curl $URL -o ${ACCESSION}.fa
done

TAB=$'\t'
cat <<EOF > replace.tsv
NC_014697${TAB}Pt
EOF

cat NC_014697.fa |
    hnsm filter -s stdin |
    hnsm replace stdin replace.tsv -o genome.fa

```

### Illumina

* Download `Metadata` from NCBI SRA Run Selector via a web browser
    * <https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA310042>
    * Save it to `SraRunTable.csv`

* Extract "Supplementary Table 18" from P67 of <https://www.nature.com/articles/ng.2586>
* Extract "Supplementary Data 2" of <https://www.nature.com/articles/s41467-018-07744-3>

```shell script
mkdir -p ~/data/plastid/peach/ena
cd ~/data/plastid/peach/ena

cat SraRunTable.csv |
    mlr --icsv --otsv cat \
    > SraRunTable.tsv

cat SraRunTable.tsv |
    tsv-select -H -f Experiment,"Sample\ Name",Organism,Cultivar |
    sed '1 s/^/#/' |
    sed 's/\b*#$//g' |
    keep-header -- sort -k2,2 |
    mlr --itsv --ocsv cat \
    > source.csv

cat << EOF >> source.csv
SRX150233,PL07,P_ferganensis,Fergana Valley
SRX150247,PL08,Sahua_Hong_Pantao,Southern China
SRX150251,PL09,Shen_Zhou_Mitao,Northern China
SRX150239,PMC05,Oro_A,Brazil
SRX150234,PMC06,GF305,France
SRX150226,PMC07,Bolero,Italy
SRX150230,PMC08,F1_Contender_x_Ambra,Italy
SRX150237,PMC09,IF7310828,Italy
SRX150253,PMC10,Yumyeong,Korea
SRX150243,PMC11,Quetta,Pakistan
SRX150229,PMC12,Earligold,USA
SRX150254,PMC13,PLov2-2N,USA
SRX150255,WK03,P_kansuensis,China
SRX150227,WD05,P_davidiana,China
#P_mira,China
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml --ascp

rgr md ena_info.tsv --fmt

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 2 "{}"

md5sum --check ena_info.md5.txt

```

| name  | srx        | platform | layout | ilength | srr        |       spots | bases  |
|-------|------------|----------|--------|---------|------------|------------:|--------|
| DB13  | SRX1558324 | ILLUMINA | PAIRED |         | SRR3141229 |  48,807,851 | 11.36G |
| DB14  | SRX1558325 | ILLUMINA | PAIRED |         | SRR3141238 |  82,279,917 | 15.33G |
| DB15  | SRX1558326 | ILLUMINA | PAIRED |         | SRR3141248 |  51,281,845 | 11.94G |
| DL01  | SRX1558273 | ILLUMINA | PAIRED |         | SRR3141032 |  57,561,297 | 13.4G  |
| DL02  | SRX1558274 | ILLUMINA | PAIRED |         | SRR3141040 |  51,777,704 | 12.06G |
| DL03  | SRX1558275 | ILLUMINA | PAIRED |         | SRR3141049 |  55,244,740 | 12.86G |
| DL04  | SRX1558276 | ILLUMINA | PAIRED |         | SRR3141057 |  53,428,518 | 12.44G |
| DL06  | SRX1558278 | ILLUMINA | PAIRED |         | SRR3141073 |  58,914,803 | 13.72G |
| DL07  | SRX1558279 | ILLUMINA | PAIRED |         | SRR3141083 |  53,940,267 | 12.56G |
| DL08  | SRX1558280 | ILLUMINA | PAIRED |         | SRR3141098 |  54,502,565 | 12.69G |
| DL09  | SRX1558281 | ILLUMINA | PAIRED |         | SRR3141113 |  56,317,357 | 13.11G |
| DL11  | SRX1558316 | ILLUMINA | PAIRED |         | SRR3141192 |  55,528,290 | 12.93G |
| DL12  | SRX1558317 | ILLUMINA | PAIRED |         | SRR3141204 |  52,518,349 | 12.23G |
| PB07  | SRX1556745 | ILLUMINA | PAIRED |         | SRR3138139 |  65,810,224 | 15.32G |
| PB08  | SRX1556747 | ILLUMINA | PAIRED |         | SRR3138145 |  76,276,852 | 14.21G |
| PB09  | SRX1556748 | ILLUMINA | PAIRED |         | SRR3138146 |  65,786,015 | 15.32G |
| PB10  | SRX1556749 | ILLUMINA | PAIRED |         | SRR3138147 |  57,785,720 | 13.45G |
| PL01  | SRX1556729 | ILLUMINA | PAIRED |         | SRR3138115 |  67,499,522 | 12.57G |
| PL02  | SRX1556731 | ILLUMINA | PAIRED |         | SRR3138117 |  81,150,229 | 15.12G |
| PL03  | SRX1556732 | ILLUMINA | PAIRED |         | SRR3138121 |  60,038,656 | 11.18G |
| PL04  | SRX1556733 | ILLUMINA | PAIRED |         | SRR3138123 |  63,382,341 | 11.81G |
| PL05  | SRX1556734 | ILLUMINA | PAIRED |         | SRR3138129 |  66,399,614 | 12.37G |
| PL06  | SRX1556744 | ILLUMINA | PAIRED |         | SRR3138132 |  73,430,514 | 13.68G |
| PL07  | SRX150233  | ILLUMINA | PAIRED | 225     | SRR502999  |  47,741,287 | 6.67G  |
| PL08  | SRX150247  | ILLUMINA | PAIRED | 400     | SRR502991  |  17,642,272 | 3.48G  |
| PL09  | SRX150251  | ILLUMINA | PAIRED | 400     | SRR502993  |  11,832,973 | 2.34G  |
| PMC05 | SRX150239  | ILLUMINA | PAIRED | 400     | SRR502986  |  37,303,732 | 7.57G  |
| PMC06 | SRX150234  | ILLUMINA | PAIRED | 400     | SRR502983  |  27,123,880 | 5.1G   |
| PMC07 | SRX150226  | ILLUMINA | PAIRED | 400     | SRR501836  |  34,819,509 | 7.07G  |
| PMC08 | SRX150230  | ILLUMINA | PAIRED | 225     | SRR502997  |  43,384,092 | 6.06G  |
| PMC09 | SRX150237  | ILLUMINA | PAIRED | 225     | SRR503001  |  37,195,543 | 5.2G   |
| PMC10 | SRX150253  | ILLUMINA | PAIRED | 400     | SRR502994  |  34,636,238 | 6.45G  |
| PMC11 | SRX150243  | ILLUMINA | PAIRED | 400     | SRR502989  |  17,051,143 | 3.37G  |
| PMC12 | SRX150229  | ILLUMINA | PAIRED | 400     | SRR502996  |  41,263,866 | 7.69G  |
| PMC13 | SRX150254  | ILLUMINA | PAIRED | 400     | SRR502985  | 123,590,441 | 23.25G |
| PW01  | SRX1556767 | ILLUMINA | PAIRED |         | SRR3138168 |  72,248,792 | 13.46G |
| PW02  | SRX1556769 | ILLUMINA | PAIRED |         | SRR3138169 |  67,947,885 | 12.66G |
| PW03  | SRX3181624 | ILLUMINA | PAIRED |         | SRR6031401 |  94,737,277 | 17.65G |
| PW05  | SRX1556771 | ILLUMINA | PAIRED |         | SRR3138171 |  58,754,980 | 10.94G |
| PW06  | SRX1558270 | ILLUMINA | PAIRED |         | SRR3141016 |  61,355,322 | 11.43G |
| PW07  | SRX1558271 | ILLUMINA | PAIRED |         | SRR3141018 |  60,013,031 | 11.18G |
| PW08  | SRX1558272 | ILLUMINA | PAIRED |         | SRR3141019 |  80,461,753 | 14.99G |
| PW09  | SRX1554823 | ILLUMINA | PAIRED |         | SRR3136174 |  68,658,961 | 19.18G |
| PW10  | SRX1554827 | ILLUMINA | PAIRED |         | SRR3136179 |  79,108,844 | 22.1G  |
| PW11  | SRX1554829 | ILLUMINA | PAIRED |         | SRR3136181 |  60,634,898 | 16.94G |
| RW01  | SRX3181586 | ILLUMINA | PAIRED |         | SRR6031356 |  56,988,341 | 15.92G |
| RW02  | SRX3181587 | ILLUMINA | PAIRED |         | SRR6031357 |  54,426,942 | 15.21G |
| RW03  | SRX3181589 | ILLUMINA | PAIRED |         | SRR6031362 |  58,368,636 | 16.31G |
| RW04  | SRX3181622 | ILLUMINA | PAIRED |         | SRR6031396 |  53,274,030 | 14.88G |
| RW05  | SRX3181623 | ILLUMINA | PAIRED |         | SRR6031398 |  72,906,337 | 20.37G |
| RW06  | SRX3181626 | ILLUMINA | PAIRED |         | SRR6031407 |  95,825,292 | 17.85G |
| RW07  | SRX3181627 | ILLUMINA | PAIRED |         | SRR6031411 |  86,782,710 | 16.16G |
| WD05  | SRX150227  | ILLUMINA | PAIRED | 400     | SRR502982  |  31,904,395 | 5.94G  |

## 采用的倍数因子值

* `2`

## Symlink

```shell script
cd ~/data/plastid/peach/

export FOLD=2
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/lovell/chr.sizes |
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
    ~/data/plastid/peach/ \
    wangq@202.119.37.251:data/plastid/peach

# rsync -avP wangq@202.119.37.251:data/plastid/peach/ ~/data/plastid/peach

```

```shell script
cd ~/data/plastid/peach/

cat opts.tsv | # head -n 30 | tail -n 5 |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        if [ ! -d {1} ]; then
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
cd ~/data/plastid/peach/

cat opts.tsv | head -n 30 |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -d {1} ]; then
            exit;
        fi

        cd {1}

        if [ -d 4_down_sampling ]; then
            rm -fr 4_down_sampling
            rm -fr 6_down_sampling
        fi
    '

cat opts.tsv | #head -n 15 | tail -n 5 |
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

cat opts.tsv | head -n 58 |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            if [ -d {1} ]; then
                echo "==> Remove {1}/"
                rm -fr {1}
            fi
        fi
    '

```

