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


## 基本信息

* Genome: GCF_000001735.3, TAIR10, 119.668 Mb
* Chloroplast: [NC_000932](https://www.ncbi.nlm.nih.gov/nuccore/NC_000932), **Columbia**, 154478 bp
* Chloroplast: [KX551970](https://www.ncbi.nlm.nih.gov/nuccore/KX551970), **Landsberg erecta**,
  154515 bp
* Mitochondrion: [Y08501](https://www.ncbi.nlm.nih.gov/nuccore/Y08501), 366924 bp


## 项目信息

* PRJNA178613

* <https://www.pnas.org/content/109/51/20992.long>

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

```shell script
mkdir -p ~/data/plastid/Atha_cross/ena
cd ~/data/plastid/Atha_cross/ena

cat SraRunTable.txt |
    mlr --icsv --otsv cat \
    > SraRunTable.tsv

cat SraRunTable.tsv |
    tsv-select -H -f Experiment,"Sample\ Name",Bases |
    tsv-filter -H \
        --ge Bases:1000000000 |
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

| name            | srx       | platform | layout | ilength | srr       | spots    | bases  |
|:----------------|:----------|:---------|:-------|:--------|:----------|:---------|:-------|
| Sample_14       | SRX202214 | ILLUMINA | PAIRED | 450     | SRR611088 | 25212780 | 4.7G   |
| Sample_18       | SRX202197 | ILLUMINA | PAIRED | 450     | SRR611072 | 12000000 | 2.23G  |
| Sample_19       | SRX202198 | ILLUMINA | PAIRED | 450     | SRR611073 | 12000000 | 2.23G  |
| Sample_20       | SRX202201 | ILLUMINA | PAIRED | 450     | SRR611074 | 12329694 | 2.3G   |
| Sample_21       | SRX202202 | ILLUMINA | PAIRED | 450     | SRR611075 | 12000000 | 2.23G  |
| Sample_4        | SRX202195 | ILLUMINA | PAIRED | 450     | SRR611076 | 12000000 | 2.23G  |
| Sample_5        | SRX202211 | ILLUMINA | PAIRED | 450     | SRR611089 | 12000000 | 2.23G  |
| Sample_5        | SRX202211 | ILLUMINA | PAIRED | 450     | SRR616963 | 13935870 | 2.6G   |
| Sample_6        | SRX202212 | ILLUMINA | PAIRED | 450     | SRR611090 | 23590718 | 4.39G  |
| Sample_7        | SRX202213 | ILLUMINA | PAIRED | 450     | SRR611091 | 25939674 | 4.83G  |
| Sample_8        | SRX202196 | ILLUMINA | PAIRED | 450     | SRR611077 | 12000000 | 2.23G  |
| Sample_Col_G    | SRX202246 | ILLUMINA | PAIRED | 450     | SRR611086 | 49891349 | 9.29G  |
| Sample_Col_G    | SRX202246 | ILLUMINA | PAIRED | 450     | SRR616966 | 24851796 | 4.63G  |
| Sample_Ler_XL_4 | SRX202247 | ILLUMINA | PAIRED | 450     | SRR611087 | 50791450 | 9.46G  |
| Sample_Ler_XL_4 | SRX202247 | ILLUMINA | PAIRED | 450     | SRR616965 | 25436255 | 4.74G  |
| Sample_c1c2     | SRX202204 | ILLUMINA | PAIRED | 450     | SRR611078 | 13983882 | 2.6G   |
| Sample_c41      | SRX202215 | ILLUMINA | PAIRED | 450     | SRR611092 | 25604873 | 4.77G  |
| Sample_c42      | SRX202216 | ILLUMINA | PAIRED | 450     | SRR611093 | 25505766 | 4.75G  |
| Sample_c45      | SRX202217 | ILLUMINA | PAIRED | 450     | SRR611094 | 25572903 | 4.76G  |
| Sample_c47      | SRX202203 | ILLUMINA | PAIRED | 450     | SRR611079 | 11642813 | 2.17G  |
| Sample_c48      | SRX202218 | ILLUMINA | PAIRED | 450     | SRR611095 | 25722535 | 4.79G  |
| Sample_c51      | SRX202219 | ILLUMINA | PAIRED | 450     | SRR611096 | 25759432 | 4.8G   |
| Sample_c52      | SRX202220 | ILLUMINA | PAIRED | 450     | SRR611097 | 25769891 | 4.8G   |
| Sample_c54      | SRX202221 | ILLUMINA | PAIRED | 450     | SRR611098 | 25762016 | 4.8G   |
| Sample_c57      | SRX202222 | ILLUMINA | PAIRED | 450     | SRR611099 | 25534424 | 4.76G  |
| Sample_c61      | SRX202223 | ILLUMINA | PAIRED | 450     | SRR611100 | 25601360 | 4.77G  |
| Sample_c62      | SRX202224 | ILLUMINA | PAIRED | 450     | SRR611101 | 25682651 | 4.78G  |
| Sample_c63      | SRX202225 | ILLUMINA | PAIRED | 450     | SRR611102 | 25535887 | 4.76G  |
| Sample_c64      | SRX202226 | ILLUMINA | PAIRED | 450     | SRR611103 | 25536189 | 4.76G  |
| Sample_c65      | SRX202227 | ILLUMINA | PAIRED | 450     | SRR611104 | 25493386 | 4.75G  |
| Sample_c66      | SRX202228 | ILLUMINA | PAIRED | 450     | SRR611105 | 25586286 | 4.77G  |
| Sample_c73      | SRX202229 | ILLUMINA | PAIRED | 450     | SRR611106 | 23348564 | 4.35G  |
| Sample_c81      | SRX202230 | ILLUMINA | PAIRED | 450     | SRR611107 | 25760277 | 4.8G   |
| Sample_c82      | SRX202231 | ILLUMINA | PAIRED | 450     | SRR611108 | 25625314 | 4.77G  |
| Sample_c83      | SRX202232 | ILLUMINA | PAIRED | 450     | SRR611109 | 25800250 | 4.81G  |
| Sample_c84      | SRX202233 | ILLUMINA | PAIRED | 450     | SRR611110 | 25880633 | 4.82G  |
| Sample_c85      | SRX202234 | ILLUMINA | PAIRED | 450     | SRR611111 | 25794536 | 4.8G   |
| Sample_c87      | SRX202235 | ILLUMINA | PAIRED | 450     | SRR611112 | 24917455 | 4.64G  |
| Sample_c88      | SRX202236 | ILLUMINA | PAIRED | 450     | SRR611113 | 24818024 | 4.62G  |
| Sample_c89      | SRX202237 | ILLUMINA | PAIRED | 450     | SRR611114 | 25551606 | 4.76G  |
| Sample_c90      | SRX202238 | ILLUMINA | PAIRED | 450     | SRR611115 | 25161363 | 4.69G  |
| Sample_c91      | SRX202240 | ILLUMINA | PAIRED | 450     | SRR611116 | 25624037 | 4.77G  |
| Sample_c92      | SRX202241 | ILLUMINA | PAIRED | 450     | SRR611117 | 25688992 | 4.78G  |
| Sample_c93      | SRX202242 | ILLUMINA | PAIRED | 450     | SRR611118 | 13983461 | 2.6G   |
| Sample_c93      | SRX202242 | ILLUMINA | PAIRED | 450     | SRR616982 | 11739658 | 2.19G  |
| Sample_c94      | SRX202243 | ILLUMINA | PAIRED | 450     | SRR611084 | 75625099 | 14.09G |
| Sample_c95      | SRX202245 | ILLUMINA | PAIRED | 450     | SRR611085 | 74530836 | 13.88G |
| Sample_l2c2     | SRX202205 | ILLUMINA | PAIRED | 450     | SRR611080 | 13984182 | 2.6G   |
| Sample_l2l3     | SRX202206 | ILLUMINA | PAIRED | 450     | SRR611081 | 13980212 | 2.6G   |
| Sample_l4c1     | SRX202207 | ILLUMINA | PAIRED | 450     | SRR611082 | 13984884 | 2.6G   |
| Sample_l4l3     | SRX202210 | ILLUMINA | PAIRED | 450     | SRR611083 | 13964774 | 2.6G   |

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
cd ~/data/plastid/Atha_cross/

cat opts.tsv | # head -n 170 | #tail -n 10 |
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

