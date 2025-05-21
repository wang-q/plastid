# *Medicago truncatula* Hapmap Project

<!-- TOC -->
* [*Medicago truncatula* Hapmap Project](#medicago-truncatula-hapmap-project)
  * [基本信息](#基本信息)
  * [项目信息](#项目信息)
  * [其他可能可用的项目](#其他可能可用的项目)
  * [数据下载](#数据下载)
    * [Reference](#reference)
    * [Illumina](#illumina)
  * [Symlink](#symlink)
  * [Run](#run)
  * [Pack and clean](#pack-and-clean)
  * [VCF](#vcf)
<!-- TOC -->

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
    hnsm filter -s stdin |
    hnsm replace stdin replace.tsv |
    hnsm order stdin <(echo Pt; echo Mt) -o genome.fa

```

### Illumina

* Download `Metadata` from NCBI SRA Run Selector via a web browser
    * https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP001874
    * Save it to `SraRunTable.csv`

* https://medicago.legumeinfo.org/tools/germplasm/
    * Save the page via browser
    * Extract table via `pup`
    * Convert xls to csv via `excel`

```shell script
mkdir -p ~/data/plastid/Mtru_384/ena
cd ~/data/plastid/Mtru_384/ena

cat SraRunTable.csv |
    mlr --icsv --otsv cat \
    > SraRunTable.tsv

 cat 'Germplasm – Medicago Analysis Portal.html' |
    pup 'table#germplasm-datatable' \
    > germplasm.xls

# Convert xls(html) to csv via `excel`
cat germplasm.csv | wc -l
#339

cat germplasm.csv |
    head -n 11 |
    mlr --icsv --otsv cat |
    rgr md stdin --num

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
    tsv-select -H -f Experiment,ID,"Country" |
    mlr --itsv --ocsv cat |
    sed 1d \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml --ascp

rgr md ena_info.tsv --fmt |
    head -n 20

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 2 "{}"

```

| ID    | Line | Accession | Country | Seeds Origin     | Latitude | Longitude |
|-------|------|-----------|---------|------------------|---------:|----------:|
| HM001 | CC8  | SA22322   | Syria   | INRA-Montpellier |   35.017 |      37.1 |
| HM002 | CC8  | SA28064   | Cyprus  | INRA-Montpellier |   34.783 |    33.167 |
| HM003 | CC8  | ESP105-L  | Spain   | INRA-Montpellier |   38.076 |    -3.816 |
| HM004 | CC8  | DZA045-6  | Algeria | INRA-Montpellier |   36.923 |     7.736 |
| HM005 | CC8  | DZA315-16 | Algeria | INRA-Montpellier |   34.716 |     0.158 |
| HM006 | CC8  | F83005-5  | France  | INRA-Montpellier |   43.571 |     6.224 |
| HM007 | CC8  | Salses71B | France  | INRA-Montpellier |    42.82 |     2.945 |
| HM008 | CC8  | DZA012-J  | Algeria | INRA-Montpellier |   36.549 |     3.183 |
| HM009 | CC16 | GRC020-B  | Greece  | INRA-Montpellier |   38.122 |    21.543 |
| HM010 | CC16 | SA24714   | Italy   | INRA-Montpellier |   37.533 |    14.517 |

| name  | srx       | platform | layout | ilength | srr        |      spots | bases |
|-------|-----------|----------|--------|--------:|------------|-----------:|-------|
| HM001 | SRX375894 | ILLUMINA | PAIRED |     257 | SRR1034054 | 18,738,482 | 3.14G |
| HM002 | SRX375896 | ILLUMINA | PAIRED |     219 | SRR1034056 | 20,086,982 | 3.37G |
| HM003 | SRX375917 | ILLUMINA | PAIRED |     283 | SRR1034077 | 18,599,526 | 3.12G |
| HM004 | SRX375905 | ILLUMINA | PAIRED |     273 | SRR1034065 | 19,031,676 | 3.19G |
| HM005 | SRX375923 | ILLUMINA | PAIRED |     247 | SRR1034083 | 18,903,805 | 3.17G |
| HM007 | SRX375930 | ILLUMINA | PAIRED |     229 | SRR1034090 | 22,213,827 | 3.72G |
| HM008 | SRX375937 | ILLUMINA | PAIRED |     214 | SRR1034097 | 21,766,454 | 3.65G |
| HM011 | SRX375943 | ILLUMINA | PAIRED |     228 | SRR1034103 | 15,624,226 | 2.62G |
| HM012 | SRX375962 | ILLUMINA | PAIRED |     351 | SRR1034122 | 16,190,747 | 2.71G |
| HM014 | SRX375910 | ILLUMINA | PAIRED |     293 | SRR1034070 | 29,010,895 | 4.86G |
| HM015 | SRX375948 | ILLUMINA | PAIRED |     298 | SRR1034108 | 20,171,851 | 3.38G |
| HM016 | SRX375953 | ILLUMINA | PAIRED |     299 | SRR1034113 | 22,543,648 | 3.78G |
| HM017 | SRX376069 | ILLUMINA | PAIRED |     256 | SRR1034229 | 18,409,555 | 3.09G |
| HM019 | SRX375978 | ILLUMINA | PAIRED |     341 | SRR1034138 | 17,367,340 | 2.91G |
| HM020 | SRX376071 | ILLUMINA | PAIRED |     241 | SRR1034231 | 17,181,113 | 2.88G |
| HM022 | SRX376075 | ILLUMINA | PAIRED |     266 | SRR1034235 | 18,096,959 | 3.03G |
| HM026 | SRX375990 | ILLUMINA | PAIRED |     315 | SRR1034150 | 12,887,454 | 2.16G |
| HM027 | SRX376004 | ILLUMINA | PAIRED |     292 | SRR1034164 | 20,996,668 | 3.52G |

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

