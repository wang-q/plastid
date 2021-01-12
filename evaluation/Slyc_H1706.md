# *Solanum lycopersicum* Heinz 1706 倍数因子测试


[TOC levels=1-3]: # ""

- [*Solanum lycopersicum* Heinz 1706 倍数因子测试](#solanum-lycopersicum-heinz-1706-倍数因子测试)
  - [基本信息](#基本信息)
  - [Symlink](#symlink)
  - [Trim and cutoff](#trim-and-cutoff)
  - [`kat hist` and `kat gcp`](#kat-hist-and-kat-gcp)
  - [Depth](#depth)
  - [Merge all results](#merge-all-results)
  - [Remove intermediate files](#remove-intermediate-files)

## 基本信息

`cutoff = 倍数因子 * 覆盖深度`

* 因子值 0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64
* 测试文件 [SRR1572628](https://www.ncbi.nlm.nih.gov/sra/SRX698770)
  * 覆盖度 5

```shell script
BASES=$(
    faops n50 -N 0 -S \
        ~/data/plastid/ena/SRR1572628_1.fastq.gz \
        ~/data/plastid/ena/SRR1572628_2.fastq.gz |
    cut -f 2
)

GENOME=$(
    faops size ~/data/plastid/genome/h1706/genome.fa |
        cut -f 2 |
        perl -nl -e '
            $sum += $_;
            END {print $sum}
        '
)

echo ${BASES}
# 4839669000

echo ${GENOME}
# 807826382

bc <<< "${BASES} / ${GENOME}"
# 5

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/h1706
cd ~/data/plastid/evaluation/h1706

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1572628_${FOLD}

    mkdir -p ${BASE_NAME}/1_genome
    pushd ${BASE_NAME}/1_genome

    ln -fs ../../../../genome/h1706/genome.fa genome.fa
    popd

    mkdir -p ${BASE_NAME}/2_illumina
    pushd ${BASE_NAME}/2_illumina

    ln -fs ../../../../ena/SRR1572628_1.fastq.gz R1.fq.gz
    ln -fs ../../../../ena/SRR1572628_2.fastq.gz R2.fq.gz
    popd

done

```

## Trim and cutoff

```shell script
cd ~/data/plastid/evaluation/h1706

# 倍数因子::cutoff
ARRAY=()
DEPTH=5
for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    CUTOFF=$(bc <<< "(${DEPTH} * ${FOLD}) / 1")
    ARRAY+=("${FOLD}::${CUTOFF}")
done
echo "${ARRAY[@]}"
#0::0 0.25::1 0.5::2 1::5 2::10 4::20 8::40 16::80 32::160 64::320

for item in "${ARRAY[@]}" ; do
    FOLD="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"

    BASE_NAME=SRR1572628_${FOLD}
    pushd ${BASE_NAME}

    rm *.sh

    if [[ "${FOLD}" == "0" ]]; then
        anchr template \
            --genome 807826382 \
            --parallel 24 \
            --xmx 80g \
            \
            --fastqc \
            --insertsize \
            --kat \
            \
            --trim "--dedupe" \
            --qual "25" \
            --len "60" \
            --filter "adapter artifact" \
            \
            --bwa Q25L60

        bsub -q mpi -n 24 -J "${BASE_NAME}" "
            bash 2_fastqc.sh
            bash 2_insert_size.sh
            bash 2_kat.sh
            bash 2_trim.sh
            bash 9_stat_reads.sh
            bash 3_bwa.sh
        "
    else
        anchr template \
            --genome 807826382 \
            --parallel 24 \
            --xmx 80g \
            \
            --trim "--dedupe --cutoff ${CUTOFF} --cutk 31" \
            --qual "25" \
            --len "60" \
            --filter "adapter artifact" \
            \
            --bwa Q25L60

        bsub -q mpi -n 24 -J "${BASE_NAME}" "
            bash 2_trim.sh
            bash 9_stat_reads.sh
            bash 3_bwa.sh
        "
    fi

    popd

done

```

## `kat hist` and `kat gcp`

```shell script
cd ~/data/plastid/evaluation/h1706

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1572628_${FOLD}

    mkdir -p ${BASE_NAME}/kat
    pushd ${BASE_NAME}/kat

    bsub -q mpi -n 24 -J "${BASE_NAME}-kat" "
        kat hist \
            -t 24 -m 31 \
            ../2_illumina/trim/Q25L60/R1.fq.gz \
            ../2_illumina/trim/Q25L60/R2.fq.gz \
            ../2_illumina/trim/Q25L60/Rs.fq.gz \
            -o R-hist-31

        kat gcp \
            -t 24 -m 31 \
            ../2_illumina/trim/Q25L60/R1.fq.gz \
            ../2_illumina/trim/Q25L60/R2.fq.gz \
            ../2_illumina/trim/Q25L60/Rs.fq.gz \
            -o R-gcp-31
    "

    popd

done

```

## Depth


```shell script
cd ~/data/plastid/evaluation/h1706

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1572628_${FOLD}

    echo 1>&2 "==> ${BASE_NAME}"

    pushd ${BASE_NAME}/3_bwa

    cat join.tsv |
        grep -v "^Mt" |
        grep -v "^Pt" |
        tsv-summarize -H --sum chrLength,covLength,bases --min min --max max |
        sed '1d' |
        perl -e '
            my $line = <>;
            chomp $line;
            my ($chrLength, $covLength, $bases, $min, $max, ) = split qq(\t), $line;
            my $covRate = sprintf qq(%.4f), $covLength / $chrLength;
            my $mean = sprintf qq(%.2f), $bases / $chrLength;
            print join qq(\t), (
                "Nc", $chrLength, $covLength, $covRate, $bases, $mean, $min, $max,
            );
            print qq(\n);
        ' |
        (cat join.tsv | sed '2,13d' && cat) \
        > combine.tsv

    popd

done

```

## Merge all results

```shell script
cd ~/data/plastid/evaluation/h1706

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1572628_${FOLD}

    echo 1>&2 "==> ${BASE_NAME}"

    pushd ${BASE_NAME}/3_bwa > /dev/null

    echo -e "Fold\tchrom\n${FOLD}\tNc\n${FOLD}\tMt\n${FOLD}\tPt" |
        tsv-join -H --filter-file combine.tsv --key-fields chrom --append-fields 2-8

    popd > /dev/null

done |
    tsv-uniq \
    > SRR1572628_folds.tsv

for PART in Nc Mt Pt; do
    cat SRR1572628_folds.tsv |
        tsv-filter -H --str-eq chrom:${PART} |
        mlr --itsv --omd cat
    echo
    echo
done

```

| Fold | chrom | chrLength | covLength | covRate | bases      | mean | min | max   |
|:-----|:------|:----------|:----------|:--------|:-----------|:-----|:----|:------|
| 0    | Nc    | 807224664 | 687361649 | 0.8515  | 3180562052 | 3.94 | 0   | 20886 |
| 0.25 | Nc    | 807224664 | 687367973 | 0.8515  | 3180568663 | 3.94 | 0   | 20936 |
| 0.5  | Nc    | 807224664 | 678841856 | 0.8410  | 3165122034 | 3.92 | 0   | 20924 |
| 1    | Nc    | 807224664 | 512507255 | 0.6349  | 2438099474 | 3.02 | 0   | 21022 |
| 2    | Nc    | 807224664 | 204747604 | 0.2536  | 1110973726 | 1.38 | 0   | 21151 |
| 4    | Nc    | 807224664 | 137292726 | 0.1701  | 819204674  | 1.01 | 0   | 20919 |
| 8    | Nc    | 807224664 | 111626374 | 0.1383  | 699778951  | 0.87 | 0   | 20716 |
| 16   | Nc    | 807224664 | 92796197  | 0.1150  | 605505425  | 0.75 | 0   | 20469 |
| 32   | Nc    | 807224664 | 75908434  | 0.0940  | 515050993  | 0.64 | 0   | 20284 |
| 64   | Nc    | 807224664 | 61797609  | 0.0766  | 436774399  | 0.54 | 0   | 20279 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean  | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:------|:----|:----|
| 0    | Mt    | 446257    | 446257    | 1.0000  | 43828430 | 98.21 | 2   | 828 |
| 0.25 | Mt    | 446257    | 446257    | 1.0000  | 43807272 | 98.17 | 1   | 799 |
| 0.5  | Mt    | 446257    | 446255    | 1.0000  | 43837766 | 98.23 | 0   | 858 |
| 1    | Mt    | 446257    | 446257    | 1.0000  | 43778626 | 98.10 | 2   | 846 |
| 2    | Mt    | 446257    | 446254    | 1.0000  | 43750446 | 98.04 | 0   | 902 |
| 4    | Mt    | 446257    | 446257    | 1.0000  | 43804934 | 98.16 | 1   | 875 |
| 8    | Mt    | 446257    | 446257    | 1.0000  | 43789320 | 98.13 | 3   | 846 |
| 16   | Mt    | 446257    | 436794    | 0.9788  | 39605474 | 88.75 | 0   | 861 |
| 32   | Mt    | 446257    | 209356    | 0.4691  | 13305819 | 29.82 | 0   | 840 |
| 64   | Mt    | 446257    | 31133     | 0.0698  | 1332256  | 2.99  | 0   | 863 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:-------|:----|:-----|
| 0    | Pt    | 155461    | 155461    | 1.0000  | 132760191 | 853.98 | 3   | 2245 |
| 0.25 | Pt    | 155461    | 155461    | 1.0000  | 132746965 | 853.89 | 3   | 2245 |
| 0.5  | Pt    | 155461    | 155461    | 1.0000  | 132756549 | 853.95 | 5   | 2245 |
| 1    | Pt    | 155461    | 155461    | 1.0000  | 132674798 | 853.43 | 2   | 2243 |
| 2    | Pt    | 155461    | 155461    | 1.0000  | 132594962 | 852.91 | 4   | 2241 |
| 4    | Pt    | 155461    | 155461    | 1.0000  | 132593950 | 852.91 | 9   | 2241 |
| 8    | Pt    | 155461    | 155461    | 1.0000  | 132612435 | 853.03 | 8   | 2241 |
| 16   | Pt    | 155461    | 155461    | 1.0000  | 132611693 | 853.02 | 4   | 2241 |
| 32   | Pt    | 155461    | 155461    | 1.0000  | 132639793 | 853.20 | 5   | 2241 |
| 64   | Pt    | 155461    | 155461    | 1.0000  | 132292659 | 850.97 | 6   | 2241 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/h1706

find . -type d -name "trim" | xargs rm -fr
find . -type f -path "*3_bwa/genome.fa*" | xargs rm
find . -type f -name "*.ba[mi]" | xargs rm
find . -type f -name "*.per-base.bed.gz" | xargs rm

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

