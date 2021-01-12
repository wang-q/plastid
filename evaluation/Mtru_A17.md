# *Medicago truncatula* A17 倍数因子测试


[TOC levels=1-3]: # ""

- [*Medicago truncatula* A17 倍数因子测试](#medicago-truncatula-a17-倍数因子测试)
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
* 测试文件 [SRR1542423](https://www.ncbi.nlm.nih.gov/sra/SRX673852)
  * 覆盖度 23

```shell script
BASES=$(
    faops n50 -N 0 -S \
        ~/data/plastid/ena/SRR1542423_1.fastq.gz \
        ~/data/plastid/ena/SRR1542423_2.fastq.gz |
    cut -f 2
)

GENOME=$(
    faops size ~/data/plastid/genome/a17/genome.fa |
        cut -f 2 |
        perl -nl -e '
            $sum += $_;
            END {print $sum}
        '
)

echo ${BASES}
# 8958357672

echo ${GENOME}
# 384862644

bc <<< "${BASES} / ${GENOME}"
# 23

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/a17
cd ~/data/plastid/evaluation/a17

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${FOLD}

    mkdir -p ${BASE_NAME}/1_genome
    pushd ${BASE_NAME}/1_genome

    ln -fs ../../../../genome/a17/genome.fa genome.fa
    popd

    mkdir -p ${BASE_NAME}/2_illumina
    pushd ${BASE_NAME}/2_illumina

    ln -fs ../../../../ena/SRR1542423_1.fastq.gz R1.fq.gz
    ln -fs ../../../../ena/SRR1542423_2.fastq.gz R2.fq.gz
    popd

done

```

## Trim and cutoff

```shell script
cd ~/data/plastid/evaluation/a17

# 倍数因子::cutoff
ARRAY=()
DEPTH=23
for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    CUTOFF=$(bc <<< "(${DEPTH} * ${FOLD}) / 1")
    ARRAY+=("${FOLD}::${CUTOFF}")
done
echo "${ARRAY[@]}"
#0::0 0.25::5 0.5::11 1::23 2::46 4::92 8::184 16::368 32::736 64::1472

for item in "${ARRAY[@]}" ; do
    FOLD="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"

    BASE_NAME=SRR1542423_${FOLD}
    pushd ${BASE_NAME}

    rm *.sh

    if [[ "${FOLD}" == "0" ]]; then
        anchr template \
            --genome 384862644 \
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
            --genome 384862644 \
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
cd ~/data/plastid/evaluation/a17

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${FOLD}

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
cd ~/data/plastid/evaluation/a17

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${FOLD}

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
        (cat join.tsv | sed '2,9d' && cat) \
        > combine.tsv

    popd

done

```

## Merge all results

```shell script
cd ~/data/plastid/evaluation/a17

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${FOLD}

    echo 1>&2 "==> ${BASE_NAME}"

    pushd ${BASE_NAME}/3_bwa > /dev/null

    echo -e "Fold\tchrom\n${FOLD}\tNc\n${FOLD}\tMt\n${FOLD}\tPt" |
        tsv-join -H --filter-file combine.tsv --key-fields chrom --append-fields 2-8

    popd > /dev/null

done |
    tsv-uniq \
    > SRR1542423_folds.tsv

for PART in Nc Mt Pt; do
    cat SRR1542423_folds.tsv |
        tsv-filter -H --str-eq chrom:${PART} |
        mlr --itsv --omd cat
    echo
    echo
done

```

| Fold | chrom | chrLength | covLength | covRate | bases      | mean | min | max    |
|:-----|:------|:----------|:----------|:--------|:-----------|:-----|:----|:-------|
| 0    | Nc    | 384466993 | 348072517 | 0.9053  | 3237657237 | 8.42 | 0   | 107159 |
| 0.25 | Nc    | 384466993 | 319243960 | 0.8304  | 3077659704 | 8.01 | 0   | 92341  |
| 0.5  | Nc    | 384466993 | 197683490 | 0.5142  | 2210392114 | 5.75 | 0   | 92229  |
| 1    | Nc    | 384466993 | 76576135  | 0.1992  | 1292591729 | 3.36 | 0   | 91636  |
| 2    | Nc    | 384466993 | 51550888  | 0.1341  | 1079035637 | 2.81 | 0   | 108024 |
| 4    | Nc    | 384466993 | 39755248  | 0.1034  | 942208106  | 2.45 | 0   | 89569  |
| 8    | Nc    | 384466993 | 31023999  | 0.0807  | 848499067  | 2.21 | 0   | 108386 |
| 16   | Nc    | 384466993 | 23226239  | 0.0604  | 720018690  | 1.87 | 0   | 87693  |
| 32   | Nc    | 384466993 | 16235103  | 0.0422  | 618721163  | 1.61 | 0   | 103831 |
| 64   | Nc    | 384466993 | 11134978  | 0.0290  | 507222214  | 1.32 | 0   | 105054 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:----|
| 0    | Mt    | 271618    | 271618    | 1.0000  | 49932049 | 183.83 | 1   | 726 |
| 0.25 | Mt    | 271618    | 271618    | 1.0000  | 49142557 | 180.93 | 1   | 721 |
| 0.5  | Mt    | 271618    | 271618    | 1.0000  | 49020706 | 180.48 | 1   | 721 |
| 1    | Mt    | 271618    | 271618    | 1.0000  | 48878964 | 179.95 | 1   | 720 |
| 2    | Mt    | 271618    | 271618    | 1.0000  | 48861022 | 179.89 | 1   | 718 |
| 4    | Mt    | 271618    | 271495    | 0.9995  | 48383733 | 178.13 | 0   | 717 |
| 8    | Mt    | 271618    | 253902    | 0.9348  | 38889573 | 143.18 | 0   | 708 |
| 16   | Mt    | 271618    | 95498     | 0.3516  | 9127266  | 33.60  | 0   | 706 |
| 32   | Mt    | 271618    | 21853     | 0.0805  | 578709   | 2.13   | 0   | 571 |
| 64   | Mt    | 271618    | 14959     | 0.0551  | 40418    | 0.15   | 0   | 127 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Pt    | 124033    | 124031    | 1.0000  | 80757838 | 651.10 | 0   | 3314 |
| 0.25 | Pt    | 124033    | 124032    | 1.0000  | 80570481 | 649.59 | 0   | 3313 |
| 0.5  | Pt    | 124033    | 124029    | 1.0000  | 80085169 | 645.68 | 0   | 3292 |
| 1    | Pt    | 124033    | 124032    | 1.0000  | 79931815 | 644.44 | 0   | 3285 |
| 2    | Pt    | 124033    | 124032    | 1.0000  | 79714333 | 642.69 | 0   | 3285 |
| 4    | Pt    | 124033    | 124024    | 0.9999  | 79818760 | 643.53 | 0   | 3281 |
| 8    | Pt    | 124033    | 123963    | 0.9994  | 79719423 | 642.73 | 0   | 3280 |
| 16   | Pt    | 124033    | 123771    | 0.9979  | 79182159 | 638.40 | 0   | 3281 |
| 32   | Pt    | 124033    | 117909    | 0.9506  | 75454137 | 608.34 | 0   | 3284 |
| 64   | Pt    | 124033    | 82647     | 0.6663  | 54141806 | 436.51 | 0   | 3278 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/a17

find . -type d -name "trim" | xargs rm -fr
find . -type f -path "*3_bwa/genome.fa*" | xargs rm
find . -type f -name "*.ba[mi]" | xargs rm
find . -type f -name "*.per-base.bed.gz" | xargs rm

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

