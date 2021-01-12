# *Arabidopsis thaliana* Ler-0 倍数因子测试


[TOC levels=1-3]: # ""

- [*Arabidopsis thaliana* Ler-0 倍数因子测试](#arabidopsis-thaliana-ler-0-倍数因子测试)
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
* 测试文件 SRR616965
  * 覆盖度 42

```shell script
BASES=$(
    faops n50 -N 0 -S \
        ~/data/plastid/ena/SRR616965_1.fastq.gz \
        ~/data/plastid/ena/SRR616965_2.fastq.gz |
    cut -f 2
)

GENOME=$(
    faops size ~/data/plastid/genome/col_0/genome.fa |
        cut -f 2 |
        perl -nl -e '
            $sum += $_;
            END {print $sum}
        '
)

echo ${BASES}
# 5087251000

echo ${GENOME}
# 119667750

bc <<< "${BASES} / ${GENOME}"
# 42

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/ler_0
cd ~/data/plastid/evaluation/ler_0

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616965_${FOLD}

    mkdir -p ${BASE_NAME}/1_genome
    pushd ${BASE_NAME}/1_genome

    ln -fs ../../../../genome/col_0/genome.fa genome.fa
    popd

    mkdir -p ${BASE_NAME}/2_illumina
    pushd ${BASE_NAME}/2_illumina

    ln -fs ../../../../ena/SRR616965_1.fastq.gz R1.fq.gz
    ln -fs ../../../../ena/SRR616965_2.fastq.gz R2.fq.gz
    popd

done

```

## Trim and cutoff

```shell script
cd ~/data/plastid/evaluation/ler_0

# 倍数因子::cutoff
ARRAY=()
DEPTH=42
for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    CUTOFF=$(bc <<< "(${DEPTH} * ${FOLD}) / 1")
    ARRAY+=("${FOLD}::${CUTOFF}")
done
echo "${ARRAY[@]}"
#0::0 0.25::10 0.5::21 1::42 2::84 4::168 8::336 16::672 32::1344 64::2688

for item in "${ARRAY[@]}" ; do
    FOLD="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"

    BASE_NAME=SRR616965_${FOLD}
    pushd ${BASE_NAME}

    rm *.sh

    if [[ "${FOLD}" == "0" ]]; then
        anchr template \
            --genome 119667750 \
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
            --genome 119667750 \
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
cd ~/data/plastid/evaluation/ler_0

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616965_${FOLD}

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
cd ~/data/plastid/evaluation/ler_0

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616965_${FOLD}

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
        (cat join.tsv | sed '2,6d' && cat) \
        > combine.tsv

    popd

done

```

## Merge all results

```shell script
cd ~/data/plastid/evaluation/ler_0

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616965_${FOLD}

    echo 1>&2 "==> ${BASE_NAME}"

    pushd ${BASE_NAME}/3_bwa > /dev/null

    echo -e "Fold\tchrom\n${FOLD}\tNc\n${FOLD}\tMt\n${FOLD}\tPt" |
        tsv-join -H --filter-file combine.tsv --key-fields chrom --append-fields 2-8

    popd > /dev/null

done |
    tsv-uniq \
    > SRR616965_folds.tsv

for PART in Nc Mt Pt; do
    cat SRR616965_folds.tsv |
        tsv-filter -H --str-eq chrom:${PART} |
        mlr --itsv --omd cat
    echo
    echo
done

```

| Fold | chrom | chrLength | covLength | covRate | bases      | mean  | min | max   |
|:-----|:------|:----------|:----------|:--------|:-----------|:------|:----|:------|
| 0    | Nc    | 119146348 | 112843627 | 0.9471  | 2393195441 | 20.09 | 0   | 23605 |
| 0.25 | Nc    | 119146348 | 112773888 | 0.9465  | 2380158986 | 19.98 | 0   | 23694 |
| 0.5  | Nc    | 119146348 | 50167841  | 0.4211  | 817795168  | 6.86  | 0   | 23632 |
| 1    | Nc    | 119146348 | 11881091  | 0.0997  | 383657385  | 3.22  | 0   | 23653 |
| 2    | Nc    | 119146348 | 8406746   | 0.0706  | 328370381  | 2.76  | 0   | 23656 |
| 4    | Nc    | 119146348 | 6090497   | 0.0511  | 290987442  | 2.44  | 0   | 23694 |
| 8    | Nc    | 119146348 | 4449987   | 0.0373  | 260503330  | 2.19  | 0   | 23550 |
| 16   | Nc    | 119146348 | 3065811   | 0.0257  | 182678949  | 1.53  | 0   | 23378 |
| 32   | Nc    | 119146348 | 2325855   | 0.0195  | 160846975  | 1.35  | 0   | 23116 |
| 64   | Nc    | 119146348 | 1876951   | 0.0158  | 152689486  | 1.28  | 0   | 22457 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:----------|:-------|:----|:----|
| 0    | Mt    | 366924    | 351472    | 0.9579  | 109519085 | 298.48 | 0   | 822 |
| 0.25 | Mt    | 366924    | 351476    | 0.9579  | 109367942 | 298.07 | 0   | 822 |
| 0.5  | Mt    | 366924    | 351297    | 0.9574  | 109310192 | 297.91 | 0   | 820 |
| 1    | Mt    | 366924    | 350967    | 0.9565  | 109281709 | 297.83 | 0   | 820 |
| 2    | Mt    | 366924    | 350815    | 0.9561  | 109308035 | 297.90 | 0   | 820 |
| 4    | Mt    | 366924    | 350948    | 0.9565  | 109197324 | 297.60 | 0   | 820 |
| 8    | Mt    | 366924    | 348567    | 0.9500  | 103864413 | 283.07 | 0   | 820 |
| 16   | Mt    | 366924    | 114147    | 0.3111  | 8703600   | 23.72  | 0   | 622 |
| 32   | Mt    | 366924    | 79795     | 0.2175  | 1378446   | 3.76   | 0   | 619 |
| 64   | Mt    | 366924    | 76285     | 0.2079  | 1374260   | 3.75   | 0   | 616 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean    | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:--------|:----|:-----|
| 0    | Pt    | 154478    | 154478    | 1.0000  | 767923755 | 4971.09 | 244 | 6898 |
| 0.25 | Pt    | 154478    | 154478    | 1.0000  | 767928954 | 4971.12 | 244 | 6927 |
| 0.5  | Pt    | 154478    | 154478    | 1.0000  | 767914496 | 4971.03 | 244 | 6968 |
| 1    | Pt    | 154478    | 154478    | 1.0000  | 767908001 | 4970.99 | 244 | 7018 |
| 2    | Pt    | 154478    | 154478    | 1.0000  | 767817464 | 4970.40 | 244 | 6957 |
| 4    | Pt    | 154478    | 154478    | 1.0000  | 767871785 | 4970.75 | 244 | 6949 |
| 8    | Pt    | 154478    | 154478    | 1.0000  | 767821009 | 4970.42 | 243 | 6945 |
| 16   | Pt    | 154478    | 154478    | 1.0000  | 767834842 | 4970.51 | 243 | 6970 |
| 32   | Pt    | 154478    | 154478    | 1.0000  | 767808163 | 4970.34 | 243 | 6924 |
| 64   | Pt    | 154478    | 154478    | 1.0000  | 767857310 | 4970.66 | 177 | 6968 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/ler_0

find . -type d -name "trim" | xargs rm -fr
find . -type f -path "*3_bwa/genome.fa*" | xargs rm
find . -type f -name "*.ba[mi]" | xargs rm
find . -type f -name "*.per-base.bed.gz" | xargs rm

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

