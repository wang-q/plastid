# *Arabidopsis thaliana* Col-0 倍数因子测试


[TOC levels=1-3]: # ""

- [*Arabidopsis thaliana* Col-0 倍数因子测试](#arabidopsis-thaliana-col-0-倍数因子测试)
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
* 测试文件 SRR616966
  * 覆盖度 4,970,359,200 / 119,667,750 = 41

```shell script
BASES=$(
    faops n50 -N 0 -S \
        ~/data/plastid/ena/SRR616966_1.fastq.gz \
        ~/data/plastid/ena/SRR616966_2.fastq.gz |
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
# 4970359200

echo ${GENOME}
# 119667750

bc <<< "${BASES} / ${GENOME}"
# 41

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/col_0
cd ~/data/plastid/evaluation/col_0

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${FOLD}

    mkdir -p ${BASE_NAME}/1_genome
    pushd ${BASE_NAME}/1_genome

    ln -fs ../../../../genome/col_0/genome.fa genome.fa
    popd

    mkdir -p ${BASE_NAME}/2_illumina
    pushd ${BASE_NAME}/2_illumina

    ln -fs ../../../../ena/SRR616966_1.fastq.gz R1.fq.gz
    ln -fs ../../../../ena/SRR616966_2.fastq.gz R2.fq.gz
    popd

done

```

## Trim and cutoff

```shell script
cd ~/data/plastid/evaluation/col_0

# 倍数因子::cutoff
ARRAY=()
DEPTH=41
for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    CUTOFF=$(bc <<< "(${DEPTH} * ${FOLD}) / 1")
    ARRAY+=("${FOLD}::${CUTOFF}")
done
echo "${ARRAY[@]}"
#0::0 0.25::10 0.5::20 1::41 2::82 4::164 8::328 16::656 32::1312 64::2624

for item in "${ARRAY[@]}" ; do
    FOLD="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"

    BASE_NAME=SRR616966_${FOLD}
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
cd ~/data/plastid/evaluation/col_0

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${FOLD}

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
cd ~/data/plastid/evaluation/col_0

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${FOLD}

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
cd ~/data/plastid/evaluation/col_0

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${FOLD}

    echo 1>&2 "==> ${BASE_NAME}"

    pushd ${BASE_NAME}/3_bwa > /dev/null

    echo -e "Fold\tchrom\n${FOLD}\tNc\n${FOLD}\tMt\n${FOLD}\tPt" |
        tsv-join -H --filter-file combine.tsv --key-fields chrom --append-fields 2-8

    popd > /dev/null

done |
    tsv-uniq \
    > SRR616966_folds.tsv

for PART in Nc Mt Pt; do
    cat SRR616966_folds.tsv |
        tsv-filter -H --str-eq chrom:${PART} |
        mlr --itsv --omd cat
    echo
    echo
done

```

| Fold | chrom | chrLength | covLength | covRate | bases      | mean  | min | max   |
|:-----|:------|:----------|:----------|:--------|:-----------|:------|:----|:------|
| 0    | Nc    | 119146348 | 118877472 | 0.9977  | 2503236815 | 21.01 | 0   | 22491 |
| 0.25 | Nc    | 119146348 | 118729184 | 0.9965  | 2487450084 | 20.88 | 0   | 22583 |
| 0.5  | Nc    | 119146348 | 78370964  | 0.6578  | 1156585513 | 9.71  | 0   | 22647 |
| 1    | Nc    | 119146348 | 13891866  | 0.1166  | 374063728  | 3.14  | 0   | 22624 |
| 2    | Nc    | 119146348 | 9239981   | 0.0776  | 316271886  | 2.65  | 0   | 22640 |
| 4    | Nc    | 119146348 | 6530203   | 0.0548  | 277059862  | 2.33  | 0   | 22676 |
| 8    | Nc    | 119146348 | 3868129   | 0.0325  | 219204159  | 1.84  | 0   | 22694 |
| 16   | Nc    | 119146348 | 2627606   | 0.0221  | 196228829  | 1.65  | 0   | 22597 |
| 32   | Nc    | 119146348 | 1661074   | 0.0139  | 177680732  | 1.49  | 0   | 22486 |
| 64   | Nc    | 119146348 | 1079260   | 0.0091  | 167425566  | 1.41  | 0   | 22137 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:----|
| 0    | Mt    | 366924    | 363702    | 0.9912  | 53109668 | 144.74 | 0   | 361 |
| 0.25 | Mt    | 366924    | 363588    | 0.9909  | 52973136 | 144.37 | 0   | 360 |
| 0.5  | Mt    | 366924    | 363544    | 0.9908  | 53028165 | 144.52 | 0   | 360 |
| 1    | Mt    | 366924    | 363651    | 0.9911  | 52978802 | 144.39 | 0   | 360 |
| 2    | Mt    | 366924    | 363283    | 0.9901  | 52967591 | 144.36 | 0   | 360 |
| 4    | Mt    | 366924    | 356203    | 0.9708  | 49636499 | 135.28 | 0   | 360 |
| 8    | Mt    | 366924    | 94504     | 0.2576  | 6758709  | 18.42  | 0   | 306 |
| 16   | Mt    | 366924    | 29041     | 0.0791  | 634438   | 1.73   | 0   | 250 |
| 32   | Mt    | 366924    | 25094     | 0.0684  | 547565   | 1.49   | 0   | 252 |
| 64   | Mt    | 366924    | 18043     | 0.0492  | 367465   | 1.00   | 0   | 251 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean    | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:--------|:----|:-----|
| 0    | Pt    | 154478    | 154478    | 1.0000  | 390507326 | 2527.92 | 4   | 3444 |
| 0.25 | Pt    | 154478    | 154478    | 1.0000  | 389932578 | 2524.19 | 4   | 3442 |
| 0.5  | Pt    | 154478    | 154478    | 1.0000  | 389906371 | 2524.03 | 4   | 3442 |
| 1    | Pt    | 154478    | 154478    | 1.0000  | 389937339 | 2524.23 | 4   | 3442 |
| 2    | Pt    | 154478    | 154478    | 1.0000  | 389843785 | 2523.62 | 4   | 3442 |
| 4    | Pt    | 154478    | 154478    | 1.0000  | 389872488 | 2523.81 | 4   | 3442 |
| 8    | Pt    | 154478    | 154478    | 1.0000  | 389847713 | 2523.65 | 4   | 3442 |
| 16   | Pt    | 154478    | 154478    | 1.0000  | 389820409 | 2523.47 | 4   | 3442 |
| 32   | Pt    | 154478    | 154478    | 1.0000  | 389713923 | 2522.78 | 2   | 3442 |
| 64   | Pt    | 154478    | 148548    | 0.9616  | 258962605 | 1676.37 | 0   | 3437 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/col_0

find . -type d -name "trim" | xargs rm -fr
find . -type f -path "*3_bwa/genome.fa*" | xargs rm
find . -type f -name "*.ba[mi]" | xargs rm
find . -type f -name "*.per-base.bed.gz" | xargs rm

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

