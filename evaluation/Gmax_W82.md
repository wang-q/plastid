# *Glycine max* Williams 82 倍数因子测试


[TOC levels=1-3]: # ""

- [*Glycine max* Williams 82 倍数因子测试](#glycine-max-williams-82-倍数因子测试)
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
* 测试文件 [SRR10296600](https://www.ncbi.nlm.nih.gov/sra/SRX7009428)
  * 覆盖度 51

```shell script
BASES=$(
    faops n50 -N 0 -S \
        ~/data/plastid/ena/SRR10296600_1.fastq.gz \
        ~/data/plastid/ena/SRR10296600_2.fastq.gz |
    cut -f 2
)

GENOME=$(
    faops size ~/data/plastid/genome/w82/genome.fa |
        cut -f 2 |
        perl -nl -e '
            $sum += $_;
            END {print $sum}
        '
)

echo ${BASES}
# 48633106500

echo ${GENOME}
# 949738161

bc <<< "${BASES} / ${GENOME}"
# 51

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/w82
cd ~/data/plastid/evaluation/w82

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${FOLD}

    mkdir -p ${BASE_NAME}/1_genome
    pushd ${BASE_NAME}/1_genome

    ln -fs ../../../../genome/w82/genome.fa genome.fa
    popd

    mkdir -p ${BASE_NAME}/2_illumina
    pushd ${BASE_NAME}/2_illumina

    ln -fs ../../../../ena/SRR10296600_1.fastq.gz R1.fq.gz
    ln -fs ../../../../ena/SRR10296600_2.fastq.gz R2.fq.gz
    popd

done

```

## Trim and cutoff

```shell script
cd ~/data/plastid/evaluation/w82

# 倍数因子::cutoff
ARRAY=()
DEPTH=51
for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    CUTOFF=$(bc <<< "(${DEPTH} * ${FOLD}) / 1")
    ARRAY+=("${FOLD}::${CUTOFF}")
done
echo "${ARRAY[@]}"
#0::0 0.25::12 0.5::25 1::51 2::102 4::204 8::408 16::816 32::1632 64::3264

for item in "${ARRAY[@]}" ; do
    FOLD="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"

    BASE_NAME=SRR10296600_${FOLD}
    pushd ${BASE_NAME}

    rm *.sh

    if [[ "${FOLD}" == "0" ]]; then
        anchr template \
            --genome 949738161 \
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
            --genome 949738161 \
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
cd ~/data/plastid/evaluation/w82

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${FOLD}

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
cd ~/data/plastid/evaluation/w82

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${FOLD}

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
        (cat join.tsv | sed '2,21d' && cat) \
        > combine.tsv

    popd

done

```

## Merge all results

```shell script
cd ~/data/plastid/evaluation/w82

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${FOLD}

    echo 1>&2 "==> ${BASE_NAME}"

    pushd ${BASE_NAME}/3_bwa > /dev/null

    echo -e "Fold\tchrom\n${FOLD}\tNc\n${FOLD}\tMt\n${FOLD}\tPt" |
        tsv-join -H --filter-file combine.tsv --key-fields chrom --append-fields 2-8

    popd > /dev/null

done |
    tsv-uniq \
    > SRR10296600_folds.tsv

for PART in Nc Mt Pt; do
    cat SRR10296600_folds.tsv |
        tsv-filter -H --str-eq chrom:${PART} |
        mlr --itsv --omd cat
    echo
    echo
done

```

| Fold | chrom | chrLength | covLength | covRate | bases       | mean  | min | max   |
|:-----|:------|:----------|:----------|:--------|:------------|:------|:----|:------|
| 0    | Nc    | 949183385 | 929332756 | 0.9791  | 31776930265 | 33.48 | 0   | 51488 |
| 0.25 | Nc    | 949183385 | 929081317 | 0.9788  | 31745893291 | 33.45 | 0   | 51456 |
| 0.5  | Nc    | 949183385 | 919354021 | 0.9686  | 30562383338 | 32.20 | 0   | 51544 |
| 1    | Nc    | 949183385 | 492738118 | 0.5191  | 14345206570 | 15.11 | 0   | 51455 |
| 2    | Nc    | 949183385 | 374329757 | 0.3944  | 11672184468 | 12.30 | 0   | 51233 |
| 4    | Nc    | 949183385 | 330601251 | 0.3483  | 10280230761 | 10.83 | 0   | 50771 |
| 8    | Nc    | 949183385 | 297401067 | 0.3133  | 9086488145  | 9.57  | 0   | 50382 |
| 16   | Nc    | 949183385 | 265816257 | 0.2800  | 7900059739  | 8.32  | 0   | 49239 |
| 32   | Nc    | 949183385 | 231520544 | 0.2439  | 6652082721  | 7.01  | 0   | 47417 |
| 64   | Nc    | 949183385 | 195752318 | 0.2062  | 5455695707  | 5.75  | 0   | 44196 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:-------|:----|:-----|
| 0    | Mt    | 402558    | 388550    | 0.9652  | 344439689 | 855.63 | 0   | 2498 |
| 0.25 | Mt    | 402558    | 388553    | 0.9652  | 344086398 | 854.75 | 0   | 2494 |
| 0.5  | Mt    | 402558    | 388555    | 0.9652  | 344053953 | 854.67 | 0   | 2494 |
| 1    | Mt    | 402558    | 388453    | 0.9650  | 344072197 | 854.71 | 0   | 2494 |
| 2    | Mt    | 402558    | 388451    | 0.9650  | 344109580 | 854.81 | 0   | 2493 |
| 4    | Mt    | 402558    | 388553    | 0.9652  | 344052540 | 854.67 | 0   | 2495 |
| 8    | Mt    | 402558    | 388455    | 0.9650  | 344016392 | 854.58 | 0   | 2496 |
| 16   | Mt    | 402558    | 360127    | 0.8946  | 285101256 | 708.22 | 0   | 2494 |
| 32   | Mt    | 402558    | 164183    | 0.4078  | 52219150  | 129.72 | 0   | 2199 |
| 64   | Mt    | 402558    | 96723     | 0.2403  | 6420182   | 15.95  | 0   | 1480 |


| Fold | chrom | chrLength | covLength | covRate | bases      | mean     | min | max   |
|:-----|:------|:----------|:----------|:--------|:-----------|:---------|:----|:------|
| 0    | Pt    | 152218    | 152218    | 1.0000  | 2061751537 | 13544.73 | 18  | 17319 |
| 0.25 | Pt    | 152218    | 152218    | 1.0000  | 2061427561 | 13542.60 | 21  | 17318 |
| 0.5  | Pt    | 152218    | 152218    | 1.0000  | 2060905315 | 13539.17 | 23  | 17307 |
| 1    | Pt    | 152218    | 152218    | 1.0000  | 2060498502 | 13536.50 | 19  | 17305 |
| 2    | Pt    | 152218    | 152218    | 1.0000  | 2060143096 | 13534.16 | 18  | 17305 |
| 4    | Pt    | 152218    | 152218    | 1.0000  | 2060378268 | 13535.71 | 17  | 17304 |
| 8    | Pt    | 152218    | 152218    | 1.0000  | 2060293624 | 13535.15 | 18  | 17305 |
| 16   | Pt    | 152218    | 152218    | 1.0000  | 2060104592 | 13533.91 | 20  | 17304 |
| 32   | Pt    | 152218    | 152218    | 1.0000  | 2060315413 | 13535.29 | 19  | 17304 |
| 64   | Pt    | 152218    | 152218    | 1.0000  | 2060234061 | 13534.76 | 15  | 17305 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/w82

find . -type d -name "trim" | xargs rm -fr
find . -type f -path "*3_bwa/genome.fa*" | xargs rm
find . -type f -name "*.ba[mi]" | xargs rm
find . -type f -name "*.per-base.bed.gz" | xargs rm

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

