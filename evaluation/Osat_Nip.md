# *Oryza sativa* Nipponbare 倍数因子测试


[TOC levels=1-3]: # ""

- [*Oryza sativa* Nipponbare 倍数因子测试](#oryza-sativa-nipponbare-倍数因子测试)
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
* 测试文件 [SRR545231](https://www.ncbi.nlm.nih.gov/sra/SRX179254)
  * 覆盖度 46

```shell script
BASES=$(
    faops n50 -N 0 -S \
        ~/data/plastid/ena/SRR545231_1.fastq.gz \
        ~/data/plastid/ena/SRR545231_2.fastq.gz |
    cut -f 2
)

GENOME=$(
    faops size ~/data/plastid/genome/nip/genome.fa |
        cut -f 2 |
        perl -nl -e '
            $sum += $_;
            END {print $sum}
        '
)

echo ${BASES}
# 17220721594

echo ${GENOME}
# 373870564

bc <<< "${BASES} / ${GENOME}"
# 46

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/nip
cd ~/data/plastid/evaluation/nip

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${FOLD}

    mkdir -p ${BASE_NAME}/1_genome
    pushd ${BASE_NAME}/1_genome

    ln -fs ../../../../genome/nip/genome.fa genome.fa
    popd

    mkdir -p ${BASE_NAME}/2_illumina
    pushd ${BASE_NAME}/2_illumina

    ln -fs ../../../../ena/SRR545231_1.fastq.gz R1.fq.gz
    ln -fs ../../../../ena/SRR545231_2.fastq.gz R2.fq.gz
    popd

done

```

## Trim and cutoff

```shell script
cd ~/data/plastid/evaluation/nip

# 倍数因子::cutoff
ARRAY=()
DEPTH=46
for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    CUTOFF=$(bc <<< "(${DEPTH} * ${FOLD}) / 1")
    ARRAY+=("${FOLD}::${CUTOFF}")
done
echo "${ARRAY[@]}"
#0::0 0.25::11 0.5::23 1::46 2::92 4::184 8::368 16::736 32::1472 64::2944

for item in "${ARRAY[@]}" ; do
    FOLD="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"

    BASE_NAME=SRR545231_${FOLD}
    pushd ${BASE_NAME}

    rm *.sh

    if [[ "${FOLD}" == "0" ]]; then
        anchr template \
            --genome 373870564 \
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
            --genome 373870564 \
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
cd ~/data/plastid/evaluation/nip

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${FOLD}

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
cd ~/data/plastid/evaluation/nip

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${FOLD}

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
cd ~/data/plastid/evaluation/nip

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${FOLD}

    echo 1>&2 "==> ${BASE_NAME}"

    pushd ${BASE_NAME}/3_bwa > /dev/null

    echo -e "Fold\tchrom\n${FOLD}\tNc\n${FOLD}\tMt\n${FOLD}\tPt" |
        tsv-join -H --filter-file combine.tsv --key-fields chrom --append-fields 2-8

    popd > /dev/null

done |
    tsv-uniq \
    > SRR545231_folds.tsv

for PART in Nc Mt Pt; do
    cat SRR545231_folds.tsv |
        tsv-filter -H --str-eq chrom:${PART} |
        mlr --itsv --omd cat
    echo
    echo
done

```

| Fold | chrom | chrLength | covLength | covRate | bases       | mean  | min | max  |
|:-----|:------|:----------|:----------|:--------|:------------|:------|:----|:-----|
| 0    | Nc    | 373245519 | 368402032 | 0.9870  | 12983797798 | 34.79 | 0   | 8186 |
| 0.25 | Nc    | 373245519 | 366583918 | 0.9822  | 12964540235 | 34.73 | 0   | 8170 |
| 0.5  | Nc    | 373245519 | 356144812 | 0.9542  | 12665609206 | 33.93 | 0   | 8173 |
| 1    | Nc    | 373245519 | 174199288 | 0.4667  | 4466878641  | 11.97 | 0   | 8192 |
| 2    | Nc    | 373245519 | 120942207 | 0.3240  | 3151703014  | 8.44  | 0   | 8172 |
| 4    | Nc    | 373245519 | 98773643  | 0.2646  | 2556632934  | 6.85  | 0   | 8140 |
| 8    | Nc    | 373245519 | 82490186  | 0.2210  | 2129813140  | 5.71  | 0   | 8140 |
| 16   | Nc    | 373245519 | 67192387  | 0.1800  | 1740850801  | 4.66  | 0   | 7965 |
| 32   | Nc    | 373245519 | 50568041  | 0.1355  | 1272065380  | 3.41  | 0   | 7797 |
| 64   | Nc    | 373245519 | 33468069  | 0.0897  | 798576696   | 2.14  | 0   | 7212 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:----|
| 0    | Mt    | 490520    | 488459    | 0.9958  | 62088079 | 126.58 | 0   | 465 |
| 0.25 | Mt    | 490520    | 488638    | 0.9962  | 62070637 | 126.54 | 0   | 530 |
| 0.5  | Mt    | 490520    | 488594    | 0.9961  | 62057759 | 126.51 | 0   | 489 |
| 1    | Mt    | 490520    | 488642    | 0.9962  | 62097034 | 126.59 | 0   | 466 |
| 2    | Mt    | 490520    | 488581    | 0.9960  | 61976607 | 126.35 | 0   | 484 |
| 4    | Mt    | 490520    | 432170    | 0.8810  | 46324406 | 94.44  | 0   | 501 |
| 8    | Mt    | 490520    | 88134     | 0.1797  | 7604591  | 15.50  | 0   | 458 |
| 16   | Mt    | 490520    | 43163     | 0.0880  | 4466508  | 9.11   | 0   | 498 |
| 32   | Mt    | 490520    | 33508     | 0.0683  | 3760523  | 7.67   | 0   | 447 |
| 64   | Mt    | 490520    | 6881      | 0.0140  | 486583   | 0.99   | 0   | 405 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Pt    | 134525    | 132943    | 0.9882  | 56062147 | 416.74 | 0   | 1746 |
| 0.25 | Pt    | 134525    | 132903    | 0.9879  | 56268862 | 418.28 | 0   | 1746 |
| 0.5  | Pt    | 134525    | 132771    | 0.9870  | 56165888 | 417.51 | 0   | 1746 |
| 1    | Pt    | 134525    | 132892    | 0.9879  | 56113124 | 417.12 | 0   | 1746 |
| 2    | Pt    | 134525    | 132957    | 0.9883  | 56090073 | 416.95 | 0   | 1746 |
| 4    | Pt    | 134525    | 132841    | 0.9875  | 56059193 | 416.72 | 0   | 1746 |
| 8    | Pt    | 134525    | 132867    | 0.9877  | 56121372 | 417.18 | 0   | 1746 |
| 16   | Pt    | 134525    | 132857    | 0.9876  | 55970896 | 416.06 | 0   | 1746 |
| 32   | Pt    | 134525    | 115373    | 0.8576  | 29714063 | 220.88 | 0   | 1743 |
| 64   | Pt    | 134525    | 18374     | 0.1366  | 2141712  | 15.92  | 0   | 450  |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/nip

find . -type d -name "trim" | xargs rm -fr
find . -type f -path "*3_bwa/genome.fa*" | xargs rm
find . -type f -name "*.ba[mi]" | xargs rm
find . -type f -name "*.per-base.bed.gz" | xargs rm

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

