# *Prunus persica* Lovell 倍数因子测试


[TOC levels=1-3]: # ""

- [*Prunus persica* Lovell 倍数因子测试](#prunus-persica-lovell-倍数因子测试)
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
* 测试文件 [SRR502985](https://www.ncbi.nlm.nih.gov/sra/SRX150254)
  * 覆盖度 110

```shell script
BASES=$(
    faops n50 -N 0 -S \
        ~/data/plastid/ena/SRR502985_1.fastq.gz \
        ~/data/plastid/ena/SRR502985_2.fastq.gz |
    cut -f 2
)

GENOME=$(
    faops size ~/data/plastid/genome/lovell/genome.fa |
        cut -f 2 |
        perl -nl -e '
            $sum += $_;
            END {print $sum}
        '
)

echo ${BASES}
# 24965269082

echo ${GENOME}
# 225852601

bc <<< "${BASES} / ${GENOME}"
# 110

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/lovell
cd ~/data/plastid/evaluation/lovell

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR502985_${FOLD}

    mkdir -p ${BASE_NAME}/1_genome
    pushd ${BASE_NAME}/1_genome

    ln -fs ../../../../genome/lovell/genome.fa genome.fa
    popd

    mkdir -p ${BASE_NAME}/2_illumina
    pushd ${BASE_NAME}/2_illumina

    ln -fs ../../../../ena/SRR502985_1.fastq.gz R1.fq.gz
    ln -fs ../../../../ena/SRR502985_2.fastq.gz R2.fq.gz
    popd

done

```

## Trim and cutoff

```shell script
cd ~/data/plastid/evaluation/lovell

# 倍数因子::cutoff
ARRAY=()
DEPTH=110
for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    CUTOFF=$(bc <<< "(${DEPTH} * ${FOLD}) / 1")
    ARRAY+=("${FOLD}::${CUTOFF}")
done
echo "${ARRAY[@]}"
#0::0 0.25::27 0.5::55 1::110 2::220 4::440 8::880 16::1760 32::3520 64::7040

for item in "${ARRAY[@]}" ; do
    FOLD="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"

    BASE_NAME=SRR502985_${FOLD}
    pushd ${BASE_NAME}

    rm *.sh

    if [[ "${FOLD}" == "0" ]]; then
        anchr template \
            --genome 225852601 \
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
            --genome 225852601 \
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
cd ~/data/plastid/evaluation/lovell

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR502985_${FOLD}

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
cd ~/data/plastid/evaluation/lovell

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR502985_${FOLD}

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
cd ~/data/plastid/evaluation/lovell

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR502985_${FOLD}

    echo 1>&2 "==> ${BASE_NAME}"

    pushd ${BASE_NAME}/3_bwa > /dev/null

    echo -e "Fold\tchrom\n${FOLD}\tNc\n${FOLD}\tMt\n${FOLD}\tPt" |
        tsv-join -H --filter-file combine.tsv --key-fields chrom --append-fields 2-8

    popd > /dev/null

done |
    tsv-uniq \
    > SRR502985_folds.tsv

for PART in Nc Mt Pt; do
    cat SRR502985_folds.tsv |
        tsv-filter -H --str-eq chrom:${PART} |
        mlr --itsv --omd cat
    echo
    echo
done

```

| Fold | chrom | chrLength | covLength | covRate | bases       | mean  | min | max   |
|:-----|:------|:----------|:----------|:--------|:------------|:------|:----|:------|
| 0    | Nc    | 225694811 | 222793155 | 0.9871  | 16134784059 | 71.49 | 0   | 13656 |
| 0.25 | Nc    | 225694811 | 222526077 | 0.9860  | 16069824244 | 71.20 | 0   | 13612 |
| 0.5  | Nc    | 225694811 | 217489368 | 0.9636  | 14380700791 | 63.72 | 0   | 13505 |
| 1    | Nc    | 225694811 | 120284718 | 0.5330  | 5382692445  | 23.85 | 0   | 13760 |
| 2    | Nc    | 225694811 | 94097073  | 0.4169  | 4188422371  | 18.56 | 0   | 14003 |
| 4    | Nc    | 225694811 | 78282996  | 0.3469  | 3526307401  | 15.62 | 0   | 14036 |
| 8    | Nc    | 225694811 | 65389460  | 0.2897  | 2950822870  | 13.07 | 0   | 13879 |
| 16   | Nc    | 225694811 | 52490204  | 0.2326  | 2371883298  | 10.51 | 0   | 13386 |
| 32   | Nc    | 225694811 | 39754554  | 0.1761  | 1893829293  | 8.39  | 0   | 12877 |
| 64   | Nc    | 225694811 | 28562638  | 0.1266  | 1569281213  | 6.95  | 0   | 12182 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean    | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:--------|:----|:-----|
| 0    | Pt    | 157790    | 157790    | 1.0000  | 663531327 | 4205.15 | 13  | 8249 |
| 0.25 | Pt    | 157790    | 157790    | 1.0000  | 662940802 | 4201.41 | 6   | 8244 |
| 0.5  | Pt    | 157790    | 157790    | 1.0000  | 662841137 | 4200.78 | 4   | 8243 |
| 1    | Pt    | 157790    | 157790    | 1.0000  | 662659333 | 4199.63 | 3   | 8242 |
| 2    | Pt    | 157790    | 157790    | 1.0000  | 662550325 | 4198.94 | 2   | 8242 |
| 4    | Pt    | 157790    | 157659    | 0.9992  | 662245038 | 4197.00 | 0   | 8241 |
| 8    | Pt    | 157790    | 157480    | 0.9980  | 660744133 | 4187.49 | 0   | 8242 |
| 16   | Pt    | 157790    | 157276    | 0.9967  | 654323214 | 4146.80 | 0   | 8243 |
| 32   | Pt    | 157790    | 155986    | 0.9886  | 571454194 | 3621.61 | 0   | 8241 |
| 64   | Pt    | 157790    | 145257    | 0.9206  | 273302539 | 1732.07 | 0   | 8189 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/lovell

find . -type d -name "trim" | xargs rm -fr
find . -type f -path "*3_bwa/genome.fa*" | xargs rm
find . -type f -name "*.ba[mi]" | xargs rm
find . -type f -name "*.per-base.bed.gz" | xargs rm

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

