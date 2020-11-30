# *Medicago truncatula* A17 倍数因子测试


[TOC levels=1-3]: # ""

- [*Medicago truncatula* A17 倍数因子测试](#medicago-truncatula-a17-倍数因子测试)
  - [基本信息](#基本信息)
  - [Symlink](#symlink)
  - [Trim and cutoff](#trim-and-cutoff)
  - [`kat hist` and `kat gcp`](#kat-hist-and-kat-gcp)
  - [Mapping](#mapping)
  - [Depth](#depth)
  - [Merge all results](#merge-all-results)
  - [Remove intermediate files](#remove-intermediate-files)

## 基本信息

`cutoff = 倍数因子 * 覆盖深度`

+ 因子值 0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64
+ 测试文件 [SRR1542423](https://www.ncbi.nlm.nih.gov/sra/SRX673852)
+ 覆盖度 8,958,357,672 / 384,862,644 = 23

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

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/a17
cd ~/data/plastid/evaluation/a17

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${NAME}
    
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
ARRAY=(
    '0::0'
    '0.25::5'
    '0.5::11'
    '1::23'
    '2::46'
    '4::92'
    '8::184'
    '16::368'
    '32::736'
    '64::1472'
)

for item in "${ARRAY[@]}" ; do
    NAME="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"
    
    BASE_NAME=SRR1542423_${NAME}
    pushd ${BASE_NAME}
    
    rm *.sh
        
    if [[ "${NAME}" == "0" ]]; then
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
            --filter "adapter artifact"

        bsub -q mpi -n 24 -J "${BASE_NAME}" "
            bash 2_fastqc.sh
            bash 2_insert_size.sh
            bash 2_kat.sh
            bash 2_trim.sh
            bash 9_stat_reads.sh
        "        
    else
        anchr template \
            --genome 384862644 \
            --parallel 24 \
            --xmx 80g \
            \
            --fastqc \
            --insertsize \
            --kat \
            \
            --trim "--dedupe --cutoff ${CUTOFF} --cutk 31" \
            --qual "25" \
            --len "60" \
            --filter "adapter artifact"

        bsub -q mpi -n 24 -J "${BASE_NAME}" "
            bash 2_trim.sh
            bash 9_stat_reads.sh
        "
    fi

    popd

done

```

## `kat hist` and `kat gcp`

```shell script
cd ~/data/plastid/evaluation/a17

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${NAME}
    
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

## Mapping

```shell script
cd ~/data/plastid/evaluation/a17

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${NAME}
    
    mkdir -p ${BASE_NAME}/mapping
    pushd ${BASE_NAME}/mapping
    
    bsub -q mpi -n 24 -J "${BASE_NAME}-mapping" '

        # Pipe all reads together as we do not need mate info
        gzip -dcf \
            ../2_illumina/trim/Q25L60/R1.fq.gz \
            ../2_illumina/trim/Q25L60/R2.fq.gz \
            ../2_illumina/trim/Q25L60/Rs.fq.gz |
            bwa mem -M -t 20 \
                ../../../../genome/a17/genome.fa \
                /dev/stdin |
            pigz -p 4 \
            > R.sam.gz
        
        picard CleanSam \
            --INPUT R.sam.gz \
            --OUTPUT R.clean.bam \
            --VALIDATION_STRINGENCY LENIENT
    
        picard SortSam \
            --INPUT R.clean.bam \
            --OUTPUT R.sort.bam \
            --SORT_ORDER coordinate \
            --VALIDATION_STRINGENCY LENIENT
    
        picard BuildBamIndex \
            --INPUT R.sort.bam \
            --VALIDATION_STRINGENCY LENIENT
            
        find . -name "R.sam.gz" | xargs rm
        find . -name "R.clean.bam" | xargs rm
    '
    
    popd

done

```


## Depth

* Depth via `mosdepth`

```shell script
cd ~/data/plastid/evaluation/a17

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth

    mosdepth R ../mapping/R.sort.bam
    
    popd

done

```

* Covered regions via `spanr`

```shell script
cd ~/data/plastid/evaluation/a17

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth

    gzip -dcf R.per-base.bed.gz |
        perl -nla -F"\t" -e '
            $F[3] == 0 and next;
            $start = $F[1] + 1;
            if ($start == $F[2]) {
                print qq($F[0]:$start);
            }
            else {
                print qq($F[0]:$start-$F[2]);
            }
        ' |
        spanr cover stdin -o covered.yml
        
    spanr stat ../../../../genome/a17/chr.sizes covered.yml -o stdout |
        grep -v "^all" |
        sed 's/^chr/chrom/' |
        sed 's/,size/,covLength/' |
        sed 's/,coverage/,covRate/' |
        sed 's/,/\t/g' \
        > coverage.tsv
    
    popd

done

```

* Combine results

```shell script
cd ~/data/plastid/evaluation/a17

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth

    cat coverage.tsv |
        tsv-join -H --filter-file R.mosdepth.summary.txt \
            --key-fields chrom --append-fields 3-6 \
        > join.tsv

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

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth > /dev/null

    echo -e "Fold\tchrom\n${NAME}\tNc\n${NAME}\tMt\n${NAME}\tPt" |
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

| Fold | chrom | chrLength | covLength | covRate | bases      | mean  | min | max     |
|:-----|:------|:----------|:----------|:--------|:-----------|:------|:----|:--------|
| 0    | Nc    | 384466993 | 349624192 | 0.9094  | 4510879694 | 11.73 | 0   | 1407943 |
| 0.25 | Nc    | 384466993 | 320317797 | 0.8331  | 4322435031 | 11.24 | 0   | 1406474 |
| 0.5  | Nc    | 384466993 | 198580741 | 0.5165  | 3306410912 | 8.60  | 0   | 1406500 |
| 1    | Nc    | 384466993 | 77234923  | 0.2009  | 2228530088 | 5.80  | 0   | 1405043 |
| 2    | Nc    | 384466993 | 52136316  | 0.1356  | 1965633146 | 5.11  | 0   | 1404882 |
| 4    | Nc    | 384466993 | 40294864  | 0.1048  | 1820778536 | 4.74  | 0   | 1402544 |
| 8    | Nc    | 384466993 | 31554654  | 0.0821  | 1694721993 | 4.41  | 0   | 1397985 |
| 16   | Nc    | 384466993 | 23746431  | 0.0618  | 1561019339 | 4.06  | 0   | 1393249 |
| 32   | Nc    | 384466993 | 16716014  | 0.0435  | 1429771844 | 3.72  | 0   | 1385978 |
| 64   | Nc    | 384466993 | 11597103  | 0.0302  | 1286133231 | 3.35  | 0   | 1376485 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Mt    | 271618    | 271618    | 1.0000  | 69136307 | 254.54 | 1   | 1294 |
| 0.25 | Mt    | 271618    | 271618    | 1.0000  | 67672293 | 249.15 | 1   | 1286 |
| 0.5  | Mt    | 271618    | 271618    | 1.0000  | 67441873 | 248.30 | 1   | 1285 |
| 1    | Mt    | 271618    | 271618    | 1.0000  | 67266172 | 247.65 | 1   | 1284 |
| 2    | Mt    | 271618    | 271618    | 1.0000  | 67167647 | 247.29 | 1   | 1279 |
| 4    | Mt    | 271618    | 271495    | 0.9995  | 66472305 | 244.73 | 0   | 1277 |
| 8    | Mt    | 271618    | 254216    | 0.9359  | 53346996 | 196.40 | 0   | 1268 |
| 16   | Mt    | 271618    | 97444     | 0.3588  | 13137668 | 48.37  | 0   | 1267 |
| 32   | Mt    | 271618    | 25529     | 0.0940  | 973637   | 3.58   | 0   | 1059 |
| 64   | Mt    | 271618    | 17850     | 0.0657  | 55453    | 0.20   | 0   | 239  |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Pt    | 124033    | 124033    | 1.0000  | 89734517 | 723.47 | 1   | 4610 |
| 0.25 | Pt    | 124033    | 124033    | 1.0000  | 89252242 | 719.58 | 1   | 4588 |
| 0.5  | Pt    | 124033    | 124032    | 1.0000  | 88523820 | 713.71 | 0   | 4558 |
| 1    | Pt    | 124033    | 124032    | 1.0000  | 88370854 | 712.48 | 0   | 4549 |
| 2    | Pt    | 124033    | 124032    | 1.0000  | 88192394 | 711.04 | 0   | 4533 |
| 4    | Pt    | 124033    | 124032    | 1.0000  | 88210157 | 711.18 | 0   | 4541 |
| 8    | Pt    | 124033    | 123939    | 0.9992  | 88073301 | 710.08 | 0   | 4530 |
| 16   | Pt    | 124033    | 123791    | 0.9980  | 87458644 | 705.12 | 0   | 4536 |
| 32   | Pt    | 124033    | 117773    | 0.9495  | 83404438 | 672.44 | 0   | 4545 |
| 64   | Pt    | 124033    | 83533     | 0.6735  | 61285085 | 494.10 | 0   | 4528 |

## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/a17

find . -type d -name "trim" | xargs rm -fr
find . -type d -name "mapping" | xargs rm -fr

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

