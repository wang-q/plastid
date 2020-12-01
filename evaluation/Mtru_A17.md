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

## Mapping

```shell script
cd ~/data/plastid/evaluation/a17

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${FOLD}
    
    mkdir -p ${BASE_NAME}/mapping
    pushd ${BASE_NAME}/mapping
    
    if [ -f R.sort.bai ]; then
        echo >&2 '    R.sort.bai already presents'
        popd;
        continue;
    fi

    bsub -q mpi -n 24 -J "${BASE_NAME}-mapping" '
        # export JAVA_OPTS="-Xmx20G"

        # Pipe all reads together as we do not need mate info
        gzip -dcf \
            ../2_illumina/trim/Q25L60/R1.fq.gz \
            ../2_illumina/trim/Q25L60/R2.fq.gz \
            ../2_illumina/trim/Q25L60/Rs.fq.gz |
            faops filter -l 0 stdin stdout | # ignore QUAL
            bowtie2 -p 20 --very-fast -t \
                -x ../../../../genome/a17/genome.fa \
                -f -U /dev/stdin |
            picard CleanSam \
                --INPUT /dev/stdin \
                --OUTPUT /dev/stdout \
                --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 0 |
            picard SortSam \
                --INPUT /dev/stdin \
                --OUTPUT R.sort.bam \
                --SORT_ORDER coordinate \
                --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 1
    
        picard BuildBamIndex \
            --INPUT R.sort.bam \
            --VALIDATION_STRINGENCY LENIENT
    '
    
    popd

done

```

* Stats of mapping

```shell script
cd ~/data/plastid/evaluation/a17

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${FOLD}
    
    echo "==> ${BASE_NAME}"

    cat ${BASE_NAME}/mapping/output.* |
        grep -Pzo '(?s)Multiseed full.+overall alignment rate'

done

```

```text
==> SRR1542423_0
Multiseed full-index search: 00:08:30
39795514 reads; of these:
  39795514 (100.00%) were unpaired; of these:
    4093232 (10.29%) aligned 0 times
    16047972 (40.33%) aligned exactly 1 time
    19654310 (49.39%) aligned >1 times
89.71% overall alignment rate
==> SRR1542423_0.25
Multiseed full-index search: 00:08:08
37826383 reads; of these:
  37826383 (100.00%) were unpaired; of these:
    3673779 (9.71%) aligned 0 times
    14752532 (39.00%) aligned exactly 1 time
    19400072 (51.29%) aligned >1 times
90.29% overall alignment rate
==> SRR1542423_0.5
Multiseed full-index search: 00:06:50
29669403 reads; of these:
  29669403 (100.00%) were unpaired; of these:
    3276378 (11.04%) aligned 0 times
    7976711 (26.89%) aligned exactly 1 time
    18416314 (62.07%) aligned >1 times
88.96% overall alignment rate
==> SRR1542423_1
Multiseed full-index search: 00:05:27
20991325 reads; of these:
  20991325 (100.00%) were unpaired; of these:
    2825125 (13.46%) aligned 0 times
    1782829 (8.49%) aligned exactly 1 time
    16383371 (78.05%) aligned >1 times
86.54% overall alignment rate
==> SRR1542423_2
Multiseed full-index search: 00:04:53
18816136 reads; of these:
  18816136 (100.00%) were unpaired; of these:
    2629155 (13.97%) aligned 0 times
    1161418 (6.17%) aligned exactly 1 time
    15025563 (79.85%) aligned >1 times
86.03% overall alignment rate
==> SRR1542423_4
Multiseed full-index search: 00:04:36
17602656 reads; of these:
  17602656 (100.00%) were unpaired; of these:
    2508267 (14.25%) aligned 0 times
    1020641 (5.80%) aligned exactly 1 time
    14073748 (79.95%) aligned >1 times
85.75% overall alignment rate
==> SRR1542423_8
Multiseed full-index search: 00:04:21
16448954 reads; of these:
  16448954 (100.00%) were unpaired; of these:
    2422824 (14.73%) aligned 0 times
    848697 (5.16%) aligned exactly 1 time
    13177433 (80.11%) aligned >1 times
85.27% overall alignment rate
==> SRR1542423_16
Multiseed full-index search: 00:04:02
15020170 reads; of these:
  15020170 (100.00%) were unpaired; of these:
    2361314 (15.72%) aligned 0 times
    525547 (3.50%) aligned exactly 1 time
    12133309 (80.78%) aligned >1 times
84.28% overall alignment rate
==> SRR1542423_32
Multiseed full-index search: 00:03:41
13814223 reads; of these:
  13814223 (100.00%) were unpaired; of these:
    2294616 (16.61%) aligned 0 times
    408210 (2.95%) aligned exactly 1 time
    11111397 (80.43%) aligned >1 times
83.39% overall alignment rate
==> SRR1542423_64
Multiseed full-index search: 00:03:20
12453870 reads; of these:
  12453870 (100.00%) were unpaired; of these:
    2214213 (17.78%) aligned 0 times
    372792 (2.99%) aligned exactly 1 time
    9866865 (79.23%) aligned >1 times
82.22% overall alignment rate

```

## Depth

* Depth via `mosdepth`

```shell script
cd ~/data/plastid/evaluation/a17

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${FOLD}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth

    if [ -f R.mosdepth.summary.txt ]; then
        echo >&2 '    R.mosdepth.summary.txt already presents'
        popd;
        continue;
    fi

    mosdepth R ../mapping/R.sort.bam
    
    popd

done

```

* Covered regions via `spanr`

```shell script
cd ~/data/plastid/evaluation/a17

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${FOLD}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth

    gzip -dcf R.per-base.bed.gz |
        perl -nla -F"\t" -e '
            $F[3] == 0 and next;
            $start = $F[1] + 1;
            $end = $F[2];
            if ($start == $F[2]) {
                print qq($F[0]:$start);
            }
            else {
                print qq($F[0]:$start-$end);
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

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${FOLD}
    
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

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1542423_${FOLD}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth > /dev/null

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

| Fold | chrom | chrLength | covLength | covRate | bases      | mean  | min | max    |
|:-----|:------|:----------|:----------|:--------|:-----------|:------|:----|:-------|
| 0    | Nc    | 384466993 | 349437629 | 0.9089  | 4327426276 | 11.26 | 0   | 858687 |
| 0.25 | Nc    | 384466993 | 319962957 | 0.8322  | 4143247285 | 10.78 | 0   | 858229 |
| 0.5  | Nc    | 384466993 | 197843727 | 0.5146  | 3141188884 | 8.17  | 0   | 858023 |
| 1    | Nc    | 384466993 | 76653164  | 0.1994  | 2086980331 | 5.43  | 0   | 857746 |
| 2    | Nc    | 384466993 | 51894129  | 0.1350  | 1836298182 | 4.78  | 0   | 857145 |
| 4    | Nc    | 384466993 | 40208484  | 0.1046  | 1698972847 | 4.42  | 0   | 856196 |
| 8    | Nc    | 384466993 | 31558696  | 0.0821  | 1578702338 | 4.11  | 0   | 854628 |
| 16   | Nc    | 384466993 | 23828317  | 0.0620  | 1449524262 | 3.77  | 0   | 852300 |
| 32   | Nc    | 384466993 | 16793009  | 0.0437  | 1322816896 | 3.44  | 0   | 848676 |
| 64   | Nc    | 384466993 | 11639440  | 0.0303  | 1184656229 | 3.08  | 0   | 843430 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Mt    | 271618    | 271618    | 1.0000  | 68924581 | 253.76 | 1   | 1294 |
| 0.25 | Mt    | 271618    | 271618    | 1.0000  | 67529569 | 248.62 | 1   | 1285 |
| 0.5  | Mt    | 271618    | 271618    | 1.0000  | 67342935 | 247.93 | 1   | 1285 |
| 1    | Mt    | 271618    | 271618    | 1.0000  | 67186006 | 247.35 | 1   | 1283 |
| 2    | Mt    | 271618    | 271618    | 1.0000  | 67111050 | 247.08 | 1   | 1279 |
| 4    | Mt    | 271618    | 271495    | 0.9995  | 66427833 | 244.56 | 0   | 1277 |
| 8    | Mt    | 271618    | 254074    | 0.9354  | 53310986 | 196.27 | 0   | 1268 |
| 16   | Mt    | 271618    | 94380     | 0.3475  | 13137540 | 48.37  | 0   | 1266 |
| 32   | Mt    | 271618    | 20746     | 0.0764  | 967564   | 3.56   | 0   | 1057 |
| 64   | Mt    | 271618    | 13964     | 0.0514  | 51110    | 0.19   | 0   | 244  |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Pt    | 124033    | 124033    | 1.0000  | 89460546 | 721.26 | 1   | 4626 |
| 0.25 | Pt    | 124033    | 124033    | 1.0000  | 89004547 | 717.59 | 1   | 4617 |
| 0.5  | Pt    | 124033    | 124030    | 1.0000  | 88424794 | 712.91 | 0   | 4592 |
| 1    | Pt    | 124033    | 124024    | 0.9999  | 88170772 | 710.87 | 0   | 4577 |
| 2    | Pt    | 124033    | 124024    | 0.9999  | 88091785 | 710.23 | 0   | 4567 |
| 4    | Pt    | 124033    | 124024    | 0.9999  | 88034561 | 709.77 | 0   | 4567 |
| 8    | Pt    | 124033    | 124009    | 0.9998  | 87919456 | 708.84 | 0   | 4565 |
| 16   | Pt    | 124033    | 123783    | 0.9980  | 87330660 | 704.09 | 0   | 4565 |
| 32   | Pt    | 124033    | 117673    | 0.9487  | 83399317 | 672.40 | 0   | 4559 |
| 64   | Pt    | 124033    | 82755     | 0.6672  | 61217216 | 493.56 | 0   | 4554 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/a17

find . -type d -name "trim" | xargs rm -fr
find . -type d -name "mapping" | xargs rm -fr
# find . -type d -name "depth" | xargs rm -fr

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

