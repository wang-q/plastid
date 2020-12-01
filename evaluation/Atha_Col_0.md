# *Arabidopsis thaliana* Col-0 倍数因子测试


[TOC levels=1-3]: # ""

- [*Arabidopsis thaliana* Col-0 倍数因子测试](#arabidopsis-thaliana-col-0-倍数因子测试)
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
+ 测试文件 SRR616966
+ 覆盖度 4,970,359,200 / 119,667,750 = 41

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

```

* 末被包含在 `anchr` 里的工具

```shell script
# 将 reads 比对到参考序列上
brew install bowtie2

# 计算深度
brew install mosdepth

# 转换得到 bigwig, 使用 IGV 查看比对图
pip install deeptools

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/col_0
cd ~/data/plastid/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${NAME}
    
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
ARRAY=(
    '0::0'
    '0.25::10'
    '0.5::20'
    '1::40'
    '2::80'
    '4::160'
    '8::320'
    '16::640'
    '32::1280'
    '64::2560'
)

for item in "${ARRAY[@]}" ; do
    NAME="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"
    
    BASE_NAME=SRR616966_${NAME}
    pushd ${BASE_NAME}
    
    rm *.sh
        
    if [[ "${NAME}" == "0" ]]; then
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
            --genome 119667750 \
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
cd ~/data/plastid/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${NAME}
    
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
cd ~/data/plastid/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${NAME}
    
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
                -x ../../../../genome/col_0/genome.fa \
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
cd ~/data/plastid/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${NAME}
    
    echo "==> ${BASE_NAME}"

    cat ${BASE_NAME}/mapping/output.* |
        grep -Pzo '(?s)Multiseed full.+overall alignment rate'

done

```

```text
==> SRR616966_0
Multiseed full-index search: 00:05:38
42753121 reads; of these:
  42753121 (100.00%) were unpaired; of these:
    8115126 (18.98%) aligned 0 times
    24636878 (57.63%) aligned exactly 1 time
    10001117 (23.39%) aligned >1 times
81.02% overall alignment rate
==> SRR616966_0.25
Multiseed full-index search: 00:04:59
34873004 reads; of these:
  34873004 (100.00%) were unpaired; of these:
    431365 (1.24%) aligned 0 times
    24462648 (70.15%) aligned exactly 1 time
    9978991 (28.62%) aligned >1 times
98.76% overall alignment rate
==> SRR616966_0.5
Multiseed full-index search: 00:02:56
20354104 reads; of these:
  20354104 (100.00%) were unpaired; of these:
    311278 (1.53%) aligned 0 times
    10796915 (53.05%) aligned exactly 1 time
    9245911 (45.43%) aligned >1 times
98.47% overall alignment rate
==> SRR616966_1
Multiseed full-index search: 00:01:48
11994918 reads; of these:
  11994918 (100.00%) were unpaired; of these:
    261222 (2.18%) aligned 0 times
    3631401 (30.27%) aligned exactly 1 time
    8102295 (67.55%) aligned >1 times
97.82% overall alignment rate
==> SRR616966_2
Multiseed full-index search: 00:01:50
11321357 reads; of these:
  11321357 (100.00%) were unpaired; of these:
    222341 (1.96%) aligned 0 times
    3546897 (31.33%) aligned exactly 1 time
    7552119 (66.71%) aligned >1 times
98.04% overall alignment rate
==> SRR616966_4
Multiseed full-index search: 00:01:37
10840616 reads; of these:
  10840616 (100.00%) were unpaired; of these:
    195607 (1.80%) aligned 0 times
    3493232 (32.22%) aligned exactly 1 time
    7151777 (65.97%) aligned >1 times
98.20% overall alignment rate
==> SRR616966_8
Multiseed full-index search: 00:01:32
9696276 reads; of these:
  9696276 (100.00%) were unpaired; of these:
    178257 (1.84%) aligned 0 times
    3248561 (33.50%) aligned exactly 1 time
    6269458 (64.66%) aligned >1 times
98.16% overall alignment rate
==> SRR616966_16
Multiseed full-index search: 00:01:31
9362828 reads; of these:
  9362828 (100.00%) were unpaired; of these:
    166340 (1.78%) aligned 0 times
    3225383 (34.45%) aligned exactly 1 time
    5971105 (63.77%) aligned >1 times
98.22% overall alignment rate
==> SRR616966_32
Multiseed full-index search: 00:01:27
9129676 reads; of these:
  9129676 (100.00%) were unpaired; of these:
    156036 (1.71%) aligned 0 times
    3215393 (35.22%) aligned exactly 1 time
    5758247 (63.07%) aligned >1 times
98.29% overall alignment rate
==> SRR616966_64
Multiseed full-index search: 00:01:16
7634180 reads; of these:
  7634180 (100.00%) were unpaired; of these:
    141731 (1.86%) aligned 0 times
    2003031 (26.24%) aligned exactly 1 time
    5489418 (71.91%) aligned >1 times
98.14% overall alignment rate
```

## Depth

* Depth via `mosdepth`

```shell script
cd ~/data/plastid/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${NAME}
    
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
cd ~/data/plastid/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${NAME}
    
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
        
    spanr stat ../../../../genome/col_0/chr.sizes covered.yml -o stdout |
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
cd ~/data/plastid/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${NAME}
    
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
        (cat join.tsv | sed '2,6d' && cat) \
        > combine.tsv

    popd

done

```

## Merge all results

```shell script
cd ~/data/plastid/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth > /dev/null

    echo -e "Fold\tchrom\n${NAME}\tNc\n${NAME}\tMt\n${NAME}\tPt" |
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

| Fold | chrom | chrLength | covLength | covRate | bases      | mean  | min | max    |
|:-----|:------|:----------|:----------|:--------|:-----------|:------|:----|:-------|
| 0    | Nc    | 119146348 | 118888104 | 0.9978  | 2849796891 | 23.92 | 0   | 106651 |
| 0.25 | Nc    | 119146348 | 118743926 | 0.9966  | 2832037578 | 23.77 | 0   | 106647 |
| 0.5  | Nc    | 119146348 | 78534610  | 0.6591  | 1422224302 | 11.94 | 0   | 106645 |
| 1    | Nc    | 119146348 | 14286243  | 0.1199  | 604519043  | 5.07  | 0   | 106633 |
| 2    | Nc    | 119146348 | 9511521   | 0.0798  | 542266145  | 4.55  | 0   | 106620 |
| 4    | Nc    | 119146348 | 6793685   | 0.0570  | 501052722  | 4.21  | 0   | 106587 |
| 8    | Nc    | 119146348 | 4094381   | 0.0344  | 438937996  | 3.68  | 0   | 106477 |
| 16   | Nc    | 119146348 | 2853174   | 0.0239  | 414090706  | 3.48  | 0   | 106210 |
| 32   | Nc    | 119146348 | 1854590   | 0.0156  | 392469946  | 3.29  | 0   | 105590 |
| 64   | Nc    | 119146348 | 1282693   | 0.0108  | 376843211  | 3.16  | 0   | 104014 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Mt    | 366924    | 363784    | 0.9914  | 58791874 | 160.23 | 0   | 1924 |
| 0.25 | Mt    | 366924    | 363781    | 0.9914  | 58694616 | 159.96 | 0   | 1921 |
| 0.5  | Mt    | 366924    | 363781    | 0.9914  | 58692222 | 159.96 | 0   | 1921 |
| 1    | Mt    | 366924    | 363781    | 0.9914  | 58686680 | 159.94 | 0   | 1920 |
| 2    | Mt    | 366924    | 363490    | 0.9906  | 58667344 | 159.89 | 0   | 1920 |
| 4    | Mt    | 366924    | 358221    | 0.9763  | 55571688 | 151.45 | 0   | 1920 |
| 8    | Mt    | 366924    | 97136     | 0.2647  | 7959947  | 21.69  | 0   | 1919 |
| 16   | Mt    | 366924    | 30005     | 0.0818  | 1377573  | 3.75   | 0   | 1919 |
| 32   | Mt    | 366924    | 27615     | 0.0753  | 1298721  | 3.54   | 0   | 1919 |
| 64   | Mt    | 366924    | 21074     | 0.0574  | 1173507  | 3.20   | 0   | 1919 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean    | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:--------|:----|:-----|
| 0    | Pt    | 154478    | 154478    | 1.0000  | 482866004 | 3125.79 | 4   | 4507 |
| 0.25 | Pt    | 154478    | 154478    | 1.0000  | 481938779 | 3119.79 | 4   | 4504 |
| 0.5  | Pt    | 154478    | 154478    | 1.0000  | 481866034 | 3119.32 | 4   | 4504 |
| 1    | Pt    | 154478    | 154478    | 1.0000  | 481850340 | 3119.22 | 4   | 4504 |
| 2    | Pt    | 154478    | 154478    | 1.0000  | 481833933 | 3119.11 | 4   | 4504 |
| 4    | Pt    | 154478    | 154478    | 1.0000  | 481823604 | 3119.04 | 4   | 4504 |
| 8    | Pt    | 154478    | 154478    | 1.0000  | 481815378 | 3118.99 | 4   | 4504 |
| 16   | Pt    | 154478    | 154478    | 1.0000  | 481810317 | 3118.96 | 4   | 4504 |
| 32   | Pt    | 154478    | 154478    | 1.0000  | 481645944 | 3117.89 | 2   | 4504 |
| 64   | Pt    | 154478    | 150475    | 0.9741  | 352739147 | 2283.43 | 0   | 4498 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/col_0

find . -type d -name "trim" | xargs rm -fr
find . -type d -name "mapping" | xargs rm -fr
# find . -type d -name "depth" | xargs rm -fr

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

