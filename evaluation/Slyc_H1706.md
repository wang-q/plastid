# *Solanum lycopersicum* Heinz 1706 倍数因子测试


[TOC levels=1-3]: # ""

- [*Solanum lycopersicum* Heinz 1706 倍数因子测试](#solanum-lycopersicum-heinz-1706-倍数因子测试)
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
+ 测试文件 [SRR1572628](https://www.ncbi.nlm.nih.gov/sra/SRX698770)
+ 覆盖度 5

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
            --genome 807826382 \
            --parallel 24 \
            --xmx 80g \
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

# find . -type f -wholename "*trim/R*.fq.gz"

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

## Mapping

```shell script
cd ~/data/plastid/evaluation/h1706

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1572628_${FOLD}
    
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
                -x ../../../../genome/h1706/genome.fa \
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
cd ~/data/plastid/evaluation/h1706

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1572628_${FOLD}
    
    echo "==> ${BASE_NAME}"

    cat ${BASE_NAME}/mapping/output.* |
        grep -Pzo '(?s)Multiseed full.+overall alignment rate'

done

```

```text
==> SRR1572628_0
Multiseed full-index search: 00:06:46
39077164 reads; of these:
  39077164 (100.00%) were unpaired; of these:
    1525053 (3.90%) aligned 0 times
    21316709 (54.55%) aligned exactly 1 time
    16235402 (41.55%) aligned >1 times
96.10% overall alignment rate
==> SRR1572628_0.25
Multiseed full-index search: 00:06:38
39076276 reads; of these:
  39076276 (100.00%) were unpaired; of these:
    1524954 (3.90%) aligned 0 times
    21316562 (54.55%) aligned exactly 1 time
    16234760 (41.55%) aligned >1 times
96.10% overall alignment rate
==> SRR1572628_0.5
Multiseed full-index search: 00:06:40
38848839 reads; of these:
  38848839 (100.00%) were unpaired; of these:
    1481947 (3.81%) aligned 0 times
    21152295 (54.45%) aligned exactly 1 time
    16214597 (41.74%) aligned >1 times
96.19% overall alignment rate
==> SRR1572628_1
Multiseed full-index search: 00:05:31
30765650 reads; of these:
  30765650 (100.00%) were unpaired; of these:
    1376252 (4.47%) aligned 0 times
    14282164 (46.42%) aligned exactly 1 time
    15107234 (49.10%) aligned >1 times
95.53% overall alignment rate
==> SRR1572628_2
Multiseed full-index search: 00:03:12
15790006 reads; of these:
  15790006 (100.00%) were unpaired; of these:
    1161835 (7.36%) aligned 0 times
    2761722 (17.49%) aligned exactly 1 time
    11866449 (75.15%) aligned >1 times
92.64% overall alignment rate
==> SRR1572628_4
Multiseed full-index search: 00:02:33
12434140 reads; of these:
  12434140 (100.00%) were unpaired; of these:
    999765 (8.04%) aligned 0 times
    1301411 (10.47%) aligned exactly 1 time
    10132964 (81.49%) aligned >1 times
91.96% overall alignment rate
==> SRR1572628_8
Multiseed full-index search: 00:02:20
11022472 reads; of these:
  11022472 (100.00%) were unpaired; of these:
    826018 (7.49%) aligned 0 times
    1023931 (9.29%) aligned exactly 1 time
    9172523 (83.22%) aligned >1 times
92.51% overall alignment rate
==> SRR1572628_16
Multiseed full-index search: 00:02:02
9863824 reads; of these:
  9863824 (100.00%) were unpaired; of these:
    670451 (6.80%) aligned 0 times
    852379 (8.64%) aligned exactly 1 time
    8340994 (84.56%) aligned >1 times
93.20% overall alignment rate
==> SRR1572628_32
Multiseed full-index search: 00:01:45
8495816 reads; of these:
  8495816 (100.00%) were unpaired; of these:
    507027 (5.97%) aligned 0 times
    672550 (7.92%) aligned exactly 1 time
    7316239 (86.12%) aligned >1 times
94.03% overall alignment rate
==> SRR1572628_64
Multiseed full-index search: 00:01:34
7435598 reads; of these:
  7435598 (100.00%) were unpaired; of these:
    389120 (5.23%) aligned 0 times
    570460 (7.67%) aligned exactly 1 time
    6476018 (87.09%) aligned >1 times
94.77% overall alignment rate

```

## Depth

* Depth via `mosdepth`

```shell script
cd ~/data/plastid/evaluation/h1706

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1572628_${FOLD}
    
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
cd ~/data/plastid/evaluation/h1706

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1572628_${FOLD}
    
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
        
    spanr stat ../../../../genome/h1706/chr.sizes covered.yml -o stdout |
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
cd ~/data/plastid/evaluation/h1706

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR1572628_${FOLD}
    
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
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth > /dev/null

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
| 0    | Nc    | 807224664 | 690315363 | 0.8552  | 3481804917 | 4.31 | 0   | 75317 |
| 0.25 | Nc    | 807224664 | 690314439 | 0.8552  | 3481735850 | 4.31 | 0   | 75305 |
| 0.5  | Nc    | 807224664 | 681202580 | 0.8439  | 3464064210 | 4.29 | 0   | 75305 |
| 1    | Nc    | 807224664 | 514441861 | 0.6373  | 2681096142 | 3.32 | 0   | 75305 |
| 2    | Nc    | 807224664 | 205815919 | 0.2550  | 1230703754 | 1.52 | 0   | 75298 |
| 4    | Nc    | 807224664 | 138500990 | 0.1716  | 917562580  | 1.14 | 0   | 75286 |
| 8    | Nc    | 807224664 | 112994900 | 0.1400  | 796799014  | 0.99 | 0   | 75265 |
| 16   | Nc    | 807224664 | 94282042  | 0.1168  | 703686906  | 0.87 | 0   | 75232 |
| 32   | Nc    | 807224664 | 77461765  | 0.0960  | 613603053  | 0.76 | 0   | 75173 |
| 64   | Nc    | 807224664 | 63406970  | 0.0785  | 534451456  | 0.66 | 0   | 74928 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Mt    | 446257    | 446257    | 1.0000  | 44833383 | 100.47 | 1   | 1076 |
| 0.25 | Mt    | 446257    | 446257    | 1.0000  | 44833295 | 100.47 | 1   | 1076 |
| 0.5  | Mt    | 446257    | 446257    | 1.0000  | 44818599 | 100.43 | 1   | 1076 |
| 1    | Mt    | 446257    | 446257    | 1.0000  | 44806115 | 100.40 | 1   | 1076 |
| 2    | Mt    | 446257    | 446257    | 1.0000  | 44794921 | 100.38 | 1   | 1076 |
| 4    | Mt    | 446257    | 446257    | 1.0000  | 44793552 | 100.38 | 1   | 1076 |
| 8    | Mt    | 446257    | 446257    | 1.0000  | 44777908 | 100.34 | 1   | 1076 |
| 16   | Mt    | 446257    | 437039    | 0.9793  | 40369362 | 90.46  | 0   | 1076 |
| 32   | Mt    | 446257    | 210209    | 0.4710  | 13698654 | 30.70  | 0   | 1076 |
| 64   | Mt    | 446257    | 31785     | 0.0712  | 2024847  | 4.54   | 0   | 1076 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:-------|:----|:-----|
| 0    | Pt    | 155461    | 155461    | 1.0000  | 139601678 | 897.99 | 6   | 2546 |
| 0.25 | Pt    | 155461    | 155461    | 1.0000  | 139597788 | 897.96 | 6   | 2546 |
| 0.5  | Pt    | 155461    | 155461    | 1.0000  | 139582061 | 897.86 | 6   | 2546 |
| 1    | Pt    | 155461    | 155461    | 1.0000  | 139516109 | 897.43 | 6   | 2544 |
| 2    | Pt    | 155461    | 155461    | 1.0000  | 139460102 | 897.07 | 6   | 2543 |
| 4    | Pt    | 155461    | 155461    | 1.0000  | 139446273 | 896.99 | 6   | 2543 |
| 8    | Pt    | 155461    | 155461    | 1.0000  | 139443737 | 896.97 | 6   | 2543 |
| 16   | Pt    | 155461    | 155461    | 1.0000  | 139442218 | 896.96 | 6   | 2543 |
| 32   | Pt    | 155461    | 155461    | 1.0000  | 139438208 | 896.93 | 6   | 2543 |
| 64   | Pt    | 155461    | 155461    | 1.0000  | 138990512 | 894.05 | 6   | 2543 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/h1706

find . -type d -name "trim" | xargs rm -fr
find . -type d -name "mapping" | xargs rm -fr
# find . -type d -name "depth" | xargs rm -fr

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

