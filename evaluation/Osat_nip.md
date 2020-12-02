# *Oryza sativa* Nipponbare 倍数因子测试


[TOC levels=1-3]: # ""

- [*Oryza sativa* Nipponbare 倍数因子测试](#oryza-sativa-nipponbare-倍数因子测试)
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
+ 测试文件 [SRR545231](https://www.ncbi.nlm.nih.gov/sra/SRX179254)
+ 覆盖度 17,220,721,594 / 373,870,564 = 46

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
            --genome 373870564 \
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

# find . -type f -wholename "*trim/R*.fq.gz"

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

## Mapping

```shell script
cd ~/data/plastid/evaluation/nip

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${FOLD}
    
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
                -x ../../../../genome/nip/genome.fa \
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
cd ~/data/plastid/evaluation/nip

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${FOLD}
    
    echo "==> ${BASE_NAME}"

    cat ${BASE_NAME}/mapping/output.* |
        grep -Pzo '(?s)Multiseed full.+overall alignment rate'

done

```

```text
==> SRR545231_0
Multiseed full-index search: 00:19:07
142307290 reads; of these:
  142307290 (100.00%) were unpaired; of these:
    493206 (0.35%) aligned 0 times
    91963290 (64.62%) aligned exactly 1 time
    49850794 (35.03%) aligned >1 times
99.65% overall alignment rate
==> SRR545231_0.25
Multiseed full-index search: 00:18:14
141947630 reads; of these:
  141947630 (100.00%) were unpaired; of these:
    434089 (0.31%) aligned 0 times
    91709630 (64.61%) aligned exactly 1 time
    49803911 (35.09%) aligned >1 times
99.69% overall alignment rate
==> SRR545231_0.5
Multiseed full-index search: 00:18:49
138070556 reads; of these:
  138070556 (100.00%) were unpaired; of these:
    419977 (0.30%) aligned 0 times
    88339753 (63.98%) aligned exactly 1 time
    49310826 (35.71%) aligned >1 times
99.70% overall alignment rate
==> SRR545231_1
Multiseed full-index search: 00:07:43
49909647 reads; of these:
  49909647 (100.00%) were unpaired; of these:
    138820 (0.28%) aligned 0 times
    8920577 (17.87%) aligned exactly 1 time
    40850250 (81.85%) aligned >1 times
99.72% overall alignment rate
==> SRR545231_2
Multiseed full-index search: 00:05:48
35846086 reads; of these:
  35846086 (100.00%) were unpaired; of these:
    106888 (0.30%) aligned 0 times
    3191488 (8.90%) aligned exactly 1 time
    32547710 (90.80%) aligned >1 times
99.70% overall alignment rate
==> SRR545231_4
Multiseed full-index search: 00:04:50
29281184 reads; of these:
  29281184 (100.00%) were unpaired; of these:
    97082 (0.33%) aligned 0 times
    2099977 (7.17%) aligned exactly 1 time
    27084125 (92.50%) aligned >1 times
99.67% overall alignment rate
==> SRR545231_8
Multiseed full-index search: 00:03:51
24250163 reads; of these:
  24250163 (100.00%) were unpaired; of these:
    87053 (0.36%) aligned 0 times
    1408351 (5.81%) aligned exactly 1 time
    22754759 (93.83%) aligned >1 times
99.64% overall alignment rate
==> SRR545231_16
Multiseed full-index search: 00:03:25
19986642 reads; of these:
  19986642 (100.00%) were unpaired; of these:
    78712 (0.39%) aligned 0 times
    940042 (4.70%) aligned exactly 1 time
    18967888 (94.90%) aligned >1 times
99.61% overall alignment rate
==> SRR545231_32
Multiseed full-index search: 00:02:40
14609798 reads; of these:
  14609798 (100.00%) were unpaired; of these:
    70443 (0.48%) aligned 0 times
    591576 (4.05%) aligned exactly 1 time
    13947779 (95.47%) aligned >1 times
99.52% overall alignment rate
==> SRR545231_64
Multiseed full-index search: 00:01:52
9137019 reads; of these:
  9137019 (100.00%) were unpaired; of these:
    60855 (0.67%) aligned 0 times
    335178 (3.67%) aligned exactly 1 time
    8740986 (95.67%) aligned >1 times
99.33% overall alignment rate

```

## Depth

* Depth via `mosdepth`

```shell script
cd ~/data/plastid/evaluation/nip

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${FOLD}
    
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
cd ~/data/plastid/evaluation/nip

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${FOLD}
    
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
        
    spanr stat ../../../../genome/nip/chr.sizes covered.yml -o stdout |
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
cd ~/data/plastid/evaluation/nip

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${FOLD}
    
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
cd ~/data/plastid/evaluation/nip

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${FOLD}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth > /dev/null

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

| Fold | chrom | chrLength | covLength | covRate | bases       | mean  | min | max   |
|:-----|:------|:----------|:----------|:--------|:------------|:------|:----|:------|
| 0    | Nc    | 373245519 | 368880124 | 0.9883  | 13560236257 | 36.33 | 0   | 26925 |
| 0.25 | Nc    | 373245519 | 366974601 | 0.9832  | 13537010581 | 36.27 | 0   | 26925 |
| 0.5  | Nc    | 373245519 | 356485614 | 0.9551  | 13202521828 | 35.37 | 0   | 26925 |
| 1    | Nc    | 373245519 | 174583986 | 0.4677  | 4686493369  | 12.56 | 0   | 26925 |
| 2    | Nc    | 373245519 | 121424118 | 0.3253  | 3323311688  | 8.90  | 0   | 26922 |
| 4    | Nc    | 373245519 | 99390325  | 0.2663  | 2705233257  | 7.25  | 0   | 26913 |
| 8    | Nc    | 373245519 | 83242987  | 0.2230  | 2259422275  | 6.05  | 0   | 26893 |
| 16   | Nc    | 373245519 | 68162616  | 0.1826  | 1853810568  | 4.97  | 0   | 26860 |
| 32   | Nc    | 373245519 | 51786399  | 0.1387  | 1360926235  | 3.65  | 0   | 26769 |
| 64   | Nc    | 373245519 | 34821148  | 0.0933  | 863762280   | 2.31  | 0   | 26473 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:----|
| 0    | Mt    | 490520    | 488275    | 0.9954  | 62173965 | 126.75 | 0   | 504 |
| 0.25 | Mt    | 490520    | 488275    | 0.9954  | 62170129 | 126.74 | 0   | 504 |
| 0.5  | Mt    | 490520    | 488275    | 0.9954  | 62170062 | 126.74 | 0   | 504 |
| 1    | Mt    | 490520    | 488275    | 0.9954  | 62167182 | 126.74 | 0   | 504 |
| 2    | Mt    | 490520    | 488275    | 0.9954  | 62143743 | 126.69 | 0   | 504 |
| 4    | Mt    | 490520    | 434688    | 0.8862  | 46167804 | 94.12  | 0   | 504 |
| 8    | Mt    | 490520    | 106603    | 0.2173  | 8021224  | 16.35  | 0   | 504 |
| 16   | Mt    | 490520    | 59121     | 0.1205  | 5024311  | 10.24  | 0   | 504 |
| 32   | Mt    | 490520    | 43100     | 0.0879  | 4254130  | 8.67   | 0   | 489 |
| 64   | Mt    | 490520    | 9695      | 0.0198  | 532392   | 1.09   | 0   | 373 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Pt    | 134525    | 133022    | 0.9888  | 49789479 | 370.11 | 0   | 1885 |
| 0.25 | Pt    | 134525    | 133022    | 0.9888  | 49787918 | 370.10 | 0   | 1885 |
| 0.5  | Pt    | 134525    | 133022    | 0.9888  | 49787918 | 370.10 | 0   | 1885 |
| 1    | Pt    | 134525    | 133022    | 0.9888  | 49786736 | 370.09 | 0   | 1885 |
| 2    | Pt    | 134525    | 133022    | 0.9888  | 49786736 | 370.09 | 0   | 1885 |
| 4    | Pt    | 134525    | 133022    | 0.9888  | 49786674 | 370.09 | 0   | 1885 |
| 8    | Pt    | 134525    | 133022    | 0.9888  | 49786173 | 370.09 | 0   | 1885 |
| 16   | Pt    | 134525    | 133022    | 0.9888  | 49785068 | 370.08 | 0   | 1885 |
| 32   | Pt    | 134525    | 116368    | 0.8650  | 26754146 | 198.88 | 0   | 1882 |
| 64   | Pt    | 134525    | 19414     | 0.1443  | 1993578  | 14.82  | 0   | 423  |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/nip

find . -type d -name "trim" | xargs rm -fr
find . -type d -name "mapping" | xargs rm -fr
# find . -type d -name "depth" | xargs rm -fr

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

