# *Glycine max* Williams 82 倍数因子测试


[TOC levels=1-3]: # ""

- [*Glycine max* Williams 82 倍数因子测试](#glycine-max-williams-82-倍数因子测试)
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

* 因子值 0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64
* 测试文件 [SRR10296600](https://www.ncbi.nlm.nih.gov/sra/SRX7009428)
  * 覆盖度 48,633,106,500 / 949,738,161 = 51

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
DEPTH=50
for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    CUTOFF=$(bc <<< "(${DEPTH} * ${FOLD}) / 1")
    ARRAY+=("${FOLD}::${CUTOFF}")
done
echo "${ARRAY[@]}"
#0::0 0.25::12 0.5::25 1::50 2::100 4::200 8::400 16::800 32::1600 64::3200

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
            --filter "adapter artifact"

        bsub -q mpi -n 24 -J "${BASE_NAME}" '
            bash 2_trim.sh
        '        
    else
        anchr template \
            --genome 949738161 \
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

        bsub -q mpi -n 24 -J "${BASE_NAME}" '
            bash 2_trim.sh
        '
    fi

    popd

done

# find . -type f -wholename "*trim/R*.fq.gz"

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

## Mapping

```shell script
cd ~/data/plastid/evaluation/w82

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${FOLD}
    
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
                -x ../../../../genome/w82/genome.fa \
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
cd ~/data/plastid/evaluation/w82

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${FOLD}
    
    echo "==> ${BASE_NAME}"

    cat ${BASE_NAME}/mapping/output.* |
        grep -Pzo '(?s)Multiseed full.+overall alignment rate'

done

```

```text
==> SRR10296600_0.5
Multiseed full-index search: 00:59:56
266230345 reads; of these:
  266230345 (100.00%) were unpaired; of these:
    5588122 (2.10%) aligned 0 times
    113971603 (42.81%) aligned exactly 1 time
    146670620 (55.09%) aligned >1 times
97.90% overall alignment rate
==> SRR10296600_1
Multiseed full-index search: 00:38:44
138424613 reads; of these:
  138424613 (100.00%) were unpaired; of these:
    4007832 (2.90%) aligned 0 times
    17178829 (12.41%) aligned exactly 1 time
    117237952 (84.69%) aligned >1 times
97.10% overall alignment rate
==> SRR10296600_2
Multiseed full-index search: 00:33:16
117315224 reads; of these:
  117315224 (100.00%) were unpaired; of these:
    3776062 (3.22%) aligned 0 times
    12605868 (10.75%) aligned exactly 1 time
    100933294 (86.04%) aligned >1 times
96.78% overall alignment rate

```

## Depth

* Depth via `mosdepth`

```shell script
cd ~/data/plastid/evaluation/w82

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${FOLD}
    
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
cd ~/data/plastid/evaluation/w82

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${FOLD}
    
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
        
    spanr stat ../../../../genome/w82/chr.sizes covered.yml -o stdout |
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
cd ~/data/plastid/evaluation/w82

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${FOLD}
    
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
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth > /dev/null

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
| 0    | Nc    | 949183385 | 929976190 | 0.9798  | 35678243492 | 37.59 | 0   | 47362 |
| 0.25 | Nc    | 949183385 | 929709716 | 0.9795  | 35616400666 | 37.52 | 0   | 47289 |
| 0.5  | Nc    | 949183385 | 919870612 | 0.9691  | 34228418977 | 36.06 | 0   | 47279 |
| 1    | Nc    | 949183385 | 496767108 | 0.5234  | 16068387605 | 16.93 | 0   | 47250 |
| 2    | Nc    | 949183385 | 376680778 | 0.3968  | 13059515516 | 13.76 | 0   | 47180 |
| 4    | Nc    | 949183385 | 332662176 | 0.3505  | 11489481456 | 12.10 | 0   | 47071 |
| 8    | Nc    | 949183385 | 299776763 | 0.3158  | 10162316738 | 10.71 | 0   | 46779 |
| 16   | Nc    | 949183385 | 268631500 | 0.2830  | 8852676948  | 9.33  | 0   | 45931 |
| 32   | Nc    | 949183385 | 234774194 | 0.2473  | 7471267400  | 7.87  | 0   | 44314 |
| 64   | Nc    | 949183385 | 199579698 | 0.2103  | 6154626907  | 6.48  | 0   | 41854 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:-------|:----|:-----|
| 0    | Mt    | 402558    | 389025    | 0.9664  | 387717396 | 963.13 | 0   | 2293 |
| 0.25 | Mt    | 402558    | 389025    | 0.9664  | 387105433 | 961.61 | 0   | 2291 |
| 0.5  | Mt    | 402558    | 389025    | 0.9664  | 387099329 | 961.60 | 0   | 2291 |
| 1    | Mt    | 402558    | 388914    | 0.9661  | 387051379 | 961.48 | 0   | 2291 |
| 2    | Mt    | 402558    | 388914    | 0.9661  | 387038288 | 961.45 | 0   | 2291 |
| 4    | Mt    | 402558    | 388914    | 0.9661  | 387038090 | 961.45 | 0   | 2291 |
| 8    | Mt    | 402558    | 388914    | 0.9661  | 387036160 | 961.44 | 0   | 2291 |
| 16   | Mt    | 402558    | 367645    | 0.9133  | 331422533 | 823.29 | 0   | 2291 |
| 32   | Mt    | 402558    | 167533    | 0.4162  | 59693089  | 148.28 | 0   | 2270 |
| 64   | Mt    | 402558    | 100878    | 0.2506  | 7817654   | 19.42  | 0   | 1977 |


| Fold | chrom | chrLength | covLength | covRate | bases      | mean     | min | max   |
|:-----|:------|:----------|:----------|:--------|:-----------|:---------|:----|:------|
| 0    | Pt    | 152218    | 152218    | 1.0000  | 2730112733 | 17935.54 | 6   | 23677 |
| 0.25 | Pt    | 152218    | 152218    | 1.0000  | 2728573539 | 17925.43 | 6   | 23660 |
| 0.5  | Pt    | 152218    | 152218    | 1.0000  | 2726619641 | 17912.60 | 6   | 23640 |
| 1    | Pt    | 152218    | 152218    | 1.0000  | 2725438049 | 17904.83 | 6   | 23635 |
| 2    | Pt    | 152218    | 152218    | 1.0000  | 2725048751 | 17902.28 | 6   | 23635 |
| 4    | Pt    | 152218    | 152218    | 1.0000  | 2724946237 | 17901.60 | 6   | 23635 |
| 8    | Pt    | 152218    | 152218    | 1.0000  | 2724905760 | 17901.34 | 6   | 23635 |
| 16   | Pt    | 152218    | 152218    | 1.0000  | 2724887039 | 17901.21 | 6   | 23635 |
| 32   | Pt    | 152218    | 152218    | 1.0000  | 2724870977 | 17901.11 | 6   | 23635 |
| 64   | Pt    | 152218    | 152218    | 1.0000  | 2724865865 | 17901.08 | 6   | 23635 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/w82

find . -type d -name "trim" | xargs rm -fr
find . -type d -name "mapping" | xargs rm -fr
# find . -type d -name "depth" | xargs rm -fr

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

