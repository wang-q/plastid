# *Prunus persica* Lovell 倍数因子测试


[TOC levels=1-3]: # ""

- [*Prunus persica* Lovell 倍数因子测试](#prunus-persica-lovell-倍数因子测试)
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
+ 测试文件 [SRR502985](https://www.ncbi.nlm.nih.gov/sra/SRX150254)
+ 覆盖度 110

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
            --genome 225852601 \
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

## Mapping

```shell script
cd ~/data/plastid/evaluation/lovell

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR502985_${FOLD}
    
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
                -x ../../../../genome/lovell/genome.fa \
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
cd ~/data/plastid/evaluation/lovell

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR502985_${FOLD}
    
    echo "==> ${BASE_NAME}"

    cat ${BASE_NAME}/mapping/output.* |
        grep -Pzo '(?s)Multiseed full.+overall alignment rate'

done

```

```text
==> SRR502985_0
Multiseed full-index search: 00:27:47
195725737 reads; of these:
  195725737 (100.00%) were unpaired; of these:
    6375330 (3.26%) aligned 0 times
    112834060 (57.65%) aligned exactly 1 time
    76516347 (39.09%) aligned >1 times
96.74% overall alignment rate
==> SRR502985_0.25
Multiseed full-index search: 00:27:44
194375766 reads; of these:
  194375766 (100.00%) were unpaired; of these:
    5862073 (3.02%) aligned 0 times
    112136516 (57.69%) aligned exactly 1 time
    76377177 (39.29%) aligned >1 times
96.98% overall alignment rate
==> SRR502985_0.5
Multiseed full-index search: 00:24:55
174597854 reads; of these:
  174597854 (100.00%) were unpaired; of these:
    5496363 (3.15%) aligned 0 times
    94688160 (54.23%) aligned exactly 1 time
    74413331 (42.62%) aligned >1 times
96.85% overall alignment rate
==> SRR502985_1
Multiseed full-index search: 00:12:15
73213603 reads; of these:
  73213603 (100.00%) were unpaired; of these:
    4704514 (6.43%) aligned 0 times
    8344017 (11.40%) aligned exactly 1 time
    60165072 (82.18%) aligned >1 times
93.57% overall alignment rate
==> SRR502985_2
Multiseed full-index search: 00:10:13
59790266 reads; of these:
  59790266 (100.00%) were unpaired; of these:
    4561251 (7.63%) aligned 0 times
    5776775 (9.66%) aligned exactly 1 time
    49452240 (82.71%) aligned >1 times
92.37% overall alignment rate
==> SRR502985_4
Multiseed full-index search: 00:09:01
52321041 reads; of these:
  52321041 (100.00%) were unpaired; of these:
    4480459 (8.56%) aligned 0 times
    4932088 (9.43%) aligned exactly 1 time
    42908494 (82.01%) aligned >1 times
91.44% overall alignment rate
==> SRR502985_8
Multiseed full-index search: 00:08:00
44767296 reads; of these:
  44767296 (100.00%) were unpaired; of these:
    3342570 (7.47%) aligned 0 times
    4262806 (9.52%) aligned exactly 1 time
    37161920 (83.01%) aligned >1 times
92.53% overall alignment rate
==> SRR502985_16
Multiseed full-index search: 00:06:34
35650781 reads; of these:
  35650781 (100.00%) were unpaired; of these:
    723812 (2.03%) aligned 0 times
    3392122 (9.51%) aligned exactly 1 time
    31534847 (88.45%) aligned >1 times
97.97% overall alignment rate
==> SRR502985_32
Multiseed full-index search: 00:05:55
28895551 reads; of these:
  28895551 (100.00%) were unpaired; of these:
    261684 (0.91%) aligned 0 times
    2543038 (8.80%) aligned exactly 1 time
    26090829 (90.29%) aligned >1 times
99.09% overall alignment rate
==> SRR502985_64
Multiseed full-index search: 00:04:51
21676327 reads; of these:
  21676327 (100.00%) were unpaired; of these:
    181675 (0.84%) aligned 0 times
    364129 (1.68%) aligned exactly 1 time
    21130523 (97.48%) aligned >1 times
99.16% overall alignment rate

```

## Depth

* Depth via `mosdepth`

```shell script
cd ~/data/plastid/evaluation/lovell

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR502985_${FOLD}
    
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
cd ~/data/plastid/evaluation/lovell

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR502985_${FOLD}
    
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
        
    spanr stat ../../../../genome/lovell/chr.sizes covered.yml -o stdout |
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
cd ~/data/plastid/evaluation/lovell

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR502985_${FOLD}
    
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
cd ~/data/plastid/evaluation/lovell

for FOLD in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR502985_${FOLD}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth > /dev/null

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
| 0    | Nc    | 225694811 | 222818605 | 0.9873  | 17829586428 | 79.00 | 0   | 16993 |
| 0.25 | Nc    | 225694811 | 222551398 | 0.9861  | 17750937077 | 78.65 | 0   | 16981 |
| 0.5  | Nc    | 225694811 | 217501787 | 0.9637  | 15854226925 | 70.25 | 0   | 16972 |
| 1    | Nc    | 225694811 | 120139806 | 0.5323  | 5973794768  | 26.47 | 0   | 16971 |
| 2    | Nc    | 225694811 | 94106754  | 0.4170  | 4670093957  | 20.69 | 0   | 16971 |
| 4    | Nc    | 225694811 | 78436242  | 0.3475  | 3945727928  | 17.48 | 0   | 16969 |
| 8    | Nc    | 225694811 | 65717918  | 0.2912  | 3318907037  | 14.71 | 0   | 16968 |
| 16   | Nc    | 225694811 | 53064727  | 0.2351  | 2690363089  | 11.92 | 0   | 16968 |
| 32   | Nc    | 225694811 | 40578290  | 0.1798  | 2164801488  | 9.59  | 0   | 16968 |
| 64   | Nc    | 225694811 | 29442059  | 0.1305  | 1794716029  | 7.95  | 0   | 16968 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean    | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:--------|:----|:-----|
| 0    | Pt    | 157790    | 157790    | 1.0000  | 717395489 | 4546.52 | 14  | 9629 |
| 0.25 | Pt    | 157790    | 157790    | 1.0000  | 716514879 | 4540.94 | 6   | 9618 |
| 0.5  | Pt    | 157790    | 157790    | 1.0000  | 716441195 | 4540.47 | 5   | 9618 |
| 1    | Pt    | 157790    | 157790    | 1.0000  | 716350026 | 4539.89 | 3   | 9618 |
| 2    | Pt    | 157790    | 157790    | 1.0000  | 716294603 | 4539.54 | 2   | 9618 |
| 4    | Pt    | 157790    | 157647    | 0.9991  | 716072755 | 4538.14 | 0   | 9618 |
| 8    | Pt    | 157790    | 157414    | 0.9976  | 714571960 | 4528.63 | 0   | 9618 |
| 16   | Pt    | 157790    | 157180    | 0.9961  | 707413558 | 4483.26 | 0   | 9618 |
| 32   | Pt    | 157790    | 155947    | 0.9883  | 616861680 | 3909.38 | 0   | 9618 |
| 64   | Pt    | 157790    | 145972    | 0.9251  | 287832288 | 1824.15 | 0   | 9555 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/lovell

find . -type d -name "trim" | xargs rm -fr
find . -type d -name "mapping" | xargs rm -fr
# find . -type d -name "depth" | xargs rm -fr

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

