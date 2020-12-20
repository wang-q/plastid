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

