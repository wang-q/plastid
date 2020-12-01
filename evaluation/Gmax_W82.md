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

+ 因子值 0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64
+ 测试文件 [SRR10296600](https://www.ncbi.nlm.nih.gov/sra/SRX7009428)
+ 覆盖度 48,633,106,500 / 949,738,161 = 51

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

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/w82
cd ~/data/plastid/evaluation/w82

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${NAME}
    
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
ARRAY=(
    '0::0'
    '0.25::12'
    '0.5::25'
    '1::50'
    '2::100'
    '4::200'
    '8::400'
    '16::800'
    '32::1600'
    '64::3200'
)

for item in "${ARRAY[@]}" ; do
    NAME="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"
    
    BASE_NAME=SRR10296600_${NAME}
    pushd ${BASE_NAME}
    
    rm *.sh
        
    if [[ "${NAME}" == "0" ]]; then
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

        bsub -q mpi -n 24 -J "${BASE_NAME}" "
            bash 2_trim.sh
        "        
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

        bsub -q mpi -n 24 -J "${BASE_NAME}" "
            bash 2_trim.sh
        "
    fi

    popd

done

# find . -type f -wholename "*trim/R*.fq.gz"

```

## `kat hist` and `kat gcp`

```shell script
cd ~/data/plastid/evaluation/w82

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${NAME}
    
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

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${NAME}
    
    mkdir -p ${BASE_NAME}/mapping
    pushd ${BASE_NAME}/mapping
    
    if [ -f R.sort.bai ]; then
        echo >&2 '    R.sort.bai already presents'
        popd;
        continue;
    fi

    bsub -q mpi -n 24 -J "${BASE_NAME}-mapping" '
        export JAVA_OPTS="-Xmx20G"

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


## Depth

* Depth via `mosdepth`

```shell script
cd ~/data/plastid/evaluation/w82

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${NAME}
    
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

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${NAME}
    
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

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${NAME}
    
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
cd ~/data/plastid/evaluation/w82

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth > /dev/null

    echo -e "Fold\tchrom\n${NAME}\tNc\n${NAME}\tMt\n${NAME}\tPt" |
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

## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/w82

find . -type d -name "trim" | xargs rm -fr
find . -type d -name "mapping" | xargs rm -fr

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

