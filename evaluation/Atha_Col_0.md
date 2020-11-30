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
brew install bwa

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
    
    bsub -q mpi -n 24 -J "${BASE_NAME}-mapping" '

        # Pipe all reads together as we do not need mate info
        gzip -dcf \
            ../2_illumina/trim/Q25L60/R1.fq.gz \
            ../2_illumina/trim/Q25L60/R2.fq.gz \
            ../2_illumina/trim/Q25L60/Rs.fq.gz |
            bwa mem -M -t 20 \
                ../../../../genome/col_0/genome.fa \
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
cd ~/data/plastid/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR616966_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth

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
            if ($start == $F[2]) {
                print qq($F[0]:$start);
            }
            else {
                print qq($F[0]:$start-$F[2]);
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
| 0    | Nc    | 119146348 | 118886986 | 0.9978  | 2866772031 | 24.06 | 0   | 110517 |
| 0.25 | Nc    | 119146348 | 118741324 | 0.9966  | 2848469903 | 23.91 | 0   | 110633 |
| 0.5  | Nc    | 119146348 | 78550490  | 0.6593  | 1438110229 | 12.07 | 0   | 110731 |
| 1    | Nc    | 119146348 | 14305813  | 0.1201  | 619823713  | 5.20  | 0   | 110584 |
| 2    | Nc    | 119146348 | 9518949   | 0.0799  | 557173606  | 4.68  | 0   | 110505 |
| 4    | Nc    | 119146348 | 6785818   | 0.0570  | 515540352  | 4.33  | 0   | 110641 |
| 8    | Nc    | 119146348 | 4074559   | 0.0342  | 452938381  | 3.80  | 0   | 110311 |
| 16   | Nc    | 119146348 | 2821962   | 0.0237  | 427477095  | 3.59  | 0   | 110344 |
| 32   | Nc    | 119146348 | 1813743   | 0.0152  | 405003054  | 3.40  | 0   | 109429 |
| 64   | Nc    | 119146348 | 1233227   | 0.0104  | 388257853  | 3.26  | 0   | 107774 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Mt    | 366924    | 363942    | 0.9919  | 58950023 | 160.66 | 0   | 1916 |
| 0.25 | Mt    | 366924    | 363819    | 0.9915  | 58793180 | 160.23 | 0   | 1920 |
| 0.5  | Mt    | 366924    | 364006    | 0.9920  | 58751358 | 160.12 | 0   | 1905 |
| 1    | Mt    | 366924    | 363815    | 0.9915  | 58760072 | 160.14 | 0   | 1948 |
| 2    | Mt    | 366924    | 363524    | 0.9907  | 58719881 | 160.03 | 0   | 1929 |
| 4    | Mt    | 366924    | 358250    | 0.9764  | 55570271 | 151.45 | 0   | 1917 |
| 8    | Mt    | 366924    | 98518     | 0.2685  | 7948565  | 21.66  | 0   | 1930 |
| 16   | Mt    | 366924    | 30839     | 0.0840  | 1369424  | 3.73   | 0   | 1904 |
| 32   | Mt    | 366924    | 27902     | 0.0760  | 1291105  | 3.52   | 0   | 1928 |
| 64   | Mt    | 366924    | 20563     | 0.0560  | 1174510  | 3.20   | 0   | 1945 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean    | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:--------|:----|:-----|
| 0    | Pt    | 154478    | 154478    | 1.0000  | 483357613 | 3128.97 | 4   | 4507 |
| 0.25 | Pt    | 154478    | 154478    | 1.0000  | 482266768 | 3121.91 | 4   | 4504 |
| 0.5  | Pt    | 154478    | 154478    | 1.0000  | 482118136 | 3120.95 | 4   | 4504 |
| 1    | Pt    | 154478    | 154478    | 1.0000  | 482067384 | 3120.62 | 4   | 4504 |
| 2    | Pt    | 154478    | 154478    | 1.0000  | 481988086 | 3120.11 | 4   | 4504 |
| 4    | Pt    | 154478    | 154478    | 1.0000  | 481952699 | 3119.88 | 4   | 4504 |
| 8    | Pt    | 154478    | 154478    | 1.0000  | 481915822 | 3119.64 | 4   | 4504 |
| 16   | Pt    | 154478    | 154478    | 1.0000  | 481919945 | 3119.67 | 4   | 4504 |
| 32   | Pt    | 154478    | 154478    | 1.0000  | 481755638 | 3118.60 | 2   | 4504 |
| 64   | Pt    | 154478    | 150494    | 0.9742  | 352816898 | 2283.93 | 0   | 4498 |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/col_0

find . -type d -name "trim" | xargs rm -fr
find . -type d -name "mapping" | xargs rm -fr

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

